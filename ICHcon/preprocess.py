import os
import re
import subprocess
import numpy as np
import ants
import SimpleITK as sitk
import nibabel as nib
import skimage
from scipy import ndimage
from fsl.data.image import Image
from fsl.utils.image.resample import resample
from fsl.wrappers import bet
from fsl.wrappers import fslmaths


def bash_in_python(cmd):
    out, err = unix_cmd(cmd)
    show_error(err)
    return out, err


def show_error(err):
    if len(err) > 0:
        print(err)


def unix_cmd(cmd):
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         shell=True)
    out, err = p.communicate()
    return out, err


def read_ants(img_path):
    img = sitk.ReadImage(img_path)
    direction = np.array(img.GetDirection()).reshape(3, 3)
    spacing = img.GetSpacing()
    origin = img.GetOrigin()
    np_img = sitk.GetArrayFromImage(img).astype(float)
    img_ant = ants.from_numpy(np.transpose(np_img, (2, 1, 0)),
                              origin=origin,
                              spacing=spacing, direction=direction)
    return img_ant


def select_depth(D, depths):
    if D > depths:
        start_depth = (D - depths)//2
        return start_depth, D - (start_depth + depths)
    else:
        return 0, 0


def sel_central(img, depths):
    if img.shape[2] > depths:
        start_depth = (img.shape[2] - depths)//2
        return img[..., start_depth:(start_depth + depths)]
    else:
        return img


def skull_strip(img_path, save_path, window=[0, 100], sigma=0):
    ori_ob = nib.load(img_path)
    ori_img = ori_ob.get_fdata()
    if len(ori_img.shape) == 4:
        ori_img = ori_img[..., 0]
    if sigma != 0:
        img = skimage.filters.gaussian(ori_img, sigma=sigma)
    else:
        img = ori_img

    lab_list = []
    for i in range(img.shape[2]):
        # thresholding and erosion
        mask = (img[:, :, i] > window[0]).astype(float) + \
                (img[:, :, i] < window[1]).astype(float)
        mask = skimage.morphology.binary_erosion(mask == 2)

        # largest connected component
        labels = skimage.measure.label(mask)
        lab_count = np.bincount(labels.flat)
        if len(lab_count) == 1:
            LCC = np.zeros(labels.shape)
        else:
            LCC = (labels == np.argmax(lab_count[1:]) + 1)

        # dilation and fill holes
        LCC = skimage.morphology.binary_dilation(LCC)
        lab_list.append(ndimage.binary_fill_holes(LCC))

    masks = (np.stack(lab_list, axis=-1)).astype(int)
    save_img = np.round(ori_img*masks)
    save_ob = nib.Nifti1Image(save_img.clip(window[0], window[1]),
                              header=ori_ob.header,
                              affine=ori_ob.affine, dtype=np.int16)
    save_ob.to_filename(save_path)


def preprocess_ct(img_path, out_path, target_res=5):
    img = Image(img_path)
    res = img.pixdim[2]

    # sometimes the pixdim is inaccurate
    # it is not possible for head CT to be smaller than 50 slices yet having
    # 1mm resolution
    abort = img.shape[2] < 50 and res < 1
    if not abort:
        new_shape = [*img.shape[:2], int(img.shape[2]*res/target_res)]
        new_img, affine = resample(img, new_shape)
        new_img = remove_blank(new_img)
        new_img = nib.Nifti1Image(new_img, affine=affine)
        new_img.to_filename(out_path)


def empty_index(vec):
    indices = np.where(vec > 1)[0]
    return [min(indices), max(indices)]


def apply_square_bbox(img, square=True):
    if len(img.shape) == 4:
        img = img.sum(3)
    empty_x = empty_index(img.sum(2).sum(1))
    empty_y = empty_index(img.sum(2).sum(0))
    empty_z = empty_index(img.sum(1).sum(0))
    start_xy = min(empty_x[0], empty_y[0])
    end_xy = max(empty_x[1], empty_y[1])
    if square:
        return img[start_xy:end_xy, start_xy:end_xy, empty_z[0]:empty_z[1]]
    else:
        return img[empty_x[0]:empty_x[1], empty_y[0]:empty_y[1]]


def remove_blank(img, HU_thres=1, thres=0.1, return_index=False):
    '''
    Args:
        `HU_thres`: Hounsfield unit above which a pixel is counted as being
        inside brain tissues
        `thres`: percentage of non-zero pixel below which a particular slice is
        removed
    '''
    red_img = apply_square_bbox(np.round(img), square=False)
    H, W = red_img.shape[:2]
    z_perc = (red_img > HU_thres).sum(axis=(0, 1))/(H*W)
    z_perc = (z_perc > thres).astype(float)
    # reduce the threshold if the image is too small
    # if z_perc.max() == 0:
    #     z_perc = (img > HU_thres).sum(axis=(0, 1))/(H*W)
    #     z_perc = (z_perc > thres/10).astype(float)

    N = len(z_perc)//2
    neck = np.where(z_perc[:N] == 0)[0]
    start = max(neck) if len(neck) > 0 else 0
    vertex = np.where(z_perc[N:] == 0)[0]
    end = min(vertex) if len(vertex) > 0 else len(z_perc[N:])

    if return_index:
        return [start+1, end+N]
    else:
        return img[..., (start+1):(end+N)]


def append_str(str_name, append_name):
    filename = os.path.basename(str_name)
    first_name = filename.split(".")[0]
    return re.sub("/"+first_name+"\\.",
                  "/"+first_name+"_"+append_name+".", str_name)


def brain(img_path, save_path):
    """
    Brain Extraction with FSL
    Params:
    - image: nifti object, scan to brain extract
    Output:
    - brain_image: nifti object, extracted brain
    """
    image = nib.load(img_path)
    affine = image.affine
    header = image.header
    tmpfile = append_str(img_path, "tmp")
    image.to_filename(tmpfile)

    # FSL calls
    ID = ''.join([str(i) for i in np.random.choice(9, 10)])
    tmp_dir = os.path.dirname(save_path)+"/"+ID
    os.mkdir(tmp_dir)

    # bash_in_python(
    #     "fslmaths {} -thr 0 -uthr 100".format(img_path) +
    #     " -bin -fillh {}/mask.nii.gz;".format(tmp_dir) +
    #     "fslmaths {} -mas {}/mask.nii.gz".format(img_path, tmp_dir) +
    #     " {}/tmp.nii.gz;".format(tmp_dir) +
    #     "bet {}/tmp.nii.gz {}/tmp.nii.gz -f 0.01;".format(tmp_dir, tmp_dir) +
    #     "fslmaths {}/tmp.nii.gz -bin".format(tmp_dir) +
    #     " -fillh {}/img.nii.gz;".format(tmp_dir) +
    #     "fslmaths {}"
    #     )
    mask = fslmaths(image).thr("0.000000").uthr("100.000000").bin().fillh().run()
    fslmaths(image).mas(mask).run(tmpfile)
    bet(tmpfile, tmpfile, fracintensity=0.01)
    mask = fslmaths(tmpfile).bin().fillh().run()
    image = fslmaths(image).mas(mask).run()
    image = nib.Nifti1Image(image.get_fdata(), affine, header, dtype=np.int16)
    os.remove(tmpfile)
    image.to_filename(save_path)


def preprocess_imgs(ct_path, save_path, skullstrip=True, preprocess=True):
    if not os.path.exists(save_path):
        ID = ''.join([str(i) for i in np.random.choice(9, 10)])
        tmp_path = re.sub(".nii.gz", ID+".nii.gz", save_path)

        skull_strip(ct_path, tmp_path)
        preprocess_ct(tmp_path, tmp_path)
        bash_in_python("mri_synthstrip -i {}".format(tmp_path) +
                       " -o {}".format(save_path))
        os.remove(tmp_path)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--ct_dir', type=str)
    parser.add_argument('--save_dir', type=str)

    parser.add_argument("--num", type=int, default=1)
    parser.add_argument("--arrayID", type=int, default=0)
    args = parser.parse_args()

    for index, i in enumerate(sorted(os.listdir(args.ct_dir))):
        if index % args.num == args.arrayID:
            preprocess_imgs(args.ct_dir+"/"+i, args.save_dir+"/"+i)
