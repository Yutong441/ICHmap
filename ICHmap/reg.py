# register lesion mask to MNI space
import os
import shutil
import numpy as np
import ants
import nibabel as nib
from fsl.wrappers import fslmaths
from fsl.data.image import Image
from fsl.utils.image.resample import resample
import SimpleITK as sitk
import skimage
from scipy import ndimage


def read_ants(img_path):
    '''
    Read ants using simple itk

    In reading nifti images from ants, this error may be encountered:
    "If antspy throws in error: ITK ERROR: ITK only supports orthonormal
    direction cosines. No orthonormal definition found!"

    The solution is to install SimpleITK version 2.0.2. However, ants does not
    seem to benefit from this version change and keeps throwing off this error.
    This function uses SimpleITK to load nifti image and convert into ants
    '''
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


def register(ct_path, lesion_path, MNI, save_path):
    CT = ants.image_read(ct_path)
    lesion = ants.image_read(lesion_path)

    # no_les = ants.mask_image(CT, 1 - lesion)
    tx = ants.registration(MNI, CT, type_of_transform="SyN")
    lesion_MNI = ants.apply_transforms(MNI, lesion, tx["fwdtransforms"],
                                       interpolator="genericLabel")

    ants.image_write(lesion_MNI, save_path)
    # ants.image_write(tx["warpedmovout"],
    #                  re.sub(".nii.gz$", "_img.nii.gz", save_path))
    # for i in tx["fwdtransforms"]:
    #     os.remove(i)


def skull_strip(img_path, save_path, window=[0, 100], sigma=0):
    ori_ob = nib.load(img_path)
    ori_img = ori_ob.get_fdata()
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
    save_ob = nib.Nifti1Image(np.round(ori_img*masks), header=ori_ob.header,
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
        new_img = nib.Nifti1Image(new_img, affine=affine)
        new_img.to_filename(out_path)


def register_lesion(ct_path, lesion_path, MNI, MNI_mask, save_path,
                    sel_depth=20, skullstrip=True, preprocess=True):
    '''
    Register ICH volume mask to MNI512
    decrease the effect of hemorrhage on registration:
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4846589/
    original paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2938189/

    If antspy throws in error: ITK ERROR: ITK only supports orthonormal
    direction cosines. No orthonormal definition found! Make sure to install
    SimpleITK version 2.0.2
    '''
    ID = ''.join([str(i) for i in np.random.choice(9, 10)])
    tmp_dir = os.path.dirname(save_path)+"/"+ID
    os.mkdir(tmp_dir)

    if skullstrip:
        skull_strip(ct_path, tmp_dir+"/CT_stripped.nii.gz")
        ct_path = tmp_dir+"/CT_stripped.nii.gz"

    if preprocess:
        preprocess_ct(ct_path, tmp_dir+"/CT_preproc.nii.gz")
        ct_path = tmp_dir+"/CT_preproc.nii.gz"

    if os.path.exists(ct_path):
        CT_nib = nib.load(ct_path)
        CT = read_ants(ct_path)
        lesion = read_ants(lesion_path)

        if lesion.shape[2] == sel_depth:
            start, finis = select_depth(CT.shape[2], sel_depth)
            H, W, D = lesion.shape
            padded_les = np.concatenate([np.zeros([H, W, start]),
                                         lesion.numpy(),
                                        np.zeros([H, W, finis])], axis=2)
        else:
            padded_les = lesion.numpy()

        # binary CT, fill hole to get brain mask
        output = fslmaths(CT_nib).bin().fillh().run()
        # remove the ICH region
        no_les = output.get_fdata()*(1 - padded_les)

        ori, direct, spacing = CT.origin, CT.direction, CT.spacing
        no_les = ants.from_numpy(no_les, origin=ori,
                                 direction=direct, spacing=spacing)
        padded_les = ants.from_numpy(padded_les, origin=ori,
                                     direction=direct, spacing=spacing)

        # tell ants to ignore the hemorrhage
        tx = ants.registration(MNI, CT, type_of_transform="SyN", mask=MNI_mask,
                               moving_mask=no_les,
                               outprefix=tmp_dir+"/ants")
        lesion_MNI = ants.apply_transforms(
            MNI, padded_les, tx["fwdtransforms"],
            interpolator="genericLabel")

        ants.image_write(lesion_MNI, save_path)
    shutil.rmtree(tmp_dir)


def pipeline(ct_dir, mask_dir, save_dir, atlas_dir,
             num=1, arrayID=0, **kwargs):
    MNI_path = atlas_dir+"/template/MNI.nii.gz"
    MNI_mask_path = atlas_dir+"/template/MNI_mask.nii.gz"

    MNI = ants.image_read(MNI_path)
    MNI_mask = ants.image_read(MNI_mask_path)

    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

    for index, i in enumerate(os.listdir(ct_dir)):
        if index % num == arrayID:
            ct_path = ct_dir+"/"+i
            lesion_path = mask_dir+"/"+i
            save_path = save_dir+"/"+i

            if os.path.exists(lesion_path) and os.path.exists(ct_path):
                if not os.path.exists(save_path):
                    register_lesion(ct_path, lesion_path, MNI, MNI_mask,
                                    save_path, **kwargs)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--ct_dir', type=str)
    parser.add_argument('--mask_dir', type=str)
    parser.add_argument('--save_dir', type=str)
    parser.add_argument('--atlas_dir', type=str, default="None")
    parser.add_argument('--preprocess', action="store_true")
    parser.add_argument('--skullstrip', action="store_true")

    parser.add_argument("--num", type=int, default=1)
    parser.add_argument("--arrayID", type=int, default=0)
    args = parser.parse_args()

    if args.atlas_dir == "None":
        atlas_dir = os.path.dirname(os.path.realpath(__file__)) + \
            "/../atlas/"
    else:
        atlas_dir = args.atlas_dir

    pipeline(args.ct_dir, args.mask_dir, args.save_dir,
             atlas_dir, args.num, args.arrayID,
             preprocess=args.preprocess,
             skullstrip=args.skullstrip)
