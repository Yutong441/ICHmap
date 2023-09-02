import os
import shutil
import numpy as np
import nibabel as nib
import ants
from fsl.wrappers import fslmaths
from preprocess import read_ants
from preprocess import sel_central
from preprocess import select_depth
from preprocess import bash_in_python


def register_lesion(ct_path, lesion_path, MNI, MNI_mask, save_path):
    '''
    Method:
    1. skullstripping (optional)
    2. downsampling to 5mm, select central slices
    3. lesion-aware registration

    Args:
        `ct_path`: CT image path
        `lesion_path`: ICH mask path
        `MNI`: MNI brain template path
        `MNI_mask`: MNI brain mask path
        `save_path`: where to save the registered lesion
        `sel_depth`: only select the central slices
    '''
    ID = ''.join([str(i) for i in np.random.choice(9, 10)])
    tmp_dir = os.path.dirname(save_path)+"/"+ID
    os.mkdir(tmp_dir)

    if os.path.exists(ct_path):
        CT_nib = nib.load(ct_path)
        CT = read_ants(ct_path)
        lesion = read_ants(lesion_path)

        if lesion.shape[2] != CT.shape[2]:
            sel_depth = lesion.shape[2]
            lesion = sel_central(lesion.numpy(), sel_depth)
            start, finis = select_depth(CT.shape[2], sel_depth)
            H, W, D = lesion.shape
            padded_les = np.concatenate([np.zeros([H, W, start]),
                                         lesion,
                                        np.zeros([H, W, finis])], axis=2)
        else:
            padded_les = lesion.numpy()

        # binary CT, fill hole to get brain mask
        if "FSLDIR" in os.environ and os.path.exists(os.environ["FSLDIR"]):
            print("Using fslmaths")
            output = fslmaths(CT_nib).bin().fillh().run()
        else:
            bash_in_python(
                "fslmaths {} -bin -fillh {}/output.nii.gz".format(
                    ct_path, tmp_dir)
            )
            output = nib.load(tmp_dir+"/output.nii.gz")

        no_les = output.get_fdata()*(1 - padded_les)

        ori, direct, spacing = CT.origin, CT.direction, CT.spacing
        no_les = ants.from_numpy(no_les, origin=ori,
                                 direction=direct, spacing=spacing)
        padded_les = ants.from_numpy(padded_les, origin=ori,
                                     direction=direct, spacing=spacing)

        # tell ants to ignore the hemorrhage
        os.environ["ANTS_RANDOM_SEED"] = "1"
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

    parser.add_argument("--num", type=int, default=1)
    parser.add_argument("--arrayID", type=int, default=0)
    args = parser.parse_args()

    if args.atlas_dir == "None":
        atlas_dir = os.path.dirname(os.path.realpath(__file__)) + \
            "/../atlas/"
    else:
        atlas_dir = args.atlas_dir

    pipeline(args.ct_dir, args.mask_dir, args.save_dir,
             atlas_dir, args.num, args.arrayID)
