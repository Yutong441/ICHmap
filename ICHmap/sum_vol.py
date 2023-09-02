import re
import os
import pandas as pd
import nibabel as nib


def sum_vol(img_path):
    img_nib = nib.load(img_path)
    img = img_nib.get_fdata()
    vol = (img > 0).sum()

    pixdim = img_nib.header["pixdim"]
    voxdim = pixdim[1]*pixdim[2]*pixdim[3]
    return voxdim*vol


def sum_all(img_dir, save_path):
    out_dict = {}
    for i in sorted(os.listdir(img_dir)):
        one_vol = sum_vol(img_dir+"/"+i)
        out_dict[re.sub(".nii.gz", "", i)] = one_vol/1000

    out_dict = pd.DataFrame.from_dict(out_dict, orient="index")
    out_dict.to_csv(save_path)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--img_dir', type=str)
    parser.add_argument('--save_path', type=str)
    args = parser.parse_args()
    sum_all(args.img_dir, args.save_path)
