import re
import os
import ants
import pandas as pd
import conn_metric as ME
from preprocess import read_ants
from register import register_lesion


def lesion_metric(final_dir, atlas_dir):
    all_me, all_nodes = [], []
    if not os.path.exists(final_dir+"/metric"):
        os.mkdir(final_dir+"/metric")

    if os.path.exists(final_dir+"/ICH_reg.nii.gz"):
        metrics, nodes = ME.get_all_metrics(
            final_dir+"/ICH_reg.nii.gz", final_dir+"/metric/", atlas_dir)
        all_me.append(metrics)
        all_nodes.append(nodes)

    all_me = pd.concat(all_me, axis=0)
    all_nodes = pd.concat(all_nodes, axis=0)
    all_me.to_csv(final_dir+"/graph_metrics.csv")
    all_nodes.to_csv(final_dir+"/node_metrics.csv")


def lesion_disconn(ct_dir, mask_dir, save_dir, atlas_dir,
                   num=1, arrayID=0, **kwargs):
    '''
    Lesion connectome

    Args:
        `ct_dir`: directory containing CT data
        `mask_dir`: directory containing ICH mask
        `save_dir`: save registered ICH mask + connectivity matrix + metrics
        `atlas_dir`: path to the directory containing the HCP-MMP atlas
        `num`: number of jobs to submit in parallel
        `arrayID`: job ID. The `num` and `arrayID` arguments are only suitable
        in HPC setting.

    Returns:
        the folder structure in `save_dir`:
        patient1
        ---ICH_seg.nii.gz
        ---graph_metrics.csv
        ---node_metrics.csv
        ---metric/

        patient2
        ---ICH_seg.nii.gz
        ---graph_metrics.csv
        ---node_metrics.csv
        ---metric/
    '''
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

    MNI = read_ants(atlas_dir+"/template/MNI_header.nii.gz")
    MNI_mask = ants.threshold_image(MNI, low_thresh=0.5, binary=True)

    for index, i in enumerate(os.listdir(ct_dir)):
        if index % num == arrayID:
            ID = re.sub(".nii.gz", "", i)
            if not os.path.exists(save_dir+"/"+ID):
                os.mkdir(save_dir+"/"+ID)

            lesion_path = mask_dir+"/"+i
            ct_path = ct_dir+"/"+i
            save_path = save_dir+"/"+ID+"/ICH_reg.nii.gz"

            if not os.path.exists(save_path):
                register_lesion(ct_path, lesion_path, MNI, MNI_mask,
                                save_path, **kwargs)

            if os.path.exists(save_path):
                if not os.path.exists(save_dir+"/"+ID+"/metric/tract_num.csv"):
                    lesion_metric(save_dir+"/"+ID, atlas_dir)


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

    lesion_disconn(
        args.ct_dir, args.mask_dir, args.save_dir,
        atlas_dir, args.num, args.arrayID)
