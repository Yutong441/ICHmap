import os
import pandas as pd
from conn_matrix import get_conn_mat_all
import conn_metric_for as CMF


def normalize_metrics(W, k, method="fortran", *args, **kwargs):
    if method == "fortran":
        return CMF.normalize_metrics_for(W, k, *args, **kwargs)


def get_all_metrics(lesion_path, save_M_prefix, atlas_dir, *args, **kwargs):
    atlas = ["HCP-MMP", "visual", "somatomotor", "dorsal_attention",
             "ventral_attention", "limbic", "frontoparietal", "default"]

    all_mets, all_nodes = [], []
    for i in atlas:
        parcel_path = atlas_dir+"/"+i+".nii.gz"
        # first, obtain the paths passing through the lesion area
        M, nM, T = get_conn_mat_all(atlas_dir+"/fiber/", parcel_path,
                                    lesion_path=lesion_path,
                                    ref_path=atlas_dir+"/conn_mat/")

        # mat_dict = {"count": M, "ncount": nM}
        mat_dict = {"count": M}

        for key, val in mat_dict.items():
            # threshold the normalized count matrix
            if key == "ncount":
                val *= val > 1
            val.to_csv(save_M_prefix+"/"+os.path.basename(i)+"_"+key+".csv")
            norm_met, node_met = normalize_metrics(val, 15,
                                                   *args, **kwargs)
            norm_met["atlas"] = os.path.basename(i)
            norm_met["measure"] = key
            all_mets.append(norm_met)

            node_met["nodes"] = node_met.index
            node_met = pd.melt(node_met, id_vars="nodes")
            node_met["atlas"] = os.path.basename(i)
            node_met["measure"] = key
            all_nodes.append(node_met)

        if i == "HCP-MMP":
            T.to_csv(save_M_prefix+"tract_num.csv")

    all_mets = pd.concat(all_mets, axis=0)
    all_nodes = pd.concat(all_nodes, axis=0)
    return all_mets, all_nodes
