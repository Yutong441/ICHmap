# Spatial correlates of ICH symptoms
## Installation
Download the singualrity container from [here](https://drive.google.com/drive/folders/16IiDczcBSkzWHiUPlcKM6NROUzeT7D62?usp=sharing).
Place the container file to wherever you want.

Please install singularity according to this [page](https://docs.sylabs.io/guides/3.0/user-guide/installation.html).
In HPC, singularity is usually pre-installed.
However, you need to load the relevant module such as `module load singularity`.

## Input data
Create a result folder called `$result`.
In this folder, create a subfolder called `ct`.
Put all the **original** CT scans are in that folder. 

Put the outcome data in the `label` directory.
For VLSM, the outcome data must contain a column called "ID" which is the same
as the names of the CT images.
It contains two columns: one denoting the time to event, the other the event
status. 
For example, if the outcome data is called `Disability.csv`, then the two
columns must be named "Disability_time" and "Disability_bool".
It also needs to contain the "Age" and "Sex" column.

For SCCAN, the outcome data should be prepared differently.
It must contain a column called "ID".
It can contain any number of columns, which will be analysed individually.
The file must be called `SCCAN.csv`.

Overall, the directory should look like:

```
ct
---MGH001.nii.gz
---MGH002.nii.gz
label
---Disability.csv
---SCCAN.csv
```

## Analysis
### Preprocessing
Run the following command
```bash
singularity exec --bind /opt/data:$result ichmap.sif /bin/bash /opt/preprocess.sh 1 0 3
```
Unless running multiple parallel jobs on HPC, the first two arguments should be
`1` and `0`.
The last number `3` refers to the number of cpu allocated.
Please see this [script](scripts/container_preprocess.sh) for an example to run
on HPC to take advantage of array job.

The output in the `$result` folder would look like, for instance:

```
skullstrip
---MGH001.nii.gz
---MGH002.nii.gz
mask
---MGH001.nii.gz
---MGH002.nii.gz
reg_lesion
---MGH001.nii.gz
---MGH002.nii.gz

connectome
---MGH001/
---------ICH_reg.nii.gz
---------graph_metrics.csv
---------node_metrics.csv
---------metric/
---MGH002/

stats
---all_graphs.csv
---all_tracts.csv
---all_vols.csv
```

### VLSM
Running VLSM is very simple:
```bash
singularity exec --bind $result:/opt/data ichmap.sif /bin/bash /opt/VLSM.sh Disability
```

The output would be located in the `$result/VLSM/Disability` folder.

### SCCAN
```bash
singularity exec --bind $result:/opt/data ichmap.sif /bin/bash /opt/SCCAN.sh
```
The output would be located in the `$result/SCCAN` folder.
