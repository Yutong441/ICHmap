# Spatial correlates of ICH symptoms
## Installation
* step 0: obtain ICH mask
    - Please run the deepbleed algorithm **without removing fragments**
    - i.e., Don't put the `--remove_frag` argument
    - I found removing fragments tend to remove the segmentation in the IVH
    region.
    - Without removing the fragments, it may mis-segment the choroid plexus and
    skull calcification, neither of which are in the way of any major tracts.
    - Thus, there are more benefits keeping the fragments

* step 1: set up python virtual environment
```bash
python3.9 -m venv venv
source venv/bin/activate
pip install -r requirement.txt
deactivate
```
 
* step 2: compile code
    - skip this step if you are using python 3.9 in Ubuntu 22.04
    - install gcc, gfortran and liblapack (usually pre-installed in Linux)
    
```bash
source venv/bin/activate
cd ICHcon/cdipy
python setup.py build_ext --inplace
cd -

cd ICHcon/bct_for
python -m numpy.f2py -llapack -c -m fort *.f90
cd -
deactivate
```
 
* step 3: download atlas (for lesion disconnectome) 
    - download atlas [here](https://drive.google.com/drive/folders/10I3bosYYpxn0fwjvneHlj4AY3YcdXr3R?usp=sharing)
    - put the files into the `atlas/` folder
    - NB: I need to approve your access before you can download the atlas
 
* step 4: set up R environment (for SCCAN)
    - Open an R interactive session 
    - `install.packages("data.table")`
    - `devtools::install_github('dorianps/LESYMAP')`
    - If it does not work, follow the instruction [here](https://github.com/dorianps/LESYMAP)
 
* step 5: set up Julia environment (for Cox-VLSM)
    - Install Julia if it is not done already, follow [here](https://julialang.org/downloads/)
    - Type `julia` in the terminal to open an interactive shell, type `]` to
    enter package environment
    - `add Survival StatsModels CSV DataFrames GLM NIfTI ArgParse MultipleTesting HypothesisTests`

## Lesion connectome
* prepare data
Make sure all the **original** CT scans are in one folder (`ct_dir`), and
all the ICH mask are in another folder (`mask_dir`). 
An example folder is:

```
ct_dir
---MGH001.nii.gz
---MGH002.nii.gz
mask_dir
---MGH001.nii.gz
---MGH002.nii.gz
```

* run analysis
```bash
python ICHcon/lesion.py \
    --ct_dir=$ct_dir \
    --mask_dir=$mask_dir \
    --save_dir=$save_dir \
    --preprocess --skullstrip
```

Processing one patient takes roughly 15 minutes.

## SCCAN
It is very hard to automate this analysis because I did not write the SCCAN
function. Here are the steps:

* step 1: register images to MNI template. NB: VLSM uses a **different**
template from lesion disconnectome

```bash
python ICHmap/reg.py \
    --ct_dir=$ct_dir \
    --mask_dir=$mask_dir \
    --save_dir=$reg_mask_dir \
    --preprocess --skullstrip
```

* step 2: prepare two text files
    - First, make a new directory `$SCCAN`, which should contain two
    subdirectories: `label` and `results`
    - Within the `$SCCAN/label` directory, save 2 files:
    - The first file called `Disability_img.txt`, contains the **absolute path**
    of all registered ICH masks
    - The second file is called `Disability_behavior.txt`, which contains the
    mRS of patients at 6 month or 12 month. NB: each row in the
    `Disability_behavior.txt` file must **correspond** to each row in the
    `Disability_img.txt` file
    - NB: **No missing** values are allowed
 
* step 3: run analysis
```bash
Rscript ICHmap/lesy.R $SCCAN Disability
```

This analysis takes 2~6 hours depending on the sample size.

## Cox-VLSM
* step 1: register images to MNI template if not already done
 
* step 2: Prepare a dataframe called `Disability.csv` with the following attributes:
    - `Disability_time`: time to develop disability (mRS $>2$) or censor
    - `Disability_bool`: whether disability was developed (mRS $>2$)
    - `Age`
    - `Sex`
    - `ICH_vol`
    - `paths`: image paths
    - NB: I have prepared this dataframe for the Yale dataset, you just need to
    add the `ICH_vol` and `paths` columns.

* step 3: run the analysis
```bash
julia ICHmap/VLSMsurv.jl \
    --root=Disability.csv \
    --outcome=Disability \
    --save_dir=$VLSM_output \
    --covar=Age,Sex,ICH_vol
```

This analysis takes 1~5 minutes depending on the sample size.
