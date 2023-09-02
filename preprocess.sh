#! /bin/bash
# perform skullstripping, ICH segmentation, and registration
# input: directory contains one directory named `ct`
# output:

# --ct
# --skullstripped
# --mask
# --connectome
# --reg_lesion
# --stats

num_jobs=${1:-1}
arrayID=${2:-0}
num_cpu=${3:-3}

cd /opt/
mkdir -p data/skullstripped
mkdir -p data/mask
mkdir -p data/connectome
mkdir -p data/reg_lesion
mkdir -p data/stats

source ~/.bashrc
source /opt/.miniconda3/etc/profile.d/conda.sh
export FREESURFER_HOME=/opt/.synthstrip/
export PATH=$PATH:$FREESURFER_HOME/bin
export PATH=$PATH:/opt/.fsl
export PATH=$PATH:/opt/.julia/bin

# step 1: skullstripping
conda activate ich
python ICHcon/preprocess.py --ct_dir=data/ct \
    --save_dir=data/skullstripped \
    --num=$num_jobs \
    --arrayID=$arrayID
conda deactivate

# step 2: ICH segmentation
conda activate seg
# bash deepbleed/run.sh data/skullstripped data/mask $num_cpu
python deepbleed/predict2.py --indir data/skullstripped \
    --outdir data/mask \
    --weights deepbleed/weights/weights \
    --cpus $num_cpu \
    --verbose \
    --brain \
    --remove_frag \
    --num=$num_jobs \
    --arrayID=$arrayID
conda deactivate

# step 3: registration
conda activate ich
python ICHcon/register.py --ct_dir=data/skullstripped \
    --mask_dir=data/mask \
    --save_dir=data/reg_lesion \
    --num=$num_jobs \
    --arrayID=$arrayID

# step 4: lesion disconnectome
python ICHcon/lesion.py --ct_dir=data/skullstripped \
    --mask_dir=data/mask \
    --save_dir=data/connectome \
    --num=$num_jobs \
    --arrayID=$arrayID

# step 5: summary statistics
if [ $arrayID == 0 ]
then
    Rscript ICHmap/summary.R data/connectome data/stats
    python ICHmap/sum_vol.py --img_dir=data/mask \
        --save_path=data/stats/vol.csv
fi
