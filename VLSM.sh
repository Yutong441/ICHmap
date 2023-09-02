#! /bin/bash
cd /opt
outcome=${1:-"Disability"}
mkdir data/tmp -p
mkdir data/VLSM -p
mkdir data/VLSM/$outcome -p

export FREESURFER_HOME=/opt/.synthstrip/
export PATH=$PATH:$FREESURFER_HOME/bin
export PATH=$PATH:/opt/.fsl
export PATH=$PATH:/opt/.julia/bin

Rscript ICHmap/format_df.R data/labels/$outcome".csv" \
    data/tmp/$outcome".csv"

num=$(cat data/tmp/$outcome".csv" | wc -l)
num=$((num - 1))
echo "n="$num > data/VLSM/$outcome/"sample.txt"

JULIA_DEPOT_PATH=/opt/.julia_lib julia ICHmap/VLSMsurv.jl --root=data/tmp \
    --outcome=$outcome \
    --save_dir="data/VLSM/"$outcome

rm data/tmp -r
