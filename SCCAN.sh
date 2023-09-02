# prepare dataset
cd /opt
mkdir -p data/SCCAN
mkdir -p data/tmp
Rscript ICHmap/split_data.R

# run SCCAN
for i in $(ls data/tmp/*_img.txt)
do
    outcome=${i%"_img.txt"}
    outcome=${outcome#"data/tmp/"}
    Rscript ICHmap/lesy.R data/tmp data/SCCAN $outcome
    num=$(cat $i | wc -l)
    echo "n="$num > data/SCCAN/$outcome/"sample_size.txt"
done

rm data/tmp -r
