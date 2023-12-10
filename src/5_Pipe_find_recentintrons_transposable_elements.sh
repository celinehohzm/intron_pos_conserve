conda create -n repeatmasker_env
conda activate repeatmasker_env
conda config --add channels bioconda
conda config --add channels conda-forge
conda install repeatmasker
RepeatMasker -species human -pa 4 recentintrons_multidnaseq.fa
less recentintrons_multidnaseq.fa.out
python recentintrons_filter_repeatmasker.py
less recentintrons_filtered_repeats.csv
