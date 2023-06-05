echo "start=$(date)"

source /home/eanderson/Virtual_Envs/SnakeMake/bin/activate
export PATH=/home/eanderson/bin:$PATH
export PATH=/usr/bin:$PATH

#smk=/home/ycai/hm/dev/pipeline/dev/CGI_WGS_nonstandard_Pipeline/Snakefiles/run.all.smk
smk=./split_only_runall.smk
snakemake -p -j 32 -k -s ${smk} 2>> snakemake.err.txt

end=$(date)
echo "end=$(date)"
