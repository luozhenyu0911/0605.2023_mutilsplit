configfile: "config.yaml"
# from fastq/*.fq.gz, merge lane first

def run_all_input(wildcards):
    # These will always be returned and include the summary report, GC tables and graphs
    # Coverage, flagstat, and other metrics
    run_all_files = ["data/split_read.2.fq.gz", "data/split_read.1.fq.gz", "data/read_2_fastqc.html", "split_stat_read1.log"]
    part_id=["{0:03}".format(i) for i in range(config['params']['splitNfq']+1)][1:]

    # for id in part_id:
    #     run_all_files.append("data/split_read_{id}.2.fq.gz".format(id=id))
    
    return run_all_files

rule run_all:
    input:
        run_all_input
        

src_dir = config['params']['src_dir']
#include: "/home/ycai/hm/test/CGI_WGS_nonstandard_Pipeline/test/pool/split_py/full_data/splitreads.smk"
include: "./split_only.smk"
#include: "./split_only.smk"

# /home/ycai/hm/dev/pipeline/splitreads_only_from_fq.smk
# /home/ycai/hm/dev/pipeline/splitreads_only_from_combined.reads.smk
