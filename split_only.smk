# import pathlib
from pathlib import Path
import re
import glob
# We've gotten some inconsistent fastq names in the past
# This code scrapes the directories for a file matching "_1.fq.gz"
# This becomes the input for aggregating the reads
def get_fastqs_one(wildcards):
    # get fq base path from config
    fq_path = Path(config['samples']['fq_path'])
    # get lanes to add to the path
    lanes = config['samples']['lanes']
    fq_files = []
    for lane in lanes:
        # set aggregate path for each lane
        loc_path = fq_path / lane
        # iterate through the files in the filepath
        path = "{}/{}_read_1*.fq.gz".format(loc_path, lane)
        for fp in glob.glob(path):
            if "discard" not in str(fp):
                print(fp)
                fq_files.append(str(fp))
    return fq_files            


def get_fastqs_two(wildcards):
    # get fq base path from config
    fq_path = Path(config['samples']['fq_path'])
    # get lanes to add to the path
    lanes = config['samples']['lanes']
    fq_files= []
    for lane in lanes:
        # set aggregate path for each lane
        loc_path = fq_path / lane
        # iterate through the files in the filepath
        
        if config['params']['sequence_type']=='pe':
            path = "{}/{}_read_2*.fq.gz".format(loc_path, lane)
            for fp in glob.glob(path):
                if "discard" not in str(fp):
                    print(fp)
                    fq_files.append(str(fp))
        elif config['params']['sequence_type']=='se':
            path = "{}/{}_read*.fq.gz".format(loc_path, lane)
            for fp in glob.glob(path):
                if "discard" not in str(fp):
                    print(fp)
                    fq_files.append(str(fp))
    return fq_files            


# aggregate all fq1 files in data/
rule cat_read_one:
    input:
        get_fastqs_one
    output:
        "data/read_1.fq.gz"
    run:
        if config['params']['sequence_type']=='pe':
            # if their are multiple files concatenate them
            if len(config['samples']['lanes']) > 1:
                shell("cat {input} > {output}")
            # otherwise just link the one path
            else:
                shell("ln -s ../{input} {output}")
        elif config['params']['sequence_type']=='se':
            shell("touch data/read_1.fq.gz")


# aggregate all fq2 files in data/
rule cat_read_two:
    input:
        get_fastqs_two
    output:
        "data/read_2.fq.gz"
    run:
        if config['params']['sequence_type']=='se':
            shell("ln -s ../{input} {output}")
        elif config['params']['sequence_type']=='pe':
            # if their are multiple files concatenate them
            if len(config['samples']['lanes']) > 1:
                shell("cat {input} > {output}")
            # otherwise just link the one path
            else:
                shell("ln -s ../{input} {output}")
        

# barcodes can also get messy
# sometimes they're from the barcode list and sometimes they're from the RC list
# This function samples the barcodes and attempts to determine the appropriate list
def determine_barcode_list(n_samples):
    import gzip
    import sys
    gc_swap = config['params']['gc_swap']
    if gc_swap:
        # GC dye swapped
        barcode_file = "/home/ycai/hm/dev/pipeline/barcode.GCswap.list"
        barcode_rc_file = "/home/ycai/hm/dev/pipeline/barcode_RC.GCswap.list"
        
    else:
        # set potential barcode files based on the config file, also set the fastq path
        barcode_file = config['params']['toolsdir'] + config['params']['barcode']
        barcode_rc_file = config['params']['toolsdir'] + config['params']['barcode_RC']

    fastq_path = "data/read_2.fq.gz"

    # Read in barcodes fro, the barcodes file
    # We allow for one mismatch for every barcode
    # so we add those as well
    def get_barcodes(barcodes_file):
        print(f"reading in {barcodes_file}", file=sys.stderr)
        nucs = ['A', 'C', 'G', 'T']
        barcodes = []
        with open(barcodes_file, "r") as bcs_file:
            for line in bcs_file:
                bc = line.strip().split()[0]
                # iterate through each nucleotide in the barcode
                for i in range(0, len(bc)):
                    # iterate through potential mismatch nucleotides
                    for nuc in nucs:
                        bc_alt = list(bc)
                        bc_alt[i] = nuc
                        # add mismatched barcode to the barcodes list
                        barcodes.append("".join(bc_alt))
            
            # return the list as a set, their shouldn't be any duplicates
            # this is mostly to make the search faster
            barcodes_set = set(barcodes)
            return barcodes_set

    
    # scrape the fastq file for barcodes
    def get_fq_barcodes(fastq_path, n_samples):
        
        with gzip.open(fastq_path, "rb") as fastq:
            fq_bcs = []
            counter = 1
            print("Getting fastq barcodes", file=sys.stderr)
            for line in fastq:
                # keep track of number of barcodes we've looked at
                counter += 1
                # break after we've sampled enough barcodes
                if counter == n_samples + 1:
                    break
                if counter % 400000 == 0:
                    print(f"Read in {counter/4} lines", file=sys.stderr)
                # append barcodes for each read
                if counter % 4 == 3:
                    seq = line.decode('utf-8').strip()
                    bc_start_idx = config['params']['bc_start']-1
                    _bc_len = 10
                    _bc_gap = 6 
                    fq_bcs.append(seq[bc_start_idx: bc_start_idx+_bc_len])
                    fq_bcs.append(seq[bc_start_idx+_bc_len+_bc_gap: bc_start_idx+_bc_gap+_bc_len*2])
                    fq_bcs.append(seq[bc_start_idx+_bc_gap*2+_bc_len*2: bc_start_idx+_bc_gap*2+_bc_len*3])

            return fq_bcs

    # get barcodes, rc_barcodes, and fq_barcodes
    barcodes = get_barcodes(barcode_file)
    barcodes_rc = get_barcodes(barcode_rc_file)
    fq_bcs = get_fq_barcodes(fastq_path, n_samples)

    bcs_found = 0
    rc_bcs_found = 0
    # check fq barcodes against the barcodes and rc_barcodes
    for bc in fq_bcs:
        if bc in barcodes:
            bcs_found += 1
        if bc in barcodes_rc:
            rc_bcs_found += 1

    print(f"Barcodes: {bcs_found}\nRC_Barcodes: {rc_bcs_found}\n",file=sys.stderr)
    # if more barcodes were found than RC barcodes return the barcodes path
    if bcs_found > rc_bcs_found:
        print("Barcode_list: {}".format(config['params']['barcode']), file=sys.stderr)
        return config['params']['barcode']
    # if more RC barcodes found retunr RC barcodes path
    elif rc_bcs_found > bcs_found:
        print("Barcode_list: {}".format(config['params']['barcode_RC']), file=sys.stderr)
        return config['params']['barcode_RC']
    # if the same number of barcodes were found try again with more samples (unlikely)
    # if we sample more than 10 mil reads, exit - something is wrong
    else:
        if n_samples > 40000000:
            print("Can't determine correct barcode file for sample", file=sys.stderr)
            sys.exit(1)
        else:
            return determine_barcode_list(n_samples * 2)
            

# We can also get barcodes in two formats
# BGI started omitting the 6bp spacers in between barcodes
# We want to check ti see which format we're getting and run the correct command accordingly
def get_read_diff(read1, read2):
    import sys
    import gzip
    with gzip.open(read1, "r") as r1, gzip.open(read2, "r") as r2:
        # r1.readline()
        r2.readline()
        # seq1 = r1.readline().strip()
        seq2 = r2.readline().strip()
        # check to see that our supplied read length is the same as the determined barcode
        # exit if they're not
        # if len(seq1) != config['params']['read_len']:
        #     print(f"interpreted read length ({len(seq1)}) does not match supplied read length ({config['params']['read_len']})", file=sys.stderr)
        #     sys.exit(1)
        if len(seq2) == 63:
            print(f"21+42BP barcode detected", file=sys.stderr)
            return True
        elif len(seq2) == 142:
            print(f"100+42BP barcode detected", file=sys.stderr)
            return True
        elif len(seq2) == 92:
            print(f"50+42BP barcode detected", file=sys.stderr)
            return True
        elif len(seq2) == 200:
            print(f"158+42BP barcode detected", file=sys.stderr)
            return True
        else:
            print(f"Unknown barcode length detected", file=sys.stderr)
            sys.exit(1)

part_ids=["{0:03}".format(i) for i in range(config['params']['splitNfq']+1)][1:]
rule splitNfq:
    input:
        r1="data/read_1.fq.gz",
        r2="data/read_2.fq.gz"
    output:
        expand("out/read_2.part_{part_id}.fq.gz", part_id=part_ids)
    params:
        splitN = config['params']['splitNfq']
    run:
        if config['params']['sequence_type']=="se":
            shell("~/anaconda3/bin/seqkit split2 {input.r2} -p {params.splitN} -O ./out -f ")
        else:
            shell("~/anaconda3/bin/seqkit split2 -1 {input.r1} -2 {input.r2} -p {params.splitN} -O ./out -f ")

# rule for splitting reads
rule split_reads:
    input:
        "out/read_2.part_{part_id}.fq.gz"
    output:
        "data/split_read_{part_id}.2.fq.gz"
        # "split_stat_read1.log"
    params:
        r2_len = config['params']['read_len'],
        r1_len = config['params']['read_len_r1'], 
        toolsdir = config['params']['toolsdir'],
        src_dir = config['params']['src_dir'],
        randomBC_dir = config['params']['randomBC_dir'],
        general_python = config['params']['general_python'],
        bc_len = config['params']['sbc_len'],
        cbc_len = config['params']['cbc_len'],
        bc_len_redundant = cbc_len
        bc_start = config['params']['bc_start'],
        gdna_start = config['params']['gdna_start'],
        additional_bc_start = config['params']['additional_bc_start'],
        additional_bc_len = config['params']['additional_bc_len'],
        gdna_start_r1 = config['params']['gdna_start_r1'], 
        additional_bc_len_r1 = config['params']['additional_bc_len_r1'],
        adapter_len = config['params']['adapter_len'], 
    run:
        if config['modules']['BCsplit']==True:
            if config['params']['bc_condition'] == 'random_bc':
                params.barcode = '/home/ycai/hm/dev/pipeline/CGI_WGS_nonstandard_Pipeline/barcode.list'
                shell("perl /home/ycai/hm/test/new_bc_random/test/split_py/full_data/split_randomBC_index_dualBC.pool.pl "
                "{params.toolsdir}{params.barcode} "
                "out/read_1.part_{wildcards.part_id}.fq.gz out/read_2.part_{wildcards.part_id}.fq.gz {params.r2_len} data/split_read_{wildcards.part_id} {params.cbc_len} {params.bc_len_redundant} {params.bc_start} {params.gdna_start} {params.additional_bc_start} {params.additional_bc_len} {params.gdna_start_r1} {params.r1_len} {params.additional_bc_len_r1} "
                "2> data/split_read_{wildcards.part_id}.err ")
                # shell("{params.general_python} {params.randomBC_dir}/src/split_randomBC.py "
                # "--read1_fq_gz {input.r1} --read2_fq_gz {input.r2} --read_len {params.r2_len} --output_suffix data/split_read --cbc_len {params.cbc_len} --bc_len_redundant {params.bc_len_redundant} --bc_start {params.bc_start} --gdna_start {params.gdna_start} --additional_bc_start {params.additional_bc_start} --additional_bc_len {params.additional_bc_len} "
                # "2> data/split_stat_read.err")
            elif config['params']['bc_condition'] == 'standard':
                # determine which barcode list to use
                if config['params']['enforced_bc_list']=='barcode':
                    params.barcode = '/barcode.list'
                elif config['params']['enforced_bc_list']=='barcode_RC':
                    params.barcode = '/barcode_RC.list'
                else:
                    params.barcode = determine_barcode_list(400000)
                    # check the difference in length then use the appropriate command
                    
                # if get_read_diff(input[0], input[1]):
                if params.bc_len==42:
                    shell("perl {params.src_dir}/src/split_barcode_PEXXX_42_reads_index.pl "
                        "{params.toolsdir}{params.barcode} "
                        "{input} {params.r2_len} data/split_read {params.bc_start} {params.gdna_start} {params.additional_bc_start} {params.additional_bc_len} {params.gdna_start_r1} {params.r1_len} "
                        "2> data/split_stat_read.err")
                    # shell("{params.general_python} {params.src_dir}/split_barcode_PEXXX_42_reads.py "
                    # "--barcode_list {params.toolsdir}{params.barcode} "
                    # "--read1_fq_gz {input.r1} --read2_fq_gz {input.r2} --read_len {params.r2_len} --output_suffix data/split_read --bc_start {params.bc_start} --gdna_start {params.gdna_start} "
                    # "2> data/split_stat_read.err")

                elif params.bc_len==30:
                    # this uses the 30 bp split script
                    shell("perl {params.toolsdir}/tools/split_barcode_PEXXX_30_reads.pl "
                        "{params.toolsdir}{params.barcode} "
                        "{input} {params.len} data/split_read "
                        "2> data/split_stat_read.err")
            elif config['params']['bc_condition'] == 'BCgDNA':
                params.barcode = determine_barcode_list(400000)
                shell("perl {params.src_dir}/src/split_barcode_PEXXX_42_reads_BCgDNA.pl "
                        "{params.toolsdir}{params.barcode} "
                        "{input} {params.r2_len} data/split_read {params.bc_start} {params.gdna_start} {params.additional_bc_start} {params.additional_bc_len} {params.gdna_start_r1} {params.r1_len} {params.adapter_len} "
                        "2> data/split_stat_read.err")
            else:
                print("unknown type")
        else:
            shell("ln -s ./read_1.fq.gz data/split_read.1.fq.gz && ln -s ./read_2.fq.gz data/split_read.2.fq.gz")

rule cat_fq2:
    input:
        expand("data/split_read_{part_id}.2.fq.gz", part_id=part_id)
    output:
        "data/split_read.2.fq.gz"
    run:
        shell("cat data/split_read_*.2.fq.gz > data/split_read.2.fq.gz && rm data/split_read_*.2.fq.gz && rm out/read_2.part_*.fq.gz ")

rule cat_fq1:
    input:
        expand("data/split_read_{part_id}.2.fq.gz", part_id=part_id)
    output:
        "data/split_read.1.fq.gz"
    run:
        if config['params']['sequence_type']=="se":
            shell("touch data/split_read.1.fq.gz && rm data/split_read_*.1.fq.gz && rm out/read_1.part_*.fq.gz ")
        else:
            shell("cat data/split_read_*.1.fq.gz > data/split_read.1.fq.gz && rm data/split_read_*.1.fq.gz && rm out/read_1.part_*.fq.gz ")

### merge split_stat_read1.log in each thread into one file
rule cat_log:
    input:
        "data/split_read.2.fq.gz"
    output:
        "split_stat_read1.log"
    params: 
        splitN = config['params']['splitNfq']
    run:
        shell("python /home/ycai/hm/test/CGI_WGS_nonstandard_Pipeline/test/pool/split_py/full_data/combine_split.py --splitNfq {params.splitN} && rm data/*log && rm data/*err ")

rule fastqc:
    input:
        "data/read_2.fq.gz"
    output:
        "data/read_2_fastqc.html"
    run:
        if config['params']['sequence_type']=="se": 
            shell("/home/ycai/hm/tools/FastQC/fastqc {input} -o ./data/")
        else:
            shell("/home/ycai/hm/tools/FastQC/fastqc data/read_1.fq.gz data/read_2.fq.gz -o ./data/")
