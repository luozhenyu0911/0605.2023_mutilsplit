import os, errno
import sys
# from extract import extract_basecalls_np
# from extract import calReader_old as calReader
import numpy as np
import gzip
# from multiprocessing import Pool
from datetime import datetime
from glob import glob
import re
import shutil
import subprocess
from collections import defaultdict
# import json
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--splitNfq", type=int)
args = parser.parse_args() 
N = args.splitNfq

# def combine_fastqs(paths,save_path):
    # # with gzip.open(save_path,'wb') as outfile:
    # #     for path in paths:
    # #         with gzip.open(path,'rb') as infile:
    # #            shutil.copyfileobj(infile,outfile)

    # ## use Linux cat command to concatenate fastqs fastest
    # command = ['cat'] + paths
    # with gzip.open(save_path,'w') as outfile:
        # processor = subprocess.run(command, stdin=None, stderr=None, stdout=outfile)

class combine_split(object):
    def __init__(self):
        self._dict = defaultdict(int)
        self.str_line0 = ''
        self.barcode_types = 0
        # self.split_barcode_num = 0
        self.reads_num = 0
        self.split_reads_num = 0
        self.bc1_cnt = 0
        self.bc2_cnt = 0
        self.bc3_cnt = 0


    def combine_split_log(self, n):
        path_log = 'data/split_read_{}.log'.format(str(n))
        with open(path_log, 'r') as f:
            tmp = [next(f) for x in range(4)]
            self.str_line0 = tmp[0].strip()
            self.barcode_types = int(tmp[0].strip().split('=')[2])
            self.reads_num+=int(tmp[2].strip().split('=')[1])
            self.split_reads_num+=int(tmp[3].strip().split('=')[1].split('(')[0])
            try:
                for line in f:
                    info = line.strip().split('\t')
                    bc = info[2]
                    bc_cnt = int(info[1])
                    self._dict[bc] +=bc_cnt
            except Exception as e: print(e)
            
    def combine_bcDiversity_log(self, n):
        path_log = 'data/split_read_{}.bc.diversity.count.log'.format(str(n))
        with open(path_log, 'r') as f:
            tmp = [next(f) for x in range(5)]
            # self.reads_num+=int(tmp[0])
            self.bc1_cnt += int(tmp[1].strip().split('=')[1])
            self.bc2_cnt += int(tmp[2].strip().split('=')[1])
            self.bc3_cnt += int(tmp[3].strip().split('=')[1])

if __name__ == "__main__":
    # N=2
    ins = combine_split()
    for n in ["{0:03}".format(i) for i in range(N+1)][1:]:
        ins.combine_split_log(n)
        ins.combine_bcDiversity_log(n)

    with open('split_stat_read1.log', 'w') as f:
        n=1
        split_barcode_num = len(ins._dict)
        r=100*split_barcode_num/ins.barcode_types
        f.write("{}\n".format(ins.str_line0))
        f.write("Real_Barcode_types = {split_barcode_num} ({r} %)\n".format(split_barcode_num=str(split_barcode_num), r=str(r)))
        f.write("Reads_pair_num  = {reads_num} \n".format(reads_num=str(ins.reads_num)))
        r2 = 100*ins.split_reads_num/ins.reads_num
        f.write("Reads_pair_num(after split) = {split_reads_num} ({r} %)\n".format(split_reads_num=str(ins.split_reads_num), r=str(r2)))
        # dict_sorted_by_val= {k: v for k, v in sorted(ins._dict.items(), key=lambda item: item[1], reverse=True)}
        
        for k,v in ins._dict.items():
            f.write('\t'.join([str(n), str(v), k])+'\n')
            n+=1
    with open('bc.individual.count.log', 'w') as bcf:
        print(ins.reads_num, file=bcf)
        print("bc1_cnt = "+ str(ins.bc1_cnt), file=bcf)
        print("bc2_cnt = "+ str(ins.bc2_cnt), file=bcf)
        print("bc3_cnt = "+ str(ins.bc3_cnt), file=bcf)
        print(ins.bc1_cnt, ins.bc2_cnt, ins.bc3_cnt, sep='\t', file=bcf)
        print(ins.bc1_cnt/ins.reads_num, ins.bc2_cnt/ins.reads_num, ins.bc3_cnt/ins.reads_num, sep='\t', file=bcf)
    # json.dump(ins._dict, open('split_read.json', 'w'))
