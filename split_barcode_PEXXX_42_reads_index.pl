#!/usr/bin/env perl
use strict;

if(@ARGV != 11)
{
        print "Example: perl split_barcode_PEXXX_42_reads.pl barcode.list read_1.fq.gz read_2.fq.gz 100 split_read bc_start gdna_start additional_bc_start additional_bc_len gdna_start_r1 read_len_r1 \n";
        exit(0);
}

my $read_len = $ARGV[3]; # Get read length
my $bc_start_idx = $ARGV[5]-1; # on seq start from 1, index start from 0
my $gdna_start_idx = $ARGV[6]-1;
my $additional_bc_start_idx = $ARGV[7]-1;
my $additional_bc_len = $ARGV[8];
my $gdna_start_idx_r1 = $ARGV[9]-1;
my $read_len_r1 = $ARGV[10];

my ($n1, $n2, $n3, $n4, $n5) = (10, 6, 10, 6, 10); # get barcode and intermediate lengths
my %barcode_hash; # hash map of barcodes
my ($bc1_cnt, $bc2_cnt, $bc3_cnt) = (0, 0, 0);
open IN,"$ARGV[0]" or die "can't open barcode.list"; # open barcode list
my $n = 0; # initialize counter
while(<IN>){ # for line in file
  $n ++; # counter plus one
  my @line = split; # split line into barcode and ID
  my @barcode = split(//,$line[0]); # split barcode into array of chars
  my $barcode_ID = $line[1]; # save ID num
  for(my $num = 0; $num <= 9; $num++){ # add a hash key for each barcode with 1 base mismatch
    my @barcode_mis = @barcode; # duplicate original barcode for barcode mismatch
    $barcode_mis[$num] = "A"; # sub first base with A
    my $barcode_mis = join("",@barcode_mis); # join into string
    $barcode_hash{$barcode_mis} = $barcode_ID; # add barcode as key and the ID as value
    @barcode_mis = @barcode; # repeat for other nucleotides for all bases
    $barcode_mis[$num] = "G";
    my $barcode_mis = join("",@barcode_mis);
    $barcode_hash{$barcode_mis} = $barcode_ID;
    @barcode_mis = @barcode;
    $barcode_mis[$num] = "C";
    my $barcode_mis = join("",@barcode_mis);
    $barcode_hash{$barcode_mis} = $barcode_ID;
    @barcode_mis = @barcode;
    $barcode_mis[$num] = "T";
    my $barcode_mis = join("",@barcode_mis);
    $barcode_hash{$barcode_mis} = $barcode_ID;
  }
}
close IN; # close file
my $barcode_types = $n * $n * $n; # num potential combinations
my $barcode_each = $n; # num barcodes

open IN1,"gzip -dc $ARGV[1] |" or die "cannot open read1"; # open read1
open IN2,"gzip -dc $ARGV[2] |" or die "cannot open read2"; # open read2
open OUT1, "| gzip > $ARGV[4].1.fq.gz" or die "Can't write file";
open OUT2, "| gzip > $ARGV[4].2.fq.gz" or die "Can't write file";
$n = 0; # initialize couner
my $reads_num; # number of reads
my $progress; # n million reads 
my %index_hash; # key is c_barcode, val is index
my %index_hash_reverse; # key is index, val is c_barcode
my $split_barcode_num; # num of represented c_barcodes
my $T; 
my $id;
my $reads_num;
my @line;
my @Read_num;
$Read_num[0] = 0; 
my $split_reads_num;
my $additional_bc;
my $seq_r1;
my $qual_r1;

while(<IN2>){ # read lines from Read2
  chomp; # get line
  @line = split; # split line (this seems unnecessary)
  $n ++; # counter += 1
  if($n % 4 == 1){ # if counter modulo 4 = 1
    $reads_num ++; # add 1 to the number of reads processed
    my @A  = split(/\//,$line[0]); # split line by forward slash, save as A
         $id = $A[0]; # id is the readname
         if($reads_num % 1000000 == 1) # if 1mil reads 
         {
              print "reads processed $progress (M) reads ...\n"; # print message
              $progress ++; # add 1 mil reads for next time
         }

  } # now we have the readname as $id
  if($n % 4 == 2){ # if we have the second line (sequence)
    my $read = substr($line[0], $gdna_start_idx, $read_len);
    my $b1 = substr($line[0], $bc_start_idx, $n1); # get barcode seqs from R2
    my $b2 = substr($line[0], $bc_start_idx+$n1+$n2, $n3);
    my $b3 = substr($line[0], $bc_start_idx+$n1+$n2+$n3+$n4, $n5);

    if($additional_bc_len>0){
      $additional_bc = substr($line[0], $additional_bc_start_idx, $additional_bc_len);
      }
    if(exists $barcode_hash{$b1}) {
      $bc1_cnt ++;
    }
    if(exists $barcode_hash{$b2}) {
      $bc2_cnt ++;
    }
    if(exists $barcode_hash{$b3}) {
      $bc3_cnt ++;
    }

    if((exists $barcode_hash{$b1}) && (exists $barcode_hash{$b2}) && (exists $barcode_hash{$b3})){ # check the barcodes are in the hash map
      my $hash = $barcode_hash{$b1}."_".$barcode_hash{$b2}."_".$barcode_hash{$b3}; # key is the three barcode ids
      if(!(exists $index_hash{$hash})){ # if the combined id isn't in the combined hash map
        $split_barcode_num ++; # add another represented barcode
        $index_hash{$hash} = $split_barcode_num; # add the combined id to the hash map
        $index_hash_reverse{$split_barcode_num} = $hash; # add the reverse (key is n combined barcodes, value is c_bc)
        $Read_num[$index_hash{$hash}] = 0; # this is an array of num reads for each represented barcode
                                           # because the barcodes represented start at 0 they can act as the index
      }
      $split_reads_num ++; # add a split read
      $Read_num[$index_hash{$hash}] ++; # add a another read for the barcode
      if($read_len_r1>0){
      $T = <IN1>; chomp($T); # remove new line
      print OUT1 "$id\#$hash\#\/1\tBX:Z:$hash\n";
      $T = <IN1>; chomp($T);
      $seq_r1 = substr($T, $gdna_start_idx_r1, $read_len_r1);
      print OUT1 "$seq_r1\n"; # JUST PRINT
      $T = <IN1>; chomp($T);
      print OUT1 "$T\n"; # JUST PRINT
      $T = <IN1>; chomp($T);
      $qual_r1 = substr($T, $gdna_start_idx_r1, $read_len_r1);
      print OUT1 "$qual_r1\n"; # JUST PRINT
        }
      # print OUT2 "$id\#$hash\#\/2\tBX:Z:$hash\n";
      if($additional_bc_len>0){
        print OUT2 "$id\#$hash\_$additional_bc\#\/2\tBX:Z:$hash\n";
      } else {
        print OUT2 "$id\#$hash\#\/2\tBX:Z:$hash\n";
      }
      print OUT2 "$read\n";
      $T = <IN2>; $n++;chomp($T);
      print OUT2 "$T\n";
      $T = <IN2>; $n++;chomp($T);
      my $qual = substr($T,$gdna_start_idx,$read_len);
      print OUT2 "$qual\n";
    }
    else{
      $Read_num[0] ++;
      if($read_len_r1>0) {
      $T = <IN1>; chomp($T);
      print OUT1 "$id\#0_0_0\#\/1\tBX:Z:0_0_0\n";
      $T = <IN1>; chomp($T);
      $seq_r1 = substr($T, $gdna_start_idx_r1, $read_len_r1);
      print OUT1 "$seq_r1\n";
      $T = <IN1>; chomp($T);
      print OUT1 "$T\n";
      $T = <IN1>; chomp($T);
      $qual_r1 = substr($T, $gdna_start_idx_r1, $read_len_r1);
      print OUT1 "$qual_r1\n";
      }
      # print OUT2 "$id\#0_0_0\#\/2\tBX:Z:0_0_0\n";
      if($additional_bc_len>0){
        print OUT2 "$id\#0_0_0\_$additional_bc\#\/2\tBX:Z:0_0_0\n";
      } else {
        print OUT2 "$id\#0_0_0\#\/2\tBX:Z:0_0_0\n";
      }
      print OUT2 "$read\n";
      $T = <IN2>; $n++;chomp($T);
      print OUT2 "$T\n";
      $T = <IN2>; $n++;chomp($T);
      my $qual = substr($T,$gdna_start_idx,$read_len);
      print OUT2 "$qual\n";
    }

  }
}
close IN1;
close IN2;
close OUT1;
close OUT2;


open OUT3, ">split_stat_read1.log" or die "Can't write file";
print OUT3 "Barcode_types = $barcode_each * $barcode_each * $barcode_each = $barcode_types\n";
my $r;
$r = 100 *  $split_barcode_num/$barcode_types;
print OUT3 "Real_Barcode_types = $split_barcode_num ($r %)\n";
$r = 100 *  $split_reads_num/$reads_num;
print OUT3 "Reads_pair_num  = $reads_num \n";
print OUT3 "Reads_pair_num(after split) = $split_reads_num ($r %)\n";
for(my $i=1;$i<=$split_barcode_num;$i++){
  print OUT3 "$i\t$Read_num[$i]\t$index_hash_reverse{$i}\n";
}
close OUT3;

### check split rate of each individual bc
open OUT4, ">bc.individual.count.log" or die "Can't write file";
if($read_len_r1==0) {
  print OUT4 "$reads_num\n";
} else{
  my $reads_num_pe = $reads_num*2;
  print OUT4 "$reads_num_pe\n";
}
print OUT4 "bc1_cnt=$bc1_cnt \n";
print OUT4 "bc2_cnt=$bc2_cnt \n";
print OUT4 "bc3_cnt=$bc3_cnt \n";
print OUT4 "$bc1_cnt\t$bc2_cnt\t$bc3_cnt\n";
my $bc1_cnt_pct = $bc1_cnt/$reads_num;
my $bc2_cnt_pct = $bc2_cnt/$reads_num;
my $bc3_cnt_pct = $bc3_cnt/$reads_num;
print OUT4 "$bc1_cnt_pct\t$bc2_cnt_pct\t$bc3_cnt_pct\n";
close OUT4;

print "all done!\n";

