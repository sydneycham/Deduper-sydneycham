#!/usr/bin/env python

import argparse
import re 

def get_args():
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='This code is used to identify duplicates in SAM file reads. On single end reads, and is assuming the UMIs have been added to the read.THis code takes into consideration all possile cigar strings (inclduing soft clipping), strandedness, single-end reads, and known UMI')
    parser.add_argument("-f", "--filename", help="Input filename", required=True)
    parser.add_argument("-o", "--outfile", help="output file name", type=str)
    parser.add_argument("-o2", "--outfile2", help="output file name for duplicates", type=str)
    parser.add_argument("-u", "--umi", help="umi file name", type=str)
    return parser.parse_args()

args=get_args()
f=args.filename
o=args.outfile
o2=args.outfile
u=args.umi

umis = [] #empty umis list intialize

def umis_adder(u):
    with open (u, "r") as fh:
        while True:
            line = fh.readline().split()
            if line == []:
                break
            umis.append(line)
        return(u)
umis_adder(u)

def strand_dir(bit_flag):
    '''this function will determine the direction of the strand, returns True(Forward) or False(reverse) '''
    return (bit_flag & 0x10) !=0


def soft_clip(line):
    '''This function calculates and returns the adjusted position, taking into account the original position, the number of matches (M), right soft clipping (S), deletions (D), and skips (N) in the CIGAR string.'''
    line = line.split()
    cigar = line[5]
    pos = int(line[3])
    bit_flag = int(line[1])
    reverse_strand = strand_dir(bit_flag) #determines direction of strand
    cigar_fields = re.findall(r'(\d+)([MDNS])', cigar)
        #if there is soft clipping and the strand is reverse, referred to as "reverse_strand", we will get an integer to add later to position
    adjusted_pos = 0
    if reverse_strand:
        if cigar_fields[0][1] == "S":
            cigar_fields = cigar_fields[1::]
        for i in range(len(cigar_fields)):
            adjusted_pos += int(cigar_fields[i][0])
        return adjusted_pos + pos
        #if there is soft clipping and the strand is FORWARD, referred to as "not reverse_strand", we will get an integer to add later to position
    else:
        if cigar_fields[0][1] == "S":
            adjusted_pos = int(cigar_fields[0][0])
            return pos - adjusted_pos
        else:
            return pos


def main(f):
    prev_chrome = None
    dupcount = 0
    with open (f, "r") as inSam, open(o, "w") as outfile, open(o2, "w") as dupfile, open(u, "w") as unknown_umis:
        main_set = set()
        counter=0
        while True:
            line = inSam.readline().strip()
            if line.startswith("@"):
                outfile.write(f'{line}\n')
                continue
            if(line == ""):
                break
            parts = line.split('\t')
            startpos = soft_clip(line)
            chrom = parts[2]
            umi = parts[0].split(":")[-1]
            #if umi not in umis:
                #unknown_umis.write(str(line))
                #continue
            bitflag = int(parts[1])
            stranded = strand_dir(bitflag)
            #create a set of tuples with startpos,strand,umi,chrom. search in set for matching tuples, if tuple present send to output duplicate file, if not present write to main file and save in set
            if (umi, stranded, chrom, startpos) not in main_set:
                main_set.add(tuple((umi, stranded, chrom, startpos)))
                outfile.write(f'{line}\n')
            else:
                dupcount+=1
        # if parts[2] != prev_chrom:     
        #     membership_set.clear()
        #     prev_chrom = chromosome
    print(dupcount)
        
main(f)

  
# #plus strand example, should adjust the position from 101 to 100
#line="ReadID    163    ReferenceName    101    30    1S1M5S    *    0    0    ATCGG...    ~~~~    NM:i:0"

# #minus strand exmaple, should adjust the position from 101 to 107
#line = "ReadID    16    ReferenceName    101    30    1S1M5N5D5S    *    0    0    ATCGG...    ~~~~    NM:i:0"  

# cigar = line.split()[5]
# adjusted_pos = soft_clip(line, cigar)
# print(f"The original positon is: {int(line.split()[3])}")
# print(f"Adjusted Position: {adjusted_pos}")            