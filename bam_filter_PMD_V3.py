#!/usr/bin/env python

import re
import argparse
from tqdm import tqdm
import pysam
import io
import gzip

parser = argparse.ArgumentParser(prog='Parse bams to fix c-t before ngsLCA', description='')
parser.add_argument('--input', "-in",  help='input bam', required=True)
parser.add_argument('--min_distance', "-md",  help='maximum distance from reference ignoring C->T, default 0.95', type=float, default=0.95)
parser.add_argument('--max_distance', "-Md",  help='maximum distance from reference ignoring C->T, default 1', type=float, default=1.0)
parser.add_argument('--output', "-o",  help='Output file', required=False)
parser.add_argument('--distant_assignments', "-d",help='Detecting taxa with no exact matches', action="store_true")
parser.add_argument('--taxon_filer', "-tf",help='taxon filter', action="store_true")
parser.add_argument('--acc2taxid', "-t",help='acc2taxid')
parser.add_argument('--list_of_taxids_to_keep', "-l",help='list of taxids to keep')
args=parser.parse_args()

#ignore C->T and G->A

def read_taxid():
    #progress=tqdm(total=None, desc="Number of accessions parsed ")
    acc2taxid = io.TextIOWrapper(gzip.open(args.acc2taxid))
    acc2taxid_dict={}
    for item in acc2taxid:
        accession, accession_version, taxid, number = item.split("\t")
        taxid = taxid.strip()
        acc2taxid_dict[accession_version] = taxid
        #progress.update()
    #progress.close()

    return(acc2taxid_dict)

def parsing_bam_file():
    i=0
    if args.taxon_filer == True:
        taxids = open(args.list_of_taxids_to_keep)
        taxid_list=set()
        for item in taxids:
            taxid_list.add(item.strip())
        
    #progress = tqdm(total=None, desc="total reads parsed")
    #progress_2 = tqdm(total=None, desc="reads with :" + str(args.min_distance) + "-"+ str(args.max_distance) + "%")
    samfile = pysam.AlignmentFile(args.input, "rb")
    selected_reads = pysam.AlignmentFile(args.output, "wb", template=samfile)
    for read in samfile :
        if args.taxon_filer == True:
            taxid = acc2taxid_dict[read.reference_name]
            if taxid not in taxid_list:
                continue    

        if "I" in read.cigarstring or "D" in read.cigarstring:
            #progress.update()
            continue

        if read.flag == 0 or read.flag == 256 :
            if re.match('.*[C].*', read.get_tag("MD")) :
                MD = str(read.get_tag("MD"))
            else:
                distance=1-(read.get_tag("NM")/len(read.query))
                if (args.min_distance <= distance <= args.max_distance):
                    #print(distance)
                    selected_reads.write(read)
                #progress.update()
                continue
        elif read.flag == 16 or read.flag == 272 :
            if re.match('.*[G].*', read.get_tag("MD")) :
                MD = str(read.get_tag("MD"))    
            else:
                distance=1-(read.get_tag("NM")/len(read.query))
                if (args.min_distance <= distance <= args.max_distance):
                    selected_reads.write(read)
                #progress.update()
                continue

        read_bases = read.query_sequence[read.query_alignment_start:read.query_alignment_end]
        ref_bases = read.get_reference_sequence()
        count=0
        if read.flag == 0 or read.flag == 256 :
            for read_base, ref_base in zip(read_bases, ref_bases):
                if read_base != ref_base:
                    if ref_base.upper() == "C" and read_base.upper() == "T":
                            count+=1
        if read.flag == 16 or read.flag == 272 :
            for read_base, ref_base in zip(read_bases, ref_bases):
                if read_base != ref_base:
                    if ref_base.upper() == "G" and read_base.upper() == "A":
                        count+=1

        
        distance=1-((read.get_tag("NM")-count)/len(read.query))
        if (args.min_distance <= distance <= args.max_distance):
            selected_reads.write(read)
        
        #progress.update()
        #progress_2.update()
        #if i == 1000:
        #    quit()
        #i+=1

def parsing_bam_file_larger_distance():
    i=0
    #progress = tqdm(total=None, desc="total reads parsed")
    #progress_2 = tqdm(total=None, desc="reads with :" + str(args.min_distance) + "-"+ str(args.max_distance) + "%")
    samfile = pysam.AlignmentFile(args.input, "rb")
    selected_reads = pysam.AlignmentFile(args.output, "wb", template=samfile)
    
    last_reads = ""
    first_read = True
    for read in samfile :
        #print("\n")
        
        if "I" in read.cigarstring or "D" in read.cigarstring:
            #progress.update()
            continue


        if last_reads != read.query_name :
            first_read = True
            last_reads = read.query_name

        elif last_reads == read.query_name :
            first_read = False

        #print(first_read, read.query_name)
        

        #if read.query_name == "K00233:327:HVL3FBBXY:8:1101:1164:17016":
        #    quit()

        
        #print("\n", first_read,
         #     read.query_name)
        
        if read.flag == 0 or read.flag == 256 :
            if re.match('.*[T].*', read.get_tag("MD")) :
                MD = str(read.get_tag("MD"))
            else:
                distance=1-(read.get_tag("NM")/len(read.query))
                if first_read == True :
                    first_distance = distance
                #print(first_distance, distance)
                if (args.min_distance <= distance < args.max_distance) and first_distance<args.max_distance:
                    #progress_2.update()
                    selected_reads.write(read)
                    continue
                #progress.update()
                continue
        elif read.flag == 16 or read.flag == 272 :
            if re.match('.*[A].*', read.get_tag("MD")) :
                MD = str(read.get_tag("MD"))    
            else:
                distance=1-(read.get_tag("NM")/len(read.query))
                if first_read == True :
                    first_distance = distance
                    
                #print(first_distance, distance)
                if (args.min_distance <= distance < args.max_distance) and first_distance<args.max_distance:
                    selected_reads.write(read)
                    #progress_2.update()
                    continue
                #progress.update()
                continue

        read_bases = read.query_sequence[read.query_alignment_start:read.query_alignment_end]
        ref_bases = read.get_reference_sequence()
        count=0
        if read.flag == 0 or read.flag == 256 :
            for read_base, ref_base in zip(read_bases, ref_bases):
                if read_base != ref_base:
                    if ref_base.upper() == "C" and read_base.upper() == "T":
                            count+=1
        if read.flag == 16 or read.flag == 272 :
            for read_base, ref_base in zip(read_bases, ref_bases):
                if read_base != ref_base:
                    if ref_base.upper() == "G" and read_base.upper() == "A":
                        count+=1
        distance=1-((read.get_tag("NM")-count)/len(read.query))
        if first_read == True :
            first_distance = distance
    
        #print(first_distance, distance)

        if (args.min_distance <= distance < args.max_distance) and first_distance<args.max_distance:
            #progress_2.update()
            selected_reads.write(read)
            continue

        
        progress.update()
        #if i == 1000000:
        #    quit()
        #i+=1

if __name__ == '__main__':
    print("start main\n")
    if args.distant_assignments != True:
        if args.taxon_filer == True:
            acc2taxid_dict = read_taxid()
        parsing_bam_file()
    if args.distant_assignments == True:
        print("start looking for distant reads")
        parsing_bam_file_larger_distance()