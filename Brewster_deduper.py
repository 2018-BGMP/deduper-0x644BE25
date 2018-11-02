from subprocess import run
import argparse
import re

r = re.compile("[0-9]+[A-Z]")
shiftR = {'N','M','=','X','S','D'}

def getArgs():
    """ use argparse to return arguments from invocation """
    parser=argparse.ArgumentParser(description="Remove PCR duplicates from SAM file")
    parser.add_argument("-f", "--file", help="absolute path to SAM file", required=True, type=str)
    parser.add_argument("-p", '--paired', help="(optional) indicates paired-end reads", action="store_true")
    parser.add_argument("-u", "--umi", help="path to UMI file", type=str)
    
    return(parser.parse_args())
    
def sortSAM(unsorted_file_name):
    """ take name of SAM file, returns the name of a sorted SAM file """
    temp = unsorted_file_name.split("/")
    short=temp[len(temp)-1][:-4]
    sorted_file_name = short+"_sorted.sam"
    run("/usr/bin/time samtools sort -o "+sorted_file_name+" "+unsorted_file_name, shell=True)
    return(sorted_file_name)

def buildDict(umi_list):
    """ take a list of UMIs, return a dict of (UMI: empty set) pairs """
    """ TESTING: with UMI file """
    umi_dict={}
    for umi in umi_list:
        umi_dict[umi.strip('\n')]=set([0])
    
    return(umi_dict)

def getStrand(flag):
    """ take a SAM bitwise flag, and return int for validity/strand 
       (0 = not mapped, 1 or -1 = strand)
       TESTING: 165 -> 0, 20 -> 0, 769 -> 1, 784 -> -1"""
    if 4 & flag > 0: # checks for 'unmapped'
        return(0)
    if 16 & flag == 0: # checks for 'reverse strand'
        return(1) # return forward
    return(-1) # return reverse

def getForPos(bits,pos):
    """ correct for leftmost soft clipping in forward reads """
    for b in bits:
        c = b[-1:]
        if c == 'M':
            return(pos)
        if c == 'S':
            pos = pos-int(b[:-1])
    return(pos)

def getRevPos(bits,pos):
    """ find true 5' position for reverse reads """
    mSeen = False
    for b in bits:
        c = b[-1:]
        if c == 'M': 
            mSeen = True
        if c in shiftR and mSeen:
            pos += int(b[:-1])
    return(-1*pos)

def parseRead(raw_read):
    """ parses , returns chromosome, strand, corrected start pos, UMI """
    """ TESTING: use reads from unit test input file """
    read = raw_read.strip('\n').split('\t')
    return(read)
    
def main():
    """ DO IT ALL """
    # get command line arguments
    args = getArgs()
    
    # if unsupported flag given, quit with error.
    if args.paired:
        print("Error: paired-end reads not supported")
        return()
    if args.umi == None:
        print("Error: randomers not supported. Please provide UMI file")
        return()
    
    # sort our SAM file
    short_name = sortSAM(args.file)
    sorted_file = open(sortSAM(args.file))
    umi_file = open(args.umi)
    umi_list = umi_file.readlines()
    umi_file.close()
    
    results = open(short_name[:-10]+'deduped.sam', 'w+') 
    curr_chrom = ''
    
    # iterate over lines
    for line in sorted_file:
        # write all header lines
        if line.startswith('@'):
            results.write(line)
        else:
            read = parseRead(line)
            umi = read[0][-8:]
            
            # switch to new chromosome, if necessary
            if read[2] != curr_chrom:
                umi_dict = buildDict(umi_list)
                curr_chrom = read[2]
            
            # process line, write if not PCR dup
            if umi in umi_dict:
                strand = getStrand(int(read[1]))
                pos=0
                bits = re.findall(r, read[5])
                if strand == 1: pos = getForPos(bits, int(read[3]))
                if strand == -1: pos = getRevPos(bits,int(read[3]))
                if pos not in umi_dict[umi]:
                    umi_dict[umi].add(pos)
                    results.write(line)
    # close files
    sorted_file.close()
    results.close()
    
if __name__ == "__main__": 
    main()