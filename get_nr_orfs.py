#!/usr/bin/python3
# -*- coding: utf-8 -*-
###############################>GENERAL-INFORMATIONS<###############################
"""
Build in Python 3.6

Author:
Filipe Dezordi
zimmer.filipe@gmail.com
https://dezordi.github.io/

Script repository:
https://github.com/dezordi/PEVEI

Dependencies:
EMBOSS:6.6.0.0+
CD-HIT version 4.7+
"""
###############################>LIBRARIES<###############################

import pandas as pd
import numpy as np
import argparse, csv, os, subprocess, shlex, sys, time, re, glob
from Bio import SeqIO

###############################>ARGUMENTS<###############################
parser = argparse.ArgumentParser(description = 'This script receives a fasta file and execute a getorf analysis and remove redundant ORFs from getorf output based on clusterization analysis',formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-in", "--input", help="Fasta file", required=True)
parser.add_argument("-ms", "--minsize", help="minsize getorf argument 'Minimum nucleotide size of ORF to report.'", default='100',type=str)
parser.add_argument("-p","--threads",help="Threads for cd-hit-est analyze.", default = 1, type=str)
parser.add_argument("-m","--memory",help="Memory in MB for cd-hit-est analyze",default = 2000, type=str)
#Storing argument on variables
args = parser.parse_args()
input_file = args.input
min_size = args.minsize
var_threads = args.threads
var_memory = args.memory
###############################>FUNCTIONS<###############################
def format_name(line):
    """
    This function remove some special characters that getorf recognize as pattern of end fasta header.
    """
    if '>' in line:
        try:
           line = re.sub(r'\/','--',line)
           line = re.sub(r':','__',line)
        except:
            pass
    return line

def reformat_name(line):
    """
    This function back the sequence modified names to original names.
    """
    if '>' in line:
        try:
           line = re.sub(r'--','/',line)
           line = re.sub(r'__',':',line)
        except:
            pass
    return line

def format_getorf(line):
    """
    This function formats getorf fasta names output.
    """
    if '>' in line:
        try:
           line = re.sub(r' \[','/ORF:',line)
           line = re.sub(r' - ','-',line)
           line = re.sub(r']','',line)
           line = re.sub(r' \(REVERSE SENSE\)','/REVERSE',line)
        except:
            pass
    return line
###############################>EXECUTION<###############################
#Formating names to remove special characters

fasta_file =open(input_file,'r')
fasta_reader = fasta_file.readlines()
fasta_file.close()
with open(input_file,'w') as fasta_out:
    for line in fasta_reader:
        line = format_name(line)
        fasta_out.write(line)
fasta_out.close()
#Running getorf
with open(input_file+'.getorf.log','w') as getorf_log:
    getorf_out = input_file+'.getorf'
    getorf_est_cmd = 'getorf -seq '+input_file+' -outseq '+getorf_out+' -table 1 -minsize '+min_size+' -find 1'
    getorf_est_cmd = shlex.split(getorf_est_cmd)
    getorf_est_cmd_process = subprocess.Popen(getorf_est_cmd,stdout = getorf_log)
    getorf_est_cmd_process.wait()
#Reformating names from input file and getorf output
fasta_file =open(input_file,'r')
fasta_reader = fasta_file.readlines()
fasta_file.close()
with open(input_file,'w') as fasta_out:
    for line in fasta_reader:
        line = reformat_name(line)
        fasta_out.write(line)
fasta_out.close()
getorf_file =open(getorf_out,'r')
getorf_reader = getorf_file.readlines()
getorf_file.close()
with open(getorf_out,'w') as getorf_out_renamed:
    for line in getorf_reader:
        line = reformat_name(line)
        line = format_getorf(line)
        getorf_out_renamed.write(line)
getorf_out_renamed.close()
#Excluding redundant
getorf_csv = open(getorf_out+'.csv','w')
getorf_csv_writer = csv.writer(getorf_csv,delimiter=',')
getorf_csv_writer.writerow(['ORF-ID','Element-ID','ORF','ORF-Start','ORF-End','ORF-Length','Sense'])
getorf_info_list = []
with open(getorf_out,'r') as getorf_file:
    getorf_reader = getorf_file.readlines()
    for line in getorf_reader:
        if '>' in line:
            ORF_id = re.sub(r'>','',line).rstrip('\n').strip()
            Element_id = re.sub(r'_[0-9]*\/ORF.*','',ORF_id).rstrip('\n')
            ORF = re.sub(r'\/ORF.*','',ORF_id).rstrip('\n')
            ORF = re.sub(r'.*_','',ORF)
            ORF_start = re.sub(r'.*ORF:','',ORF_id).rstrip('\n')
            ORF_start = re.sub(r'-.*','',ORF_start)
            ORF_end = re.sub(r'.*ORF:','',ORF_id).rstrip('\n')
            ORF_end = re.sub(r'.*-','',ORF_end)
            ORF_end = re.sub(r'\/REVERSE','',ORF_end)
            if 'REVERSE' in ORF_id:
                ORF_length = int(ORF_start) - int(ORF_end)
                sense = 'neg'
                getorf_info_list.append([ORF_id,Element_id,ORF,ORF_end,ORF_start,ORF_length,sense])
            else:
                ORF_length = int(ORF_end) - int(ORF_start)
                sense = 'pos'
                getorf_info_list.append([ORF_id,Element_id,ORF,ORF_start,ORF_end,ORF_length,sense])
                
getorf_csv_writer.writerows(getorf_info_list)
getorf_csv.close()
Element_id_list = []
with open(getorf_out+'.csv','r') as getorf_csv:
    getorf_csv_reader = csv.reader(getorf_csv,delimiter=',')
    for line in getorf_csv_reader:
        if line[1] not in Element_id_list:
            Element_id_list.append(line[1])
for element in Element_id_list:
    orf_ids = []
    with open(getorf_out+'.csv','r') as getorf_csv:
        getorf_csv_reader = csv.reader(getorf_csv,delimiter=',')
        for line in getorf_csv_reader:
            if element == line[1]:
                orf_ids.append(line[0])
    all_ORF_sequences = SeqIO.parse(open(getorf_out),'fasta')
    element_output_orfs = str(element+'.ORFs')
    try:
        element_output_orfs = re.sub(r'.*\/','',element_output_orfs)
    except:
        pass
    with open(element_output_orfs, "w") as element_output:
        for seq in all_ORF_sequences:
            if seq.id in orf_ids:
                SeqIO.write([seq], element_output, "fasta")
    with open(element_output_orfs+'.cd-hit.log','w') as cd_hit_log:
        cd_hit_out = element_output_orfs+'.orf_cd'
        cd_hit_est_cmd = 'cd-hit -i '+element_output_orfs+' -o '+cd_hit_out+' -c 0.9 -aL 0.1 -M '+str(var_memory)+' -T '+str(var_threads)+' -d 200'
        cd_hit_est_cmd = shlex.split(cd_hit_est_cmd)
        cd_hit_est_cmd_process = subprocess.Popen(cd_hit_est_cmd,stdout = cd_hit_log)
        cd_hit_est_cmd_process.wait()

nr_orfs = "cat *orf_cd > non_redudant_orfs.fasta"
subprocess.call(nr_orfs,shell=True)
clstr = "cat *orf_cd.clstr > cd_hit_clstrs.txt"
subprocess.call(clstr,shell=True)
cd_hit_log = "cat *.cd-hit.log > cd_hit_logs.txt"
subprocess.call(cd_hit_log,shell=True)

fileList = glob.glob('./*.ORFs*', recursive=True)
for filePath in fileList:
    try:
        os.remove(filePath)
    except OSError:
        print("Error while deleting file")

fileList = glob.glob('./*.orf_cd.cd*', recursive=True)
for filePath in fileList:
    try:
        os.remove(filePath)
    except OSError:
        print("Error while deleting file")

os.rename('non_redudant_orfs.fasta',input_file+'.getorf.nr')