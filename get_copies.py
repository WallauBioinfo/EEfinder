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

Test files:
https://github.com/dezordi/PEVEI/test_files/get_copies/

Dependencies:
CD-HIT version 4.7+
"""
###############################>LIBRARIES<###############################

import pandas as pd
import numpy as np
import argparse, csv, os, subprocess, shlex, sys, time, re

###############################>ARGUMENTS<###############################
parser = argparse.ArgumentParser(description = 'This scripts receives a fasta file and execute a cd-hit-est analysis and formats the cd-hit output',formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-in", "--input", help="Fasta file, e.g. example.fasta", required=True)
parser.add_argument("-p","--threads",help="Threads for cd-hit-est analyze.", default = 1, type=str)
parser.add_argument("-m","--memory",help="Memory in MB for cd-hit-est analyze",default = 2000, type=str)
#Storing argument on variables
args = parser.parse_args()
input_file = args.input
var_threads = args.threads
var_memory = args.memory
###############################>EXECUTION<###############################
#Creating lists to store genome assemblies and viral taxonomy
with open(input_file+'.cd-hit.log','w') as cd_hit_log:
    cd_hit_out = input_file+'.cd'
    cd_hit_est_cmd = 'cd-hit-est -i '+input_file+' -o '+cd_hit_out+' -c 0.8 -aL 0.8 -g 1 -n 8 -M '+var_memory+' -T '+var_threads+' -d 200'
    cd_hit_est_cmd = shlex.split(cd_hit_est_cmd)
    cd_hit_est_cmd_process = subprocess.Popen(cd_hit_est_cmd,stdout = cd_hit_log)
    cd_hit_est_cmd_process.wait()
with open(input_file+'.cd.clstr','r') as cluster_file, open(input_file+'.cd.clstr.tsv','w') as cluster_formated:
    cluster_file_reader = cluster_file.readlines()
    cluster_formated_writer = csv.writer(cluster_formated,delimiter='\t')
    cluster_formated_writer.writerow(['Sequence','Assembly-Code','Specie-Code','Cluster','Representative'])
    cluster_list = []
    for line in cluster_file_reader:
        line = line.rstrip('\n')
        if 'Cluster' in line:
            cluster_number = re.sub(r'>Cluster ','',line)
        else:
            if 'at ' in line:
                representative = 'FALSE'
            else:
                representative = 'TRUE'
            sequence_name = re.sub(r'.*>','',line)
            sequence_name = re.sub(r'\.\.\..*','',sequence_name)
            prefix_assembly = re.sub(r'\/.*','',sequence_name)
            prefix_spp = re.sub(r'(_[a-zA-Z0-9]*$){1}?','',prefix_assembly)
            cluster_list.append([sequence_name,prefix_assembly,prefix_spp,cluster_number,representative])
    cluster_formated_writer.writerows(cluster_list)
df = pd.read_csv(input_file+'.cd.clstr.tsv',sep='\t')
count_clusters = df.pivot_table(index=['Cluster'], aggfunc='size')
df_count_cluster = count_clusters.to_frame()
df_count_cluster = df_count_cluster.reset_index()
df_count_cluster.columns = ['Cluster','Copies']
df_inter_assembly = df.groupby('Cluster')['Assembly-Code'].nunique()
df_inter_assembly = df_inter_assembly.to_frame()
df_inter_assembly = df_inter_assembly.reset_index()
df_inter_assembly.columns = ['Cluster','Assemblies']
df_inter_specie = df.groupby('Cluster')['Specie-Code'].nunique()
df_inter_specie = df_inter_specie.to_frame()
df_inter_specie = df_inter_specie.reset_index()
df_inter_specie.columns = ['Cluster','Species']
df2 = pd.merge(df, df_count_cluster, on=['Cluster'])
df2 = pd.merge(df2, df_inter_assembly, on=['Cluster'])
df2 = pd.merge(df2, df_inter_specie, on=['Cluster'])
df2 = df2.rename(columns={"Species": "inter-Species", "Assemblies": "inter-Assemblies"})
df2['inter-Assemblies'] = np.where(df2['inter-Assemblies']==1, 'FALSE', 'TRUE')
df2['inter-Species'] = np.where(df2['inter-Species']==1, 'FALSE', 'TRUE')
df2.drop(['Assembly-Code', 'Specie-Code'], axis='columns', inplace=True)
df2.to_csv(input_file+'.cd.clstr.tsv',index=False,sep='\t')
