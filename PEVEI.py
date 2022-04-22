#!/usr/bin/python3
# -*- coding: utf-8 -*-
###############################>GENERAL-INFORMATIONS<###############################
"""
Build in Python 3.6

Author:
Filipe Dezordi
zimmer.filipe@gmail.com
https://github.com/dezordi

Dependencies:
BLAST 2.9.0+
bedtools v2.26.0
"""

###############################>LIBRARIES<###############################

import pandas as pd
import numpy as np
import argparse, csv, os, subprocess, shlex, sys, time, re
from Bio import SeqIO
from Bio import Entrez
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast.Applications import NcbitblastxCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline

###############################>ARGUMENTS<###############################
parser = argparse.ArgumentParser(description = 'This script receives a fasta file, and predict regions of Endogenous Viral Elements.',usage='\n>Running with default parameters:\n    PEVEI.py -in genome.fasta -db viral_prot.fasta -db1 host_genes.fasta -db2 host_tranposons.fasta\n>Using a prefix:\n    PEVEI.py -in genome.fasta -db viral_prot.fasta -db1 host_genes.fasta -db2 host_tranposons.fasta -pr teste_1\n>Retriving only EVEs in contigs greater than 50,000bp, with 15,000bp of each flanking region, considering as same EVE element in a range of 500bp and removing EVEs with more than 10 percent of lowercase letters:\n    PEVEI.py -in genome.fasta -db viral_prot.fasta -db1 host_genes.fasta -db2 host_tranposons.fasta -ln 50000 -fl 15000 -lm 500 -mp 10',formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-in", "--input", help="Fasta file (nucleotides).", required=True)
parser.add_argument("-db", "--database", help="Viral proteins database file.", required=True)
parser.add_argument("-db1",'--databasefilter', help="Host and heatshock proteins database file.")
parser.add_argument("-db2","--databaseTE",help="Transposon database file.")
parser.add_argument("-bt","--blasttype",help="Choose the blast strategy against tranposons database, default = blastn.",default='blastn',choices=['blastn','blastx','tblastx'])
parser.add_argument("-ln", "--length", help="Minimum length of contigs used for BLAST, default = 10,000.",type=int, default=10000)
parser.add_argument("-fl","--flank", help="Length of flanking regions of EVEs to be extracted, default = 10,000.", type=int,default=10000)
parser.add_argument("-lm", "--limit", help="Limit of bases used to merge regions on bedtools merge, default = 1.",type=str, default=str(1))
parser.add_argument("-mp","--mask_per",help="Limit of lowercase letters in percentage to consider a putative EVE as a repetitive region, default = 50.",type=str,default=50)
parser.add_argument("-p", "--threads", help="Threads for multi-thread analysis, default = 1.",type=int, default=1)
parser.add_argument("-rm","--removetmp",help="Remove temporary files generated through analysis? default = True.", default = False, choices=['True','False'])
parser.add_argument("-mk","--makeblastdb",help="Make blast database?, default = True.", default = True, choices=['True','False'])
parser.add_argument("-st","--step",help="Select step of analysis, default = All.", default = 'All', choices=['All','EVEfinder','Tax','TEs','Flank','Flank-blast'])
parser.add_argument("-pr","--prefix",help="Write the prefix name for output files. We strongly recommend the use of genome name assembly, this prefix will be used to create the EVEname (The EVE name will be formated as PREFIX|CONTIG/SCAFFOLD:START-END) default = input file name.")
#Storing argument on variables
args = parser.parse_args()
input_file = args.input
database_file = args.database
database_filter = args.databasefilter
database_TE = args.databaseTE
blast_type = args.blasttype
length_cutoff = args.length
flank_region = args.flank
limit_merge = args.limit
mask_per = args.mask_per
threads = args.threads
temp_remove = args.removetmp
make_db = args.makeblastdb
step = args.step
prefix = args.prefix
if prefix == None: #format prefix name, if the -pr argument was not used
    prefix = input_file
    prefix = re.sub("\..*","",prefix)
    prefix = re.sub('.*/','',prefix).rstrip('\n')

###############################>BIO.ENTREZ<###############################
"""
If you prefer put your ncbi information direct on code, put # in lines 74 and 75, remove # from lines 76 and 77
and put your ncbi email and ncbi api-key inside the ''.
"""
#Entrez.email = str(input('Digite your e-mail to logging on NCBI: ')) #add this line after...
#Entrez.api_key = str(input('Digite your api key from NCBI: ')) #add this line after...
Entrez.email = 'zimmer.filipe@gmail.com'
Entrez.api_key = '29a8e7fd2add3e853ce3d3c18f140399a809'

###############################>FUNCTIONS<###############################
def makeblastdb(data,in_db_type):
    """
    This function creates blast databases.

    Keyword arguments:
    data - the database_file, parsed with -db argument
    in_db_type - 'nucl' or 'prot' strings
    """

    clinedb = NcbimakeblastdbCommandline(dbtype=in_db_type, input_file = data)
    stdout, stderr = clinedb()
    return(print(f'Makeblastdb for {data}: DONE'))

def runblastx(query_file,database_file):
    """
    This function runs blastx with genome and viruses.

    Keyword arguments:
    query_file - input_file - parsed with -in argument
    database - pre-formated blast database - parsed with -db argument
    
    Blast arguments:
    max_intron_length = "Length of the largest intron allowed in a translated nucleotide sequence when linking multiple distinct alignments (a negative value disables linking)."
    soft_masking = "Apply filtering locations as soft masks (i.e., only for finding initial matches)."
    """

    cline = NcbiblastxCommandline(query = query_file, db = database_file, out = query_file+'.blastx', 
    outfmt = 6, word_size = 3, evalue = 0.00001, num_threads = threads, matrix = 'BLOSUM45', max_intron_length = 100, soft_masking = 'true')
    stdout, stderr = cline()
    return(print(f'BLASTx with {query_file} against {database_file}: DONE'))

def blastx_filter(blast_result,tag):
    """
    This function receives a blastx result and filter based on query ID and ranges of qstart and qend.

    Keyword arguments:
    blast_result = blastx result outfmt 6
    """

    header_outfmt6 = ['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore'] #creates a blast header output in format = 6
    df = pd.read_csv(blast_result, sep='\t',header = None,names = header_outfmt6).sort_values(by='bitscore',ascending = False) #convert the csv blast file in a dataframe
    df['sense'] = '' #creates a new column for sense of hit
    df['bed_name'] = '' #creates a new column with bedtools formated name
    df['tag'] = '' #creates a new columns for tag of blast VIR or HOST
    df.to_csv(blast_result+'.csv',sep='\t') #convert de dataframe in a csv file
    rows = list()
    with open(blast_result+'.csv','r') as csv_file: 
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader: #this loop verify the sense of match
            if row[0] == '':
                pass
            else:
                if float(row[7]) > float(row[8]):
                    row[13] = 'neg'
                    a = row[7]
                    b = row[8]
                    row[7] = b
                    row[8] = a
                else:
                    row[13] = 'pos'
                if tag == 'VIR':
                    row[14] = row[1]+":"+row[7]+"-"+row[8]
                    row[15] = 'VIR'
                else:
                    row[14] = row[1]
                    row[15] = 'HOST'
            rows.append(row)
        with open(blast_result+'.csv.mod', 'w') as writeFile:
            writer = csv.writer(writeFile, delimiter='\t')
            writer.writerows(rows)
    csv_file.close()
    writeFile.close()
    pd.options.display.float_format = "{:,.2f}".format
    df = pd.read_csv(blast_result+'.csv.mod', sep='\t')
    df["evalue"] = pd.to_numeric(df["evalue"], downcast="float") #this line format the evalue as float, to avoid a representation by a large nuber pd.dataframe creates, for a limitation of ndarray, numbers fewer than -9223372036854775808 (np.iinfo(np.int64).min) are converted to 0.0
    '''
    The next three line is a trick used to remove redundant hits () in in 3 decimal places, in this case, as we are
    using a blast do recovery the viral signature in queries, that represent our genome, the filter is applied by
    query name and query start and end ranges, an example:
    INPUT:
    qseqid	sseqid	pident	length	mismatch	gapopen	qstart	qend	sstart	send	evalue	bitscore
    aag2_ctg_162	AAC97621	30.636	173	108	3	130612	130100	132	294	2.43e-08	69.7
    aag2_ctg_162	AAU10897	23.611	216	163	2	130717	130073	134	348	2.52e-10	75.3
    aag2_ctg_162	AOC55195	24.535	269	197	4	130864	130073	84	351	4.49e-11	77.8
    OUTPUT:
    qseqid	sseqid	pident	length	mismatch	gapopen	qstart	qend	sstart	send	evalue	bitscore	sense
    aag2_ctg_162	AOC55195	24.535	269	197	4	130073	130864	84	351	4.49e-11	77.8	-
    '''
    df['qstart_rng'] = df.qstart.floordiv(100)
    df['qend_rng'] = df.qend.floordiv(100)
    df_2 = df.drop_duplicates(subset=['qseqid','qstart_rng','sense']).drop_duplicates(subset=['qseqid','qstart_rng','sense']).sort_values(by=['qseqid'])
    df_2 = df_2[df_2.length >= 33]
    df_3 = df_2[['qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore','sense','bed_name','tag']]
    out_csv = blast_result+'.filtred'
    df_3.to_csv(out_csv, sep='\t', index = False)
    if tag == "VIR":
        df_4 = df_2[['qseqid','qstart','qend']]
        out_bed = blast_result+'.filtred.bed'
        df_4.to_csv(out_bed, sep='\t', index = False, header=False)    
    os.remove(blast_result+'.csv')
    os.remove(blast_result+'.csv.mod')
    return(print(f'{blast_result} filtred!'))

def get_fasta(in_file,bed_file,out_file):
    """
    This function execute the bedtools getfasta.

    Keyword arguments:
    in_file - input_file, parsed with -in argument.
    bed_file - bed file, genereated along the pipeline
    out_file - output file 
    """

    get_fasta = 'bedtools getfasta -fi '+in_file+' -bed '+bed_file+' -fo '+out_file
    get_fasta = shlex.split(get_fasta)
    cmd_get_fasta = subprocess.Popen(get_fasta)
    cmd_get_fasta.wait()

def runblastdb1(query_file,database_file): ##colocar opções de blastn ou blastx
    """
    This function runs blast with putative EVEs against host genes database.

    Keyword arguments:
    query_file - input_file - parsed with -in argument
    database - pre-formated blast database - parsed with -db argument
    """

    cline = NcbiblastxCommandline(query = query_file, db = database_file, out = query_file+'.db_filter.blastx',
    outfmt = 6, evalue = 0.00001, num_threads = threads, max_target_seqs = 5,  max_intron_length = 100, soft_masking = 'true')
    stdout, stderr = cline()
    return(print(f'BLASTx with {query_file} against {database_file}: DONE'))

def compare_blasts(vir_blast_filtred,filter_blast_filtred):
    """
    This function compares 2 blast results, for queries with same ID, only the one with the major bitscore is keept. In a final step
    Only queries with tag VIR are maintained.

    Keyword arguments:
    vir_blast_filtred - filtred blast against vir database
    filter_blast_filtred - viltred blast against filter database
    """

    df_vir = pd.read_csv(vir_blast_filtred, sep='\t')
    df_vir['qseqid'] = df_vir['bed_name']
    df_host = pd.read_csv(filter_blast_filtred, sep='\t')
    df_hybrid = pd.concat([df_vir, df_host], ignore_index=True)
    df_hybrid = df_hybrid.sort_values(by=['qseqid','bitscore'], ascending = False)
    df_hybrid.to_csv(filter_blast_filtred+'.concat', sep='\t', index = False)
    df_nr = df_hybrid.drop_duplicates(subset=['qseqid'])
    df_nr.to_csv(filter_blast_filtred+'.concat.nr', sep='\t', index = False)
    df_nr_vir = df_nr[df_nr.tag == 'VIR']
    df_nr_vir.to_csv(filter_blast_filtred+'.concat.nr.vir', sep='\t', index = False)
    return(print("DONE"))

def runblastdb2(query_file,database_file):
    """
    This function runs blast with putative EVEs against mosquito transposon database.

    Keyword arguments:
    query_file - input_file - parsed with -in argument
    database - pre-formated blast database - parsed with -db argument
    """

    if blast_type == 'blastn':
        cline = NcbiblastnCommandline(query = query_file, db = database_file, out = query_file+'.blast',
        outfmt = 6, evalue = 0.00001, num_threads = threads, task = 'blastn', max_target_seqs = 5)
        stdout, stderr = cline()
        return(print(f'BLASTn with {query_file} against {database_file}: DONE'))
    if blast_type == 'blastx':
        cline = NcbiblastxCommandline(query = query_file, db = database_file, out = query_file+'.blast',
        outfmt = 6, evalue = 0.00001, num_threads = threads, max_target_seqs = 5)
        stdout, stderr = cline()
        return(print(f'BLASTx with {query_file} against {database_file}: DONE'))
    if blast_type == 'tblastx':
        cline = NcbitblastxCommandline(query = query_file, db = database_file, out = query_file+'.blast',
        outfmt = 6, evalue = 0.00001, num_threads = threads, max_target_seqs = 5)
        stdout, stderr = cline()
        return(print(f'tBLASTx with {query_file} against {database_file}: DONE'))

 
def cut_seq(seq_file,cutoff):
    """
    This function remove contigs bellow the threshold.

    Keyword arguments:
    seq_file - fasta file, parsed with -in argument
    cutoff - cutoff length, parsed with -ln argument
    """

    new_sequences = []
    input_handle=open(seq_file,'r')
    output_handle=open(seq_file+'.fmt','w')
    for record in SeqIO.parse(input_handle,'fasta'):
        if len(record.seq) >= int(cutoff):
            new_sequences.append(record)
    SeqIO.write(new_sequences,output_handle,"fasta")
    return(print(f'Sequences bellow than {length_cutoff} bp are removed from {input_file}, results are stored in {prefix}.rn.fmt!'))

def mask_clean(in_file, m_per):
    """
    This function execute the cleanning step.

    Keyword arguments:
    in_file - fasta file, with putative EVEs.
    m_per - treshold masked percentage value, parsed with -mp argument.
    """

    sequences={}
    for seq_record in SeqIO.parse(in_file, "fasta"):
        sequence = str(seq_record.seq)
        if (float(sequence.count("a")+sequence.count("t")+sequence.count("c")+sequence.count("g")+sequence.count("n")+sequence.count("N"))/float(len(sequence)))*100 <= float(m_per):
            if sequence not in sequences:
                sequences[sequence] = seq_record.id
    with open(in_file+'.cl', "w+") as output_file:
        for sequence in sequences:
            output_file.write(">" + sequences[sequence] + "\n" + sequence + "\n")
    return(print(f'Sequences with {m_per} percent of lower-case letters are removed, results are stored in {in_file}.clean!'))

if __name__ == '__main__':
    '''
    Main Routine
    This block of code is executed, whenever the script
    is started from the command line.
    '''
    ###############################>EVE-FINDER<###############################
    if step == 'All' or step == 'EVEfinder':
        with open(input_file,'r') as input_file, open(prefix+'.rn','w+') as output_file:
            read_file = input_file.read()
            for line in read_file:
                if '>' in line:
                    line = re.sub(r">",">"+prefix+"/",line)
                output_file.write(line)
        if make_db == True: #making blastdb
            makeblastdb(database_file,'prot')
        #remove contigs bellow the threshold
        cut_seq(prefix+'.rn',length_cutoff)
        #running blastx against viral proteins
        runblastx(prefix+'.rn.fmt',database_file) ##FEITA
        #filtring blastx outputs and generating bed file
        blastx_filter(prefix+'.rn.fmt.blastx',"VIR")
        #extracting sequences based on query start and query end matches
        get_fasta(prefix+".rn.fmt",prefix+'.rn.fmt.blastx.filtred.bed',prefix+'.rn.fmt.blastx.filtred.bed.fasta')
        makeblastdb(database_filter,'prot')
        runblastdb1(prefix+'.rn.fmt.blastx.filtred.bed.fasta',database_filter)
        blastx_filter(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx',"HOST")
        compare_blasts(prefix+'.rn.fmt.blastx.filtred',prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred')
        
    ###############################>TAXONOMY<###############################
    if step == 'All' or step == 'Tax':
        error_tax = []
        prot_id = []
        count = 0
        #create temp files with intermediate entrez seachers
        with open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir','r') as putative_EVEs, open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tmp','w+',newline ='') as tax_out_tmp:
            putative_EVEs_reader = csv.reader(putative_EVEs, delimiter='\t')
            tax_tmp_writer = csv.writer(tax_out_tmp, delimiter='\t')
            for line in putative_EVEs_reader:
                prot_id.append(line[1])
            prot_id = list(dict.fromkeys(prot_id))
            count = 1
            tax_reference_list = []
            tax_levels_list = []
            #first round of search, get protein product and virus ID
            for var_prot_id in prot_id:
                try:
                    handle_code = Entrez.efetch(db = 'protein', id = var_prot_id, rettype = 'gb')
                    time.sleep(0.15)
                    handle_code2 = Entrez.efetch(db = 'protein', id = var_prot_id, rettype = 'gb')
                    handle_prot_var = str(re.findall(r'\/product=.*',handle_code.read()))
                    handle_code_var = str(re.findall(r'\/db_xref="taxon:.*',handle_code2.read()))
                    handle_code_var = re.sub(r"\[.*:", '', handle_code_var)
                    handle_code_var = re.sub(r'".*\]', '', handle_code_var)
                    handle_prot_var = re.sub(r'.*=', '', handle_prot_var)
                    handle_prot_var = re.sub(r"'\]", '', handle_prot_var)
                    handle_prot_var = re.sub(r'"', '', handle_prot_var)
                    tax_reference_list.append([var_prot_id,handle_code_var,handle_prot_var])
                    print(f'{count}/{len(prot_id)}')
                    count += 1
                    time.sleep(0.15)
                except:
                    print(f'Error with {var_prot_id}')
            count = 1
            #second round of search, get virus taxonomy
            for i in tax_reference_list:
                sk = ''
                od = ''
                fm = ''
                gn = ''
                sp = ''
                prot_id = i[0].rstrip('\n')
                tax_id = i[1].rstrip('\n')
                prot_prod = i[2].rstrip('\n')
                try:
                    handle_tax_data = Entrez.efetch(db="Taxonomy",id=tax_id, retmode="xml")
                    record_tax_data = Entrez.read(handle_tax_data)
                    for x in record_tax_data[0]["LineageEx"]:
                        if 'superkingdom' in x.values():
                            sk = x['ScientificName']
                        if 'order' in x.values():
                            od = x['ScientificName']
                        if 'family' in x.values():
                            fm = x['ScientificName']
                        if 'genus' in x.values():
                            gn = x['ScientificName']
                        if 'species' in x.values():
                            sp = x['ScientificName']
                except:
                    print(f'effor with: {tax_id}')
                    continue
                tax_levels_list.append([prot_id,prot_prod,tax_id,sk,od,fm,gn,sp])
                print(f'{count}/{len(tax_reference_list)}')
                count += 1
                time.sleep(0.15)
            tax_tmp_writer.writerows(tax_levels_list)
        tax_out_tmp.close()
        tax_complete_list = []
        with open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir','r') as eve_infos:
            eve_infos_reader = csv.reader(eve_infos, delimiter='\t')
            for eve_line in eve_infos_reader:
                if "qseqid" in eve_line:
                    pass
                else:
                    with open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tmp','r') as tax_out_tmp:
                        tax_tmp_reader = csv.reader(tax_out_tmp, delimiter='\t')
                        for tax_line in tax_tmp_reader:
                            if eve_line[1] == tax_line[0]:
                                tax_complete_list.append([eve_line[0],eve_line[1],tax_line[1],tax_line[2], tax_line[3], tax_line[4],tax_line[5],tax_line[6], tax_line[7],
                                eve_line[2],eve_line[3],eve_line[6],eve_line[7],eve_line[10],eve_line[11],eve_line[12]])
        df_tax = pd.DataFrame(tax_complete_list, columns=["Element-ID","Protein-ID","ProteinProduct","Virus-ID","Super-Kingdom","Order","Family","Genus","Species","Pident",
        "Length","Qstart","Qend","Evalue","Bitscore","Sense"])
        filter = df_tax['ProteinProduct'].str.contains("ypothetical")
        df_tax = df_tax[~filter]
        df_tax = df_tax.replace(r'^\s*$', np.nan, regex=True)
        df_tax.to_csv(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax', sep='\t', index = False, na_rep='Unclassified')
        df_to_bed = df_tax[['Element-ID','Order','Family',"Qstart",'Qend','Protein-ID','Sense']]
        df_to_bed = df_to_bed.sort_values(by=['Element-ID','Order','Family','Qstart'])
        df_to_bed.to_csv(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.tmp', sep ='\t',index = False, header = False,na_rep='Unclassified')
        with open (prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.tmp','r') as tmp_to_bed:
            tmp_header = ['Element','Order','Family','Start','End','Protein','Sense']
            df = pd.read_csv(tmp_to_bed, sep='\t',header = None,names = tmp_header)
            df['Element'] = df['Element'].str.replace(r":.*","")
            names_var = df[['Element', 'Order','Family','Sense']].agg('|'.join, axis=1)
            protein_sense = df[["Protein","Sense"]].agg('|'.join, axis=1)
            df_2 = df[['Start','End']].copy()
            df_2.insert(0, 'Bed_name', names_var)
            df_2.insert(3,'Protein_Sense',protein_sense)
            df_2 = df_2.sort_values(['Bed_name','Start'], ascending = (True, True))
            df_2.to_csv(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed', sep='\t', index = False, header=False)
        with open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge',"w") as bed_merge_out: #abm = annotation bed merged
            bed_merge_cmd = 'bedtools merge -d '+limit_merge+' -i '+prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed'+' -c 4 -o collapse -delim "AND"'
            bed_merge_cmd = shlex.split(bed_merge_cmd)
            bed_merge_process = subprocess.Popen(bed_merge_cmd, stdout=bed_merge_out)
            bed_merge_process.wait()
        with open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt',"w") as bed_merge_fmt:
            remove_annot = "sed -e 's/|.*ae|pos//g' -e 's/|.*ied|pos//g' -e 's/|.*ae|neg//g' -e 's/|.*ied|neg//g' "+prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge' # format bed_name
            remove_annot = shlex.split(remove_annot) # Split the string using shell-like syntax.
            cmd_remove_annot = subprocess.Popen(remove_annot, stdout=bed_merge_fmt)
            cmd_remove_annot.wait()
        get_fasta(prefix+".rn.fmt",prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt',prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta')
        #get information of merged elements
        with open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt','r') as bed_merge_file, open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.tax','w') as bed_merge_tax_out:
            bed_merge_tax_list = []
            bed_merge_file_reader = csv.reader(bed_merge_file,delimiter='\t')
            bed_merge_tax_out_writer = csv.writer(bed_merge_tax_out,delimiter='\t')
            bed_merge_tax_out_writer.writerow(["Element-ID","Protein-IDs","ProteinProducts","Order","Family","Genus","Species"])
            count = 0
            for line in bed_merge_file_reader:
                element_merged_id = line[0].rstrip('\n')+":"+line[1].strip('\n')+'-'+line[2].strip('\n')
                if 'pos' in line[3]:
                    sense = 'pos'
                    line[3] = re.sub('\|pos','',line[3]).rstrip('\n')
                elif 'neg' in line[3]:
                    sense = 'neg'
                    line[3] = re.sub('\|neg','',line[3]).rstrip('\n')
                protein_ids = line[3].rstrip('\n')
                with open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tmp','r') as prot_info:
                    prot_info_reader = csv.reader(prot_info,delimiter='\t')
                    protein_terms = ""
                    for line_prot in prot_info_reader:
                        if "AND" in protein_ids:
                            protein_ids = re.sub('AND',"|",line[3]).rstrip('\n')
                            if line_prot[0].rstrip('\n') in protein_ids:
                                if line_prot[1] not in protein_terms:
                                    protein_terms += line_prot[1]+' AND '
                                    vir_order = line_prot[4]
                                    vir_family = line_prot[5]
                                    vir_genus = line_prot[6]
                                    vir_specie = line_prot[7]
                        else:
                            if line_prot[0].rstrip('\n') in protein_ids:
                                protein_terms = line_prot[1]
                                vir_order = line_prot[4]
                                vir_family = line_prot[5]
                                vir_genus = line_prot[6]
                                vir_specie = line_prot[7]
                if vir_order == '':
                    vir_order = "Unclassified"
                if vir_family == '':
                    vir_family = "Unclassified"
                if vir_genus == '':
                    vir_genus = "Unclassified"
                if vir_specie == '':
                    vir_specie = "Unclassified"
                count+=1
                print(count)
                bed_merge_tax_list.append([element_merged_id,protein_ids,protein_terms,vir_order,vir_family,vir_genus,vir_specie])
            bed_merge_tax_out_writer.writerows(bed_merge_tax_list)
    ###############################>FILTER<###############################
    if step == 'All' or step == 'TEs':      
        #run blast againts transposon database (only EVEs that can be putative TEs)
        
        if database_TE != None:
            transposon_terms = ["recombinase","Recombinase","transposase","Transposase","Helicase","helicase","LINE","gag-like","RNA-dependent DNA polymerase","repeat"]
            putative_transposon_tax = ["Ortervirales","Phycodnaviridae","Polydnaviridae","Mimiviridae","Unclassified","Iridoviridae","Caudovirales","Poxviridae","Baculoviridae","Herpesviridae","Marseilleviridae"]
            with open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.tax','r') as tax_file, open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta','r') as EVEs_merged_seq, open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.pTEs','w') as putative_TEs_seq:
                putative_TEs_list = []
                result_records = []
                csv_reader = csv.reader(tax_file,delimiter='\t')
                for line in csv_reader:
                    if any(tax in line[3] for tax in putative_transposon_tax) or any(tax in line[4] for tax in putative_transposon_tax) or any(term in line[2] for term in transposon_terms):
                        putative_TEs_list_ID = line[0].replace('>','').rstrip('\n')
                        putative_TEs_list.append(putative_TEs_list_ID)
                record_dict = SeqIO.to_dict(SeqIO.parse(EVEs_merged_seq, 'fasta'))
                for id_ in putative_TEs_list:
                    if id_ in record_dict:
                        result_records.append(record_dict[id_])
                SeqIO.write(result_records,  putative_TEs_seq, "fasta")
            if blast_type == 'blastn' or blast_type == 'tblastx':
                makeblastdb(database_TE,'nucl')
            elif blast_type == 'blastx':
                makeblastdb(database_TE,'prot')
            runblastdb2(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.pTEs',database_TE)
            with open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta','r') as EVEs_merged_seq, open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.pTEs.blast','r') as blast_TEs_out: #creates a read file for blast output
                lines = EVEs_merged_seq.readlines() # read lines
                all_blast_TEs_lines = blast_TEs_out.read() # read lines of blast output
                no_match_TE = [] # creates a list of no-matched sequences
                for line in lines: # read each line of query file
                    if ">" in line: # for each header in query file
                        var_seq_id = line.replace('>','').rstrip('\n') # remove '>' and '\n' for comparrisons against blast output file
                        if (all_blast_TEs_lines.find(var_seq_id)) != -1: # if header of query file exist in blast output, nothing occours
                            pass
                        else: # else, append header name in 'no_match' list
                            no_match_TE.append(var_seq_id)
            for i in no_match_TE:
                print(i)
            with open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta','r') as EVEs_merged_seq, open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.fasta','w') as output_seq:
                record_dict = SeqIO.to_dict(SeqIO.parse(EVEs_merged_seq,'fasta')) # creates a dictionary with sequences of query file
                result_records = [record_dict[id_] for id_ in no_match_TE] # search for no_match sequences in query file
                SeqIO.write(result_records, output_seq, "fasta") # write an output file with queries that don't present match in previous analysis

        #remove EVEs with soft-masked,hard-masked or regions enriched with 'NNNNN' characters
        mask_clean(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.fasta',mask_per)
        with open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.fasta','r') as with_repeat_eves, open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.fasta.cl','r') as without_repeat_eves: #cl indicate 'cleaned' EVEs
            with_repeat_eves_lines = with_repeat_eves.readlines()
            without_repeat_eves_lines = without_repeat_eves.readlines()
            with_repeat_eves_IDs = []
            without_repeat_eves_IDs = []
            for line in with_repeat_eves_lines:
                if ">" in line:
                    line = line.replace('>','').rstrip('\n')
                    with_repeat_eves_IDs.append(line)
            for line in without_repeat_eves_lines:
                if ">" in line:
                    line = line.replace('>','').rstrip('\n')
                    without_repeat_eves_IDs.append(line)

        #create taxonomy files for all EVES and 'cleaned' EVEs
        with open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.tax','r') as tax_file, open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.fasta.tax','w') as with_repeat_eves_tax_out, open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.fasta.cl.tax','w') as without_repeat_eves_tax_out:
            csv_reader = csv.reader(tax_file,delimiter='\t')
            writer_with_repeat_eves_tax_out = csv.writer(with_repeat_eves_tax_out, delimiter='\t')
            writer_without_repeat_eves_tax_out = csv.writer(without_repeat_eves_tax_out, delimiter='\t')
            with_repeat_rows = []
            without_repeat_rows = []
            for line in csv_reader:
                for ID in with_repeat_eves_IDs:
                    if ID in line:
                        with_repeat_rows.append(line)
                for ID in without_repeat_eves_IDs:
                    if ID in line:
                        without_repeat_rows.append(line)
            writer_with_repeat_eves_tax_out.writerow(["Element-ID","Protein-IDs","ProteinProducts","Order","Family","Genus","Species"])
            writer_without_repeat_eves_tax_out.writerow(["Element-ID","Protein-IDs","ProteinProducts","Order","Family","Genus","Species"])
            writer_with_repeat_eves_tax_out.writerows(with_repeat_rows)
            writer_without_repeat_eves_tax_out.writerows(without_repeat_rows)
    
    ###############################>FLANKING-REGIONS<###############################
    if step == 'All' or step == 'Flank':
        with open(prefix+'.rn.fmt.length','w') as output_length:
            length_list = []
            #writer_output_length = csv.writer(output_length,delimiter='\t')
            for seq_record in SeqIO.parse(prefix+'.rn.fmt',"fasta"):
                output_length.write(seq_record.id+'\t'+str(len(seq_record))+'\n')
        
        #extrac flanking regions of both EVE files 'All EVEs' and 'Cleaned EVEs'
        #3 files are generated: '.fl.fasta' with both flanking regions and the EVE in the middle, '.fl.left.fasta' with only the left flank, 'fl.right.fasta' with only right flank
        with open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.fasta','r') as with_repeat_eves, open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.fasta.cl','r') as without_repeat_eves, open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.fasta.bed','w') as with_repeat_eves_bed_out, open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.fasta.cl.bed','w') as without_repeat_eves_bed_out:
            with_repeat_eves_lines = with_repeat_eves.readlines()
            without_repeat_eves_lines = without_repeat_eves.readlines()
            for line in with_repeat_eves_lines:
                if '>' in line:
                    line_name = line.replace('>','')
                    line_name = re.sub(':.*','',line_name).rstrip('\n')
                    line_start = re.sub('.*:','',line)
                    line_start = re.sub('-.*','',line_start).rstrip('\n')
                    line_end = re.sub('.*-','',line).rstrip('\n')
                    with_repeat_eves_bed_out.write(line_name+'\t'+line_start+'\t'+line_end+'\n')
            for line in without_repeat_eves_lines:
                if '>' in line:
                    line_name = line.replace('>','')
                    line_name = re.sub(':.*','',line_name).rstrip('\n')
                    line_start = re.sub('.*:','',line)
                    line_start = re.sub('-.*','',line_start).rstrip('\n')
                    line_end = re.sub('.*-','',line).rstrip('\n')
                    without_repeat_eves_bed_out.write(line_name+'\t'+line_start+'\t'+line_end+'\n')
        
        ###Cleaned EVEs###
        with open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.cl.fl.bed',"w") as flank_out:
            bed_flank_cmd = 'bedtools flank -i '+prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.fasta.cl.bed'+' -g '+prefix+'.rn.fmt.length'+' -b '+str(flank_region)
            bed_flank_cmd = shlex.split(bed_flank_cmd)
            bed_flank_process = subprocess.Popen(bed_flank_cmd, stdout=flank_out)
            bed_flank_process.wait()
        with open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.cl.fl.bed','r') as flank_bed, open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.cl.fl.bed.fmt','w') as flank_bed_fmt:
            reader = csv.reader(flank_bed, delimiter='\t') 
            for line in reader:
                nextline = next(reader)
                try:
                    flank_bed_fmt.write(line[0].rstrip('\n')+'\t'+line[1].rstrip('\n')+'\t'+nextline[2]+'\n')
                except:
                    pass
        with open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.cl.fl.bed','r+') as flank_bed, open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.cl.fl-left.bed.fmt','w') as flank_left_bed_fmt, open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.cl.fl-right.bed.fmt','w') as flank_right_bed_fmt:
            reader = csv.reader(flank_bed, delimiter='\t') 
            for line in reader:
                nextline = next(reader)
                try:
                    flank_left_bed_fmt.write(line[0].rstrip('\n')+'\t'+line[1].rstrip('\n')+'\t'+line[2]+'\n')
                    flank_right_bed_fmt.write(nextline[0].rstrip('\n')+'\t'+nextline[1].rstrip('\n')+'\t'+nextline[2]+'\n')
                except:
                    pass

        get_fasta(prefix+'.rn.fmt',prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.cl.fl.bed.fmt',prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.cl.fl.fasta')
        get_fasta(prefix+'.rn.fmt',prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.cl.fl-left.bed.fmt',prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.cl.fl.left.fasta')
        get_fasta(prefix+'.rn.fmt',prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.cl.fl-right.bed.fmt',prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.cl.fl.right.fasta')

        ###All EVEs###
        with open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.fasta.fl.bed',"w") as flank_out:
            bed_flank_cmd = 'bedtools flank -i '+prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.fasta.bed'+' -g '+prefix+'.rn.fmt.length'+' -b '+str(flank_region)
            bed_flank_cmd = shlex.split(bed_flank_cmd)
            bed_flank_process = subprocess.Popen(bed_flank_cmd, stdout=flank_out)
            bed_flank_process.wait()
        with open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.fasta.fl.bed','r+') as flank_bed, open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.fasta.fl.bed.fmt','w') as flank_bed_fmt, open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.fasta.fl-left.bed.fmt','w') as flank_left_bed_fmt, open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.fasta.fl-right.bed.fmt','w') as flank_right_bed_fmt:
            reader = csv.reader(flank_bed, delimiter='\t') 
            for line in reader:
                nextline = next(reader)
                try:
                    flank_bed_fmt.write(line[0].rstrip('\n')+'\t'+line[1].rstrip('\n')+'\t'+nextline[2]+'\n')
                except:
                    pass
        with open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.fasta.fl.bed','r+') as flank_bed, open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.fasta.fl-left.bed.fmt','w') as flank_left_bed_fmt, open(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.fasta.fl-right.bed.fmt','w') as flank_right_bed_fmt:
            reader = csv.reader(flank_bed, delimiter='\t') 
            for line in reader:
                nextline = next(reader)
                try:
                    flank_left_bed_fmt.write(line[0].rstrip('\n')+'\t'+line[1].rstrip('\n')+'\t'+line[2]+'\n')
                    flank_right_bed_fmt.write(nextline[0].rstrip('\n')+'\t'+nextline[1].rstrip('\n')+'\t'+nextline[2]+'\n')
                except:
                    pass
        get_fasta(prefix+'.rn.fmt',prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.fasta.fl.bed.fmt',prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.fl.fasta')
        get_fasta(prefix+'.rn.fmt',prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.fasta.fl-left.bed.fmt',prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.fl.left.fasta')
        get_fasta(prefix+'.rn.fmt',prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.fasta.fl-right.bed.fmt',prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.fl.right.fasta')

    ###############################>FLANKING-BLASTs<###############################
    if step == 'All' or step == 'Flank-blast':
        runblastdb2(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.fl.left.fasta',database_TE)
        runblastdb2(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.fl.right.fasta',database_TE)
        runblastdb2(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.cl.fl.left.fasta',database_TE)
        runblastdb2(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.cl.fl.right.fasta',database_TE)
        blastx_filter(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.fl.left.fasta.blast','HOST')
        blastx_filter(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.fl.right.fasta.blast','HOST')
        blastx_filter(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.cl.fl.left.fasta.blast','HOST')
        blastx_filter(prefix+'.rn.fmt.blastx.filtred.bed.fasta.db_filter.blastx.filtred.concat.nr.vir.tax.bed.merge.fmt.fasta.noTEs.cl.fl.right.fasta.blast','HOST')