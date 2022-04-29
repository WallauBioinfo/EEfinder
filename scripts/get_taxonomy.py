import pandas as pd
import re, csv
from Bio import SeqIO

def get_taxonomy_info(blast_file,tax_file):
    """
    This function merge the filtred blast results with taxonomy information

    Keyword arguments:
    blast_file - tsv filtred blast reslts
    tax_file - refseq table with taxonomy and other information, parsed with -mt parameter
    """
    
    df_blast_file = pd.read_csv(blast_file, sep='\t')
    df_tax_file = pd.read_csv(tax_file)
    df_tax_file.rename(columns = {"Accession":'sseqid'}, inplace = True)
    df_merged = pd.merge(df_blast_file, df_tax_file, on="sseqid", how = 'left')
    df_merged.to_csv(f"{blast_file}.tax", index = False, header = True)
    return(print(f'DONE: Get initial taxonomy info!'))

def get_final_taxonomy(bed_formated, taxonomy_info):
    with open(bed_formated,'r') as bed_merge_file, open(f"{bed_formated}.fa.tax",'w') as bed_merge_tax_out:
        bed_merge_tax_list = []
        bed_merge_file_reader = csv.reader(bed_merge_file,delimiter='\t')
        bed_merge_tax_out_writer = csv.writer(bed_merge_tax_out,delimiter='\t')
        bed_merge_tax_out_writer.writerow(["Element-ID","Sense","Protein-IDs","Protein-Products","Molecule_type","Family","Genus","Species","Host"])
        for line in bed_merge_file_reader:
            element_merged_id = line[0].rstrip('\n')+":"+line[1].strip('\n')+'-'+line[2].strip('\n')
            if 'pos' in line[3]:
                sense = 'pos'
                line[3] = re.sub('\|pos','',line[3]).rstrip('\n')
            elif 'neg' in line[3]:
                sense = 'neg'
                line[3] = re.sub('\|neg','',line[3]).rstrip('\n')
            protein_ids = line[3].rstrip('\n')
            with open(taxonomy_info,'r') as prot_info:
                prot_info_reader = csv.reader(prot_info,delimiter=',')
                protein_terms = ""
                genus = ""
                species = ""
                host = ""
                if "AND" in protein_ids:
                    for line_prot in prot_info_reader:
                            protein_ids = re.sub('AND',"|",line[3]).rstrip('\n')
                            if line_prot[1].rstrip('\n') in protein_ids:
                                if line_prot[24].rstrip('\n') not in protein_terms:
                                    protein_terms += line_prot[24]+' AND '
                                    mol_type = line_prot[19]
                                    family = line_prot[18]
                                if line_prot[17].rstrip('\n') not in genus:
                                    genus += line_prot[17]+' AND '
                                if line_prot[16].rstrip('\n') not in species:
                                    species += line_prot[16]+' AND '
                                if line_prot[27].rstrip('\n') not in host:
                                    host +=  line_prot[27]+' AND '
                else:
                    for line_prot in prot_info_reader:
                        if line_prot[1].rstrip('\n') in protein_ids:
                                protein_terms = line_prot[24]
                                mol_type = line_prot[19]
                                family = line_prot[18]
                                genus = line_prot[17]
                                species = line_prot[16]
                                host =  line_prot[27]

            protein_terms = re.sub(r" AND $","",protein_terms)
            species = re.sub(r" AND $","",species)
            host = re.sub(r" AND $","",host)
            if mol_type == '':
                vir_order = "Undefined"
            if family == '':
                family = "Unclassified"
            if genus == '':
                genus = "Unclassified"
            if species == '':
                species = "Unclassified"
            if host == '':
                host = "Undefined"

            bed_merge_tax_list.append([element_merged_id,sense,protein_ids,protein_terms,mol_type,family,genus,species,host])
        bed_merge_tax_out_writer.writerows(bed_merge_tax_list)
        return(print(f'DONE: Get Taxonomy for Merged Elements!'))

def get_cleaned_taxonomy(cleaned_file, taxonomy_file):
    output_list = [["Element-ID","Sense","Protein-IDs","Protein-Products","Molecule_type","Family","Genus","Species","Host"]]

    for seq_record in SeqIO.parse(cleaned_file, "fasta"):
        with open(taxonomy_file,'r') as tax_file:
            taxonomy_file_reader = csv.reader(tax_file, delimiter ='\t')
            for line in taxonomy_file_reader:
                if line[0] == seq_record.id:
                    output_list.append(line)

    with open(f"{cleaned_file}.tax","w") as output_file:
        output_file_writer = csv.writer(output_file,delimiter = '\t')
        output_file_writer.writerows(output_list)
    
    return(print(f'DONE: Get Taxonomy for Cleaned Merged Elements!'))

class GetTaxonomy():
    def __init__(self, blast_file, tax_file):
        self.blast_file = blast_file
        self.tax_file = tax_file
    
    def run_get_taxonomy_info(self):
        get_taxonomy_info(self.blast_file, self.tax_file)

class GetFinalTaxonomy():
    def __init__(self, bed_formated, taxonomy_info):
        self.bed_formated = bed_formated
        self.taxonomy_info = taxonomy_info

    def run_get_final_taxonomy(self):
        get_final_taxonomy(self.bed_formated, self.taxonomy_info)

class GetCleanedTaxonomy():
    def __init__(self, cleaned_file, taxonomy_file):
        self.cleaned_file = cleaned_file
        self.taxonomy_file = taxonomy_file

    def run_get_cleaned_taxonomy(self):
        get_cleaned_taxonomy(self.cleaned_file,self.taxonomy_file)