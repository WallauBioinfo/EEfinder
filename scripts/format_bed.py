import pandas as pd

def get_annotated_bed(blast_tax_info):
    """
    This function create a bed file that will be used to merge truncated EVEs of the same family in the same sense based on a limite length treshold.

    Keyword arguments:
    blast_tax_info = csv file generated in the get_taxonomy_info function on get_taxonomy.py
    """

    df_blast_tax_info = pd.read_csv(blast_tax_info, sep =',')
    df_blast_tax_info['qseqid'] = df_blast_tax_info['qseqid'].str.replace(r'\:.*','', regex = True)
    df_blast_tax_info['Family'].fillna('Unknown', inplace=True)
    df_blast_tax_info['formated_name'] = df_blast_tax_info['qseqid']+'|'+df_blast_tax_info['Family']+'|'+df_blast_tax_info['sense']
    bed_blast_info = df_blast_tax_info[['formated_name','qstart','qend','sseqid']].copy()
    bed_blast_info = bed_blast_info.sort_values(['formated_name','qstart'], ascending = (True, True))
    bed_blast_info.to_csv(f"{blast_tax_info}.bed", index = False, header = False, sep = '\t')
    return(print(f'DONE: Create annotated bed file.'))


def reformat_bed(bed_annotated_merged_file):
    """
    This function remove the annotated information generate into the get_annotated_bed function

    Keyword arguments:
    bed_annotated_merged_file = tsv file generated in the merge_bedfile function on bed_merge.py
    """

    df_merge_file = pd.read_csv(bed_annotated_merged_file, sep = '\t',header=None)
    df_merge_file.iloc[:,0] = df_merge_file.iloc[:,0].str.replace("\|.*","",regex=True)
    df_merge_file.to_csv(f"{bed_annotated_merged_file}.fmt", index = False, header = False, sep = '\t')
    return(print(f'DONE: Remove annotation of bed file sequences name.'))


def reformat_name(bed_annotated_merged_file_formated):
    cols = ['name','start','end','annotation']
    df_bed_annotated_merged_file_formated = pd.read_csv(bed_annotated_merged_file_formated, header = None, sep = '\t')
    df_bed_annotated_merged_file_formated.columns = cols
    df_bed_annotated_merged_file_formated['formated_name'] = df_bed_annotated_merged_file_formated['name']+':'+df_bed_annotated_merged_file_formated['start'].astype(str)+'-'+df_bed_annotated_merged_file_formated['end'].astype(str)
    df_bed_annotated_merged_file_formated = df_bed_annotated_merged_file_formated[['formated_name','start','end','annotation']]
    df_bed_annotated_merged_file_formated.to_csv(f"{bed_annotated_merged_file_formated}2", sep = '\t', index = False, header = False)
    return(print(f'DONE: Return pEVEs names to merged bed file'))

class GetAnnotBed():
    def __init__(self, blast_tax_info):
        self.blast_tax_info = blast_tax_info

    def run_get_annotated_bed(self):    
        get_annotated_bed(self.blast_tax_info)

class RemoveAnnotation():
    def __init__(self, bed_annotated_merged_file):
        self.bed_annotated_merged_file = bed_annotated_merged_file

    def run_reformat_bed(self):
        reformat_bed(self.bed_annotated_merged_file)

class ReformatBedName():
    def __init__(self, bed_annotated_merged_file_formated):
        self.bed_annotated_merged_file_formated = bed_annotated_merged_file_formated

    def run_reformat_name(self):
        reformat_name(self.bed_annotated_merged_file_formated)