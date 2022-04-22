import pandas as pd

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

class GetTaxonomy():
    def __init__(self, blast_file, tax_file):
        self.blast_file = blast_file
        self.tax_file = tax_file
    
    def run_get_taxonomy_info(self):
        get_taxonomy_info(self.blast_file, self.tax_file)