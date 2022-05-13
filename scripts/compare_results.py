import pandas as pd

def compare_blasts(vir_blast_filtred,filter_blast_filtred, log):
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
    df_nr_vir = df_nr[df_nr.tag == 'EE']
    df_nr_vir.to_csv(filter_blast_filtred+'.concat.nr.vir', sep='\t', index = False)
    print("DONE: Filter step", file = log)
    return(print("DONE: Filter step"))

class CompareResults:
    def __init__(self, vir_result, host_result, log):
        self.vir_result = vir_result
        self.host_result = host_result
        self.log = log
    
    def run_comparation(self):
        compare_blasts(self.vir_result, self.host_result, self.log)
        