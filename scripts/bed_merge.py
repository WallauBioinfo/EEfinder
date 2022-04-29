import shlex, subprocess

def merge_bedfile(bed_annotated_file, limit_merge, log):
    """
    This function use the bedtools merge to join regions of the same family based on a specific nucleotide range defined by the user

    Keyword arguments:
    bed_annotated_file = tsv file generated in the get_annotated_bed function on formtad_bed.py
    limit_merge = limit range number parsed with -lm parameter
    """

    bed_merge_output = open(f"{bed_annotated_file}.merge",'w')
    bed_merge_cmd = f'bedtools merge -d {limit_merge} -i {bed_annotated_file} -c 4 -o collapse -delim " AND "'
    bed_merge_cmd = shlex.split(bed_merge_cmd)
    bed_merge_process = subprocess.Popen(bed_merge_cmd, stdout=bed_merge_output)
    bed_merge_process.wait()
    print(f'DONE: Merged pEVEs of the same family within a {limit_merge} nucleotides range', file = log)
    return(print(f'DONE: Merged pEVEs of the same family within a {limit_merge} nucleotides range'))

class MergeBed():
    def __init__(self, bed_annotated_file, limit_merge, log):
        self.bed_annotated_file = bed_annotated_file
        self.limit_merge = limit_merge
        self.log = log

    def run_merge_bedfile(self):
        merge_bedfile(self.bed_annotated_file, self.limit_merge, self.log)