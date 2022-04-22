import pandas as pd
import csv
import os


def filter(blast_result, tag):
    """
    This function receives a blastx result and filter based on query ID and ranges of qstart and qend.

    Keyword arguments:
    blast_result = blastx result outfmt 6
    """

    header_outfmt6 = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart',
                      'qend', 'sstart', 'send', 'evalue', 'bitscore']  # creates a blast header output in format = 6
    df = pd.read_csv(blast_result, sep='\t', header=None, names=header_outfmt6).sort_values(
        by='bitscore', ascending=False)  # convert the csv blast file in a dataframe
    df['sense'] = ''  # creates a new column for sense of hit
    df['bed_name'] = ''  # creates a new column with bedtools formated name
    df['tag'] = ''  # creates a new columns for tag of blast VIR or HOST
    # convert de dataframe in a csv file
    df.to_csv(blast_result+'.csv', sep='\t')
    rows = list()
    with open(blast_result+'.csv', 'r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:  # this loop verify the sense of match
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
    # this line format the evalue as float, to avoid a representation by a large nuber pd.dataframe creates, for a limitation of ndarray, numbers fewer than -9223372036854775808 (np.iinfo(np.int64).min) are converted to 0.0
    df["evalue"] = pd.to_numeric(df["evalue"], downcast="float")
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
    df_2 = df.drop_duplicates(subset=['qseqid', 'qstart_rng', 'sense']).drop_duplicates(
        subset=['qseqid', 'qstart_rng', 'sense']).sort_values(by=['qseqid'])
    df_2 = df_2[df_2.length >= 33]
    df_3 = df_2[['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart',
                 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'sense', 'bed_name', 'tag']]
    out_csv = blast_result+'.filtred'
    df_3.to_csv(out_csv, sep='\t', index=False)
    if tag == "VIR":
        df_4 = df_2[['qseqid', 'qstart', 'qend']]
        out_bed = blast_result+'.filtred.bed'
        df_4.to_csv(out_bed, sep='\t', index=False, header=False)
    os.remove(blast_result+'.csv')
    os.remove(blast_result+'.csv.mod')
    return(print(f'{blast_result} filtred!'))


class FilterTable:

    def __init__(self, blastresult, tag):
        self.blastresult = blastresult
        self.tag = tag

    def run_filter(self):
        filter(self.blastresult, self.tag)
