import shlex
import subprocess
import pandas as pd


def bed_flank(input_file, length_file, flank_region, log):
    with open(f'{input_file}.flank', "w") as flank_out:
        bed_flank_cmd = f'bedtools flank -i {input_file} -g {length_file} -b {str(flank_region)}'
        bed_flank_cmd = shlex.split(bed_flank_cmd)
        bed_flank_process = subprocess.Popen(bed_flank_cmd, stdout=flank_out)
        bed_flank_process.wait()
        print(f'DONE: Get bed {str(flank_region)}pb flanks file!', file = log)
        return(print(f'DONE: Get bed {str(flank_region)}pb flanks file!'))


def flank_formater(input_file, log):
    with open(f'{input_file}.flank.fmt', 'w') as flank_bed_fmt:
        #reader = csv.reader(flank_bed, delimiter='\t')
        df = pd.read_csv(f'{input_file}.flank', sep="\t", header=None)
        print(df)
        df = df.groupby(df.index // 2).agg(lambda x: x.dropna().astype(str).str.cat(sep=','))
        df[0] = df[0].str.replace(',.*','', regex=True)
        df[1] = df[1].str.replace(',.*','', regex=True)
        df[2] = df[2].str.replace('.*,','', regex=True)
        df.to_csv(flank_bed_fmt, sep="\t", index=False, header=None)
        print('DONE: Format bed flank files!', file = log)
        return(print('DONE: Format bed flank files!'))


class BedFlank():
    def __init__(self, input_file, length_file, flank_region, log):
        self.input_file = input_file
        self.length_file = length_file
        self.flank_region = flank_region
        self.log = log

    def run_bed_flank(self):
        bed_flank(self.input_file, self.length_file, self.flank_region, self.log)
        flank_formater(self.input_file, self.log)
