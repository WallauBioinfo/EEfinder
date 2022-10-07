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
        df_up_flank = df.iloc[::2].copy()
        df_down_flank = df.drop(0).copy()
        df_down_flank = df_down_flank.iloc[::2].copy()
        df_up_flank.columns = ["prefix", "up_start", "up_end"]
        df_down_flank.columns = ["prefix", "down_start", "down_end"]
        down_end = df_down_flank["down_end"].tolist()
        df_up_flank["down_end"] = down_end
        df_bed_flanks = df_up_flank.drop("up_end", axis = 1)
        df_bed_flanks.to_csv(flank_bed_fmt, sep="\t", index=False, header=None)
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
