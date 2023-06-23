import shlex
import subprocess
import pandas as pd


def bed_flank(input_file, length_file, flank_region, log):
    with open(f"{input_file}.flank", "w") as flank_out:
        bed_flank_cmd = (
            f"bedtools slop -i {input_file} -g {length_file} -b {str(flank_region)}"
        )
        bed_flank_cmd = shlex.split(bed_flank_cmd)
        bed_flank_process = subprocess.Popen(bed_flank_cmd, stdout=flank_out)
        bed_flank_process.wait()
        print(f"DONE: Get bed {str(flank_region)}pb flanks file!", file=log)
        return print(f"DONE: Get bed {str(flank_region)}pb flanks file!")


class BedFlank:
    def __init__(self, input_file, length_file, flank_region, log):
        self.input_file = input_file
        self.length_file = length_file
        self.flank_region = flank_region
        self.log = log

    def run_bed_flank(self):
        bed_flank(self.input_file, self.length_file, self.flank_region, self.log)
