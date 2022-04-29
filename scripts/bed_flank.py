import shlex
import subprocess
import csv


def bed_flank(input_file, length_file, flank_region):
    with open(f'{input_file}.flank', "w") as flank_out:
        bed_flank_cmd = f'bedtools flank -i {input_file} -g {length_file} -b {str(flank_region)}'
        bed_flank_cmd = shlex.split(bed_flank_cmd)
        bed_flank_process = subprocess.Popen(bed_flank_cmd, stdout=flank_out)
        bed_flank_process.wait()
        print('DONE CREATING FLANKS')


def flank_formater(input_file):
    with open(f'{input_file}.flank', 'r') as flank_bed, open(f'{input_file}.flank.fmt', 'w') as flank_bed_fmt:
        reader = csv.reader(flank_bed, delimiter='\t')
        for line in reader:
            nextline = next(reader)
            try:
                flank_bed_fmt.write(line[0].rstrip(
                    '\n')+'\t'+line[1].rstrip('\n')+'\t'+nextline[2]+'\n')
            except:
                pass
        print('DONE FORMATTING FLANKS')


class BedFlank():
    def __init__(self, input_file, length_file, flank_region):
        self.input_file = input_file
        self.length_file = length_file
        self.flank_region = flank_region

    def run_bed_flank(self):
        bed_flank(self.input_file, self.length_file, self.flank_region)
        flank_formater(self.input_file)
