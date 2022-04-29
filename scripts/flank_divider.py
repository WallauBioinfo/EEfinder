import csv


def flank_divider(input_file):
    with open(f'{input_file}', 'r+') as flank_bed, open(f'{input_file}.left', 'w') as flank_left_bed_fmt, open(f'{input_file}.right', 'w') as flank_right_bed_fmt:
        reader = csv.reader(flank_bed, delimiter='\t')
        for line in reader:
            nextline = next(reader)
            try:
                flank_left_bed_fmt.write(line[0].rstrip(
                    '\n')+'\t'+line[1].rstrip('\n')+'\t'+line[2]+'\n')
                flank_right_bed_fmt.write(nextline[0].rstrip(
                    '\n')+'\t'+nextline[1].rstrip('\n')+'\t'+nextline[2]+'\n')
            except:
                pass


class FlankDivider():
    def __init__(self, input_file):
        self.input_file = input_file

    def run_flank_divider(self):
        flank_divider(self.input_file)
