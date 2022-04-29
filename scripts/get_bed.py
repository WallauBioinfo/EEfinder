import re


def get_bed(input_file, log):
    with open(f'{input_file}', 'r') as repeat_eves, open(f'{input_file}.bed', 'w') as repeat_eves_bed_out:
        repeat_eves_lines = repeat_eves.readlines()
        for line in repeat_eves_lines:
            if '>' in line:
                line_name = line.replace('>', '')
                line_name = re.sub(':.*', '', line_name).rstrip('\n')
                line_start = re.sub('.*:', '', line)
                line_start = re.sub('-.*', '', line_start).rstrip('\n')
                line_end = re.sub('.*-', '', line).rstrip('\n')
                repeat_eves_bed_out.write(
                    f'{line_name}\t{line_start}\t{line_end}\n')
        print('DONE: Get bed file!', file = log)
        return(print('DONE: Get bed file!'))


class GetBed():
    def __init__(self, input_file, log):
        self.input_file = input_file
        self.log = log

    def run_get_bed(self):
        get_bed(self.input_file, self.log)
