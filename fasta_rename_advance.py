#!/usr/bin/python3

from glob import glob
from os import mkdir
from os.path import join as join_path
from timeit import default_timer as timer
# from sys import argv


def get_format(example):
    # read line until get id
    with open(example, 'r') as raw:
        while True:
            one_line = raw.readline().strip()
            if one_line.startswith('>'):
                break
    print('This is raw id of one of sequence you give:\n')
    print(one_record.id)
    new_format = input('''Input format you want: (eg. "3461" will generate
               "family|genus|accession_id|gene"\n''')
    return new_format


def rename(file_list, new_format, out):
    for fasta in file_list:
        with open(fasta, 'r') as old, open(
                join_path(out, fasta), 'w') as new:
            for line in old:
                # skip sequence
                if line[0] != '>':
                    new.write(line)
                    continue
                else:
                    pick = list()
                    fields = line[1:].split('|')
                    for n in new_format:
                        pick.append(fields[int(n)-1].strip())
                    new.write('>{}\n'.format('|'.join(pick)))
# to be continue


def main():
    # gb = argv[1]
    file_list = glob('*.fasta')
    example = file_list[-1]
    OUT = 'renamed'
    mkdir(OUT)

    new_format = get_format(example)
    start = timer()
    rename(file_list, new_format, OUT)

    end = timer()
    print('''\nFinished with {0:.3f} s. You can find fasta file in the folder
    {1}.'''.format(end-start, OUT))


if __name__ == '__main__':
    main()
