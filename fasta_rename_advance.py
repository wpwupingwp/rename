#!/usr/bin/python3

from glob import glob
from os import mkdir
from os.path import join as join_path
from timeit import default_timer as timer
import re
# from sys import argv


def get_format(example):
    # read line until get id
    SEP = r'[\-\|/\\:;~!\?@#$%^&\*+=]'
    print(SEP)
    with open(example, 'r') as raw:
        while True:
            # omit ">" at the beginning of line
            id_example = raw.readline().strip()
            if id_example.startswith('>'):
                break
    print('This is raw id of one of sequence you give:')
    print(id_example, '\n')
    index = 1
    for match in re.split(SEP, id_example[1:]):
        # print('{}.{}|'.format(index, match, end=''))
        print('{}.{}'.format(index, match))
        index += 1
    print('''Choose fields you want by input numbers with any order you want.
    If you omit seperators, it will use "|" to seperate fields. For example,
    "3421" or "2!3!4#1"''')
    new_format = input()
    if set(' _').issubset(set(new_format)):
        raise Exception('''To avoid possible errors, space and underline is
                            prohibited.''')
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

    new_format = get_format(example)
    exit -1
    mkdir(OUT)
    start = timer()
    rename(file_list, new_format, OUT)

    end = timer()
    print('''\nFinished with {0:.3f} s. You can find fasta file in the folder
    {1}.'''.format(end-start, OUT))


if __name__ == '__main__':
    main()
