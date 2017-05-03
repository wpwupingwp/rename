#!/usr/bin/python3

from glob import glob
from os import mkdir
from os.path import join as join_path
from timeit import default_timer as timer
import re
# from sys import argv


def get_format():
    OLD = '1.gene|2.order|3.family|4.genus|5.species|6.accession_id|7.specimen'
    print(OLD)
    new_format = input('''Input format you want: (eg. "3461" will generate
               "family|genus|accession_id|gene"\n''')
    return new_format


def rename(file_list, new_format):
    pattern = new_format
    for fasta in file_list:
        with open(fasta, 'r') as old:
            re.split()
# to be continue


def main():
    start = timer()
    # gb = argv[1]
    file_list = glob('*.fasta')
    OUT = 'renamed'
    mkdir(OUT)

    new_format = get_format()
    rename(file_list, new_format)

    end = timer()
    print('''\nFinished with {0:.3f} s. You can find fasta file in the folder
    {1}.'''.format(end-start, OUT))
