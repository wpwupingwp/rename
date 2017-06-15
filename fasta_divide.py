#!/usr/bin/python3

from timeit import default_timer as timer
from Bio import SeqIO
import argparse
import os
import re


def get_choice(fasta, sep):
    # read line until get id
    with open(fasta, 'r') as raw:
        for line in raw:
            if line.startswith('>'):
                break
    print('This is raw id of one of sequence you give:')
    print(line, '\n')
    splitted = re.split(sep, line[1:])
    for index, match in enumerate(splitted, 1):
        print('{}.{}'.format(index, match))
        index += 1
    choice = input(
        'Choose seperator you want by input number of the field')
    return int(choice)


def divide(fasta, sep, choice, out):
    for record in SeqIO.parse(fasta, 'fasta'):
        item = re.split(sep, record.id)
        # remove illegal characters
        item = item[choice-1]
        with open(os.path.join(out, item+'.fasta'), 'a') as output:
            SeqIO.write(record, output, 'fasta')


def main():
    start = timer()
    args = argparse.ArgumentParser(description=main.__doc__)
    args.add_argument('input', help='input file as fasta format')
    args.add_argument('-s', type=int, dest='choice',
                      help='the field you want to use')
    args.add_argument('-o', '--output', default='out',
                            help='output directory')
    args = args.parse_args()

    SEP = re.compile(r'[\-\|/\\:;~!\?@#$%^&\*+=]')
    os.mkdir(args.out)
    if args.choice is None:
        choice = get_choice(args.input, SEP)
    start = timer()
    divide(args.input, SEP, choice, args.out)

    end = timer()
    print('''\nFinished with {0:.3f} s. You can find fasta file in the folder
    {1}.'''.format(end-start, args.out))


if __name__ == '__main__':
    main()
