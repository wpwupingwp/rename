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
    print(line)
    splitted = re.split(sep, line[1:])
    for index, match in enumerate(splitted, 1):
        print('{}.{}'.format(index, match))
        index += 1
    choice = input(
        'Choose seperator you want by input number of the field:\n')
    return int(choice)


def divide(fasta, sep, choice, out):
    for record in SeqIO.parse(fasta, 'fasta'):
        item = re.split(sep, record.id)
        try:
            item = item[choice-1]
        except:
            item = 'FAILED'
        with open(os.path.join(out, item+'.fasta'), 'a') as output:
            SeqIO.write(record, output, 'fasta')


def main():
    start = timer()
    args = argparse.ArgumentParser(description=main.__doc__)
    args.add_argument('input', help='input file as fasta format')
    args.add_argument('-o', dest='out', help='output folder')
    args.add_argument('-c', type=int, dest='choice',
                      help='the field you want to use')
    args = args.parse_args()

    SEP = re.compile(r'[\-\|/\\:;~!\?@#$%^&\*+=]')
    if args.out is None:
        out_folder = os.path.splitext(args.input)[0]
        out_folder = '_'.join([out_folder, 'out'])
    else:
        out_folder = args.out
    try:
        os.mkdir(out_folder)
    except:
        pass
    if args.choice is None:
        args.choice = get_choice(args.input, SEP)
    start = timer()
    divide(args.input, SEP, args.choice, out_folder)

    end = timer()
    print('''
Finished with {0:.3f} s. You can find fasta file in the folder {1}. Some
sequences cannot be divided because of lack of info and they were put
into failed.fasta'''.format(end-start, out_folder))


if __name__ == '__main__':
    main()
