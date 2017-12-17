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
    print(sep)
    splitted = re.split(sep, line[1:])
    for index, match in enumerate(splitted, 1):
        print('{}.{}'.format(index, match))
        index += 1
    choice = input(('Choose seperator you want by input number of the field'
                   ' (use space to seperate multiple fields):\n'))
    return choice


def divide(fasta, sep, choice, out):
    choice = [int(i)-1 for i in choice.split(' ')]
    for record in SeqIO.parse(fasta, 'fasta'):
        new_item = list()
        items = re.split(sep, record.id)
        try:
            for i in choice:
                if items[i] == '':
                    new_item.append('NOT_FOUND')
                else:
                    new_item.append(items[i])
        except IndexError:
            new_item.append('NOT_FOUND')
        with open(os.path.join(out, '-'.join(new_item)+'.fasta'),
                  'a') as output:
            SeqIO.write(record, output, 'fasta')


def main():
    start = timer()
    args = argparse.ArgumentParser(description=main.__doc__)
    args.add_argument('input', help='input file as fasta format')
    args.add_argument('-o', dest='out', help='output folder')
    args.add_argument('-s', dest='sep', help='seperator, default "|"')
    args.add_argument('-c', type=int, dest='choice',
                      help='the field you want to use')
    args = args.parse_args()

    if args.out is None:
        out_folder = os.path.splitext(args.input)[0]
        out_folder = '_'.join([out_folder, 'out'])
    else:
        out_folder = args.out
    try:
        os.mkdir(out_folder)
    except:
        pass
    if args.sep is None:
        args.sep = re.compile(r'[\-\|/\\:;~!\?@#$%^&\*+=]')
    else:
        args.sep = re.compile(r'|')
    if args.choice is None:
        args.choice = get_choice(args.input, args.sep)
    start = timer()
    divide(args.input, args.sep, args.choice, out_folder)

    end = timer()
    print('''
Finished with {0:.3f} s. You can find fasta file in the folder {1}. Some
sequences cannot be divided because of lack of info and they were put
into failed.fasta'''.format(end-start, out_folder))


if __name__ == '__main__':
    main()
