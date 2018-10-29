#!/usr/bin/python3

from Bio import SeqIO
from timeit import default_timer as timer
import argparse
import re


def parse_args():
    args = argparse.ArgumentParser(description=main.__doc__)
    args.add_argument('input', help='input file as fasta format')
    args.add_argument('-c', type=str, dest='choice',
                      help='the field you want to use')
    args.add_argument('-m', '--method', default='longest',
                      help='method to get the only sequence')
    args.print_help()
    return args.parse_args()


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
    choice = input('Choose fields you want by input number of the field and'
                   'use space as seprator :\n')
    return choice


def uniq(args, sep):
    def get_longest(records):
        info = list()
        for record in records:
            n = record.seq.count('N') + record.seq.count('n')
            length = len(record.seq) - n
            info.append([length, record])
        longest = max(info, key=lambda x: x[0])
        return longest[1]

    if args.choice is None:
        args.choice = get_choice(args.input, sep)
    choice = [int(i)-1 for i in args.choice.split(' ')]
    raw = SeqIO.parse(args.input, 'fasta')
    name_seq = dict()
    before = 0
    after = 0
    for record in raw:
        raw_name = re.split(sep, record.id)
        name = list()
        for index, item in enumerate(raw_name):
            if index in choice:
                name.append(item)
        name = '.'.join(name)
        if name in name_seq:
            name_seq[name].append(record)
        else:
            name_seq[name] = [record, ]
        before += 1

    print('Duplicated sequences:')
    log = open(args.input+'.log', 'w')
    output = open(args.input+'.uniq', 'w')
    for record in name_seq.values():
        after += 1
        longest = get_longest(record)
        SeqIO.write(longest, output, 'fasta')
        if len(record) != 1:
            id_list = [i.id for i in record]
            log.write('Longest:\t{} in ({})\n'.format(
                longest.id, '\t'.join(id_list)))
    log.write('Before\tAfter\n')
    log.write('{}\t\t{}\n'.format(before, after))
    print('Total {} sequences in {} format.'.format(before, 'fasta'))
    print('{} sequences left in the file {}.uniq.'.format(after, args.input))
    log.close()
    output.close()


def main():
    start = timer()
    args = parse_args()

    SEP = re.compile(r'[\|/\\:;~!\?@#$%^&\*+=]')
    uniq(args, SEP)
    end = timer()
    print('Cost {:.3f} seconds.'.format(end-start))


if __name__ == '__main__':
    main()
