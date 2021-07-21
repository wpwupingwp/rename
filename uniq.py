#!/usr/bin/python3

from Bio import SeqIO
import argparse
import re


def parse_arg():
    arg = argparse.ArgumentParser(description=main.__doc__)
    arg.add_argument('input', help='input file as fasta format')
    arg.add_argument('-c', type=str, dest='choice',
                     help='the field you want to use')
    arg.add_argument('-m', '--method', default='longest',
                     help='method to get the only sequence')
    arg.print_help()
    return arg.parse_args()


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


def uniq(arg, sep):
    def get_longest(records):
        info = list()
        for record in records:
            n = record.seq.count('N') + record.seq.count('n')
            length = len(record.seq) - n
            info.append([length, record])
        longest = max(info, key=lambda x: x[0])
        return longest[1]

    if arg.choice is None:
        arg.choice = get_choice(arg.input, sep)
    choice = [int(i)-1 for i in arg.choice.split(' ')]
    raw = SeqIO.parse(arg.input, 'fasta')
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

    log = open(arg.input+'.log', 'w')
    output = open(arg.input+'.uniq', 'w')
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
    print(f'{arg.input}, before {before}, after {after}')
    log.close()
    output.close()


def main():
    arg = parse_arg()
    SEP = re.compile(r'[\|/\\:;~!\?@#$%^&\*+=]')
    uniq(arg, SEP)


if __name__ == '__main__':
    main()
