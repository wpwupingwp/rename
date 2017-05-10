#!/usr/bin/python3

from glob import glob
from os import mkdir
from os.path import join as join_path
from timeit import default_timer as timer
import re
# from sys import argv


def get_format(example, SEP):
    # read line until get id
    with open(example, 'r') as raw:
        while True:
            # omit ">" at the beginning of line
            id_example = raw.readline().strip()
            if id_example.startswith('>'):
                break
    print('This is raw id of one of sequence you give:')
    print(id_example, '\n')
    splitted = re.split(SEP, id_example[1:])
    for index, match in enumerate(splitted, 1):
        # print('{}.{}|'.format(index, match, end=''))
        print('{}.{}'.format(index, match))
        index += 1
    print('''   Choose fields you want by input numbers with any order
    you want. Seperators cannot be omitted and you should not use space or
    underline. For example, "3|1|2|3!4#1"''')
    new_format = input()
    if '_' in set(new_format) or ' ' in set(new_format):
        print('To avoid possible errors, space(" ") and underline("_")'
              ' is prohibited.')
        option = input('Still continue?  [y/n]')
        if option.startswith(('y', 'Y')):
            pass
        else:
            return None

    def minus_one(letter):
        return '{{{}}}'.format((int(letter)-1))
    new_format = re.sub(r'(\d+)', lambda match: minus_one(match.group(1)),
                        new_format)
    print('Your new id will look like'
          'this:\n{}'.format(new_format.format(*splitted)))
    return new_format


def rename(file_list, new_format, out, SEP):
    for fasta in file_list:
        with open(fasta, 'r') as old, open(
                join_path(out, fasta), 'w') as new:
            for line in old:
                # skip sequence
                if line[0] != '>':
                    new.write(line)
                    continue
                else:
                    line = new_format.format(*(re.split(SEP, line[1:])))
                    new.write('>{}\n'.format(line))


def main():
    # gb = argv[1]
    file_list = glob('*.fasta')
    example = file_list[-1]
    OUT = 'renamed'
    SEP = r'[\-\|/\\:;~!\?@#$%^&\*+=]'

    new_format = get_format(example, SEP)
    if new_format is None:
        return
    mkdir(OUT)
    start = timer()
    rename(file_list, new_format, OUT, SEP)

    end = timer()
    print('''\nFinished with {0:.3f} s. You can find fasta file in the folder
    {1}.'''.format(end-start, OUT))


if __name__ == '__main__':
    main()
