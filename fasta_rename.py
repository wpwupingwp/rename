#!/usr/bin/python3

import re
from sys import argv
from timeit import default_timer as timer


def get_format_string(example, SEP):
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
    print('''Choose fields you want by input numbers with any order
you want. Seperators cannot be omitted and you should not use space or
underline. For example, "3|1|2|3!4#1"''')
    new_format = input()
    if '_' in new_format or ' ' in new_format:
        print('To avoid possible errors, space(" ") and underline("_")'
              ' is prohibited.')
        option = input('Still continue?  [y/n]')
        if option.startswith(('y', 'Y')):
            pass
        else:
            return None
    print('Your new id will look like '
          'this:\n{}'.format(get_format(new_format).format(*splitted)))
    return new_format


def get_format(string):
    def minus_one(letter):
        return '{{{}}}'.format((int(letter)-1))

    new_format = re.sub(r'(?<!\\)(\d+)(?!\d)', lambda match: minus_one(
        match.group(1)), string)
    new_format = new_format.replace('\\', '')
    return new_format


def rename(old_file, new_format, SEP):
    old = open(old_file, 'r')
    new = open('rename-'+old_file, 'w')
    for line in old:
        if line[0] != '>':
            new.write(line)
        else:
            line = line.strip()
            try:
                line = new_format.format(*(re.split(SEP, line[1:])))
            except IndexError:
                print('Skip invalid sequence id:')
                print(line)
                continue
            new.write('>{}\n'.format(line))
    old.close()
    new.close()


def main():
    fasta_file = argv[1]
    SEP = r'[\|/:;~!\?@#$%^&\*+=]'

    if len(argv) != 3:
        format_string = get_format_string(fasta_file, SEP)
    else:
        format_string = argv[2].strip('"')
    new_format = get_format(format_string)
    start = timer()
    rename(fasta_file, new_format, SEP)
    end = timer()
    print('Finished with {0:.3f} s.'.format(end-start))


if __name__ == '__main__':
    main()
