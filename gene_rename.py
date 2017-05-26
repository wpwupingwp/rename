#!/usr/bin/python3
from Bio.Seq import Seq
import re


def normalize(old_name):
    old_name = old_name.lower()
    # (trna|trn(?=[b-z]))
    s = re.compile(r'(\d+\.?\d+?)(s|rrn)')
    if old_name.startswith('trn'):
        pattern = re.compile(r'([atcgu]{3})')
        codon = Seq(re.search(pattern, old_name).group(1))
        new_name = 'trn{}{}'.format(codon.translate(), codon)
        type = 'tRNA'
    elif old_name.startswith('rrn'):
        pattern = re.compile(r'(\d+)|rrn(.+)')
        number = re.search(pattern, old_name).group(1)
        new_name = 'rrn{}'.format(number)
        type = 'rRNA'
    elif re.search(s, old_name) is not None:
        new_name = 'rrn{}'.format(re.search(s, old_name).group(1))
        type = 'rRNA'
    else:
        pattern = re.compile(r'[^a-z]*'
                             '(?P<gene>[a-z]+)'
                             '[^a-z0-9]*'
                             '(?P<suffix>[a-z]|[0-9]+)')
        match = re.search(pattern, old_name)
        gene = match.group('gene')
        suffix = match.group('suffix')
        new_name = '{}{}'.format(gene, suffix.upper())
        type = 'normal'
    return new_name, type
