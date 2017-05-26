#!/usr/bin/python3
from Bio.Seq import Seq
import re


def normalize(old_name):
    old_name = old_name.lower()
    # (trna|trn(?=[b-z]))
    s = re.compile(r'(\d+\.?\d+?)(s|rrn|rdna)')
    if old_name.startswith('trn'):
        pattern = re.compile(r'([atcgu]{3})')
        try:
            codon = Seq(re.search(pattern, old_name).group(1))
        except:
            new_name = old_name
            gene_type = 'bad_name'
            return new_name, gene_type
        new_name = 'trn{}{}'.format(codon.reverse_complement().translate(),
                                    codon.transcribe())
        # if 'u' in codon:
        #     new_name = 'trn{}{}'.format(codon.translate(), codon)
        # else:
        #     new_name = 'trn{}{}'.format(codon.transcribe().translate(), codon)
        gene_type = 'tRNA'
    elif old_name.startswith('rrn'):
        pattern = re.compile(r'(\d+)|rrn(.+)')
        number = re.search(pattern, old_name).group(1)
        new_name = 'rrn{}'.format(number)
        gene_type = 'rRNA'
    elif re.search(s, old_name) is not None:
        new_name = 'rrn{}'.format(re.search(s, old_name).group(1))
        gene_type = 'rRNA'
    else:
        pattern = re.compile(r'[^a-z]*'
                             '(?P<gene>[a-z]+)'
                             '[^a-z0-9]*'
                             '(?P<suffix>[a-z]|[0-9]+)')
        match = re.search(pattern, old_name)
        gene = match.group('gene')
        suffix = match.group('suffix')
        new_name = '{}{}'.format(gene, suffix.upper())
        gene_type = 'normal'
    return new_name, gene_type
