#!/usr/bin/python3

from Bio import Entrez


EMAIL = 'wpwupingwp@outlook.com'
Entrez.email = EMAIL
# set 'noexp' to fetch only this level
query = input('Query string:\n')
print(query)
handle = Entrez.read(Entrez.esearch(db='nuccore', term=query,
                                    usehistory='y',))
count = handle['Count']
count = int(count)
print(count)
Retstart = 1
output_file = open('sequence.gb', 'w')
while Retstart <= count:
    print(Retstart)
    try:
        genome_content = Entrez.efetch(db='nuccore',
                                       webenv=handle['WebEnv'],
                                       query_key=handle['QueryKey'],
                                       rettype='gb',
                                       retmode='text',
                                       retstart=Retstart,
                                       retmax=2000)
        output_file.write(genome_content.read())
    except:
        continue
    Retstart = Retstart + 2000
