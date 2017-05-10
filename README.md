# gb2fasta.py

## Name format

Sequence id follows this format:

>  1.gene|2.order|3.family|4.genus|5.species|6.accession_id|7.specimen_voucher

Covert genbank format to fasta format with reformated id.

## Usage

Put this program in same folder with genbank file, double click to run. Then input file name.
- Put this program in same folder with gb file you want to convert.
- Double click to run. 
-  Input file name.
-  Find output files in "file_name_out"

# fasta_rename.py

## Name format

Sequence id follows this format:

>  1.gene|2.order|3.family|4.genus|5.species|6.accession_id|7.specimen_voucher

## Usage

- Put this program in same folder with fasta files you want to change id. 
- Double click to run. 
-  Input number of fields you want in the order you want. Notice that every field could only be used once.
-  Find output files in "renamed" folder.

# fasta_rename_advance.py

## Name format

*This program can be used in any name format.*

## Usage

- Put this program in same folder with fasta files you want to change id. 
- Double click to run. 
-  Input number of fields you want in the order you want and seperators you want to use. Notice that every
field could only be used once.
-  Find output files in "renamed" folder.

# Requirement

1. [python3](https://www.python.org/downloads/)

Be sure to install python3 rather than python 2.7. Besides, to use subprocess.run, you would better install python **3.5** or above

2. [biopython](http://biopython.org/wiki/Download)
 
Be sure to download same version of python3
