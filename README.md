# useful-unix-commands
## Description 
This repository serves as a database of useful command-line statemnts with a Bioinformatics approach. This includes vanilla Unix commands and commands from popular publicly available Bioinformatics tools. Most of these commands will be useful solving fairly basic questions about various topics such as data-wrangling and omics anaylysis.

## Vanilla unix (with AWK/SED)
- Reading a gzipped file: `zcat SRR1039508_1.fastq.gz`
- Returns all sequences that both start and end with C: `grep "^C" data.fa | grep "C$"`
- Counts the number of genes on the reverse strand from a BED file: `cut -f6 grcz10_refseq.bed | grep -c "-"`
- Counts the number of mRNA grcz10_refseq.bed that are not on chr14: `grep -wv "^chr14" grcz10_refseq.bed | cut -f4 | grep -c "^NM_"`
- Replacing matching entries in one column of a file by its matching second column from a different file (matches column 1 in input file to column 1 in second input file. Then replaces the name of the **first** column in the first file with the **second** column in the second file): `awk 'NR==FNR{a[$1]=$2; next}{$1=a[$1]; print}' file2 file1`
- Gets the unique count of each entry in a specified column: `awk -F ',' '{print $7}' text_file | sort | uniq -c`
- Removes blank lines from file, excludes spaces: `sed -i '/^$/d' file.txt`

## BEDTools
- Deduces how many genes overlap a different gene in the same BED file: `bedtools intersect -wa -a genes.bed -b genes.bed | wc -l`
- Produces a BED12 file of genes that do not overlap any repeat elements: `bedtools intersect -v -a genes.bed -b repeats.bed`
- Finds genes in genes.bed that overlap antisense to another gene (by any amount): `bedtools intersect -wa -wb -S -a genes.bed -b genes.bed`
- Create a BED file of promoter regions (-1000/+100 nt window around the transcription start site (TSS) of a gene). We would want to include some parts downstream of the stat coordinate as transcription factors will often bind there: ` bedtools flank -s -l 1 -r 0 -i genes.bed -g chr_sizes.txt | bedtools slop -s -l 1000 -r 100 -i - -g chr_sizes.txt > promoters.bed`
- Produces a BED6 file called representing a +/- 1kb window around every transcription termination site for genes in my_gene.bed: `bedtools flank -s -l 0 -r 1 -i my_gene.bed -g chr_sizes.txt | bedtools slop -s -l 1000 -r 1000 -i - -g chr_sizes.txt | cut -f -6 > output.bed`
- Makes a FASTA file called q6 that contains the sequences of all Kolobok repeat elements in my_repeat.bed: 
  -`sort -k1,1 -k2,2n my_repeat.bed | grep "Kolobok" > output1.bed`
  -`bedtools getfasta -name -s -fi genome.fa -bed output1.bed -fo > output1.fa`
- Makes a BED4 file (i.e., chr, start, end, geneID) of all introns for the genes in my_gene.bed: `bed12ToBed6 -i my_gene.bed | bedtools subtract -a my_gene.bed -b - | cut -f -4 > output.bed`

**Will be continually updated to include a wider range of bioinformatics tools**
