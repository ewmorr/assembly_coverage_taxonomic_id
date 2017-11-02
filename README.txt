Workflow for getting assembly coverage and assigning taxonomic information to abyss assembly or spades assembly
Final collated output can be used as input to R to create "blob plots" (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3843372/) and other plots to visually examine assembly quality.
Coverage file from samtools is useful for other downstream applications (e.g. estimating coverage of predicted mRNAs)


get lengths of contigs or scaffolds
```

get_seq_lens.pl contigs.fa > genome.txt

```
index with bwa
```

bwa index contigs.fa

```
map reads - bwa
```

bwa mem -M -t 24 contigs.fa reads1.fa reads2.fa > raw_mapped.sam

```
extract mapped reads and get coverage - samtools
```

samtools view -@ 24 -Sb -F 4 -o mapped.bam raw_mapped.sam 
samtools sort -@ 24 -T sorted_mapped -o sorted_mapped.bam mapped.bam
samtools index sorted_mapped.bam
samtools depth -a sorted_mapped.bam > coverage.txt

```
Calculate average coverage per contig (or scaffold)
```

coverage_per_contig_samtools.pl coverage.txt > coverage_by_sequence.txt

```
blast against NCBI for taxonmic ID
nt should be path to nt db
```

blastn -task megablast -query contigs.fa -db nt -outfmt '6 qseqid staxids bitscore std stitle' -culling_limit 5 -num_threads 24 -evalue 1e-25 -out contigs.blast.txt

```
Parse blast for taxonomic info using NCBI taxdump
taxdump databases at https://www.ncbi.nlm.nih.gov/guide/taxonomy/
```

taxdump_blast_parse.pl contigs.blast.txt taxdump/nodes.dmp taxdump/names.dmp scaffolds.blast_parsed.txt

```
Collate taxonomy, coverage, length, GC content for input to R
If is spades assembly instead of abyss assembly remove --abyss flag at end of command
This flag is required to process sequence IDs of abyss scaffolds which are not processed properly by bwa/samtools workflow (i.e. part of sequence header is removed)
```

combine_taxonomy_coverage_len_GC.pl contigs.blast_parsed.txt coverage_by_sequence.txt contigs.fa contigs_tax_cov_len_gc.txt --abyss

```
remove mapping files if storage is limited
```

rm raw_mapped.sam
rm mapped.bam

```

See R script for plotting

assmb_visualization.r

