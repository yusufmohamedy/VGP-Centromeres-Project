# required inputs for the pipeline
you will need the following inputs for the pipeline to run: chromosome_lengths.tsv, scaffold_lengths.tsv, consensus.bed

i quickly go through what you need to do to obtain them



## chromosome and scaffold lengths
so you need a file from chromosome 

needs minimum chr lengths file but if species has any scaffold its highly recommed you also make a scaffold file

both chr length and scaf lenghts should be tsv like below

can be extracted from fa.fai you might need to make one using samtools faidx
then you can separate scaffolds from chromosomes the majority of the time by looking at the first two characters of the prefix as they should be differnt for chr and scaffolds
(be careful though as some scaffold only genomes have MT prefix that is different to scaffolds)


```bash
NC_069499.1	1601664463
NC_069500.1	1450603269
NC_069501.1	1296922325
NC_069502.1	1241340601
NC_069503.1	1167632593
NC_069504.1	1160268632
NC_069505.1	643478974
NC_069506.1	354776354
NC_069507.1	291274852
NC_069508.1	220236004
```



## how to make the consensus bed file

made from running fastan, tancons, satmatch and nhmmer on a fasta files

### 1ano file


```bash
FAtoGDB Unicorn.fa

FasTAN -mp Unicorn.1gdb
```



### need to obtain a bed file made from consensus sequence of satellites
use tancons, satmatch from https://github.com/richarddurbin/alntools
use nhmmer https://github.com/EddyRivasLab/hmmer

```bash
tancons -o Unicorn-tancons.fa -s Unicorn.fa.gz Unicorn.1ano

satmatch -fa Unicorn-sat_clusters.fa -loop Unicorn-tancons.fa

nhmmer --cpu=8 --dna --qformat fasta  --tblout Unicorn-sat_positions.tbl -o Unicorn-nhmmer.out Unicorn-sat_clusters.fa Unicorn.fa

```



### to make a bed file using nhmmer output (Unicorn-sat_positions.tbl)
```bash

awk -v OFS='\t' '!/^#/ && $1 != "" {
    start = ($7 < $8) ? $7 : $8
    end   = ($7 < $8) ? $8 : $7
    print $1, start - 1, end, $3
}' Unicorn-sat_positions.tbl \
| sort -k1,1V -k4,4 -k2,2n \
| awk -v OFS='\t' '
{
    split($4, a, "-")
    unit_size = a[length(a)] + 0
    if (unit_size <= 0) unit_size = 100

    if (NR == 1) {
        chr = $1; start = $2; end = $3; sat = $4; max_gap = unit_size
        next
    }

    # Merge if same chromosome, same satellite, and gap <= unit_size
    if ($1 == chr && $4 == sat && ($2 - end) <= max_gap) {
        if ($3 > end) end = $3
    } else {
        print chr, start, end, sat
        chr = $1; start = $2; end = $3; sat = $4; max_gap = unit_size
    }
}
END {
    if (NR > 0) print chr, start, end, sat
}' \
| sort -k1,1V -k2,2n \
> Unicorn-consensus_sat.bed
```


## i have added some example slurm scripts that should be helpful for you to run this on the HPC cluster


the scripts are: example_slurm_array_task_list.tsv, example_consensus_run_slurm.sh, example_centromeric_pipeline_slurm.sh

