# Basseto

**Bas**ic **se**quence **to**ol is a tool which allows multiple basic sequence manipulations.

## Requirements

 - cmake 
 - rust

## Convert

This is the most crucial aspect of the tool and contains 2 key aspects:

 - conversion in-between formats
 - subsetting of entries
 - combination of both

The allowed formats are 

 1. uBAM  = unaligned BAM files, currently PacBio targeted but might work with others
 2. fastq = sequences with quality information added
 3. fasta = sequences with no quality added

Both, fasta and fastq can take as well compressed input and based on valid gz extensions extract on the fly.
Similarly it will automatically compress output if a gz extension is provided in the output file name.

An important point of this tool is that input fasta/q is verified to contain valid header,
sequences and qualities. Subsetting can be chosen either via a black- or a white-list.
Upon execution one has to define the input and output format, which can be the same format as well, e.g.
in the case of subsetting.

Evidently conversion between some format loose information, e.g. fastq --> fa or uBAM --> fq.
The other way around one has to provide dummy quality values if information is absent.
E.g. fa --> uBAM or fasta --> fastq



## Stats

This recent addition to the tool provides the possibility to get statistics of sequencing data.
In the most simple fashion, it will generate a small report of observed statstics such as NG50, median, max,...
With the exhaustive flag it will provide in a separate file statistics for every single sequence in a tab-separated form.

Example:

```bash
# program:	basetto
# version:	0.3.0
# author:	Emanuel Schmid-Siegert <emanuel.schmid-siegert@selexis.com>
# command:	target/release/basseto stats --input ../FusioneR/test/GRCh37.p13.genome_CHR7.fa --input-format fasta --output test
gc:	0.407
min:	159138663
max:	159138663
mean:	159138663
median:	159138663
total_bp:	159138663
total_seq:	1
total_nbp:	3785000
nx:	159138663,159138663,159138663,159138663,159138663,159138663,159138663,159138663,159138663,159138663,159138663
lx:	1,1,1,1,1,1,1,1,1,1,1
ngx:	NA
lgs:	NA


# program:	basetto
# version:	0.3.0
# author:	Emanuel Schmid-Siegert <emanuel.schmid-siegert@selexis.com>
# command:	target/release/basseto stats --input test.fa --input-format fasta --output test --exhaustive
# fields:	ID	Lengths	GC	Qual	Gap
10A	10	  0	0	0
10G	10	  1	0	0
10T10C10A10G	40	0.5	0	0
5each	20	0.5	0	0
```

N(x),L(x),NG(x) and LG(x) values are computed for 0.0-1.0 in 0.1 steps in a comma-separated list.
Providing a genome size in Mbp (e.g. human = 3100) allows to calculate NG(x) and LG(x) values, otherwise NA will be returned.

## Exhaustive mode

In exhaustive mode, there is an additional file produced which contains for each sequence a 

- ID
- Lengths
- GC content
- mean Quality
- number of gaps


### PacBio specifics for exhaustive mode

The exchaustive mode returns the quality of each read and it's name, which is very handy for qc aspects.

**Important:** Please be aware that in case of uBAM, qualities per sequence are first inquired via the rq-tag.
This can be forced to be calculated instead from nucleotide base qualities with the parameter `-r`.

Pacbio is actually cheating and caps the rq  at 1 which is not correct. As this returns otherwise just  the max of u32 I will instead cap to phred 60 for the moment,which is  99.9999% quality. If this is absent, qualities of the entire sequence will be read and a median value reported.

They describe [here](https**://ccs.how/how-does-ccs-work#9-qv-calculation) here the generation of the scores.

> QV calculation
> The log-likelihood ratio between the most likely template sequence and all of its mutated counterparts is used to calculate a quality for each base in the final 
> consensus. The average of the per-base qualities is the read accuracy, rq.

## Comparison with other tools 

### Subsetting fastq with seqtk 

[Seqtk](https://github.com/lh3/seqtk) is from Heng Li

Extracting reads from a data-set.

```bash
 time  seqtk subseq SLX1021.1P.fastq SLX1021_vectorReads.txt > seqtk.fq

real	1m35.674s
user	0m0.090s
sys	0m0.184s

time ./basseto convert --input reads.fastq --input-format fastq \
  --output test.fq --output-format fastq --whitelist Reads.txt 
INFO: Converting formats
INFO: Input format : Fastq , Output format: Fastq
INFO: ID list Reads.txt provided, reading entries...

real    1m16.341s
user    0m59.828s
sys     0m16.368s

diff test.fq seqtk.fq
```

When working on a `.fastq` file this tool is a tiny bit faster than seqtk.

### Subsetting with quality uBAM against bamtools

This is a lengthy process

```bash
 time bamtools filter -in HiFiWGS.ccs.bam \
  -tag "rq:>=0.99" -forceCompression -out  bamtools_out.bam

real    47m44.957s
user    46m59.177s
sys     0m38.162s
```

This on the other hand is blazing fast with our tool:

```bash
time ./basseto convert --input HiFiWGS.ccs.bam \
  --input-format uBAM --output basseto.ccs.bam  --output-format uBAM --ubam-quality 0.99 --threads 40
INFO: Converting formats
INFO: Input format : Ubam , Output format: Ubam

real    1m15.546s
user    50m15.573s
sys     0m59.009s
```

So single threaded that would be kind of similar but in constrast to bamtools we can go threaded here and the order remains similar.

The same but writing HIFI and LQ the same time

First lets establish number of entries:

```bash
samtools view -h  CCS_WGS.ccs.bam |  grep -oh "rq:f:[0-9]*\.[0-9]*" | awk '{FS=":"}{if($3>=0.99) HIFI+=1; if($3<0.99) LQ+=1}END{print "HIFI",HIFI,"LQ",LQ}'

HIFI 1479246 LQ 350137
```

```bash
#!/bin/bash -ue
bamtools filter -in CCS_WGS.ccs.bam -tag "rq:>=0.99"  -forceCompression -out HIFI.bam &
bamtools filter -in CCS_WGS.ccs.bam -tag "rq:<0.99"   -forceCompression -out LQ.bam &
```

```bash
#!/bin/bash -ue
basseto convert     --input CCS_WGS.bam --input-format uBAM     --output PA3029-HIFI.bam --output-format uBAM     --unfiltered PA3029-LQ.bam     --ubam-quality 0.99     --threads 40
```

This results for both of them in 1'492'758 HIFI entries with only 6 threads being used for basseto (guess due to IO throtteling) but still a difference of 45 minutes compared to 12 minutes.
We get as well similarly for both 1'917'017 in LQ results.

### Subsetting uBAM with read IDs

```bash
time ./basseto convert --input LQ.bam \
  --input-format uBAM --output basseto_out.bam --output-format uBAM --whitelist IDs.txt
INFO: Converting formats
INFO: Input format : Ubam , Output format: Ubam
INFO: ID list IDs.txt provided, reading entries...

real	7m59.816s
user	7m52.550s
sys	0m6.553s

```


With multi-threading:

```bash


time ./basseto convert --input LQ.bam \
  --input-format uBAM --output basseto_out.bam --output-format uBAM --whitelist IDs.txt -t 40
INFO: Converting formats
INFO: Input format : Ubam , Output format: Ubam
INFO: ID list PA3029-R4P9_61b7643b5b6e7a1879754762_GraftIDs.txt provided, reading entries...

real	0m22.569s
user	12m57.142s
sys	0m15.366s

```
