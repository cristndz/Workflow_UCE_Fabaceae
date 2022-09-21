##convert genomes to 2bit
faToTwoBit name_file.fasta name_file.2bit

# simulate reads from genomes using ART
> art_illumina \
    --paired \
    --in ../genomes/ariBati/ariBati.fasta \
    --out ariBati-pe100-reads \
    --len 100 \
    --fcov 2 \
    --mflen 200 \
    --sdev 150 \
    -ir 0.0 -ir2 0.0 -dr 0.0 -dr2 0.0 -qs 100 -qs2 100 -na

--La última línea apaga la simulación de error de secuenciación en cada una 
  de estas lecturas y también desactiva la creación de un archivo de alineación
  que muestra de dónde provienen las lecturas en la secuencia del genoma.--

# Combine paired reads and compress

for critter in ariBati ariCarde ariDio ariDuran ariHelo ariHoeh ariHypo ariHypo2 
ariHypo3 ariHypo4 ariHypo5 ariHypo6 ariHypo7 ariIpae ariMonti ariParagu ariStenos ariVillo;
> do
> echo "working on $critter";
> touch $critter-pe100-reads.fq;
> cat $critter-pe100-reads1.fq > $critter-pe100-reads.fq;
> cat $critter-pe100-reads2.fq >> $critter-pe100-reads.fq;
> rm $critter-pe100-reads1.fq;
> rm $critter-pe100-reads2.fq;
> gzip $critter-pe100-reads.fq;
> done;

# Prep base genome (Stampy)

/source/stampy/stampy 1.0../stampy.py --species="arachis_hypogea" --assembly="araHypo" -G araHypo araHypo.fasta
/source/stampy/stampy 1.0../stampy.py -g araHypo -H araHypo

# Align simulated reads to the base genome

export cores=4
export base=araHypo
export base_dir=/scratch_local/JMA/PROYECTO_CLOROPLASTO_CRISTIAN/uce-arachis/alignments
for critter in araBati araCarde araDio araDuran araHelo araHoeh araHypo2 araHypo3 araHypo4 araHypo5 araHypo6 araHypo7 araIpae araMonti araParagu araStenos araVillo;
    do
        export reads=$critter-pe100-reads.fq.gz;
        mkdir -p $base_dir/$critter;
        cd $base_dir/$critter;
        python2 /source/stampy/stampy-1.0.32/stampy.py --maxbasequal 93 -g ../../base/$base -h ../../base/$base \
        --substitutionrate=0.05 -t$cores --insertsize=400 -M \
        ../../reads/$reads | samtools view -Sb - > $critter-to-$base.bam;
    done;

# Remove unmapped reads from files BAM

for critter in araBati araCarde araDio araDuran araHelo araHoeh araHypo2 araHypo3 araHypo4 araHypo5 araHypo6 araHypo7 araIpae araMonti araParagu araStenos araVillo;
    do
        samtools view -h -F 4 -b $critter/$critter-to-araHypo.bam > $critter/$critter-to-araHypo-MAPPING.bam;
        rm $critter/$critter-to-araHypo.bam;
        ln -s ../$critter/$critter-to-araHypo-MAPPING.bam all/$critter-to-araHypo-MAPPING.bam;
    done;

# Now, convert our *-MAPPING.bam files to BED format
for i in ../alignments/all/*.bam; do echo $i; bedtools bamtobed -i $i -bed12 > `basename $i`.bed; done

# Sort the bed files
for i in *.bed; do echo $i; bedtools sort -i $i > ${i%.*}.sort.bed; done

# Merge overlapping or nearly overlapping intervals
for i in *.bam.sort.bed; do echo $i; bedtools merge -i $i > ${i%.*}.merge.bed; done

# To get some idea of the total number of merged, putatively conserved regions in each exemplar taxon
for i in *.bam.sort.merge.bed; do wc -l $i; done

--------------------output------------------------
415 araBati-to-araHypo-MAPPING.bam.sort.merge.bed
412 araCarde-to-araHypo-MAPPING.bam.sort.merge.bed
410 araDio-to-araHypo-MAPPING.bam.sort.merge.bed
416 araDuran-to-araHypo-MAPPING.bam.sort.merge.bed
405 araHelo-to-araHypo-MAPPING.bam.sort.merge.bed
416 araHoeh-to-araHypo-MAPPING.bam.sort.merge.bed
422 araHypo2-to-araHypo-MAPPING.bam.sort.merge.bed
405 araHypo3-to-araHypo-MAPPING.bam.sort.merge.bed
421 araHypo4-to-araHypo-MAPPING.bam.sort.merge.bed
395 araHypo5-to-araHypo-MAPPING.bam.sort.merge.bed
399 araHypo6-to-araHypo-MAPPING.bam.sort.merge.bed
402 araHypo7-to-araHypo-MAPPING.bam.sort.merge.bed
441 araIpae-to-araHypo-MAPPING.bam.sort.merge.bed
404 araMonti-to-araHypo-MAPPING.bam.sort.merge.bed
410 araParagu-to-araHypo-MAPPING.bam.sort.merge.bed
398 araStenos-to-araHypo-MAPPING.bam.sort.merge.bed
413 araVillo-to-araHypo-MAPPING.bam.sort.merge.bed

# Remove repetitive intervals

 for i in *.sort.merge.bed;
> do
> phyluce_probe_strip_masked_loci_from_set \
> --bed $i \
> --twobit ../genomes/araHypo/araHypo.2bit \
> --output ${i%.*}.strip.bed \
> --filter-mask 0.25 \
> --min-length 80
> done;

--------------------------------output-----------------------------------
Screened 415 sequences from araBati-to-araHypo-MAPPING.bam.sort.merge.bed.  Filtered 0 with > 25.0% masked bases or > 0 N-bases or < 80 length. Kept 415.
Screened 412 sequences from araCarde-to-araHypo-MAPPING.bam.sort.merge.bed.  Filtered 0 with > 25.0% masked bases or > 0 N-bases or < 80 length. Kept 412.
Screened 410 sequences from araDio-to-araHypo-MAPPING.bam.sort.merge.bed.  Filtered 0 with > 25.0% masked bases or > 0 N-bases or < 80 length. Kept 410.
Screened 416 sequences from araDuran-to-araHypo-MAPPING.bam.sort.merge.bed.  Filtered 0 with > 25.0% masked bases or > 0 N-bases or < 80 length. Kept 416.
Screened 405 sequences from araHelo-to-araHypo-MAPPING.bam.sort.merge.bed.  Filtered 0 with > 25.0% masked bases or > 0 N-bases or < 80 length. Kept 405.
Screened 416 sequences from araHoeh-to-araHypo-MAPPING.bam.sort.merge.bed.  Filtered 0 with > 25.0% masked bases or > 0 N-bases or < 80 length. Kept 416.
Screened 422 sequences from araHypo2-to-araHypo-MAPPING.bam.sort.merge.bed.  Filtered 0 with > 25.0% masked bases or > 0 N-bases or < 80 length. Kept 422.
Screened 405 sequences from araHypo3-to-araHypo-MAPPING.bam.sort.merge.bed.  Filtered 0 with > 25.0% masked bases or > 0 N-bases or < 80 length. Kept 405.
Screened 421 sequences from araHypo4-to-araHypo-MAPPING.bam.sort.merge.bed.  Filtered 0 with > 25.0% masked bases or > 0 N-bases or < 80 length. Kept 421.
Screened 395 sequences from araHypo5-to-araHypo-MAPPING.bam.sort.merge.bed.  Filtered 0 with > 25.0% masked bases or > 0 N-bases or < 80 length. Kept 395.
Screened 399 sequences from araHypo6-to-araHypo-MAPPING.bam.sort.merge.bed.  Filtered 0 with > 25.0% masked bases or > 0 N-bases or < 80 length. Kept 399.
Screened 402 sequences from araHypo7-to-araHypo-MAPPING.bam.sort.merge.bed.  Filtered 0 with > 25.0% masked bases or > 0 N-bases or < 80 length. Kept 402.
Screened 441 sequences from araIpae-to-araHypo-MAPPING.bam.sort.merge.bed.  Filtered 0 with > 25.0% masked bases or > 0 N-bases or < 80 length. Kept 441.
Screened 404 sequences from araMonti-to-araHypo-MAPPING.bam.sort.merge.bed.  Filtered 0 with > 25.0% masked bases or > 0 N-bases or < 80 length. Kept 404.
Screened 410 sequences from araParagu-to-araHypo-MAPPING.bam.sort.merge.bed.  Filtered 0 with > 25.0% masked bases or > 0 N-bases or < 80 length. Kept 410.
Screened 398 sequences from araStenos-to-araHypo-MAPPING.bam.sort.merge.bed.  Filtered 0 with > 25.0% masked bases or > 0 N-bases or < 80 length. Kept 398.
Screened 413 sequences from araVillo-to-araHypo-MAPPING.bam.sort.merge.bed.  Filtered 0 with > 25.0% masked bases or > 0 N-bases or < 80 length. Kept 413.

# Determining locus presence in multiple genomes 
----------create file configuration--------------

phyluce_probe_get_multi_merge_table \
> --conf bed-files.conf \
> --base-taxon araHypo \
> --output arachis-to-araHypo.sqlite

sqlite> select * from araHypo limit 20;
uce  chromo      start  stop  arabati  aracarde  aradio  araduran  arahelo  arahoeh  arahypo2  arahypo3  arahypo4  arahypo5  arahypo6  arahypo7  araipae  aramonti  araparagu  arastenos  aravillo
---  ----------  -----  ----  -------  --------  ------  --------  -------  -------  --------  --------  --------  --------  --------  --------  -------  --------  ---------  ---------  --------
1    MG814009.1  116    216   1        1         1       0         1        1        1         1         1         1         1         1         1        1         1          1          0
2    MG814009.1  358    1179  1        1         1       1         1        1        1         1         1         1         1         1         1        1         1          1          1
3    MG814009.1  1214   1345  1        1         1       1         1        1        1         1         1         1         1         1         1        1         1          1          1
4    MG814009.1  1388   1488  1        1         1       1         0        1        1         1         1         1         1         1         1        1         1          1          1
5    MG814009.1  1537   1715  1        1         1       1         0        1        1         1         1         1         1         1         1        1         1          1          0
6    MG814009.1  1851   2254  1        1         1       1         1        1        1         1         1         1         1         1         1        1         1          1          1
7    MG814009.1  2348   2983  1        1         1       1         1        1        1         1         1         1         1         1         1        1         1          1          1
8    MG814009.1  3049   3243  1        1         1       1         1        1        1         1         1         1         1         1         1        1         1          1          1
9    MG814009.1  3365   3518  1        1         1       1         1        1        1         1         1         0         1         1         1        1         1          1          1
10   MG814009.1  3522   4483  1        1         1       1         1        1        1         1         1         1         1         1         1        1         1          1          1

# Determining shared, conserved, loci

phyluce_probe_query_multi_merge_table \
> --db arachis-to-araHypo.sqlite \
> --base-taxon araHypo

Loci shared by araHypo + 0 taxa:        446.0
Loci shared by araHypo + 1 taxa:        446.0
Loci shared by araHypo + 2 taxa:        446.0
Loci shared by araHypo + 3 taxa:        446.0
Loci shared by araHypo + 4 taxa:        445.0
Loci shared by araHypo + 5 taxa:        444.0
Loci shared by araHypo + 6 taxa:        443.0
Loci shared by araHypo + 7 taxa:        440.0
Loci shared by araHypo + 8 taxa:        439.0
Loci shared by araHypo + 9 taxa:        436.0
Loci shared by araHypo + 10 taxa:       434.0
Loci shared by araHypo + 11 taxa:       433.0
Loci shared by araHypo + 12 taxa:       431.0
Loci shared by araHypo + 13 taxa:       427.0
Loci shared by araHypo + 14 taxa:       426.0
Loci shared by araHypo + 15 taxa:       419.0
Loci shared by araHypo + 16 taxa:       405.0
Loci shared by araHypo + 17 taxa:       350.0

# Extract those loci that match n taxa

phyluce_probe_query_multi_merge_table \
> --db arachis-to-araHypo.sqlite \
> --base-taxon araHypo \
> --output araHypo+17.bed \
> --specific-counts 17
Counter({'arahoeh': 350, 'arahypo2': 350, 'arabati': 350, 'araparagu': 350, 'arahypo4': 350, 'arastenos': 350, 'arahypo6': 350, 'arahypo7': 350, 'aradio': 350, 'araduran': 350, 'arahypo3': 350, 'arahypo5': 350, 'aracarde': 350, 'araipae': 350, 'arahelo': 350, 'aravillo': 350, 'aramonti': 350})

##################_CONSERVED_LOCUS_VALIDATION_############################

#Extract FASTA sequence from base genome for temp bait design

phyluce_probe_get_genome_sequences_from_bed \
> --bed araHypo+17.bed \
> --twobit ../genomes/araHypo/araHypo.2bit \
> --buffer-to 160 \
> --output araHypo+17.fasta
Screened 350 sequences.  Filtered 0 < 160 bp or with > 25.0% masked bases or > 0 N-bases. Kept 350.

------------------------output------------------------------
>slice_0 |MG814009.1:358-1179
CAGCTAGGTCTAGAGGGAAGTTATGAGCATTACGTTCATGCATAACTTCCATACCAAGAT
TAGCTCGGTTAATAATATCAGCCCAGGTGTTAATTACACGACCTTGACTATCAACTACAG
ATTGGTTGAAATTGAAACCATTTAAATTGAAAGCCATAGTGCTAATACCTAACGCGGTAA
ACCAGATACCTACTACAGGCCAAGCAGCTAGGAAGAAATGTAAAGAACGAGAATTGTTGA
AACTAGCATATTGGAAGATCAATCGTCCAAAATAACCATGAGCAGCCACAATATTATAGG
TTTCTTCCTCTTGACCGAATCTGTAACCTTCATTAGCAGATTCATTTTCTGTGGTTTCCC
TAATCAAACTAGAAGTTACCAAGGAACCGTGCATTGCACTGAATAGGGAGCCGCCGAATA
CACCAGCTACGCCTAACATATGAAATGGATGCATAAGAATATTGTGTTCAGCCTGAAATA
CAATCATAAAATTGAAGGTACCAGAAATTCCTAGAGGCATACCATCCGAAAAGCTTCCTT
GACCAATTGGATAGATCAAGAAAACAGCAGTAGCCGCTGCAACAGGAGCTGAATATGCAA
CAGCAATCCAAGGGCGCATACCCAGACGAAAACTCAGTTCCCACTCACGACCCATGTAAC
AAGCTACACCAAGTAAGAAGTGTAGAACAATTAGTTCATAAGGACCGCCATTGTATAACC
ATTCATCAACAGATGCCGCTTCCCATATCGGGTAAAAGTGCAAACCTATAGCCGCCGAAG
TAGGAATAATGGCACCCGAAATAATATTGTTTCCATAAAGT

# Design a temporary bait set from the base taxon

phyluce_probe_get_tiled_probes \
> --input araHypo+17.fasta \
> --probe-prefix "uce-" \
> --design arachis-V1 \
> --designer c-diaz \
> --tiling-density 3 \
> --two-probes \
> --overlap middle \
> --masking 0.25 \
> --remove-gc \
> --output araHypo+17.temp.probes
Probes removed for masking (.) / low GC % (G) / ambiguous bases (N):
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG

Conserved locus count = 288
Probe Count = 539

# Remove duplicates from our temporary bait set

phyluce_probe_easy_lastz \
> --target araHypo+17.temp.probes \
> --query araHypo+17.temp.probes \
> --identity 50 \
> --coverage 50 \
> --output araHypo+17.temp.probes-TO-SELF-PROBES.lastz
Started:  Sat Nov 14, 2020  17:53:01
Ended:  Sat Nov 14, 2020  17:53:02
Time for execution:  0.00261706908544 minutes

phyluce_probe_remove_duplicate_hits_from_probes_using_lastz \
> --fasta araHypo+17.temp.probes \
> --lastz araHypo+17.temp.probes-TO-SELF-PROBES.lastz \
> --probe-prefix=uce-
Parsing lastz file...
Screening results...
Screened 538 fasta sequences.  Filtered 60 duplicates. Kept 419.

# Align baits against exemplar genomes

phyluce_probe_run_multiple_lastzs_sqlite 
--probefile ../bed/araHypo+17.temp-DUPE-SCREENED.probes 
--scaffoldlist araBati araCarde araDio araDuran araHelo araHoeh araHypo2 araHypo3 araHypo4 araHypo5 araHypo6 araHypo7 araIpae araMonti araParagu araStenos araVillo bauBlak 
--genome-base-path ../genomes 
--identity 50 
--cores 2 
--db araHypo+17+bauBlak.sqlite 
--output arachis-genome-lastz

#Extract sequence around conserved loci from exemplar genomes

 phyluce_probe_slice_sequence_from_genomes \
> --conf arachis-genome.conf \
> --lastz arachis-genome-lastz \
> --probes 180 \
> --name-pattern "araHypo+17.temp-DUPE-SCREENED.probes_v_{}.lastz.clean" \
> --output arachis-genome-fasta
2020-11-15 17:36:27,591 - Phyluce - INFO - =================== Starting Phyluce: Slice Sequence ===================
2020-11-15 17:36:27,592 - Phyluce - INFO - ------------------- Working on bauBlak genome -------------------
2020-11-15 17:36:27,592 - Phyluce - INFO - Reading bauBlak genome
2020-11-15 17:36:28,674 - Phyluce - INFO - bauBlak: 209 uces, 0 dupes, 209 non-dupes, 62 orient drop, 2 length drop, 145 written
2020-11-15 17:36:28,675 - Phyluce - INFO - ------------------- Working on araBati genome -------------------
2020-11-15 17:36:28,675 - Phyluce - INFO - Reading araBati genome
2020-11-15 17:36:28,973 - Phyluce - INFO - araBati: 228 uces, 0 dupes, 228 non-dupes, 64 orient drop, 2 length drop, 162 written
2020-11-15 17:36:28,973 - Phyluce - INFO - ------------------- Working on araCarde genome ------------------
2020-11-15 17:36:28,974 - Phyluce - INFO - Reading araCarde genome
2020-11-15 17:36:29,277 - Phyluce - INFO - araCarde: 227 uces, 0 dupes, 227 non-dupes, 64 orient drop, 2 length drop, 161 written
2020-11-15 17:36:29,277 - Phyluce - INFO - -------------------- Working on araDio genome -------------------
2020-11-15 17:36:29,279 - Phyluce - INFO - Reading araDio genome
2020-11-15 17:36:29,574 - Phyluce - INFO - araDio: 228 uces, 0 dupes, 228 non-dupes, 64 orient drop, 2 length drop, 162 written

#Find which loci we detect consistently

phyluce_probe_get_multi_fasta_table \
> --fastas arachis-genome-fasta \
> --output multifastas.sqlite \
> --base-taxon araHypo

sqlite> select * from araHypo limit 20;
uce-86|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1
uce-87|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1
uce-278|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1
uce-279|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1
uce-272|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1
uce-273|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1
uce-271|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1
uce-277|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1
uce-274|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1
uce-275|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1
uce-104|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1
uce-105|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1
uce-106|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1

phyluce_probe_query_multi_fasta_table \
> --db multifastas.sqlite \
> --base-taxon araHypo
Loci shared by 0 taxa:  162.0
Loci shared by 1 taxa:  162.0
Loci shared by 2 taxa:  162.0
Loci shared by 3 taxa:  162.0
Loci shared by 4 taxa:  162.0
Loci shared by 5 taxa:  162.0
Loci shared by 6 taxa:  162.0
Loci shared by 7 taxa:  162.0
Loci shared by 8 taxa:  162.0
Loci shared by 9 taxa:  162.0
Loci shared by 10 taxa: 162.0
Loci shared by 11 taxa: 162.0
Loci shared by 12 taxa: 162.0
Loci shared by 13 taxa: 162.0
Loci shared by 14 taxa: 162.0
Loci shared by 15 taxa: 162.0
Loci shared by 16 taxa: 162.0
Loci shared by 17 taxa: 161.0

phyluce_probe_query_multi_fasta_table \
> --db multifastas.sqlite \
> --base-taxon araHypo \
> --output araHypo+17-back-to-17.conf \
> --specific-counts 17
Counter({'arahoeh': 161, 'arastenos': 161, 'arabati': 161, 'araparagu': 161, 'arahypo4': 161, 'arahypo6': 161, 'arahypo7': 161, 'aradio': 161, 'arahypo2': 161, 'arahypo3': 161, 'arahypo5': 161, 'aracarde': 161, 'araipae': 161, 'arahelo': 161, 'aravillo': 161, 'araduran': 161, 'aramonti': 161, 'baublak': 145})
Total loci = 161

#################_FINAL_BAIT_SET_DESIGN_#########################

# Design a bait set using all exemplar genomes (and the base)

phyluce_probe_get_tiled_probe_from_multiple_inputs \
> --fastas arachis-genome-fasta \
> --multi-fasta-output araHypo+17-back-to-17.conf \
> --probe-prefix "uce-" \
> --designer cdiaz \
> --design arachis-v1 \
> --tiling-density 3 \
> --overlap middle \
> --masking 0.25 \
> --remove-gc \
> --two-probes \
> --output arachis-v1-master-probe-list.fasta
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG


Conserved locus count = 161
Probe Count = 5413

# Remove duplicates from our bait set

phyluce_probe_easy_lastz \
> --target arachis-v1-master-probe-list.fasta \
> --query arachis-v1-master-probe-list.fasta \
> --identity 50 \
> --coverage 50 \
> --output arachis-v1-master-probe-list-TO-SELF-PROBES.lastz
Started:  Sun Nov 15, 2020  18:04:35
Ended:  Sun Nov 15, 2020  18:04:43
Time for execution:  0.130896615982 minutes

# Now, screen the alignements and filter our master probe list to remove duplicates

phyluce_probe_remove_duplicate_hits_from_probes_using_lastz 
--fasta arachis-v1-master-probe-list.fasta 
--lastz arachis-v1-master-probe-list-TO-SELF-PROBES.lastz 
--probe-prefix=uce-

Parsing lastz file...
Screening results...
Screened 5412 fasta sequences.  Filtered 0 duplicates. Kept 5413

The master probe list that has been filtered of putatively duplicate loci is now located in arachis-v1-master-probe-list-DUPE-SCREENED.fasta.
