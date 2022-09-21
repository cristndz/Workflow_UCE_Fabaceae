## genome location
/scratch_local/JMA/PROYECTO_CLOROPLASTO_CRISTIAN/FABACEAE-UCE/genomes

## put each genome in its own directory
for critter in *; do mkdir ${critter%.*}; mv $critter ${critter%.*}; done

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

# Combine paired reads and compress

for critter in ariBati ariCarde ariDio ariDuran ariHelo ariHoeh ariHypo ariHypo2 
ariHypo3 ariHypo4 ariHypo5 ariHypo6 ariHypo7 ariIpae ariMonti ariParagu ariStenos ariVillo ...;
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

/source/stampy/stampy 1.0../stampy.py --species="arachis_hypogea" --assembly="Inga_leiocalycina_KT428296" -G Inga_leiocalycina_KT428296 Inga_leiocalycina_KT428296.fasta
/source/stampy/stampy 1.0../stampy.py -g Inga_leiocalycina_KT428296 -H Inga_leiocalycina_KT428296

# Align simulated reads to the base genome

mkdir alignments

export cores=4
export base=Inga_leiocalycina_KT428296
export base_dir=/scratch_local/JMA/PROYECTO_CLOROPLASTO_CRISTIAN/FABACEAE-UCE/alignments
for critter in araBati araCarde araDio araDuran araHelo araHoeh araHypo2 araHypo3 araHypo4 araHypo5 araHypo6 araHypo7 araIpae araMonti araParagu araStenos araVillo ...;
    do
        export reads=$critter-pe100-reads.fq.gz;
        mkdir -p $base_dir/$critter;
        cd $base_dir/$critter;
        python2 /source/stampy/stampy-1.0.32/stampy.py --maxbasequal 93 -g ../../base/$base -h ../../base/$base \
        --substitutionrate=0.05 -t$cores --insertsize=400 -M \
        ../../reads/$reads | samtools view -Sb - > $critter-to-$base.bam;
    done;

# Remove unmapped reads from files BAM
cd alignments
mkdir all

for critter in araBati araCarde araDio araDuran araHelo araHoeh araHypo2 araHypo3 araHypo4 araHypo5 araHypo6 araHypo7 araIpae araMonti araParagu araStenos araVillo ...;
    do
        samtools view -h -F 4 -b $critter/$critter-to-araHypo.bam > $critter/$critter-to-araHypo-MAPPING.bam;
        rm $critter/$critter-to-araHypo.bam;
        ln -s ../$critter/$critter-to-Inga_leiocalycina_KT428296-MAPPING.bam all/$critter-to-Inga_leiocalycina_KT428296-MAPPING.bam;
    done;

------------------------------------------------------- PHYLUCE ----------------------------------------------- 

conda activate phyluce-1.7.1 
mkdir bed # in main directory
cd bed    

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


------------------------------------------------------- PHYLUCE ----------------------------------------------- 
conda activate phyluce-1.7.1 


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
