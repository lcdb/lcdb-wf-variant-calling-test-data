import gzip
import os
import tempfile
from tempfile import TemporaryDirectory

shell.prefix(' set -euo pipefail ;')

VDB_CONFIG_PRELUDE = 'export VDB_CONFIG=/usr/local/apps/ncbi/config/biowulf.kfg'
#WES test data comes from a published "reference dataset" for WGS and WES https://www.nature.com/articles/s41597-021-01077-5#Sec28
WGS_accessions = {
    'normal':'SRR7890855',
    'tumor':'SRR7890854',
}


reference = 'data/GRCh38.6.20.fa.gz'

mapped_n = 2000000
unmapped_n = 1000

known = 'data/known_variation_noiupac.vcf.gz'
dbnsfp = 'data/dbnsfp_6_20.vcf.gz'
#Final results of the test data should be a subset exon bedfile, a subset fasta reference, and subset reads
rule all:
    input:
        expand('data/{sample}_R{n}.6.20.fq.gz', n=[1,2], sample=WGS_accessions.keys()),
        'data/exons_subset.bed',
        known,
        dbnsfp


rule known_variation:
    """
    Download the known variation file. Remove iupac codes and subset.
    """
    input:
        fai = 'data/GRCh38.6.20.fa.gz.fai',
        limit = 'data/LIMIT.bed'
    resources:
        mem_mb=1024 * 16,
        disk_mb=1024 * 16,
        runtime=60
    output: 'data/known_variation_noiupac.vcf.gz'
    run:
        workdir = os.getcwd()
        with tempfile.TemporaryDirectory() as tmpdir:
            shell(
                '( cd {tmpdir}; '
                'curl -O ftp://ftp.ensembl.org/pub/release-102/variation/vcf/homo_sapiens/homo_sapiens-chr6.vcf.gz -O ftp://ftp.ensembl.org/pub/release-102/variation/vcf/homo_sapiens/homo_sapiens-chr6.vcf.gz.csi -O ftp://ftp.ensembl.org/pub/release-102/variation/vcf/homo_sapiens/homo_sapiens-chr20.vcf.gz -O ftp://ftp.ensembl.org/pub/release-102/variation/vcf/homo_sapiens/homo_sapiens-chr20.vcf.gz.csi && bcftools concat -Oz --naive homo_sapiens-chr6.vcf.gz homo_sapiens-chr20.vcf.gz > concat.vcf.gz && bcftools reheader --fai {workdir}/{input.fai} concat.vcf.gz > tmp_known_variation.vcf.gz && '
                'rbt vcf-fix-iupac-alleles < tmp_known_variation.vcf.gz | bcftools view -Oz > tmp_no_iupac.vcf.gz && '
                'zgrep ^# tmp_no_iupac.vcf.gz > {workdir}/{output} && '
                'tabix -p vcf tmp_no_iupac.vcf.gz && '
                'tabix -R {workdir}/{input.limit} tmp_no_iupac.vcf.gz >> {workdir}/{output} ) '
                # This line may not work. Do it manually
                #'bcftools view {output} -Oz > {output} '
            )


rule dbnsfp:
    """
    Download example chromosomes from dbnsfp
    """
    resources:
        mem_mb=1024*200,
        disk_mb=1024*200,
        runtime=60*8
    threads: 16
    input:
        bed = 'data/LIMIT.bed'
    output:
        vcf='data/dbnsfp_6_20.vcf.gz',
        tbi='data/dbnsfp_6_20.vcf.gz.tbi'
    run:
        workdir = os.getcwd()
        with tempfile.TemporaryDirectory() as tmpdir:
            shell(
                '''(cd {tmpdir}; wget -O- https://usf.box.com/shared/static/bvfzmkpgtphvbmmrvb2iyl2jl21o49kc > dbnsfp.zip && '''
                '''unzip dbnsfp.zip && zcat dbNSFP*_variant.chr1* | awk "NR<=1" > h && '''
                '''zgrep -v "^#" dbNSFP*_variant.chr6* > chrs && '''
                '''zgrep -v "^#" dbNSFP*_variant.chr20* >> chrs && '''
                '''sort -S 50% --parallel=8 chrs -k1,1 -k2,2n > chrs_sorted && '''
                '''cat h chrs_sorted > chrs_sorted_header && '''
                '''bgzip -c chrs_sorted_header > tmp.vcf.gz && '''
                '''tabix -p vcf tmp.vcf.gz && '''
                '''zgrep ^# tmp.vcf.gz > full.tmp.vcf && '''
                '''tabix -R {workdir}/{input.bed} tmp.vcf.gz >> full.tmp.vcf && '''
                '''bgzip -c full.tmp.vcf > {workdir}/{output.vcf}) && '''
                '''tabix -s 1 -b 2 -e 2 {output.vcf} '''
            )


rule exons:
    """
    Get the exons for chromosomes 6 and 20. The paper that published the reference dataset indicates a high coverage in chromosome 6 for Illumina HiSeq 4000 reads
    """
    resources:
        mem_mb= 1024 * 2,
        disk_mb= 1024 * 2,
        runtime=30
    output: 'data/full_exons.bed'
    shell:
        'wget https://ftp.ensembl.org/pub/release-108/gff3/homo_sapiens/Homo_sapiens.GRCh38.108.chromosome.6.gff3.gz -O- > gff ; '
        'wget https://ftp.ensembl.org/pub/release-108/gff3/homo_sapiens/Homo_sapiens.GRCh38.108.chromosome.20.gff3.gz -O- >> gff ; '
        '''zgrep 'exon' gff | awk '{{print $1, $4, $5, $6, $7}}' > tmp ; '''
        '''zgrep 'exon' gff | awk '{{print $9}}' > IDs ; '''
        'sed -i "s/Parent=transcript:\\(.*\\);Name.*/\\1/" IDs ; '
        'paste tmp IDs > {output} ; '
        'sed -i "s/[[:blank:]]/\\t/g" {output} ; '
        'rm tmp IDs gff  '


rule LIMIT:
    """
    Get a limited bed file to subset the exon bed and select reads
    """
    resources:
        mem_mb= 1024 * 2,
        disk_mb= 1024 * 2,
        runtime=20
    output: 'data/LIMIT.bed'
    shell:
        'echo "6	42900000	42970000	chr6" > {output}; '
        'echo "20	46900000	46970000	chr20" >> {output}'


rule fasta:
    """
    Get the fasta file for chromosome 6 and 20
    """
    resources:
        mem_mb= 1024 * 4,
        disk_mb= 1024 * 4,
        runtime=60
    output: 'data/GRCh38.6.20.fa.gz',
    shell:
        'curl -L  ftp://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.6.fa.gz > temp ' 
        '&& curl -L ftp://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz >> temp '
        '&& zcat temp | bgzip -c > {output} '
        '&& rm temp '


rule fasta_index:
    """
    Make index for fasta file
    """
    input:
        reference
    resources:
        mem_mb = 1024*4,
        runtime = 60
    output: reference + '.fai'
    shell:
        "samtools faidx {input} > {output} "


rule WGS_download_normal:
    """
    Get the WGS samples using the RRA accession code for R1
    """
    resources:
        mem_mb= 1024 * 64,
        disk_mb= 1024 * 2000,
        runtime=24*60
    output:
        normal_R1='data/normal_R1_full.fq.gz',
        normal_R2='data/normal_R2_full.fq.gz',
    run:
        accession = WGS_accessions['normal']
        shell('{VDB_CONFIG_PRELUDE}; fasterq-dump --split-files -t /lscratch/$SLURM_JOBID -O /data/NICHD-core0/test/$USER/sra {accession} ' )
        shell('bgzip -c /data/NICHD-core0/test/$USER/sra/{accession}_1.fastq > {output.normal_R1} ')
        shell('bgzip -c /data/NICHD-core0/test/$USER/sra/{accession}_2.fastq > {output.normal_R2} ')
        shell('rm -rf /data/NICHD-core0/test/$USER/sra ')


rule WGS_download_tumor:
    """
    Get the WGS samples using the RRA accession code for R2
    """
    resources:
        mem_mb= 1024 * 64,
        disk_mb= 1024 * 2000,
        runtime=24*60
    output:
        tumor_R2='data/tumor_R2_full.fq.gz',
        tumor_R1='data/tumor_R1_full.fq.gz'
    run:
        accession = WGS_accessions['tumor']
        shell('{VDB_CONFIG_PRELUDE}; fasterq-dump --split-files -t /lscratch/$SLURM_JOBID -O /data/NICHD-core0/test/$USER/sra {accession} ' )
        shell('bgzip -c /data/NICHD-core0/test/$USER/sra/{accession}_1.fastq > {output.tumor_R1} ')
        shell('bgzip -c /data/NICHD-core0/test/$USER/sra/{accession}_2.fastq > {output.tumor_R2} ')
        shell('rm -rf /data/NICHD-core0/test/$USER/sra ')


rule bwa_index:
    """
    Get the index for BWA
    """
    resources:
        mem_mb= 1024 * 16,
        disk_mb= 1024 * 16,
        runtime=4*60
    input:
        reference,
    output:
        multiext(reference, '.amb', '.ann', '.bwt', '.pac', '.sa')
    shell:
        'bwa index -a bwtsw {input}'
    

rule align_reads:
    """
    Align the reads to the reference, since the reference is only Chr6, and Chr20 we will only have Chr6/20 reads
    """
    resources:
        mem_mb= 1024 * 32,
        disk_mb= 1024 * 24,
        runtime=54*60
    input:
        reads=['data/{sample}_R1_full.fq.gz', 'data/{sample}_R2_full.fq.gz'],
        idx=rules.bwa_index.output,
    output: 
        bam=temp('data/{sample}.sorted.bam')
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
    threads: 32
    shell:
        "bwa mem -t {threads} {params.index} {input.reads} | samtools view -bh | samtools sort -o {output} -O BAM " 


rule small_fastq:
    """
    Extract reads that are paired and mapped in proper pair
    Extract reads that are unmapped
    Collect the tags of unmapped and mapped reads, store them separately
    Use Seqtk to subset the fastq using the read tags
    """
    resources:
        disk_mb= 1024 * 100,
        mem_mb= 1024 * 32,
        runtime=24*60
    input:
        bam=rules.align_reads.output,
        limits=rules.LIMIT.output,
        r1='data/{sample}_R1_full.fq.gz',
        r2='data/{sample}_R2_full.fq.gz',
    output:
        R1='data/{sample}_R1.6.20.fq.gz',
        R2='data/{sample}_R2.6.20.fq.gz',
        mapped_names = 'data/{sample}.names.mapped.list',
        unmapped_names= 'data/{sample}.names.unmapped.list'
    threads: 32
    run:
        N = mapped_n
        UN = unmapped_n
        shell(
            'samtools view -h -L {input.limits} {input.bam} | samtools view -f 3 - | cut -f1 | sort -u > {wildcards.sample}.intermediate.file'

        )
        shell(
            'head -n {N} {wildcards.sample}.intermediate.file  > {output.mapped_names} '
        )
        shell(
            'samtools view -f 4 {input.bam} | cut -f1 | sort -u > {wildcards.sample}.int2.file'
        )
        shell (
            'head -n {UN} {wildcards.sample}.int2.file  > {output.unmapped_names} '
        )

        shell(
            'seqtk subseq {input.r1} {output.mapped_names} '
            '> {wildcards.sample}.r1.tmp '
        )

        shell(
            'seqtk subseq {input.r1} {output.unmapped_names} '
            '>> {wildcards.sample}.r1.tmp '
        )

        shell(
            'bgzip -c {wildcards.sample}.r1.tmp > {output.R1} '
        )

        shell(
            'seqtk subseq {input.r2} {output.mapped_names} '
            '> {wildcards.sample}.r2.tmp '
        )

        shell(
            'seqtk subseq {input.r2} {output.unmapped_names} '
            '>> {wildcards.sample}.r2.tmp '
        )

        shell(
            'bgzip -c {wildcards.sample}.r2.tmp > {output.R2} '
        )

        shell(
            'rm {wildcards.sample}.intermediate.file {wildcards.sample}.int2.file {wildcards.sample}.r1.tmp {wildcards.sample}.r2.tmp'
        )


rule small_bed:
    """
    Create a subset bed file from the exons file of chr 6 and chr 20.
    The small bed file will correspond to the reads that are in the small fastqs via subsetting with the LIMITs bed file
    """
    resources:
        mem_mb= 1024 * 2,
        disk_mb= 1024 * 2,
        runtime=4*60
    input:
        limit=rules.LIMIT.output,
        bed=rules.exons.output
    output:
        'data/exons_subset.bed'
    shell:
        'bedtools intersect -a {input.bed} -b {input.limit} > {output} '
