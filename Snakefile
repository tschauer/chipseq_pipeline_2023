configfile: "configs/config.yaml"

import pandas as pd

sample_table = pd.read_table("SampleTable.txt")

SRRs = list(sample_table['Run'].unique())


################################################################################


rule all:
    input:
        "multiqc_report.html",

################################################################################




################################################################################


rule multiqc:
    input:
        expand("Output/fastqc/{SRR}_1_fastqc.html", SRR=SRRs),
        expand("Output/peaks/{SRR}_filtered.bed", SRR=SRRs),
        expand("Output/metaplots/{SRR}.TSS_plot.png", SRR=SRRs),
        expand("Output/metaplots/{SRR}.scaled_plot.png", SRR=SRRs),
        expand("Output/fragsize/{SRR}.fragmentsize.png", SRR=SRRs),
        "Output/results/heatmap_SpearmanCorr_readCounts.pdf",
        "Output/results/heatmap_Jaccard_peaks.pdf",
    output:
        "multiqc_report.html",
    threads: 4
    shell:
        """
        multiqc . --config configs/multiqc_config.yaml
        """

################################################################################


rule jaccard:
    input:
        expand("Output/peaks/{SRR}_filtered.bed", SRR=SRRs),
    output:
        "Output/results/heatmap_Jaccard_peaks.pdf",
    shell:
        """
        Rscript --vanilla scripts/jaccard.R {input} {output}
        """


################################################################################


rule peakcalling:
    input:
        removeddup="Output/BAM/{SRR}.removeddup.bam",
        bed="genome/GRCm38-blacklist.v2.bed",
    output:
        bed="Output/peaks/{SRR}_filtered.bed",
        log1="logs/{SRR}_macs2.out",
        log2="logs/{SRR}_macs2.err",
    params:
        mapqc=12
    threads: 8
    shell:
        """
        macs2 callpeak --treatment {input.removeddup} --name {wildcards.SRR} \
        --gsize 1.87e9 --format BAMPE --nomodel --nolambda \
        --cutoff-analysis --broad --broad-cutoff 0.001 --pvalue 0.001 \
        --outdir Output/peaks > {output.log1} 2> {output.log2}

        bedtools intersect -v -a "Output/peaks/{wildcards.SRR}_peaks.broadPeak" \
        -b {input.bed} | cut -f 1,2,3 > {output.bed}

        """

################################################################################


rule scaled_plot:
    input:
        mat="Output/matrix/{SRR}.scaled_matrix.gz",
    output:
        plot="Output/metaplots/{SRR}.scaled_plot.png",
    shell:
        """
        plotProfile -m {input.mat} --plotHeight 10 --plotWidth 14 -o {output.plot}
        """

rule scaled_matrix:
    input:
        bw="Output/coverage/{SRR}.coverage.bw",
        gtf="genome/genome.gtf",
    output:
        mat="Output/matrix/{SRR}.scaled_matrix.gz",
        log1="logs/{SRR}_scaled_mat.out",
        log2="logs/{SRR}_scaled_mat.err",
    threads: 8
    shell:
        """
        computeMatrix scale-regions \
        --numberOfProcessors {threads} \
        -S {input.bw} -R {input.gtf} \
        --binSize 20 --skipZeros \
        --downstream 2000 --upstream 2000 \
        --regionBodyLength 4000 \
        -o {output.mat} \
        > {output.log1} 2> {output.log2}
        """


################################################################################


rule TSSplot:
    input:
        mat="Output/matrix/{SRR}.TSS_matrix.gz",
    output:
        plot="Output/metaplots/{SRR}.TSS_plot.png",
    shell:
        """
        plotProfile -m {input.mat} --plotHeight 10 --plotWidth 12 -o {output.plot}
        """

rule TSSmatrix:
    input:
        bw="Output/coverage/{SRR}.coverage.bw",
        gtf="genome/genome.gtf",
    output:
        mat="Output/matrix/{SRR}.TSS_matrix.gz",
        log1="logs/{SRR}_TSSmat.out",
        log2="logs/{SRR}_TSSmat.err",
    threads: 8
    shell:
        """
        computeMatrix reference-point \
        --numberOfProcessors {threads} \
        -S {input.bw} -R {input.gtf} \
        --referencePoint 'TSS' \
        --binSize 20 --skipZeros \
        --downstream 2000 --upstream 2000 \
        -o {output.mat} \
        > {output.log1} 2> {output.log2}
        """


################################################################################


rule plot_corr:
    input:
        "Output/results/SpearmanCorr_readCounts.tab",
    output:
        "Output/results/heatmap_SpearmanCorr_readCounts.pdf",
    shell:
        """
        Rscript --vanilla scripts/corr_plot.R {input} {output}
        """

rule correlation:
    input:
        "Output/results/scores_per_bin.npz",
    output:
        tab="Output/results/SpearmanCorr_readCounts.tab",
        png="Output/results/heatmap_SpearmanCorr_readCounts.png",
    threads: 8
    shell:
        """
        plotCorrelation \
        --corData {input} \
        --corMethod spearman --skipZeros \
        --plotTitle "Spearman Correlation of Read Counts" \
        --whatToPlot heatmap --colorMap RdBu_r --plotNumbers \
        --plotHeight 14.25 --plotWidth 16.5 \
        -o {output.png}   \
        --outFileCorMatrix {output.tab}
        """


rule multisummary:
    input:
        bw=expand("Output/coverage/{SRR}.coverage.bw", SRR=SRRs),
        bed="genome/GRCm38-blacklist.v2.bed",
    output:
        npz="Output/results/scores_per_bin.npz",
        tab="Output/results/scores_per_bin.tab",
    threads: 8
    shell:
        """
        multiBigwigSummary bins --numberOfProcessors {threads} \
        --bwfiles {input.bw} \
        --blackListFileName {input.bed} --smartLabels \
        --outFileName  {output.npz} \
        --outRawCounts {output.tab}
        """

rule bamcoverage:
    input:
        removeddup="Output/BAM/{SRR}.removeddup.bam",
        bed="genome/GRCm38-blacklist.v2.bed",
    output:
        bw="Output/coverage/{SRR}.coverage.bw",
        log1="logs/{SRR}_bamcov.out",
        log2="logs/{SRR}_bamcov.err",
    params:
        mapqc=12
    threads: 8
    shell:
        """
        bamCoverage --bam {input.removeddup} -o {output.bw} \
        --blackListFileName {input.bed} \
        --binSize 20 --smoothLength 60 --extendReads  \
        --maxFragmentLength 700 --minMappingQuality {params.mapqc} \
        --normalizeUsing "CPM" --numberOfProcessors {threads} \
        > {output.log1} 2> {output.log2}
        """


################################################################################


rule fragsize:
    input:
        removeddup="Output/BAM/{SRR}.removeddup.bam",
    output:
        png="Output/fragsize/{SRR}.fragmentsize.png",
        log="logs/{SRR}_frags.out",
    threads: 8
    shell:
        """
        bamPEFragmentSize -hist {output.png} -p {threads} \
        -T "Fragment size" --maxFragmentLength 800 \
        --samplesLabel {wildcards.SRR} \
        -b {input.removeddup} > {output.log}
        """


################################################################################


rule picard:
    input:
        filtered="Output/BAM/{SRR}.filtered.bam",
    output:
        removeddup="Output/BAM/{SRR}.removeddup.bam",
        bai="Output/BAM/{SRR}.removeddup.bam.bai",
        log1="logs/{SRR}_picard.out",
        log2="logs/{SRR}_picard.err",
        log3="logs/{SRR}_idxstats.out",
    threads: 8
    shell:
        """
        picard MarkDuplicates -I {input.filtered} -O {output.removeddup} \
        -CREATE_INDEX TRUE -REMOVE_DUPLICATES TRUE \
        -METRICS_FILE {output.log1} 2> {output.log2}

        samtools index {output.removeddup}
        samtools idxstats {output.removeddup} | head -n 22 | sort -k 3,3nr > {output.log3}
        """

rule samtools:
    input:
        bam="Output/BAM/{SRR}.mapped.bam",
    output:
        filtered="Output/BAM/{SRR}.filtered.bam",
        bai="Output/BAM/{SRR}.filtered.bam.bai",
        log="logs/{SRR}_samtools.out",
    threads: 8
    params:
        mapqc= 12
    shell:
        """
        samtools view -bS -@ {threads} -q {params.mapqc} {input.bam} | \
        samtools sort -@ {threads} - | tee {output.filtered} | \
        samtools index - {output.bai} > {output.log}
        """


################################################################################


rule bowtie2:
    input:
        trimmed1="Output/trimmed/{SRR}_1_val_1.fq.gz",
        trimmed2="Output/trimmed/{SRR}_2_val_2.fq.gz",
        index="genome/genome.rev.2.bt2",
    output:
        bam="Output/BAM/{SRR}.mapped.bam",
        log="logs/{SRR}_bowtie2.err",
    threads: 8
    shell:
        """
        bowtie2 -x genome/genome --threads {threads} \
        --local --very-sensitive-local --no-unal \
        --no-mixed --no-discordant --dovetail -I 10 -X 700 \
        -1 {input.trimmed1} -2 {input.trimmed2} 2> {output.log} | samtools view -Sbh -o {output.bam}
        """

rule bowtie2_index:
    input:
        "genome/genome.fa"
    output:
        "genome/genome.rev.2.bt2"
    threads: 8
    shell:
        """
        bowtie2-build --threads {threads} {input} genome/genome
        """


################################################################################


rule trim_galore:
    input:
        fastq1="FastQ/{SRR}_1.fastq.gz",
        fastq2="FastQ/{SRR}_2.fastq.gz",
    output:
        trimmed1="Output/trimmed/{SRR}_1_val_1.fq.gz",
        trimmed2="Output/trimmed/{SRR}_2_val_2.fq.gz",
        log1="logs/{SRR}_trimgalore.out",
        log2="logs/{SRR}_trimgalore.err",
    params:
        dir="Output/trimmed/"
    threads: 8
    shell:
        """
        trim_galore -j {threads} --quality 28 --paired  \
        {input.fastq1} {input.fastq2} > {output.log1} 2> {output.log2}

        mv {wildcards.SRR}_1_val_1.fq.gz {params.dir}
        mv {wildcards.SRR}_2_val_2.fq.gz {params.dir}
        mv {wildcards.SRR}_1.fastq.gz_trimming_report.txt {params.dir}
        mv {wildcards.SRR}_2.fastq.gz_trimming_report.txt {params.dir}
        """


################################################################################


rule fastqc:
    input:
        fastq="FastQ/{SRR}_1.fastq.gz",
    output:
        html="Output/fastqc/{SRR}_1_fastqc.html",
        log1="logs/{SRR}_fastqc.out",
        log2="logs/{SRR}_fastqc.err",
    threads: 8
    shell:
        """
        fastqc {input.fastq} -t {threads} -o Output/fastqc > {output.log1} 2> {output.log2}
        """


################################################################################


# rule fastq_dump:
#     input:
#         "SRA/{SRR}/{SRR}.sra"
#     output:
#         "FastQ/{SRR}_1.fastq.gz",
#         "FastQ/{SRR}_2.fastq.gz",
#     shell:
#         """
#         fastq-dump --split-files --gzip {input} --outdir FastQ
#         """
#
# rule prefetch:
#     output:
#         "SRA/{SRR}/{SRR}.sra"
#     shell:
#         "prefetch --max-size 100G -O SRA {wildcards.SRR}"


################################################################################


rule get_blacklist:
    params:
        blackListLocation=config['bedBlackList'],
    output:
        mm10="genome/mm10-blacklist.v2.bed",
        ens="genome/GRCm38-blacklist.v2.bed"
    shell:
        """
        wget -O genome/mm10-blacklist.v2.bed.gz {params.blackListLocation}
        gunzip genome/mm10-blacklist.v2.bed.gz

        cat {output.mm10} | sed 's/^chr//' > {output.ens}
        """


################################################################################


rule get_genome_files:
    params:
        genomeLocation=config['genomeFTP'],
        gtfLocation=config['gtfFTP'],
    output:
        "genome/genome.fa",
        "genome/genome.fa.fai",
        "genome/genome.gtf",
    shell:
        """
        wget -O genome/genome.fa.gz {params.genomeLocation}
        gunzip genome/genome.fa.gz

        samtools faidx genome/genome.fa

        wget -O genome/genome.gtf.gz {params.gtfLocation}
        gunzip genome/genome.gtf.gz
        """


################################################################################


onsuccess:
        print("Finished!")


################################################################################
