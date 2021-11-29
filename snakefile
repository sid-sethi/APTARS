import os
import Bio
from os import path
from Bio import SeqIO
from re import search

from snakemake.utils import min_version
min_version("5.3")

# ----------------------------------------------------------------

configfile: "config.yml"
workdir: config["workdir"]

WORKDIR = config["workdir"]
SNAKEDIR = path.dirname(workflow.snakefile)

sample = config.get("sample_name", "sample")
gene_ids = config.get("genes", "")


# reading primer IDs into a list
primer_ids = []
five_p = ""
if config.get("primers"):
    records = list(SeqIO.parse(config["primers"], "fasta"))
    for seq in records:
        if search("_5p", seq.id):
            five_p = seq.id
        else:
            primer_ids.append(seq.id)

# ----------------------------------------------------------------


rule all:
    input:
        "Sqanti/" + sample + "_classification.txt"

#########################################################################


rule generate_consensus_reads:
    input:
        #bam = config["input_bam"]
        bam = config.get("input_bam", "")

    output:
        ccs_bam = "Consensus_reads/" + sample + "_ccs.bam"

    threads: config["threads"]

    shell:
        """
        ccs {input.bam} {output.ccs_bam} --minLength 10 --maxLength 50000 --minPasses 3 --minSnr 2.5 --maxPoaCoverage 0 --minPredictedAccuracy 0.99 -j {threads}
        """



rule demultiplex:
    input:
        ccs_bam = rules.generate_consensus_reads.output.ccs_bam,
        #primers = config["primers"]
        primers = config.get("primers", "")

    output:
        bams = expand("Demultiplexed/" + sample + "." + five_p + "--{ids}.bam", ids=primer_ids)

    params:
        outfile = "Demultiplexed/" + sample + ".bam"

    threads: config["threads"]

    shell:
        """
        lima {input.ccs_bam} {input.primers} {params.outfile} -j {threads} --isoseq --peek-guess
        """




rule generate_fofn:
    input:
        bams = expand("Demultiplexed/" + sample + "." + five_p + "--{ids}.bam", ids=primer_ids)

    output:
        fofn = "Demultiplexed/" + sample + ".fofn"

    run:
        textfile = open(output.fofn, "w")
        for element in primer_ids:
            textfile.write(sample + "." + five_p + "--" + element + ".bam\n")
        textfile.close()




rule refine:
    input:
        fofn = rules.generate_fofn.output.fofn,
        #primers = config["primers"]
        primers = config.get("primers", "")

    output:
        flnc = "Refine/" + sample + "_flnc.bam",
        refine_report = "Refine/" + sample + "_flnc.report.csv"

    params:
        logfile = "Refine/" + sample + "_refine.log"

    threads: config["threads"]

    shell:
        """
        isoseq3 refine --require-polya --log-level DEBUG --log-file {params.logfile} -j {threads} {input.fofn} {input.primers} {output.flnc}
        """




rule cluster:
    input:
        flnc = rules.refine.output.flnc

    output:
        clustered = "Cluster/" + sample + "_clustered.bam",
        hq_fasta = "Cluster/" + sample + "_clustered.hq.fasta.gz",
        report = "Cluster/" + sample + "_clustered.cluster_report.csv"

    params:
        logfile = "Cluster/" + sample + "_cluster.log"

    threads: config["threads"]

    shell:
        """
        isoseq3 cluster --use-qvs --verbose -j {threads} --log-file {params.logfile} {input.flnc} {output.clustered}
        """



rule minimap_mapping:
    input:
        genome = config["genome"],
        fa = rules.cluster.output.hq_fasta

    output:
        sam = "Mapping/" + sample + "_minimap.sam"

    threads: config["threads"]

    shell:
        """
        minimap2 -ax splice:hq -uf --secondary=no -t {threads} -o {output.sam} {input.genome} {input.fa}
        """


rule sort_sam:
    input:
        sam = rules.minimap_mapping.output.sam

    output:
        sortedSam = "Mapping/" + sample + "_minimap.sorted.sam"

    threads: config["threads"]

    shell:
        """
        samtools sort -O SAM -o {output.sortedSam} -@ {threads} {input.sam}
        """



rule gunzip_fa:
    input:
        fa_gz = rules.cluster.output.hq_fasta

    output:
        fa = "Cluster/" + sample + "_clustered.hq.fasta"

    shell:
        """
        gunzip -c {input.fa_gz} > {output.fa}
        """



rule get_gene_locus:
    input:
        gtf = config["gtf"]

    output:
        bed = "Mapping/" + sample + "_geneLocus.bed"

    params:
        genes = gene_ids,
        prefix = sample,
        respath = "Mapping",
        script = SNAKEDIR + "/scripts/get_gene_locus.R"

    shell:
        """
        Rscript {params.script} {input.gtf} {params.genes} {params.prefix} {params.respath}
        """




rule select_locus_from_sam:
    input:
        sam = rules.sort_sam.output.sortedSam,
        bed = rules.get_gene_locus.output.bed

    output:
        sam = "Mapping/" + sample + "_locus.sam"

    shell:
        """
        samtools view -L {input.bed} -o {output.sam} -h {input.sam}
        """



rule collapse_isoforms:
    input:
        sam = rules.sort_sam.output.sortedSam if gene_ids == "" else rules.select_locus_from_sam.output.sam,
        fa = rules.gunzip_fa.output.fa

    output:
        gff = "Collapsed_isoforms/" + sample + ".collapsed.gff",
        group = "Collapsed_isoforms/" + sample + ".collapsed.group.txt",
        mapped_fa = "Collapsed_isoforms/" + sample + ".collapsed.rep.fa"

    params:
        prefix = "Collapsed_isoforms/" + sample

    log: "Collapsed_isoforms/" + sample + "__collapse_isoforms_by_sam.log"

    shell:
        """
        (collapse_isoforms_by_sam.py --input {input.fa} -s {input.sam} -o {params.prefix} --dun-merge-5-shorter) 2> {log}
        """



rule get_abundance_all:
    input:
        group = rules.collapse_isoforms.output.group,
        cluster_report = rules.cluster.output.report

    output:
        abund = "Collapsed_isoforms/" + sample + ".collapsed.abundance.txt",
        read_stat = "Collapsed_isoforms/" + sample + ".collapsed.read_stat.txt"

    params:
        prefix = "Collapsed_isoforms/" + sample + ".collapsed"

    log: "Collapsed_isoforms/" + sample + "_get_abundance_all.log"

    shell:
        """
        (get_abundance_post_collapse.py {params.prefix} {input.cluster_report}) 2> {log}
        """



rule get_abundance_demux:
    input:
        mapped_fa = rules.collapse_isoforms.output.mapped_fa,
        read_stat = rules.get_abundance_all.output.read_stat,
        classify_csv = rules.refine.output.refine_report

    output:
        abund = "Collapsed_isoforms/" + sample + ".mapped_fl_count.txt"

    log: "Collapsed_isoforms/" + sample + "_get_abundance_demux.log"

    shell:
        """
        (demux_isoseq_with_genome.py --mapped_fafq {input.mapped_fa} --read_stat {input.read_stat} --classify_csv {input.classify_csv} -o {output.abund}) 2> {log}
        """



rule sqanti_qc:
    input:
        isoforms = rules.collapse_isoforms.output.gff,
        fl_count = rules.get_abundance_demux.output.abund if config.get("primers") else rules.get_abundance_all.output.abund,
        gtf = config["gtf"],
        genome = config["genome"],
        cage = config["cage"],
        sjs = config["intropolis"],
        poly_peak = config["polya_atlas"],
        poly_motifs = SNAKEDIR + "/data/human.polyA.list.txt"

    output:
        classification = "Sqanti/" + sample + "_classification.txt",
        juncs = "Sqanti/" + sample + "_junctions.txt"

    params:
        res_dir = 'Sqanti/',
        prefix = sample,
        sqanti_dir = config["sqanti_dir"],
        chunks = 10 if gene_ids == "" else 1

    threads: config["threads"]

    log: "Sqanti/" + sample + "_sqanti.log"

    conda: config["sqanti_dir"] + "/SQANTI3.conda_env.yml"

    shell:
        """
        (PYTHONPATH=$CONDA_PREFIX/bin python {params.sqanti_dir}/sqanti3_qc.py {input.isoforms} {input.gtf} {input.genome} --cage_peak {input.cage} --polyA_peak {input.poly_peak} --polyA_motif_list {input.poly_motifs} -c {input.sjs} -t {threads} --chunks {params.chunks} --output {params.prefix} --dir {params.res_dir} --report html -fl {input.fl_count}) 2> {log}
        """
