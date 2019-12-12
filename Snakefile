

# Manually Enter Sample Dictionary
# The order matters, group sample replicates together
sample_dict = {
    "LIB191023LC_ECL37_36_S11_R1_001":"SAL01L",
    "LIB191023LC_ECL37_25_S2_R1_001":"SAL02L",
    "LIB191023LC_ECL37_29_S6_R1_001":"SAL03L",
    "LIB191023LC_ECL37_24_S1_R1_001":"SAL04L",
    "LIB191023LC_ECL37_34_S9_R1_001":"SAL06L",
    "LIB191023LC_ECL37_26_S3_R1_001":"ITU01L",
    "LIB191023LC_ECL37_33_S8_R1_001":"ITU03L",
    "LIB191023LC_ECL37_27_S4_R1_001":"ITU04L",
    "LIB191023LC_ECL37_28_S5_R1_001":"ITU05L",
    "LIB191023LC_ECL37_32_S7_R1_001":"ITU06L",
    "LIB191023LC_ECL37_35_S10_R1_001":"ITU07L"
}

#samples = sample_dict.keys()
#suffix =


# Variables necessary for rules
GENOME = "/home/groups/hoolock/u1/genomes/Epigenetics_Core/mouse/m38"
ORGANISM = "mouse"
MY_TREAT = "0,0,0,0,0,1,1,1,1,1,1"
# coldata columns (list elements in same order as the samples appear in the sample dictionary above)
TREATMENT = "SAL,SAL,SAL,SAL,SAL,ITU,ITU,ITU,ITU,ITU,ITU"
CONDITION = "NULL"
REPLICATE = "1,2,3,4,6,1,3,4,5,6,7"
SEX = "NULL"
GENOTYPE = "NULL"
BATCH = "NULL"
PCA_COLOR_VAR = "treatment"
PCA_SHAPE_VAR = "NULL"
PCA_LABEL_VAR = "replicate"
# DMR variables
DMR_NORM = "TRUE"
DMR_NORM_METHOD = "median"
DMR_COMP_SAMPS = "SAL01L,SAL02L,SAL03L,SAL04L,SAL06L,ITU01L,ITU03L,ITU04L,ITU05L,ITU06L,ITU07L"
DMR_COMP_TREAT = "0,0,0,0,0,1,1,1,1,1,1"
DMR_MPG = "4"
DMR_WINDOW = "1000"
DMR_STEP = "1000"
DMR_COVAR1 = SEX
DMR_COVAR2 = GENOTYPE
DMR_OVERDISP = "MN"
DMR_TEST = "Chisq"
DMR_COMP_NAME = "ITU_vs_SAL"
DMR_METH_DIFF = "10"
DMR_QVAL = "0.1"

localrules: collect_trim_stats, collect_aln_stats


rule all:
    input:
        expand("data/fastqc/raw/{sample}.merged_fastqc.html", sample=sample_dict.keys()),
	"data/trimming/summary_trim",
	"data/bismark_aln/summary_aln",
	expand("data/meth_extract/{sample}.merged.canonical.cov.gz", sample=sample_dict.keys()),
	"data/ide/dim_reduct/coldata.txt",
	"data/DMR_analysis/complete.txt"

rule fastqc_raw:
    input:
        "data/fastq/{sample}.merged.fastq.gz"
    output:
        "data/fastqc/raw/{sample}.merged_fastqc.html"
    conda:
        "envs/test_BS_seq.yml"
    params:
        outdir = "data/fastqc/raw"
    shell:
        "fastqc -o {params.outdir} {input}"

rule trim_galore:
    input:
        "data/fastq/{sample}.merged.fastq.gz"
    output:
        "data/trimming/{sample}.merged_trimmed.fq.gz"
    conda:
        "envs/test_BS_seq.yml"
    params:
        outdir = "data/trimming",
	fastqc_dir = "data/fastqc/trimmed"
    shell:
        "trim_galore -o {params.outdir} --fastqc_args '--outdir {params.fastqc_dir}' "
	"--rrbs {input}"

rule bismark_aln:
    input:
        "data/trimming/{sample}.merged_trimmed.fq.gz"
    output:
        "data/bismark_aln/{sample}.merged_trimmed_bismark_bt2.bam"
    conda:
        "envs/test_BS_seq.yml"
    params:
        genome_dir = GENOME,
	outdir = "data/bismark_aln"
    shell:
        "bismark -p 4 {params.genome_dir} {input} -o {params.outdir} --bam"

rule collect_trim_stats:
    input:
        "data/trimming/*_report.txt"
    output:
        "data/trimming/summary_trim"
    conda:
        "envs/test_BS_seq.yml"
    shell:
        """
        grep "STATISTICS" *_report.txt | cut -d : -f1 > temp_name
        sed -i 's/\.fastq.gz.*\.txt/''/g; s/ /\_/g' temp_name
        grep -h --no-filename "processed in total" *_report.txt | cut -d ' ' -f1 > temp_rawReads
        grep -h --no-filename "truncated" *_report.txt | cut -f2 > temp_qualtrim
        grep -h --no-filename "shorter" *_report.txt | cut -f2 > temp_removed
        grep -h --no-filename "additional" *_report.txt | cut -f2 > temp_rrbs
        paste -d "\t" temp_name temp_rawReads temp_qualtrim temp_removed temp_rrbs > summary_trim
        sed -i '1 i\SampleName\tRawReads\tReadsTrimmed\tReadsRemoved\trrbsReadsTrimmed' summary_trim
        rm temp_*
        """

rule collect_aln_stats:
    input:
        "data/bismark_aln/*_report.txt"
    output:
        "data/bismark_aln/summary_aln"
    conda:
        "envs/test_BS_seq.yml"
    shell:
        """
        grep "Sequences analysed in total" *_report.txt | cut -d : -f1 > temp_name
	grep "Sequences analysed in total" *_report.txt | cut -f2 > temp_totReads
	grep "unique best hit" *_report.txt | cut -f2 > temp_unq
	grep "Mapping efficiency" *_report.txt | cut -f2 > temp_map
	grep "any condition" *_report.txt | cut -f2 > temp_noAln
	grep "not map uniquely" *_report.txt | cut -f2 > temp_nonUnq
	grep "extracted" *_report.txt | cut -f2 > temp_discard
	grep "C's analysed" *_report.txt | cut -f2 > temp_totCs
	grep "C methylated in CpG context" *_report.txt | cut -f2 > temp_cpg
	grep "C methylated in CHG context" *_report.txt | cut -f2 > temp_chg
	grep "C methylated in CHH context" *_report.txt | cut -f2 > temp_chh
	grep "C methylated in Unknown context" *_report.txt | cut -f2 > temp_unk
	paste -d "\t" temp_name temp_totReads temp_unq temp_map temp_noAln temp_nonUnq temp_discard temp_totCs temp_cpg temp_chg temp_chh temp_unk > summary_aln
	sed -i '1 i\SampleName\tInputReads\tUnqMapReads\tMapPer\tNoMapReads\tNonUnqMapReads\tDiscardReads\tTotalCs\tCpGcontext\tCHGcontext\tCHHcontext\tUnknownContext' summary_aln
	rm temp_*
        """

rule meth_extract:
    input:
        "data/bismark_aln/{sample}.merged_trimmed_bismark_bt2.bam"
    output:
        "data/meth_extract/{sample}.merged_trimmed_bismark_bt2.bismark.cov.gz",
	"data/meth_extract/{sample}.merged_trimmed_bismark_bt2.bedGraph.gz"
	#temp(touch("methExtract_complete"))
    conda:
        "envs/test_BS_seq.yml"
    params:
        outdir = "data/meth_extract",
	genome_dir = GENOME
    shell:
        "bismark_methylation_extractor -s --comprehensive --merge_non_CpG --multicore 3 "
	"-o {params.outdir} --gzip --bedGraph --cytosine_report "
	"--genome_folder {params.genome_dir} {input}"

rule make_can_covs:
    input:
        "data/meth_extract/{sample}.merged_trimmed_bismark_bt2.bismark.cov.gz"
    output:
        "data/meth_extract/{sample}.merged.canonical.cov.gz"
    conda:
        "envs/test_BS_seq.yml"
    params:
        organism = ORGANISM
    shell:
        "Rscript scripts/make_can_cov_files.R {input} {params.organism}"

#rule make_can_bedgraph_tracks:?

rule ide:
    input:
        expand("data/meth_extract/{sample}.merged.canonical.cov.gz", sample=sample_dict.keys())
    output:
        "data/ide/dim_reduct/coldata.txt"
    conda:
        "envs/test_BS_seq.yml"
    params:
        cov_files = ",".join(str(e) for e in ["data/meth_extract/" + s + '.merged.canonical.cov.gz' for s in sample_dict.keys()]),
        samp_names = ",".join(str(e) for e in [s for s in sample_dict.values()]),
        my_treat = MY_TREAT,
	my_treatment = TREATMENT,
	my_condition = CONDITION,
	my_rep = REPLICATE,
	my_sex = SEX,
	my_genotype = GENOTYPE,
	my_batch = BATCH,
	my_pca_color_var = PCA_COLOR_VAR,
	my_pca_shape_var = PCA_SHAPE_VAR,
	my_pca_label_var = PCA_LABEL_VAR
    shell:
        "Rscript scripts/run_ide.R {params.cov_files} {params.samp_names} {params.my_treat} "
        "{params.my_treatment} {params.my_condition} {params.my_rep} {params.my_sex} "
	"{params.my_genotype} {params.my_batch} {params.my_pca_color_var} {params.my_pca_shape_var} "
	"{params.my_pca_label_var}"

rule DMR_analysis:
    input:
        "data/ide/dim_reduct/coldata.txt"
    output:
        "data/DMR_analysis/complete.txt"
    conda:
        "envs/test_BS_seq.yml"
    params:
        cov_files = ",".join(str(e) for e in ["data/meth_extract/" + s + '.merged.canonical.cov.gz' for s in sample_dict.keys()]),
        samp_names = ",".join(str(e) for e in [s for s in sample_dict.values()]),
	my_treat = MY_TREAT,
        norm = DMR_NORM,
	norm_method = DMR_NORM_METHOD,
        comparison_samps = DMR_COMP_SAMPS,
        comparison_treat = DMR_COMP_TREAT,
        min_per_group = DMR_MPG,
        window_size = DMR_WINDOW,
        step_size = DMR_STEP,
	covariate1 = DMR_COVAR1,
        covariate2 = DMR_COVAR2,
        overdisp = DMR_OVERDISP,
        test = DMR_TEST,
        comp_name = DMR_COMP_NAME,
        meth_diff = DMR_METH_DIFF,
        qval = DMR_QVAL
    shell:
        "Rscript scripts/run_DMR.R {params.cov_files} {params.samp_names} {params.my_treat} {params.norm} {params.norm_method} "
        "{params.comparison_samps} {params.comparison_treat} {params.min_per_group} {params.window_size} "
	"{params.step_size} {params.covariate1} {params.covariate2} {params.overdisp} {params.test} "
        "{params.comp_name} {params.meth_diff} {params.qval}"
