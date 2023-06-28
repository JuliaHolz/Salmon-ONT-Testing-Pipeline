configfile: "config.yml"
tests = ["simulated_with_ir_corr", "simulated_no_ir_corr", "E2_spearman_corr", "E0_COV", "short_read_jsd", "ten_sample_jsd", "mixes_corr"]
salmon_paths = config["ont_salmons"]
transcriptome_paths = config["transcriptomes"]
ont_sra_data = list(config["ten_ONT_samples"].keys())
ont_sra_data.append(config["E0_sra"])
ont_sra_data.append(config["E2_sra"])

shortread_sra_data = list(config["four_illumina_samples"].keys())

rule plot all:
    input:
        plots = expand("plots/{test}.png", test = tests) 

rule download_gencode_transcriptome:
    input:
    output: config["transcriptomes"]["def"]
    shell: "wget " + config["transcriptome_urls"]["def"] + " -O " + config["transcriptomes"]["def"]

rule download_gencode_mouse_transcriptome:
    input:
    output: config["transcriptomes"]["mouse"]
    shell: "wget " + config["transcriptome_urls"]["mouse"] + " -O " + config["transcriptomes"]["mouse"]

rule create_SIRV_transcriptome:
    input: config["transcriptomes"]["def"]
    output: config["transcriptomes"]["sirv"]
    shell: "cat <( zcat {input})" + config["SIRV_transcripts"]  + "> {output}"

rule create_rnasequin_transcriptome:
    input: config["transcriptomes"]["def"]
    output: config["transcriptomes"]["rnasequin"]
    shell: "cat <( zcat {input})" + config["rnasequin_transcripts"]  + "> {output}"

rule download_ont_sra:
    input:
    output:
        fastq = "data/{type}/{id}.fastq"
    run: 
        shell("prefetch --max-size u -o {id}.sra {id}")
        shell("fasterq-dump -e 8 -x --split-3 {id}.sra -O data/{type}")

rule build_minimap_index: ## build minimap2 index
    input:
        genome = lambda wildcards:transcriptome_paths[wildcards.trs]
    output:
        index = "index/transcriptome_index_{trs}.mmi",
    params:
        opts = config["minimap_index_opts"]
    conda: "env.yml"
    threads: config["threads"]
    shell:"""
        minimap2 -t {threads} {params.opts} -Y -I 1000G -d {output.index} {input.genome}
    """


rule build_salmon_index: 
    input:
        salmon = config["shortread_salmon"],
        genome = config["transcriptome"]
    output:
        index = directory("index/salmon_index")
    conda: "env.yml"
    threads: config["threads"]
    shell:"""
        {input.salmon} index -t {input.genome}  -i index/salmon_index -p 8 
    """
rule map_reads: ## map reads using minimap2
    input:
       index = "index/transcriptome_index_{trs}.mmi",
       fastq = "data/ont_samples/{sample}.fastq"
    output:
       bam = "alignments/{sample}_{trs}.bam",
    log: "alignments/{sample}_{trs}.log",
    conda: "env.yml"
    threads: config["threads"]
    shell: "minimap2 -t {threads} -a -x map-ont -p 1.0 -N 100 {input.index} {input.fastq} | samtools view -Sb > {output.bam}"


rule quant_ont:
    input: 
        bam = "alignments/{sample}_{trs}.bam",
        trs =  lambda wildcards:transcriptome_paths[wildcards.trs],
        salmon = lambda wildcards: salmon_paths[wildcards.salmon_ver]
    output:
        tsv = "counts/{sample}_{salmon_ver}_{trs}/quant.sf",
    log:
        "counts/{sample}_{salmon_ver}_{trs}.log",
    benchmark:
        "counts/{sample}_{salmon_ver}_{trs}/benchmark.txt"
    conda: "env.yml"
    threads: config["threads"]
    shell: "{input.salmon} quant --ont -p 8 -t {input.trs} -l U -a {input.bam} -o counts/{wildcards.sample}_{wildcards.salmon_ver}_{wildcards.trs}" 

rule quant_shortread:
    input: 
        fastq1 = "data/four_illumina_samples/{sample}_1.fastq",
        fastq2 = "data/four_illumina_samples/{sample}_2.fastq",
        idx = "index/salmon_index",
        salmon = config["shortread_salmon"],
    output:
        tsv = "counts/{sample}_shortread/quant.sf",
    log:
        "counts/{sample}_shortread/log",
    conda: "env.yml"
    threads: config["threads"]
    shell: '''{input.salmon} quant -i {input.idx}  -la -1 {input.fastq1} -2 {input.fastq2}  --output counts/{wildcards.sample}_shortread -p 8 --gcBias --seqBias --posBias --thinningFactor 64'''
 
rule quant_no_em:
    input: 
        bam = "alignments/{sample}_{trs}.bam",
        trs =  lambda wildcards:transcriptome_paths[wildcards.trs],
        salmon = config["shortread_salmon"],
    output:
        tsv = "counts/{sample}_no_em_{trs}/quant.sf",
    log:
        "counts/{sample}_no_em_{trs}.log",
    benchmark:
        "counts/{sample}_no_em_{trs}/benchmark.txt"
    conda: "env.yml"
    threads: config["threads"]
    shell: "{input.salmon} quant --noErrorModel -p 8 -t {input.trs} -l U -a {input.bam} -o counts/{wildcards.sample}_no_em_{wildcards.trs}" 


rule calc_ten_sample_jsd:
    input:
        quantsf = expand("counts/{sample}_{sv}_def/quant.sf", sample = config["ten_ONT_samples"], sv = "{sv}"),
    output:
        jsdMat = "ten_samples/{sv}_def/JSDMat.csv",
        avgJSD = "ten_samples/{sv}_def/avgJSD.txt",
    log: "ten_samples/{sv}_def/log",
    script:
        "scripts/ten_sample_jsd.py"

rule calc_short_read_jsd:
    input:
        quantsf_ONT = expand("counts/{sample}_{sv}_def/quant.sf", sample = config["ten_ONT_samples"], sv = "{sv}"),
        quantsf_SR = expand("counts/{sample}_shortread/quant.sf", sample = config["four_illumina_samples"])
    output:
        jsdMat = "shortread/shortread_{sv}/JSDMat.csv",
        avgJSD = "shortread/shortread_{sv}/avgJSD.txt",
    log: "shortread/shortread_{sv}/log",
    script:
        "scripts/shortread_jsd.py"

rule rna_sequin:
    input:
        known_conc = "data/rnasequin_isoforms_2.4.tsv",
        quantsf_rnasequin = expand("counts/{sample}_{sv}_rnasequin/quant.sf", sample = config["rna_sequin_samples"], sv = "{sv}")
    output:
        result_file = "rna_sequin/{sv}_rnasequin/spearmancorr.txt"
    log: "rna_sequin/{sv}_rnasequin/log"
    script:
        "scripts/rna_sequin.py"

rule sirv_benchmark:
    input:
        quantsf_E0 = "counts/SRR6058584_{sv}_sirv/quant.sf",
        quantsf_E2 = "counts/SRR6058583_{sv}_sirv/quant.sf",
        sirv_conc = "data/SIRV_Set1_Norm_sequence-design-overview_20210507a.xlsx"
    output:
        result_E0 = "sirv_benchmark/{sv}_sirv/E0COV.txt",
        result_E2 = "sirv_benchmark/{sv}_sirv/E2spearmancorr.txt"
    log: "sirv_benchmark/{sv}_sirv/log"
    script:
        "scripts/sirv_benchmark.py"


rule simulated_reads:
    input:
        quantsf_no_ir = "counts/sim_mouse_no_ir_{sv}_mouse/quant.sf",
        quantsf_with_ir = "counts/sim_mouse_with_ir_{sv}_mouse/quant.sf",
        simulated_abundance = "data/sim_abundance.tsv"
    output:
        result_no_ir = "simulated_reads/{sv}_simmouse/no_ir_spearmancorr.txt",
        result_with_ir = "simulated_reads/{sv}_simmouse/with_ir_spearmancorr.txt"
    log: "simulated_reads/{sv}_simmouse/log"
    script:
        "scripts/simulated_reads.py"

rule plot_ten_ont:
    input: 
        avgJSD = expand("ten_samples/{sv}_def/avgJSD.txt", sv = config["ont_salmons"])
    output: 
        jsdplot = "plots/ten_sample_jsd.png"
    script: 
        "scripts/plot_ten_sample_jsd.py"

rule plot_short_read:
    input: 
        avgJSD = expand("shortread/shortread_{sv}/avgJSD.txt", sv = config["ont_salmons"])
    output: 
        jsdplot = "plots/short_read_jsd.png"
    script: 
        "scripts/plot_short_read_jsd.py"

rule plot_E0:
    input:
        COV = expand("sirv_benchmark/{sv}_sirv/E0COV.txt", sv = config["ont_salmons"])
    output:
        COVplot = "plots/E0_COV.png"
    script: 
        "scripts/plot_E0_COV.py"

rule plot_E2:
    input: 
        spearmancorr = expand("sirv_benchmark/{sv}_sirv/E2spearmancorr.txt", sv = config["ont_salmons"])
    output:
        corrplot = "plots/E2_spearman_corr.png"
    script: 
        "scripts/plot_E2_corr.py"

rule plot_mixes: 
    input: 
        spearmancorr = expand("rna_sequin/{sv}_rnasequin/spearmancorr.txt", sv = config["ont_salmons"])
    output:
        corrplot = "plots/mixes_corr.png"
    script: 
        "scripts/plot_mixes_corr.py"

rule plot_no_ir: 
    input: 
        spearmancorr = expand("simulated_reads/{sv}_simmouse/no_ir_spearmancorr.txt", sv = config["ont_salmons"])
    output:
        corrplot = "plots/simulated_no_ir_corr.png"
    script: 
        "scripts/plot_no_ir.py"

rule plot_with_ir: 
    input: 
        spearmancorr = expand("simulated_reads/{sv}_simmouse/with_ir_spearmancorr.txt", sv = config["ont_salmons"])
    output:
        corrplot = "plots/simulated_with_ir_corr.png"
    script: 
        "scripts/plot_with_ir.py"

