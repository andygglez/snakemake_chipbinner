
import pandas as pd
import yaml
import os

configfile: "config/config.yaml"


# default_resources: "group_components=sequential_group=1" ##### Important line for `compare` rule

if config['general']['MS_normalize']:
    with open(config['references']['MS_yaml'], "r") as file:
        MS_coefficients = yaml.safe_load(file)

def get_MS_coefficient(wildcards):
    return MS_coefficients[wildcards.sample]

def get_condition(wildcards):
    return wildcards.condition

def get_mark(wildcards):
    return wildcards.mark

def replace_suffix(sample:str):
    return sample.replace(".bam","")

def split_comparison(comparison:str):
    return comparison.split("_")

def get_treatment(wildcards):
    return wildcards.comparison.split("_")[0]

def get_parental(wildcards):
    return wildcards.comparison.split("_")[1]

sample_list = pd.read_csv("config/sample_list.csv")
conditions = list(sample_list["Condition"].unique())


condition_sample = {}
for condition in conditions:
    condition_sample[condition] = []
    rep1 = sample_list.loc[sample_list["Condition"] == condition, "Rep1"].to_list()
    if config['general']['merge_replicates']:
        rep2 = sample_list.loc[sample_list["Condition"] == condition, "Rep2"].to_list()
    condition_sample[condition].extend(rep1)
    if config['general']['merge_replicates']:
        condition_sample[condition].extend(rep2)
    # mark = sample_list.loc[sample_list["Condition"] == condition, "Mark"].to_list()[0]
    # condition_sample[condition].append(mark)


condition_mark = {}
for condition in conditions:
    if condition_mark.get(condition) is None:
        condition_mark[condition] = []
    condition_mark[condition].extend(sample_list.loc[sample_list["Condition"] == condition, "Mark"].to_list())

condition_sample_mark = {key:{mark:[] for mark in sample_list['Mark']} for key in conditions}
for condition in condition_sample_mark.keys():
    for mark in condition_sample_mark[condition]:
        condition_sample_mark[condition][mark] = list(sample_list.loc[(sample_list['Condition']==condition) & (sample_list['Mark']==mark)][['Rep1', 'Rep2']].values[0])

bigWig_folder = 'MSbigwigs' if config['general']['MS_normalize'] else 'notMSbigwigs'


# values = [v for k in dic.keys() for v in condition_sample[k]]
# merging = True if len(values)>len(conditions) else False

final_files = []
for condition in condition_sample.keys():
    if len(condition_sample[condition])==1:
        final_files.extend("outdir/notMerged/"+condition+"."+mark+".bw" for mark in condition_mark[condition])
    if len(condition_sample[condition])==2:
        final_files.extend("outdir/Merged/"+condition+"."+mark+".bw" for mark in condition_mark[condition])

final_bw_files = []
for condition in condition_sample.keys():
    final_bw_files.extend("outdir/"+bigWig_folder+"_final/"+condition+"."+mark+".bw" for mark in condition_mark[condition])


comparison_marks = {}
for comp in config['comparisons']:
    treatment = comp.split("_")[0]
    parental = comp.split("_")[1]
    common_marks = list(set(condition_mark[treatment]).intersection(set(condition_mark[parental])))
    # scatterplots.extend([f"outdir/comparisons/{comp}_"+mark+".scatterplot.pdf" for mark in common_marks])
    comparison_marks[comp] = common_marks

rule all:
    input:
        binned_bed = expand("outdir/bam2bed/{sample}.bed", sample=[samp.replace(".bam", "") for cond in condition_sample.keys() for samp in condition_sample[cond]]),
        filt_bed = expand("outdir/bam2bed/{sample}.binned_filt.bed", sample=[samp.replace(".bam", "") for cond in condition_sample.keys() for samp in condition_sample[cond]]) if config['general']['filter_low_counts'] else [],
        final_bw_files = final_bw_files,
        scatterplots = [f"outdir/comparisons/{comp}/{comp}_{mark}.scatterplot.pdf" for comp in comparison_marks for mark in comparison_marks[comp]],
        # cluster_info = [f"outdir/comparisons/{comp}/clustering/clus.{comp}_{mark}_clustered.100.100.txt" for comp in comparison_marks for mark in comparison_marks[comp]],
        # density_scatter = [f"outdir/comparisons/{comp}/clustering/{comp}_{mark}_2clusters.100.100.plots_density_scatterplot.pdf" for comp in comparison_marks for mark in comparison_marks[comp]]

rule bam2bed:
    input:
        bam = config['inputs']['sample_folder']+"/{sample}.bam"
    output:
        bed = "outdir/bam2bed/{sample}.bed",
        binned_bed = "outdir/bam2bed/{sample}.binned.bed"
    params:
        windows_ref = config['references']['binning_windows'],
        # blacklist = config['references']['blacklist']
    shell:
        """
            mkdir -p outdir/bam2bed \
        &&  bedtools bamtobed -i {input.bam} | sort -k1,1 -k2,2n > {output.bed} \
        &&  bedtools intersect -a {params.windows_ref} -b {output.bed} -sorted -c -nonamecheck > {output.binned_bed}
        """

if config['general']['filter_low_counts']:
    rule filter_low_counts:
        input:
            binned_bed = "outdir/bam2bed/{sample}.binned.bed"
        output:
            filtered_bed = "outdir/bam2bed/{sample}.binned_filt.bed"
        params:
            organism = config['general']['organism'],
            cutoff = config['general']['filter_low_counts']
        shell:
            """
                ~/environments/chipbinner/bin/Rscript scripts/filter.low.counts.R \
                                                        --binned_bed {input.binned_bed} \
                                                        --outdir outdir/bam2bed/ \
                                                        --organism {params.organism} \
                                                        --cutoff {params.cutoff}
            """

######### Optional
######### Normalization step:
# using DESeq2's median ratio method
# apply_normFactors(
#   norm_method = "DESeq2",
#   out_dir = "out_dir",
#   genome_assembly = "mm10",
#   treated_sample_bedfiles = c("Cal27_NSD1KO_H3K36me2_rep1.10kb.bed", "Cal27_NSD1KO_H3K36me2_rep2.10kb.bed"),
#   wildtype_sample_bedfiles = c("Cal27_WT_H3K36me2_rep1.10kb.bed", "Cal27_WT_H3K36me2_rep2.10kb.bed"),
#   treated_condition_label = "NSD1KO",
#   wildtype_condition_label = "WT"
# )

if config['general']['MS_normalize']:
    rule bed2bw:
        input:
            binned_bed = "outdir/bam2bed/{sample}.binned.bed" if not config['general']['filter_low_counts'] else "outdir/bam2bed/{sample}.binned_filt.bed",
            # blacklisted_bed = "outdir/bam2bed/{sample}.blacklisted.bed"
        output:
            bw = "outdir/"+bigWig_folder+"/{sample}.binned.bw" if not config['general']['filter_low_counts'] else "outdir/"+bigWig_folder+"/{sample}.binned_filt.bw"
        params:
            organism = config['general']['organism'],
            MS_factor = get_MS_coefficient
        shell:
            """
                mkdir -p outdir/MSbigwigs \
            &&  ~/environments/chipbinner/bin/Rscript \
                scripts/chipbinner.norm.bw.R --binned_bed {input.binned_bed} \
                                --outdir outdir/MSbigwigs --organism {params.organism} \
                                --MS_factor {params.MS_factor} \
            &&  echo "for {wildcards.sample} used {params.MS_factor} as MS coefficient"
            """
else:
    rule bed2bw:
        input:
            binned_bed = "outdir/bam2bed/{sample}.binned.bed" if not config['general']['filter_low_counts'] else "outdir/bam2bed/{sample}.binned_filt.bed",
            # blacklisted_bed = "outdir/bam2bed/{sample}.blacklisted.bed"
        output:
            bw = "outdir/"+bigWig_folder+"/{sample}.binned.bw" if not config['general']['filter_low_counts'] else "outdir/"+bigWig_folder+"/{sample}.binned_filt.bw"
        params:
            organism = config['general']['organism']
        shell:
            """
                mkdir -p outdir/notMSbigwigs \
            &&  ~/environments/chipbinner/bin/Rscript \
                scripts/chipbinner.norm.bw.R --binned_bed {input.binned_bed} \
                                --outdir outdir/notMSbigwigs --organism {params.organism} \
                                --MS_factor 1 \
            &&  echo "for {wildcards.sample} used 1 as MS coefficient"
            """

if config['general']['merge_replicates']:
    rule merged_bw_files:
        input:
            rep1 = lambda wildcards: f"outdir/{bigWig_folder}/{replace_suffix(condition_sample_mark[wildcards.condition][wildcards.mark][0])}.binned.bw" if not config['general']['filter_low_counts'] else lambda wildcards: f"outdir/{bigWig_folder}/{replace_suffix(condition_sample_mark[wildcards.condition][wildcards.mark][0])}.binned_filt.bw",
            rep2 = lambda wildcards: f"outdir/{bigWig_folder}/{replace_suffix(condition_sample_mark[wildcards.condition][wildcards.mark][1])}.binned.bw" if not config['general']['filter_low_counts'] else lambda wildcards: f"outdir/{bigWig_folder}/{replace_suffix(condition_sample_mark[wildcards.condition][wildcards.mark][1])}.binned_filt.bw"
        output:
            bw_merged = "outdir/"+bigWig_folder+"_final/{condition}.{mark}.bw"
        params:
            organism = config['general']['organism'],
            mark = get_mark,
            condition = get_condition,
            bw_folder=bigWig_folder
        shell:
            """
                mkdir -p outdir/{params.bw_folder}_final \
            &&  ~/environments/chipbinner/bin/Rscript scripts/merge.replicates.R \
                                                      --rep1 {input.rep1} --rep2 {input.rep2} \
                                                      --outdir outdir/{params.bw_folder}_final \
                                                      --outname {wildcards.condition}.{wildcards.mark} \
            &&  echo "merging for {input.rep1} and {input.rep2}"
            """
else:
    rule notmerged_bw_files:
        input:
            rep = lambda wildcards: f"outdir/{bigWig_folder}/{replace_suffix(condition_sample[wildcards.condition][0])}.binned.bw" if not config['general']['filter_low_counts'] else lambda wildcards: f"outdir/{bigWig_folder}/{replace_suffix(condition_sample[wildcards.condition][0])}.binned_filt.bw",
        output:
            bw_merged = "outdir/"+bigWig_folder+"_final/{condition}.{mark}.bw"
        params:
            organism = config['general']['organism'],
            mark = get_mark,
            condition = get_condition,
            bw_folder=bigWig_folder
        shell:
            """
                mkdir -p outdir/{params.bw_folder}_final \
            &&  cp {input.rep} {output.bw_merged} \
            &&  echo "not merging for {input.rep}"
            """

# def bigwigs_by_mark(wildcards, bw_folder):
#     all_bigwigs = os.listdir(bw_folder)
#     mark_bigwigs = list(filter(lambda x: wildcards.mark in x, all_bigwigs))
#     return mark_bigwigs

# rule PCA:
#     input:
#         bigwigs_mark = bigwigs_by_mark(wildcards, bw_folder="outdir/"+bigWig_folder+"_final")

if config['general']['MS_normalize']:
    rule compare:
        input:
            final_bw_files,
            bw1 = lambda wildcards: f"outdir/MSbigwigs_final/{get_treatment(wildcards)}.{get_mark(wildcards)}.bw",
            bw2 = lambda wildcards: f"outdir/MSbigwigs_final/{get_parental(wildcards)}.{get_mark(wildcards)}.bw",
        output:
            scatterplots = "outdir/comparisons/{comparison}/{comparison}_{mark}.scatterplot.pdf",
            preclust_matrix = "outdir/comparisons/{comparison}/{comparison}_{mark}_preclust_mat.csv",
            pooled_bed = "outdir/comparisons/{comparison}/{comparison}_{mark}_preclust_pooled.bed"
        resources: sequential=1
        params:
            treatment=get_treatment,
            parental=get_parental,
            organism=config['general']['organism']
        # wait:
        #     rules.merged_bw_files.output if config['general']['merge_replicates'] else rules.notmerged_bw_files.output
        shell:
            """
            mkdir -p outdir/comparisons \
            && echo "{input.bw1} // {input.bw2}" \
            && ~/environments/chipbinner/bin/Rscript scripts/genic.intergenic.plot.R \
                --bigWigWT {input.bw2} \
                --bigWigTest {input.bw1} \
                --mark {wildcards.mark} --pos_max 0.5 --organism {params.organism} \
                --wt_condition {params.parental} --test_condition {params.treatment} \
                --pow 0.75 --output_folder outdir/comparisons \
            && /home/andygg98/environments/snakemake/bin/python scripts/compare.genic.vs.intergenic.py --bigWigWT {input.bw2} \
                                                            --bigWigTest {input.bw1} \
                                                            --WTName {params.parental} \
                                                            --TreatmentName {params.treatment} \
                                                            --organism {params.organism} \
                                                            --mark {wildcards.mark} \
                                                            --output_folder outdir/comparisons
            """
else:
    rule compare:
        input:
            final_bw_files,
            bw1 = lambda wildcards: f"outdir/notMSbigwigs_final/{get_treatment(wildcards)}.{get_mark(wildcards)}.bw",
            bw2 = lambda wildcards: f"outdir/notMSbigwigs_final/{get_parental(wildcards)}.{get_mark(wildcards)}.bw"
        output:
            scatterplots = "outdir/comparisons/{comparison}/{comparison}_{mark}.scatterplot.pdf",
            preclust_matrix = "outdir/comparisons/{comparison}/{comparison}_{mark}_preclust_mat.csv",
            pooled_bed = "outdir/comparisons/{comparison}/{comparison}_{mark}_preclust_pooled.bed"
        resources: sequential=1
        params:
            treatment=get_treatment,
            parental=get_parental,
            organism=config['general']['organism']
        # wait:
        #     rules.merged_bw_files.output if config['general']['merge_replicates'] else rules.notmerged_bw_files.output
        shell:
            """
            mkdir -p outdir/comparisons \
            && echo "{input.bw1} // {input.bw2}" \
            && ~/environments/chipbinner/bin/Rscript scripts/genic.intergenic.plot.R \
                --bigWigWT {input.bw2} \
                --bigWigTest {input.bw1} \
                --mark {wildcards.mark} --pos_max 0.5 --organism {params.organism} \
                --wt_condition {params.parental} --test_condition {params.treatment} \
                --pow 0.75 --output_folder outdir/comparisons \
            && /home/andygg98/environments/snakemake/bin/python scripts/compare.genic.vs.intergenic.py --bigWigWT {input.bw2} \
                                                            --bigWigTest {input.bw1} \
                                                            --WTName {params.parental} \
                                                            --TreatmentName {params.treatment} \
                                                            --organism {params.organism} \
                                                            --mark {wildcards.mark} \
                                                            --output_folder outdir/comparisons
            """

rule clustering:
    input:
        preclust_matrix = "outdir/comparisons/{comparison}/{comparison}_{mark}_preclust_mat.csv",
        pooled_bed = "outdir/comparisons/{comparison}/{comparison}_{mark}_preclust_pooled.bed"
    output:
        cluster_info = "outdir/comparisons/{comparison}/clustering/clus.{comparison}_{mark}_clustered.100.100.txt",
        cluster2_rda = "outdir/comparisons/{comparison}/clustering/{comparison}_{mark}.2clusters.100.100.annotated_clusters.rda",
        cluster3_rda = "outdir/comparisons/{comparison}/clustering/{comparison}_{mark}.3clusters.100.100.annotated_clusters.rda"
    log:
        log_file = "outdir/comparisons/{comparison}/log.clus.{comparison}_{mark}.txt"
    shell:
        """

        for i in $(seq 100 300 1000 ); do for j in $(seq 100 300 1000 ); do python scripts/clus.py --samp {wildcards.comparison}_{wildcards.mark}_clustered --outdir outdir/comparisons/{wildcards.comparison}/clustering --csvFile  {input.preclust_matrix} --minpts $i --minsamps $j --logfile outdir/comparisons/{wildcards.comparison}/clustering/log.$i.$j.txt; ~/environments/chipbinner/bin/Rscript scripts/annotate.clusters.R --matrix {input.preclust_matrix} --pooled_bed_file {input.pooled_bed} --hdbscan_output_file outdir/comparisons/{wildcards.comparison}/clustering/clus.{wildcards.comparison}_{wildcards.mark}_clustered.$i.$j.txt --out_dir outdir/comparisons/{wildcards.comparison}/clustering; done; done

        """


if config['general']['MS_normalize']:
    rule density_plot:
        input:
            bw1 = lambda wildcards: f"outdir/MSbigwigs_final/{get_treatment(wildcards)}.{get_mark(wildcards)}.bw",
            bw2 = lambda wildcards: f"outdir/MSbigwigs_final/{get_parental(wildcards)}.{get_mark(wildcards)}.bw",

            cluster_info = "outdir/comparisons/{comparison}/clustering/clus.{comparison}_{mark}_clustered.100.100.txt",
            cluster2_rda = "outdir/comparisons/{comparison}/clustering/{comparison}_{mark}.2clusters.100.100.annotated_clusters.rda",
            cluster3_rda = "outdir/comparisons/{comparison}/clustering/{comparison}_{mark}.3clusters.100.100.annotated_clusters.rda"
        output:
            density_scatter = "outdir/comparisons/{comparison}/clustering/{comparison}_{mark}_2clusters.100.100.plots_density_scatterplot.pdf"
        params:
            organism=config['general']['organism'],
            treatment=get_treatment,
            parental=get_parental,
        shell:
            """
            for c in 2 3; do for i in outdir/comparisons/{wildcards.comparison}/clustering/*${{c}}clusters*rda; do ~/environments/chipbinner/bin/Rscript scripts/density.scatter.plot.R \
                                                        --out_name {wildcards.comparison}_{wildcards.mark}_${{c}}clusters  \
                                                        --out_dir outdir/comparisons/{wildcards.comparison}/clustering/ \
                                                        --organism {params.organism}  \
                                                        --bigWigWT {input.bw2} \
                                                        --bigWigTest {input.bw1} \
                                                        --mark {wildcards.mark} \
                                                        --annotated_clust $i \
                                                        --clust_num $c \
                                                        --wt_condition {params.parental} \
                                                        --test_condition {params.treatment} ; done; done

            """


# Run it like this:
# snakemake --cores 6 --resources sequential=1 --keep-going
# find outdir -type f -empty -delete
