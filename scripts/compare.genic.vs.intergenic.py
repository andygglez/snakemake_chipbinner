import pandas as pd
import matplotlib.pyplot as plt
import bioframe
import pybedtools
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--bigWigWT', type=str)
parser.add_argument('--bigWigTest', type=str)
parser.add_argument('--WTName', type=str)
parser.add_argument('--TreatmentName', type=str)
parser.add_argument('--organism', type=str)
parser.add_argument('--mark', type=str)
parser.add_argument('--output_folder', type=str)

args = parser.parse_args()

def read_bigwig(path:str):
    chromosomes = ["chr"+str(c) for c in range(1,23)]
    chromosomes.extend(["chrX", "chrY"])
    df = pd.DataFrame({"chrom":[],"start":[],"end":[],"value":[]}, dtype=str)
    for c in chromosomes:
        try:
            bw = bioframe.read_bigwig(path, chrom=c)
            df = pd.concat([df, bw])
        except:
            pass
    df['bin'] = df['chrom']+"."+df['start'].astype(str)+"."+df['end'].astype(str)
    return df[['bin', 'value']]

parental = read_bigwig(args.bigWigWT)
treatment = read_bigwig(args.bigWigTest)

# Merge both signal
merged = pd.merge(parental, treatment, on='bin', how='inner')
merged_bins = pybedtools.BedTool.from_dataframe(merged['bin'].str.split(".", expand=True))


def make_10kb_bins(bedtool_obj):
    new_bins = merged_bins.intersect(bedtool_obj, wa=True, f=1)
    new_bins = new_bins.to_dataframe()
    return list(new_bins['chrom']+"."+new_bins['start'].astype(str)+"."+new_bins['end'].astype(str))

if args.organism=="mm10":
    gene_coordinates = pybedtools.BedTool("/home/andygg98/links/projects/rrg-majewski-ab/andygg98/references/chipbinner/mm10.gene.body.coordinates.bed")
    intergenic_coordinates = pybedtools.BedTool("/home/andygg98/links/projects/rrg-majewski-ab/andygg98/references/mm10/mm10.intergenic.bed")
if args.organism=="hg38":
    gene_coordinates = pybedtools.BedTool("/home/andygg98/links/projects/rrg-majewski-ab/andygg98/references/chipbinner/hg38.gene.body.coordinates.bed")
    intergenic_coordinates = pybedtools.BedTool("/home/andygg98/links/projects/rrg-majewski-ab/andygg98/references/hg38/hg38.intergenic.bed")


genic_bins = make_10kb_bins(gene_coordinates)
intergenic_bins = make_10kb_bins(intergenic_coordinates)


fig, ax = plt.subplots(1,3, figsize=(14, 4))

ax[0].scatter(merged['value_x'], merged['value_y'], s=1, c='silver', alpha=0.05, label="All bins")
ax[0].plot(merged['value_x'], merged['value_x'], c='k', lw=0.8)

ax[1].scatter(merged['value_x'], merged['value_y'], s=1, c='silver', alpha=0.05, label="All bins")
ax[1].scatter(merged.loc[merged['bin'].isin(intergenic_bins)]['value_x'], merged.loc[merged['bin'].isin(intergenic_bins)]['value_y'], s=1, c='red', alpha=0.05, label="Intergenic bins")
ax[1].plot(merged['value_x'], merged['value_x'], c='k', lw=0.8)

ax[2].scatter(merged['value_x'], merged['value_y'], s=1, c='silver', alpha=0.05, label="All bins")
ax[2].scatter(merged.loc[merged['bin'].isin(genic_bins)]['value_x'], merged.loc[merged['bin'].isin(genic_bins)]['value_y'], s=1, c='blue', alpha=0.05, label="Genic bins")
ax[2].plot(merged['value_x'], merged['value_x'], c='k', lw=0.8)

ax[0].set_xlabel(args.WTName)
ax[1].set_xlabel(args.WTName)
ax[2].set_xlabel(args.WTName)

ax[0].set_ylabel(args.TreatmentName)
ax[1].set_ylabel(args.TreatmentName)
ax[2].set_ylabel(args.TreatmentName)

l1 = ax[0].legend(markerscale=2)
l2 = ax[1].legend(markerscale=2)
l3 = ax[2].legend(markerscale=2)

for i in [l1,l2, l3]:
    for handle in i.legend_handles:
        handle.set_alpha(1)


plt.savefig(f"{args.output_folder}/{args.TreatmentName}_{args.WTName}/{args.TreatmentName}_{args.WTName}_{args.mark}.genic.vs.intergenic.png")