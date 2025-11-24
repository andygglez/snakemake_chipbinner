#!/usr/bin/env python

import sys
import hdbscan
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import logging
def parse_arguments():
    parser = argparse.ArgumentParser(description="Getting clusters using HDBSCAN")
    parser.add_argument('--samp',type=str,help='prefix for output',required=True)                           
    parser.add_argument('--outdir',type=str,help='output directory',required=True)                           
    parser.add_argument('--csvFile',type=str,help='{sample}.csv derived from 10kb normalized bw files. It contains two columns of log2 scores.',required=True)
    parser.add_argument('--minpts',type=int,help='The minimum size of clusters; single linkage splits that contain fewer points than this will be considered points “falling out” of a cluster rather than a cluster splitting into two new clusters.',required=True)
    parser.add_argument('--minsamps',type=int,help='The number of samples in a neighbourhood for a point to be considered a core point.',required=True)
    parser.add_argument('--xlab',type=str,help='Label for x-axis',default='KO')
    parser.add_argument('--ylab',type=str,help='Label for y-axis',default='PA')
    parser.add_argument('--cluster_selection_method',type=str,help='A string specifying the cluster selection method. Must be one of "eom" or "leaf"',default='eom')
    parser.add_argument('--logfile',type=str,required=True,help='logfile for debugging: example.log')
    args = parser.parse_args()
    return args
    
def main():
# samp, minpts, minsamps = sys.argv[1:] 
# d = np.genfromtxt('misc/%s.csv' % samp, delimiter=',')
    args = parse_arguments()
    logging.basicConfig(filename=args.logfile, level=logging.DEBUG)
    logging.info('Starting dbscan')
    d = np.genfromtxt(args.csvFile, delimiter=',')
    #
    clusterer = hdbscan.HDBSCAN(min_cluster_size = args.minpts, min_samples = args.minsamps, core_dist_n_jobs = 1,cluster_selection_method = args.cluster_selection_method).fit(d)
    color_palette = sns.color_palette('deep', 8)
    cluster_colors = [color_palette[x] if x >= 0
                    else (0.5, 0.5, 0.5)
                    for x in clusterer.labels_]
    cluster_member_colors = [sns.desaturate(x, p) for x, p in
                            zip(cluster_colors, clusterer.probabilities_)]

    if args.outdir[-1]!="/":
        args.outdir = args.outdir+"/"

    np.savetxt(args.outdir+'clus.%s.%s.%s.txt' % (args.samp, args.minpts, args.minsamps), clusterer.labels_.astype(int), fmt = '%i')
    plt.scatter(*d.T, s=50, linewidth=0, c=cluster_member_colors, alpha=0.25)
    plt.xlabel(args.xlab)
    plt.ylabel(args.ylab)
    plt.savefig(args.outdir+'clus.%s.%s.%s.png' % (args.samp, args.minpts, args.minsamps), transparent=True)
    logging.info('Hdbscan done')

if __name__ == '__main__':
    main()
