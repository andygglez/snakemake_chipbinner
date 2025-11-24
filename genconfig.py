
import os
import ruamel.yaml
import argparse
import re
import pandas as pd
import readline

readline.set_completer_delims(' \t\n=')
readline.parse_and_bind("tab: complete")

parser = argparse.ArgumentParser()
parser.add_argument('--mode', type=str, help='Analysis. Options: [rna-seq, chip-seq, atac-seq, hi-c, cutrun, chipbinner]')
args = parser.parse_args()

def path_exists(func):
    def wrapper(*args, **kwargs):
        completed=False
        while not completed:
            file_path = func(*args, **kwargs)
            if os.path.exists(file_path):
                completed=True
                return file_path
            else:
                print('Error: Path doesn\'t exist\n')
    
    return wrapper

@path_exists
def prompt_path(prompt):
    return input('\n'+prompt+'\n')

def available_option(func):
    def wrapper(*args, **kwargs):
        available=False
        while not available:
            option = func(*args, **kwargs)
            available_options = kwargs['available_options']
            if option in available_options:
                available=True
        return option
    return wrapper

@available_option
def prompt_option(prompt, available_options: list):
    return input('\n'+prompt+f" {available_options}:\n")


######################################################   Prompting   ######################################################                 

completed1 = False
if args.mode in ['cutrun', 'hi-c'] :
    single_or_pair = 'pe'
    completed1 = True
    
if not completed1:
    single_or_pair = prompt_option(prompt='Is your data single or paired ended', available_options=['se', 'pe'])

completed2=False
while not completed2:
    sample_folder = prompt_path(prompt='Indicate sample with folders:')
    samples = os.listdir(sample_folder)
    is_gzipped = True if samples[0].endswith('fastq.gz') else False
    sample_folder = sample_folder[:-1] if sample_folder[-1]== '/' else sample_folder
    if single_or_pair=='se':
        samples = list(set(map(lambda x: x[:-6], samples))) if not is_gzipped else list(set(map(lambda x: x[:-9], samples)))
    elif single_or_pair=='pe':
        samples = list(set(map(lambda x: x[:-8], samples))) if not is_gzipped else list(set(map(lambda x: x[:-11], samples)))
    completed2 = True if len(samples)>0 else False
    samples = sorted(samples)
print(f'\nTaking samples...\n')
print(samples)

organism = prompt_option(prompt='Indidicate organism', available_options=['mm10', 'hg38'])

if organism=='hg38':
    ref_genome='/home/andygg98/projects/rrg-majewski-ab/andygg98/references/hg38/hg38.fa'
    chromSizes='/home/andygg98/projects/rrg-majewski-ab/andygg98/references/hg38/hg38.chrom.sizes'
    chrom_list='/home/andygg98/projects/rrg-majewski-ab/andygg98/references/hg38/hg38.chroms'
    blacklist='/home/andygg98/projects/rrg-majewski-ab/andygg98/references/hg38/hg38.blacklist.bed'
if organism=='mm10':
    ref_genome = '/home/andygg98/projects/rrg-majewski-ab/andygg98/references/mm10/mm10.fa'
    chromSizes = '/home/andygg98/projects/rrg-majewski-ab/andygg98/references/mm10/mm10.chrom.sizes'
    chrom_list = '/home/andygg98/projects/rrg-majewski-ab/andygg98/references/mm10/mm10.chroms'
    blacklist = '/home/andygg98/projects/rrg-majewski-ab/andygg98/references/mm10/mm10.blacklist.bed'

# if args.mode not in ['chipbinner']:
#     # ref_genome = prompt_path("Indicate location for reference genome (.fa):")
#     match = re.search(r'\/(.*?)\.', ref_genome)
#     organism = match.group(1)
# else:
#     organism = prompt_option(prompt='Indidicate organism', available_options=['mm10', 'hg38'])

#####################################################   RNA-Seq    ###################################################

if args.mode == 'rna-seq':
    
    strandedness = prompt_option(prompt="Indicate strandedness of data:", available_options=['unstranded', 'fs', 'rs'])

    alignment_add_params = '--readFilesCommand zcat' if is_gzipped else ' '

    output = {
        'general':{
            'analysis': args.mode,
            'organism': organism,
            'single/pair-ended': single_or_pair,
            'strandedness':strandedness,
            'gzipped': is_gzipped
        },
        'inputs':{
            'sample_folder':sample_folder,
            'samples': samples,
        },
        'references':{
            'reference_genome': ref_genome,
            'reference_gtf': f"references/{organism}.gtf",
        },
        'alignment':{
            'threads':10,
            'add_parameters': alignment_add_params,
        },
        'featureCounts':{
            'threads': 6,
            'rm_suffix':' '
        },
        'DESeq2':{
            'applyShrinkage': True,
            'performEnrichment': False,
            'comparisons': ["treatment_control"]
        },
        'volcano':{
            'lfc_cutoff': 1,
            'pvalue_cutoff': 0.05
        }
    }

    if not os.path.exists("output/DGEA/DESeq2"):
        os.makedirs("output/DGEA/DESeq2")
    colInfo = pd.DataFrame({'Samples':samples,'Condition':['']*len(samples)})
    colInfo.to_csv("output/DGEA/DESeq2/colInfo.csv", index=False)

    if not os.path.exists(f"output/alignment/STAR_index/{organism}"):
        os.makedirs(f"output/alignment/STAR_index/{organism}")

    if not os.path.exists(f"references/{organism}.gtf"):
        print("""Annotation gtf file hasn't been found in folder 'references'.
        Please make sure to include it before running the pipeline.""")
    
    print("""\nPlease, always remember to:
            **1- Fill the output/DGEA/DESeq2/colInfo.csv table 
            **2- Check the configfile and remove suffixes if necessary
            """)


#####################################################   CUT&RUN    ###################################################

if args.mode == 'cutrun':

    spikeIn_fasta = prompt_path('Indicate location of spikeIn reference genome (.fa):')

    # chromSize = input('Indicate location of the chromSize.txt file:\n')
    # spikeIn_fasta = input('Indicate location of spikeIn reference genome (.fa):\n')
    match = re.search(r'/([^/]+)\.fa', spikeIn_fasta)
    spikeIn_organism = match.group(1)

    output = {
        'general':{
            'analysis': args.mode,
            'organism': organism,
            'spikeIn_organism': spikeIn_organism,
            'single/pair-ended': single_or_pair
        },
        'inputs':{
            'sample_folder': sample_folder,
            'samples': samples,
        },
        'references':{
            'reference_genome': ref_genome,
            'spikeIn_fasta': spikeIn_fasta,
            'chromSize': "references/"+organism+".chrom.sizes"
        },
        'alignment_parameters':{
            'alignment_threads':16,
            "ref_alnmt_params":"--local --no-mixed --no-discordant --phred33 -I 10 -X 700",
            "spikeIn_alnmt_params":"--local --no-mixed --no-discordant --phred33 -I 10 --no-overlap --no-dovetail -X 700",
            'sam_filt_quality': 10,
        },
        'MS_norm':{
            'perform': False,
            'MS_dict': 'references/MS.norm.yaml'
        },
        'peakCalling':{
            'seacr':'/home/andygg98/softwares/WGBS/bin/SEACR_1.3.sh',
            'peakCalling_percentage': 0.01
        }
    }

#####################################################   HiC    ###################################################

if args.mode == 'hi-c' or args.mode == 'hic':

    if organism=="hg38":
        assembly_arms="/lustre06/project/6007495/andygg98/references/hg38/hg38.arms.bed"
    if organism=="mm10":
        assembly_arms="/lustre06/project/6007495/andygg98/references/mm10/mm10.arms.bed"
    
    restrict = prompt_option(prompt='Did you use restriction enzymes?', available_options=['Y', 'N'])
    restrict = True if restrict=='Y' else False
    if restrict:
        digested_genome = prompt_path(prompt='Indicate bed file with digested genome:')
    else:
        digested_genome=None

    output = {
        'general':{
            'analysis': args.mode,
            'organism': organism,
            'single/pair-ended': single_or_pair,
            'gzipped': True,
            'restrict': restrict
        },
        'inputs':{
            'sample_folder': sample_folder,
            'samples': samples,
        },
        'references':{
            'reference_genome': ref_genome,
            'chromSizes': chromSizes,
            'chrom_list': chrom_list,
            'restriction_enzyme': digested_genome,
            'blacklist': blacklist,
            'assembly_arms': assembly_arms
        },
        'alignment':{
            'threads':16,
        }
    }

#################################################   chipbinner    ###################################################

if args.mode == 'chipbinner':

    if organism=="hg38":
        ref_blacklist = "/home/andygg98/projects/rrg-majewski-ab/andygg98/references/hg38/hg38.blacklist.bed"
        binning_windows = "/lustre06/project/6007495/andygg98/references/chipbinner/hg38.100kb.windows.bed"
    if organism=="mm10":
        ref_blacklist = "/home/andygg98/projects/rrg-majewski-ab/andygg98/references/mm10/mm10.blacklist.bed" 
        binning_windows = "/home/andygg98/projects/rrg-majewski-ab/andygg98/references/chipbinner/mm10.100kb.windows.bed"

    output = {
        'general':{
            'analysis': args.mode,
            'organism': organism,
            "MS_normalize": True,
            "merge_replicates": True,
            "filter_low_counts": False
        },
        'inputs':{
            'sample_folder': sample_folder,
        },
        'references':{
            'blacklist': ref_blacklist,
            'binning_windows': binning_windows,
            'MS_yaml': 'references/MS_coefficients.yaml'
        },
        'comparisons': [
            'treatment1_control1', 
            'treatment2_control2'
        ]
    }

    if not os.path.exists("config"):
        os.mkdir("config")
    pd.DataFrame({'Rep1':[],
                  'Rep2':[],
                  'Condition':[],
                  'Mark':[]}).to_csv("config/sample_list.csv", index=False)

#### At the end ####
if not os.path.exists("config"):
    os.mkdir("config")

yaml = ruamel.yaml.YAML()
with open('config/config.yaml', 'w') as config:
    yaml.dump(output, config)