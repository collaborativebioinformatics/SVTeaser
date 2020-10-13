"""
Misc tools
"""
import os
import json

import pysam
import pandas as pd
import truvari

from acebinf import cmd_exe

from pandas.api.types import CategoricalDtype
SZBINTYPE = CategoricalDtype(categories=truvari.SZBINS, ordered=True)

def vcf_compress(fn):
    """
    Run vcftools to sort/compress/index a vcf file
    """
    ret = cmd_exe(f"vcf-sort {fn} | bgzip > {fn}.gz && tabix {fn}.gz")


def parse_truvari_dir(trudir):
    """
    creates a dataframe from all the vcfs in the Truvari directory
    loads the performance vcf
    returns vcf_dataframe, perf_dataframe
    """
    rows = []
    for state in ['fp', 'fn', 'tp-base', 'tp-call']:
        file = os.path.join(trudir, state + '.vcf.gz')
        v = pysam.VariantFile(file)
        for entry in v:  
            truscore = entry.info['TruScore'] if 'tp' in state else ''
            seq_sim = entry.info['PctSeqSimilarity'] if 'tp' in state else ''
            size_sim = entry.info['PctSizeSimilarity'] if 'tp' in state else ''
            rec_overlap = entry.info['PctRecOverlap'] if 'tp' in state else ''
            start_dist = entry.info['StartDistance'] if 'tp' in state else ''
            end_dist = entry.info['EndDistance'] if 'tp' in state else ''
            size_diff = entry.info['SizeDiff'] if 'tp' in state else ''
            num_neigh = entry.info['NumNeighbors'] if 'tp' in state else ''
            num_thresh_neigh = entry.info['NumThresholdNeighbors'] if 'tp' in state else ''
            rows.append([state,
                         truscore,
                         seq_sim,
                         size_sim,
                         rec_overlap,
                         start_dist,
                         end_dist,
                         size_diff,
                         num_neigh,
                         num_thresh_neigh,
                         truvari.entry_variant_type(entry),  
                         truvari.entry_boundaries(entry)[0],
                         truvari.entry_boundaries(entry)[1],
                         truvari.entry_size(entry)])
    df = pd.DataFrame(rows, columns=["state",
                                     "truscore",
                                     "seq_sim",
                                     "size_sim",
                                     "rec_overlap",
                                     "start_dist",
                                     "end_dist",
                                     "size_diff",
                                     "num_neigh",
                                     "num_thresh_neigh",
                                     "svtype",
                                     "start", 
                                     "end", 
                                     "svlen"])
    df['szbin'] = df['svlen'].apply(truvari.get_sizebin)
    df['szbin'] = df['szbin'].astype(SZBINTYPE)
    df["cnt"] = 1   
    perf = pd.DataFrame.from_dict(json.load(open(os.path.join(trudir, "summary.txt"))), orient='index')
    return df, perf.T

def check_gzip():
    """
    Check for presence of gzip.
    """
    ret = cmd_exe((f"gzip --help"))
    return ret.ret_code == 0

def check_samtools():
    """
    Check for presence of samtools.
    """
    ret = cmd_exe((f"samtools --help"))
    return ret.ret_code == 0

def add_fasta_entry(name, seq, fasta_fh):
    """
    Add new sequence to fasta file handle.
    """
    fasta_fh.write(">{}\n".format(name))
    fasta_fh.write("{}\n".format(seq))
    fasta_fh.flush()

