######################################
#                                    #
# Author : Joyjit Daw                #
# Email : jdaw@nvidia.com            #
#                                    #
######################################

import re
import sys
import logging
import argparse

from tempfile import NamedTemporaryFile

import pysam
import truvari

def recalibrate_vcf(ref_path, orig_vcf_path, out_vcf_path):
    """
    Re-calibrate positions of input VCF to be relative to original reference.
    """
    ref_file = pysam.FastaFile(ref_path)
    orig_vcf = pysam.VariantFile(orig_vcf_path)
    header = orig_vcf.header
 
    for chrom in ref_file.references:
        length = ref_file.get_reference_length(chrom)
        #chrom = chrom.replace("chr", "")
        header.add_line(f"##contig=<ID={chrom},length={length}>")

    with pysam.VariantFile(out_vcf_path, "w", header=header) as writer:
        for rec in orig_vcf:
            chrom = rec.chrom
            chrom, start, end = chrom.split("_")
            local_pos = rec.pos
            global_pos = int(start) + local_pos
            rec.chrom = chrom
            rec.pos = global_pos
            writer.write(rec)

def correct_survivor_vcf(in_vcf):
    """
    Correct survivor vcf mistakes so it's parsable by pysam.VariantFile
    Returns the name of the temporary file that's c
    """
    logging.debug("Correcting")
    extra_header = "\n".join(['##FILTER=<ID=LowQual,Description="Default. Manual">',
                          '##INFO=<ID=PRECISE,Number=1,Type=Flag,Description="Some type of flag">'])

    temp_file = NamedTemporaryFile(suffix=".vcf", mode='w', delete=False) 
    n_entries = 0
    with open(in_vcf, 'r') as fh:
        for line in fh:
            if line.startswith("##"):
                temp_file.write(line)
                continue
            if line.startswith("#CHROM"):
                temp_file.write(extra_header + '\n')
                line = line.strip() + "\tSAMPLE\n"
                temp_file.write(line)
                continue
            n_entries += 1
            line = re.sub(":GL:GQ:FT:RC:DR:DV:RR:RV", "", line)
            line = re.sub("LowQual", ".", line)
            temp_file.write(line)
    logging.debug("Corrected %d entries", n_entries)
    temp_file.close()
    return temp_file.name

def update_vcf(ref, insertions, survivor_vcf, out_vcf, pos_padding=0):
    """Update the SURVIVOR VCF file to have ref and alt sequences for each variant entry.

    e.g. If a variant entry has the following VCF description

    "chr1   10  INS001  N   <INS>   .   LowQual SVLEN=10"

    Then the entry will be updated with data from ref and insertions fasta to look like

    "chr1   10  INS001  A   ATTTTTTTTTTGGGGGGGGGG   .   LowQual SVLEN=10"

    Args:
        ref : Path to reference fasta file.
        insertions : Path to SURVIVOR insertions fasta file.
        survivor_vcf : Path to SURVIVOR simulated VCF file.
        out_vcf : Putput path for updated SURVIVOR VCF.
        pos_padding : Padding for start position in VCF.
    """
    survivor_vcf = correct_survivor_vcf(survivor_vcf)
    ref = pysam.FastaFile(ref)
    try:
        insertions = pysam.FastaFile(insertions)
    except OSError: # Sometimes there are no insertions?
        insertions = None
        pass
    
    vcf_reader = pysam.VariantFile(survivor_vcf)
    header = vcf_reader.header
    vcf_writer = pysam.VariantFile(out_vcf, 'w', header=header)
    n_entries = 0
    for record in vcf_reader:
        n_entries += 1
        record = truvari.copy_entry(record, header)
        chrom = record.chrom
        vcf_pos = record.pos # Position here is the VCF position, which is without padding.
        ref_pos = record.pos + pos_padding # Reference pos is VCF pos shifted by padding.
        if record.id.startswith("INS"):
            # Handle an INSERTION entry
            record.ref = ref.fetch(chrom, ref_pos, ref_pos + 1)
            survivor_insertion_key = "{}_{}".format(chrom, vcf_pos)
            record.alts = ["{}{}".format(record.ref, insertions.fetch(survivor_insertion_key))]
        elif record.id.startswith("DEL"):
            # Handle a DELETION entry
            svlen = record.info['SVLEN']
            record.ref = ref.fetch(chrom, ref_pos - 1, (ref_pos - 1) + svlen + 1)
            record.alts = [ref.fetch(chrom, ref_pos - 1, ref_pos)]
        else: # just in case inversions or something get through
            continue
        # Update the VCF position to reflect padded sequence
        record.pos = ref_pos
        vcf_writer.write(record)
    logging.info("Updated %d entries", n_entries)

def parse_args(args):
    """Build parser object with options for sample.

    Returns:
        Python argparse parsed object.
    """
    parser = argparse.ArgumentParser(
        description="A VCF editing utility which adds ref and all sequences to a SURVIVOR fasta file.")

    parser.add_argument("--reference-fasta", "-r", required=True, type=str,
                        help="Reference fasta file.")
    parser.add_argument("--survivor-insertions-fasta", "-i", required=True, type=str,
                        help="Insertions fasta file from SURVIVOR.")
    parser.add_argument("--survivor-vcf-file", "-v", required=True, type=str,
                        help="VCF file from SURVIVOR.")
    parser.add_argument("--output-vcf", "-o", required=True, type=str,
                        help="Output path of edited VCF.")
    parser.add_argument("--debug", action="store_true", 
                        help="Verbose logging")
    args = parser.parse_args(args)
    truvari.setup_logging(args.debug)
    return args


def surv_vcf_fmt_main(args):
    """
    Main entry function to the tool
    """
    args = parse_args(args)
    update_vcf(args.reference_fasta,
               args.survivor_insertions_fasta,
               args.survivor_vcf_file,
               args.output_vcf)
