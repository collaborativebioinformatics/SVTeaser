import os
import math
import inspect
import logging
import argparse
import subprocess
from collections import OrderedDict
from random import randint
import shutil

from acebinf import cmd_exe
from truvari import setup_logging
from svteaser.utils import vcf_compress
import pysam


def serialize_contigs_to_fa(contigs, fa_path):
    with open(fa_path, "w+") as fh:
        for contig, seq in contigs:
            fh.write(">{}\n".format(contig))
            fh.write("{}\n".format(seq))

def generate_altered_ref(ref_file, sv_vcf, outdir, copy_unaltered_contigs):
    """
    Generate altered ref sequence and return.
    """
    reference = pysam.FastaFile(ref_file)

    sv = pysam.VariantFile(sv_vcf)

    alt_contigs = {}

    contig = None
    contig_seq = None
    contig_pos = 0
    alt_seq = []

    for record in sv:
        if record.chrom != contig:
            # Store previously finished ref contig
            if contig is not None:
                alt_seq.append(contig_seq[contig_pos:])
                alt_contigs[contig] = "".join(alt_seq)

            # Reset variables to process new contig
            contig = record.chrom
            contig_seq = reference.fetch(contig)
            contig_pos = 0
            logging.info("Creating alt contig for {}".format(contig))

        assert(len(record.alts) == 1), "Cannot process multi allelic entries in VCF."

        ref = record.ref
        alt = record.alts[0]
        var_pos = record.pos - 1

        # For non variant positions, grab sequence from reference contig and increment ref seq pos.
        if var_pos > contig_pos:
            alt_seq.append(contig_seq[contig_pos:var_pos])
            contig_pos = var_pos

        # Add alt variant to alt sequence.
        alt_seq.append(alt)

        # Increment ref contig position by ref sequence.
        contig_pos += len(ref)

    if contig is not None:
        alt_seq.append(contig_seq[contig_pos:])
        alt_contigs[contig] = "".join(alt_seq)

    # Grab any non variant contigs in final list from original reference.
    final_contigs = []
    if copy_unaltered_contigs:
        for contig in reference.references:
            if contig not in alt_contigs:
                final_contigs.append((contig, reference.fetch(contig)))
            else:
                final_contigs.append((contig, alt_contigs[contig]))
    else:
        for contig, seq in alt_contigs.items():
            final_contigs.append((contig, seq))

    serialize_contigs_to_fa(final_contigs, os.path.join(outdir, "svteaser.altered.fa"))
    shutil.copyfile(ref_file, os.path.join(outdir, "svteaser.ref.fa"))
    if ".gz" in sv_vcf:
        shutil.copyfile(sv_vcf, os.path.join(outdir, "svteaser.sim.vcf.gz"))
    else:
        path = shutil.copyfile(sv_vcf, os.path.join(outdir, "svteaser.sim.vcf"))
        shutil.copyfile(sv_vcf, path)
        vcf_compress(path)


def known_sv_sim_main(args):
    """
    Run simulation on reference with known SVs. Output them into a directory.
    """
    args = parseArgs(args)

    logging.debug(f"Making outdir {args.output}")
    try:
        os.mkdir(args.output)
    except FileExistsError:
        logging.error(f"Output directory {args.output} already exists")
        exit(1)

    final_contigs = generate_altered_ref(args.reference,
                                         args.sv_vcf,
                                         args.output,
                                         args.copy_unaltered_contigs)

    logging.info("Finished")


def parseArgs(args):
    """
    Argument parsing
    """
    parser = argparse.ArgumentParser(prog="known_sv_sim", description=inspect.getdoc(known_sv_sim_main),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("reference", metavar="REF", type=str,
                        help="Reference file overwhich to simulate SVs")
    parser.add_argument("sv_vcf", metavar="SV_VCF", type=str,
                        help="VCF with known SVs to simulate. MUST BE SORTED.")
    parser.add_argument("output", metavar="OUT", type=str, default="output",
                        help="SVTeaser output basename (%(default)s)")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")
    parser.add_argument("--copy_unaltered_contigs", action="store_true",
                        help="Save both altered and UNaltered contigs from reference into altered fasta. Default to save only altered contigs.")
    args = parser.parse_args(args)
    args.reference = os.path.abspath(args.reference)
    args.output = args.output + ".svt"
    setup_logging(args.debug)
    return args
