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
from svteaser.utils import vcf_compress, add_fasta_entry
import pysam


def serialize_contigs_to_fa(contigs, fa_path):
    with open(fa_path, "w+") as fh:
        for contig, seq in contigs:
            fh.write(">{}\n".format(contig))
            fh.write("{}\n".format(seq))

def generate_altered_ref(ref_file, sv_vcf, outdir, copy_unaltered_contigs):
    """
    Generate altered ref sequence and return. This spikes in all the variants into a single altered seq
    and outputs all altered sequences.

    NOTE: Unused right now.
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


def generate_altered_regions(ref_file, sv_vcf, outdir, region_size, max_sv_size, padding=0):
    """
    Simulate variants from known SVs by spiking them into reference segments.
    """
    logging.info(f"Region size = {region_size}, Max SV Size = {max_sv_size}, Padding = {padding}")
    reference = pysam.FastaFile(ref_file)
    sv = pysam.VariantFile(sv_vcf)
    header = sv.header

    out_ref_path = os.path.join(outdir, "svteaser.ref.fa")
    out_ref_fh = open(out_ref_path, "w+")

    out_altered_path = os.path.join(outdir, "svteaser.altered.fa")
    out_altered_fh = open(out_altered_path, "w+")

    # Number bases to flank on either size ov variant
    flank_size = region_size // 2

    if flank_size - max_sv_size < padding:
        logging.error(f"Max SV size and required padding are not compatible. Increase region size by {2 * (padding + max_sv_size - flank_size)} or decreases max SV size to {flank_size - padding}.")
        exit(1)


    records = []

    last_chrom = ""
    chrom_seq = ""

    for record in sv:
        chrom = record.chrom
        if last_chrom != chrom:
            logging.debug("Load new chrom {}".format(chrom))
            last_chrom = chrom
            chrom_seq = reference.fetch(chrom)

        pos = record.pos - 1
        ref = record.ref
        alt = record.alts[0]

        if abs(len(alt) - len(ref)) > max_sv_size:
            logging.debug(f"Skip variations longer than {max_sv_size}")
            continue

        start_pos = max(0, pos - flank_size)
        end_pos = min(pos + flank_size, len(chrom_seq))

        ref_seq = chrom_seq[start_pos:end_pos]
        relative_pos = pos - start_pos
        alt_seq = "".join([ref_seq[:relative_pos], alt, ref_seq[relative_pos + len(ref):]])
        new_contig_name = f"{chrom}_{start_pos}_{end_pos}"

        header.add_line(f"##contig=<ID={new_contig_name},length={len(ref_seq)}>")

        add_fasta_entry(new_contig_name, ref_seq, out_ref_fh)
        add_fasta_entry(new_contig_name, alt_seq, out_altered_fh)

        # NOTE: We don't update variant record to keep variants in original coordinate frame.
        # Keeping for now in case needed, but should be removed.
        #record.chrom = new_contig_name
        #record.pos = relative_pos + 1

        records.append(record)

    out_altered_fh.close()
    out_ref_fh.close()

    out_vcf_path = os.path.join(outdir, "svteaser.sim.vcf")
    with pysam.VariantFile(out_vcf_path, "w", header=sv.header) as out_vcf_fh:
        for rec in records:
            out_vcf_fh.write(rec)

    vcf_compress(out_vcf_path)

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

    generate_altered_regions(args.reference,
                             args.sv_vcf,
                             args.output,
                             args.len_sv_region,
                             args.max_sv_size,
                             padding=args.ref_seq_padding)

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
    parser.add_argument('--len_sv_region', type=int, default=10000,
                        help='The length of regions to create.',
                        required=False)
    parser.add_argument('--max_sv_size', type=int, default=4000,
                        help='Max length of variations to spike.',
                        required=False)
    parser.add_argument('--ref_seq_padding', type=int, default=800,
                        help='Padded region around each end of reg where variation are spiked.',
                        required=False)
    args = parser.parse_args(args)
    args.reference = os.path.abspath(args.reference)
    args.output = args.output + ".svt"
    setup_logging(args.debug)
    return args
