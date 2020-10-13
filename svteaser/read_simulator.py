import os
import inspect
import logging
import argparse

from truvari import setup_logging
from acebinf import cmd_exe

from svteaser.utils import check_gzip, check_samtools

def sim_reads_art(workdir, coverage=30, readlen=150, meanfrag=400, insertsd=50, instrument="HS25"):
    """
    Run art_illumina read simulator
    """
    ret = cmd_exe("which art_illumina")
    if ret.ret_code != 0:
        logging.error("Cannot find art_illumina executable in the environment")
        exit(ret.ret_code)
    try:
        os.chdir(workdir)
    except OSError:
        logging.error(f"Cannot change into {workdir} directory")
        exit(1)
    alt_ref = 'svteaser.altered.fa'

    outdir = "sim_reads_{}_{}_{}_{}_{}".format(coverage, readlen, meanfrag, insertsd, instrument)
    os.mkdir(outdir)
    # Useful when running on same altered reference but different parameters
    out_path = os.path.join(outdir, "art_illumina.simReads")
    ret = cmd_exe((f"art_illumina -ss {instrument} -sam -na -i {alt_ref} -p "
                   f"-l {readlen} -m {meanfrag} -s {insertsd} -f {coverage} -o {out_path}"))
    if ret.ret_code != 0:
        logging.error("Problem running art_illumina")
        logging.error(ret.stderr)
        logging.error(ret.stdout)
        exit(ret.ret_code)

    # Optionally compress fq
    if check_gzip():
        ret = cmd_exe((f"gzip {out_path}1.fq"))
        if ret.ret_code != 0:
            logging.info(f"Could not compress {out_path}1.fq")
        ret = cmd_exe((f"gzip {out_path}2.fq"))
        if ret.ret_code != 0:
            logging.info(f"Could not compress {out_path}2.fq")
    if check_samtools():
        ret = cmd_exe((f"samtools view -S -b {out_path}.sam > {out_path}.bam"))
        if ret.ret_code != 0:
            logging.info(f"Could not compress {out_path}2.fq")
        os.remove(f"{out_path}.sam")

def sim_reads_main(args):
    """
    Run read simulators
    """
    args = parseArgs(args)
    # Run the commands
    sim_reads_art(args.workdir,
                  coverage=args.coverage,
                  readlen=args.read_len,
                  meanfrag=args.mean_frag,
                  insertsd=args.insert_sd,
                  instrument=args.seq_inst)
    logging.info("Finished")

def parseArgs(args):
    """
    Argument parsing
    """
    parser = argparse.ArgumentParser(prog="sim_reads", description=inspect.getdoc(sim_reads_main),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("workdir", metavar="DIR", type=str,
                        help="SVTeaser working directory")
    parser.add_argument("--coverage", type=int, default=30,
                        help="Depth of coverage to simulate (%(default)s)")
    parser.add_argument("--read-len", type=int, default=150,
                        help="Simulated read length (%(default)s)")
    parser.add_argument("--mean-frag", type=int, default=400,
                        help="Mean insert fragment length (%(default)s)")
    parser.add_argument("--insert-sd", type=int, default=50,
                        help="Insert fragment length standard deviation (%(default)s)")
    parser.add_argument("--seq-inst", type=str, default="HS25",
                        help="Sequencing instrument (%(default)s)")
    parser.add_argument("--out-dir", type=str, required=False,
                        help="Output directory to save the results to. If unspecified, \
                              will save the results at DIR")
    args = parser.parse_args(args)
    setup_logging()
    return args

