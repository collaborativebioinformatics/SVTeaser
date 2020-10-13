#!/usr/bin/env python
import sys
import argparse

from svteaser.vcfeditor import surv_vcf_fmt_main
from svteaser.surv_sim import surv_sim_main
from svteaser.read_simulator import sim_reads_main
from svteaser.known_sv_sim import known_sv_sim_main

VERSION="0.0.2"

def in_progress(args):
    """placeholder"""
    print('working on it...')

def version(args):
    """Print the version"""
    print("svteaser v%s" % VERSION)

TOOLS = {'surv_sim': surv_sim_main,
         'known_sv': known_sv_sim_main,
         'sim_reads': sim_reads_main,
        }

USAGE = """\
SVTeaser v%s - SV read simulation for rapid benchmarking

    CMDs:
        known_sv        Create genome regions from a VCF of known SVs
        surv_sim        Simulate random SVs with SURVIVOR
        sim_reads       Run read simulators
""" % VERSION

def parseArgs():
    """
    Argument parsing
    """
    parser = argparse.ArgumentParser(prog="svteaser", description=USAGE,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("cmd", metavar="CMD", choices=TOOLS.keys(), type=str,
                        help="Command to execute")
    parser.add_argument("options", metavar="OPTIONS", nargs=argparse.REMAINDER,
                        help="Options to pass to the command")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    TOOLS[args.cmd](args.options)

if __name__ == '__main__':
    parseArgs()



