"""
Misc tools
"""
from acebinf import cmd_exe

def vcf_compress(fn):
    """
    Run vcftools to sort/compress/index a vcf file
    """
    ret = cmd_exe(f"vcf-sort {fn} | bgzip > {fn}.gz && tabix {fn}.gz")

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

