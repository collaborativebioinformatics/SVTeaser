import glob
import os

from acebinf import cmd_exe

base = "hg002-chr2-reg-10kb-max-sv-len-4000.svt/HG002_2.vcf.gz"
bed = "hg002-chr2-reg-10kb-max-sv-len-4000.svt/include.bed"
ref = "/home/dnanexus/workdir/human_g1k_v37.fa"
for svcaller in ["breakseq", "lumpy", "manta", "breakdancer", "cnvnator"]:
    if svcaller != "manta": continue
    for vcf in glob.glob(f"*/*/*{svcaller}*.vcf.gz"):
        dirname = os.path.dirname(vcf)
        extra = "--pctsim 0" if svcaller != "manta" else ""
        outname = os.path.join(dirname, f"truv_{svcaller}")
        cmd = f"truvari bench -b {base} --includebed {bed} -c {vcf} -f {ref} -o {outname} {extra}"
        cmd2 = f"rm -rf {outname}"
        print(cmd)


