import joblib
import glob
from svteaser.utils import *
import pandas as pd

all_folders = []
for svcaller in ["breakseq", "lumpy", "manta", "breakdancer", "cnvnator"]:
    all_folders.extend(glob.glob(f"*/*/truv_{svcaller}"))

print(f"Found {len(all_folders)} folders")

varpieces = []
sumpieces = []
for i in all_folders:
    print(i)
    pdat = i.split('/')
    print(pdat)
    name = pdat[1]
    caller = pdat[-1].split('_')[1]
    dat = name.split('_')
    print(dat)
    cov = int(dat[2])
    readlen = int(dat[3])
    insert = int(dat[4])
    insd = int(dat[5])
    inst = dat[6]
    vdf, sdf = parse_truvari_dir(i)
    for k, v in zip(["caller", "coverage", "readlen", "insert", "insd", "inst"], [caller, cov, readlen, insert, insd, inst]):
        vdf[k] = v
        sdf[k] = v
    varpieces.append(vdf)
    sumpieces.append(sdf)

vfinal = pd.concat(varpieces)
sfinal = pd.concat(sumpieces)
joblib.dump(vfinal, "chr2.vars.final.jl")
joblib.dump(sfinal, "chr2.sums.final.jl")
