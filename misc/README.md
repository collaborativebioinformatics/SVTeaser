Generating mass benchmarking

After you have an svteaser.svt output directory (presumably with multiple read simulations),
you'll want to benchmark them against the base set of variants with truvari.

All of these tools assume you have run parliment2 and the results of it are saved 
in the directory beside the simulated reads

Step 1. sort/compress/index all of the vcf with

  $ find svteaser.svt -name "*.vcf" | python compress.py

Step 2.  Run truvari
You'll need to edit the base-vcf to have your known or simulated SVs
the ref is your svteaser.svt reference from which you inserted SVs and from which
you called your SVs
bed is an include.bed - which should be optional, but isn't currently..

Step 3. use Step1 to build sort/compress/index the truvari output vcfs

Step 4. Build your dataframe
  $ python build_df.py


