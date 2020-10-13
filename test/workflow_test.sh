svteaser -h

# SURVIVOR SIM TEST
OUTDIR="./mito_test"
rm -rf ${OUTDIR}.svt/
svteaser surv_sim chrM.fa ${OUTDIR} --num_sv_regions 2 --debug
svteaser sim_reads ${OUTDIR}.svt

# KNOWN SV SIM TEST
OUTDIR="known_sv_sim_test"
rm -rf ${OUTDIR}.svt
svteaser known_sv chrM.fa sample_known_sv_chrM.vcf ${OUTDIR} --len_sv_region 100 --max_sv_size 25 --ref_seq_padding 10 --debug
svteaser sim_reads ${OUTDIR}.svt
