
# Bash script to run parameter swap over multiple params to generate simulation reads
function print_help {
    echo "This script takes svteaser surv_sim output and runs a parameter"
    echo "sweep to generate simulated reads."
    echo "Run command : bash parameter_sweep.sh <surv_sim output dir>"
    exit
}

if [ $# == 0 ]; then
    print_help
    exit
fi

surv_sim_dir=$1

coverage=(10 20 30)
insert_frag=(400 600)

for cov in ${coverage[@]}; do
    for insert in ${insert_frag[@]}; do
        svteaser sim_reads --coverage $cov --mean-frag $insert $surv_sim_dir  
    done
done
