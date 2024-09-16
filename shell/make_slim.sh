#!/bin/bash


K_burn_in=${1}
K_small_pop=${2}
K_bottleneck=${3}
K_small_pop_recover=${4}
K_big_pop_recover=${5}
year_rescue=${6}
type_deleterious=${7}
num_deleterious=${8}
rescue_sex=${9}
rescue_time=${10}
burn_in_time=${11}
out_dir=${12}

cat > ${out_dir}/panda_${K_burn_in}K_burn_in_${K_small_pop}K_small_pop_${K_bottleneck}K_bottleneck_${K_small_pop_recover}K_small_pop_recover_${rescue_time}rescue_time_${burn_in_time}burn_in_time_${year_rescue}year_rescue_${type_deleterious}type_deleterious_${num_deleterious}num_deleterious_${rescue_sex}rescue_sex.slim << EOM
initialize() {
#Confidentiality until review
}
EOM
