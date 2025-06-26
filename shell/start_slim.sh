#!/bin/bash

# Burn-in in and first decline
${script_dir}/make_slim_burn.sh ${make_slim_args} ${slim_out_dir}
slim ${slim_out_dir}/${filename}.slim >  ${slim_out_dir}/${filename}_log_${i}.txt


# Continue the simulation
function run_slim(){
	echo "start slim :${slim_out_dir}"
	for ((i=0;i<25;i+=1));
	do
	{
		${script_dir}/make_slim_continue.sh ${make_slim_args} ${slim_burn_in_dir}/${slim_file} ${slim_out_dir}
		slim ${slim_out_dir}/${filename}.slim >  ${slim_out_dir}/${filename}_log_${i}.txt
		awk '/^gen/,0' ${slim_out_dir}/${filename}_log_${i}.txt > ${slim_out_dir}/${filename}_log_${i}_keep.txt
	}&
	done
	wait
	echo "finished"
}

run_slim

