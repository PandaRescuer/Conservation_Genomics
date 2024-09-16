#!/bin/bash

#SLiM
#only improve carrying capacity for population with 10 individuals
${script_dir}/make_slim.sh 17679 10 200 10 2000 0 none 0 Balance 50 50000 ${slim_out_dir}
${script_dir}/make_slim.sh 17679 10 200 30 2000 0 none 0 Balance 50 50000 ${slim_out_dir}
${script_dir}/make_slim.sh 17679 10 200 50 2000 0 none 0 Balance 50 50000 ${slim_out_dir}
#only improve carrying capacity for population with 30 individuals
${script_dir}/make_slim.sh 17679 30 200 30 2000 0 none 0 Balance 50 50000 ${slim_out_dir}
${script_dir}/make_slim.sh 17679 30 200 50 2000 0 none 0 Balance 50 50000 ${slim_out_dir}
${script_dir}/make_slim.sh 17679 30 200 100 2000 0 none 0 Balance 50 50000 ${slim_out_dir}
#genetic rescue with improved carrying capacity for population with 10 individuals
${script_dir}/make_slim.sh 17679 10 200 10 2000 1 none 0 Balance 50 50000 ${slim_out_dir}
${script_dir}/make_slim.sh 17679 10 200 10 2000 2 none Balance 50 50000 ${slim_out_dir}
${script_dir}/make_slim.sh 17679 10 200 10 2000 5 none Balance 50 50000 ${slim_out_dir}
${script_dir}/make_slim.sh 17679 10 200 10 2000 10 none 0 Balance 50 50000 ${slim_out_dir}
#genetic rescue with random sex and improved carrying capacity for population with 10 individuals
${script_dir}/make_slim.sh 17679 10 200 10 2000 1 none 0 random 50 50000 ${slim_out_dir}
${script_dir}/make_slim.sh 17679 10 200 10 2000 2 none 0 random 50 50000 ${slim_out_dir}
${script_dir}/make_slim.sh 17679 10 200 10 2000 5 none 0 random 50 50000 ${slim_out_dir}
${script_dir}/make_slim.sh 17679 10 200 10 2000 10 none 0 random 50 50000 ${slim_out_dir}

#genetic rescue with addtional strongly deleterious mutaions and improved carrying capacity for population with 10 individuals
# weakly deleterious mutations
${script_dir}/make_slim.sh 17679 10 200 50 2000 2 m1 100 Balance 50 50000 ${slim_out_dir}
${script_dir}/make_slim.sh 17679 10 200 50 2000 2 m1 200 Balance 50 50000 ${slim_out_dir}
${script_dir}/make_slim.sh 17679 10 200 50 2000 2 m1 500 Balance 50 50000 ${slim_out_dir}

${script_dir}/make_slim.sh 17679 10 200 50 2000 5 m1 100 Balance 50 50000 ${slim_out_dir}
${script_dir}/make_slim.sh 17679 10 200 50 2000 5 m1 200 Balance 50 50000 ${slim_out_dir}
${script_dir}/make_slim.sh 17679 10 200 50 2000 5 m1 500 Balance 50 50000 ${slim_out_dir}

# moderately deleterious mutations
${script_dir}/make_slim.sh 17679 10 200 50 2000 2 m2 10 Balance 50 50000 ${slim_out_dir}
${script_dir}/make_slim.sh 17679 10 200 50 2000 2 m2 20 Balance 50 50000 ${slim_out_dir}
${script_dir}/make_slim.sh 17679 10 200 50 2000 2 m2 50 Balance 50 50000 ${slim_out_dir}

${script_dir}/make_slim.sh 17679 10 200 50 2000 5 m2 10 Balance 50 50000 ${slim_out_dir}
${script_dir}/make_slim.sh 17679 10 200 50 2000 5 m2 20 Balance 50 50000 ${slim_out_dir}
${script_dir}/make_slim.sh 17679 10 200 50 2000 5 m2 50 Balance 50 50000 ${slim_out_dir}

#strongly deleterious mutations
${script_dir}/make_slim.sh 17679 10 200 10 2000 2 m3 5 Balance 50 50000 ${slim_out_dir}
${script_dir}/make_slim.sh 17679 10 200 10 2000 2 m3 10 Balance 50 50000 ${slim_out_dir}
${script_dir}/make_slim.sh 17679 10 200 10 2000 2 m3 50 Balance 50 50000 ${slim_out_dir}

${script_dir}/make_slim.sh 17679 10 200 10 2000 5 m3 5 Balance 50 50000 ${slim_out_dir}
${script_dir}/make_slim.sh 17679 10 200 10 2000 5 m3 10 Balance 50 50000 ${slim_out_dir}
${script_dir}/make_slim.sh 17679 10 200 10 2000 5 m3 50 Balance 50 50000 ${slim_out_dir}

#genetic rescue with addtional recessive lethal mutations and improved carrying capacity for population with 10 individuals
${script_dir}/make_slim.sh 17679 10 200 10 2000 2 m5 1 Balance 50 50000 ${slim_out_dir}
${script_dir}/make_slim.sh 17679 10 200 10 2000 2 m5 2 Balance 50 50000 ${slim_out_dir}
${script_dir}/make_slim.sh 17679 10 200 10 2000 2 m5 5 Balance 50 50000 ${slim_out_dir}
${script_dir}/make_slim.sh 17679 10 200 10 2000 2 m5 10 Balance 50 50000 ${slim_out_dir}

${script_dir}/make_slim.sh 17679 10 200 10 2000 5 m5 1 Balance 50 50000 ${slim_out_dir}
${script_dir}/make_slim.sh 17679 10 200 10 2000 5 m5 2 Balance 50 50000 ${slim_out_dir}
${script_dir}/make_slim.sh 17679 10 200 10 2000 5 m5 5 Balance 50 50000 ${slim_out_dir}
${script_dir}/make_slim.sh 17679 10 200 10 2000 5 m5 10 Balance 50 50000 ${slim_out_dir}


#Loop command
slim_list=($(ls ${slim_out_dir}/*.slim))
echo ${slim_list}
for ((j=0;j<${#slim_list[*]};j+=1));
do
{
	slim_script=${slim_list[$j]}
	for ((i=0;i<50;i+=1));
	do
	{
		slim /${slim_script} >  ${slim_script}_log_${i}.txt
		awk '/^gen/,0' ${slim_script}_log_${i}.txt > ${slim_script}_log_${i}_keep.txt
	}&
	done
	wait
	echo "finished  -----> ${slim_script}"
}
done
