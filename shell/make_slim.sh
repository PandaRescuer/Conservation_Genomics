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
	
	initializeSLiMModelType("nonWF"); 

	defineConstant("burn_in_time", ${burn_in_time});				// burn-in 
	defineConstant("decline_time", 9000);				// decline 
	defineConstant("bottleneck_time", 2014-1914);	// bottleneck 100
	defineConstant("rescue_time", ${rescue_time});				// the duration of rescue 
	defineConstant("K_burn_in", ${K_burn_in});			// burn-in, carrying capacity
	defineConstant("K_declien", 1776);			// burn-in, carrying capacity
	defineConstant("K_small_pop", ${K_small_pop}); 				// small poppulation, bottleneck, carrying capacity
	defineConstant("K_bottleneck", ${K_bottleneck}); 				// large population, bottleneck, carrying capacity
	defineConstant("K_small_pop_recover", ${K_small_pop_recover});			// small poppulation, recover, carrying capacity
	defineConstant("K_big_pop_recover", ${K_big_pop_recover}); 			// national park, carrying capacity
	
	defineConstant("year_rescue", ${year_rescue}); 				// rescue every # years
	defineConstant("type_deleterious", "${type_deleterious}"); 			// type of additional deleterious variantion from rescue	
	defineConstant("num_deleterious", ${num_deleterious}); 			// # number of additional deleterious variantion from rescue
	defineConstant("rescue_sex", "${rescue_sex}");				// sex of inds used for rescue
	
	defineConstant("sampleSize", 30); 				// sample to stats
	defineConstant("g",20000); 					// number of genes
	defineConstant("ROHcutoff", 100000);				// ROH to stats >1Mb
	defineConstant("geneLength", 1000); 				// difine gene length
	defineConstant("p_repr", 0.625); 				// annual probability of an adult female reproducing, wei,1994 
	defineConstant("seqLength", g*geneLength);			// total gene length to caculate het
	initializeSex("A"); 						// Sexual model
	defineConstant("L", c(0.6, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0,0.0,0.0,0.0,0.0, 0.1, 0.2, 0.3, 0.3, 0.5, 1.0)); // age-dependent increases in mortality
	
	initializeMutationRate(1.29e-8); 	//zhao, 2013
	
	defineConstant("h_wkDel", 0.45);	//Kyriazis, 2024
	defineConstant("h_modDel", 0.2);
	defineConstant("h_strDel", 0.05);
	defineConstant("h_semiLet", 0.0);
	defineConstant("h_let", 0.0);
	
	// set up discrete DFE with four mutation types coming from gamma DFE
	// augmented with recessive lethals
	// this approach for implementing an h-s relationship is much faster than using fitness callbacks in SLiM (see manual)
	initializeMutationType("m1", h_wkDel, "s", "do x=rgamma(1,-0.01314833,0.186); while (x < -0.001); return x;");
	initializeMutationType("m2", h_modDel, "s", "do x=rgamma(1,-0.01314833,0.186); while (x < -0.01 | x >= -0.001); return x;");
	initializeMutationType("m3", h_strDel, "s", "do x=rgamma(1,-0.01314833,0.186); while (x < -0.1 | x >= -0.01); return x;");
	initializeMutationType("m4", h_semiLet, "s", "do x=rgamma(1,-0.01314833,0.186); while (x >= -0.1); return x;");
	initializeMutationType("m5", h_let,"f", -1.0);

	// proportion of new deleterious mutations that are recessive lethal
	defineConstant("let_frac", 0.005);
	
	// set proportion of each mutation type as determined by Kim 2017 DFE augmented with lethals
	initializeGenomicElementType("g1", c(m1,m2,m3,m4,m5), c(0.491*(1-let_frac), 0.247*(1-let_frac), 0.236*(1-let_frac), 0.026*(1-let_frac), let_frac));
	
	// approach for setting up genes on different chromosomes adopted from Jacqueline's wolf scripts 
	
	// vector of # genes on 16 different crocodile lizard chromosomes - need to rescale according to number of desired genes
	gene_nums=c(1957,1838,1358,1332,1205,1211,1302,1189,954,1017,1017,752,851,981,843,840,389,351,328,285);
	
	
	for (i in 1:g){
		initializeGenomicElement(g1, ((i-1)*geneLength)+(i-1), (i*geneLength)+(i-2) );
	}
	
	rates=NULL;
	// assume no recombination within genes, a rate of 1e-3 between genes, and free recombination between chroms
	// Multiple chromosomes:
	for (i in 1:(size(gene_nums)-1)){
		rates=c(rates, 0, rep(c(1e-3, 0), asInteger(gene_nums[i-1]-1)), 0.5);
	}
	rates=c(rates, 0, rep(c(1e-3, 0), asInteger(gene_nums[size(gene_nums)-1]-1)));
	
	ends=NULL;
	for (i in 1:g){
		ends=c(ends, (i*geneLength)+(i-2), (i*geneLength)+(i-1));
	}
	ends=ends[0:(size(ends)-2)];
	
	initializeRecombinationRate(rates, ends);

}

//Kyriazis, 2024
reproduction() {
	
	for(pop in sim.subpopulations){
		
		males = pop.individuals[pop.individuals.sex=='M'];
		females = pop.individuals[pop.individuals.sex=='F'];
		
		//get reproductive age males and females
		repro_females = females[females.age >= 7 & females.age <= 20 ];
		repro_males = males[males.age >= 8 & males.age <= 20 ];
		
		//loop through females and reproduce
		//with one randomly selected male (1 offspring)
		for(mom in repro_females){
			
			// get fitness of mom and use to weight probability of reproduction
			// i.e., if fitness is 0.95, reproduction will be successful 95% of time
			// but if its 0.9 it will be successful 90% of time	
			mom_fitness = pop.cachedFitness(mom.index)/pop.individuals.tagF[mom.index];
			if(runif(1)>mom_fitness){
				next;
			}
			
			if(repro_males.size() > 0){
				// allow all repr females to mate
				//probability of actually reproducing determined by p_repr
				if(runif(1)<p_repr){
					dad = sample(repro_males,1);
					child = pop.addCrossed(mom, dad);
				}
			}
		}
	}
	self.active = 0;
}

1 early() {		
	cat("gen,popSize,meanFitness,meanHet,B_year,B_gen,FROH_1Mb,avgLethal,avgVstrDel,avgStrDel,avgModDel,avgWkDel" + "\n");
	sim.addSubpop("p1", K_burn_in);
	p1.individuals.age = rdunif(K_burn_in, min=0, max=26);
}

// genetic rescue
(burn_in_time+decline_time+bottleneck_time+1):(burn_in_time+decline_time+bottleneck_time+500) late() {
	if ((sim.cycle-burn_in_time-decline_time-bottleneck_time) <= rescue_time & year_rescue > 0){
		if((sim.cycle-burn_in_time-decline_time-bottleneck_time - 1) % year_rescue == 0 ){
			if(rescue_sex == "F"){
				rescue_individuals = p3.individuals[p3.individuals.sex=='F' & p3.individuals.age >= 2 & p3.individuals.age <= 5];
			}
			else if (rescue_sex == "M"){
				rescue_individuals = p3.individuals[p3.individuals.sex=='M' & p3.individuals.age >= 2 & p3.individuals.age <= 5];
			}
			else if (rescue_sex == "Balance"){
				females = p4.individuals[p4.individuals.sex=='F'];
				males = p4.individuals[p4.individuals.sex=='M'];
				if (females.size() > males.size()){
					rescue_individuals = p3.individuals[p3.individuals.sex=='M' & p3.individuals.age >= 2 & p3.individuals.age <= 5];
				} else {
					rescue_individuals = p3.individuals[p3.individuals.sex=='F' & p3.individuals.age >= 2 & p3.individuals.age <= 5];
				}
			}
			else {
				rescue_individuals = p3.individuals[p3.individuals.age >= 2 & p3.individuals.age <= 5];
			}
			migrants_small_pop = sample(rescue_individuals, 1);
			if (num_deleterious >0){
				for (i in 0:(num_deleterious-1))
				{
					pos = rdunif(1,0,seqLength-1);
					if (type_deleterious == 'm1') {
						sample(migrants_small_pop.genomes, 1).addNewDrawnMutation(m1, pos);
					}
					else if (type_deleterious == 'm2') {
						sample(migrants_small_pop.genomes, 1).addNewDrawnMutation(m2, pos);
					}
					else if (type_deleterious == 'm3') {
						sample(migrants_small_pop.genomes, 1).addNewDrawnMutation(m3, pos);
					}
					else if (type_deleterious == 'm5') {
						sample(migrants_small_pop.genomes, 1).addNewDrawnMutation(m5, pos);
					}
				}
			}
			p4.takeMigrants(migrants_small_pop);
		}
	}
}

//track statistics pre-bottleneck every 500 years
1:(burn_in_time+decline_time) late() {
	if (sim.cycle % 500 == 0) {
		stats = getStats(p1, sampleSize);
		cat(sim.cycle+","+p1.individuals.size() + "," + stats + "\n");
	}
}

//track statistics bottleneck every 10 years
(burn_in_time+decline_time+1):(burn_in_time+decline_time+bottleneck_time) late() {
	if (sim.cycle % 10 == 0) {
		stats = getStats(p2, sampleSize);
		cat(sim.cycle+","+p2.individuals.size() + "," + stats + "\n");
	}
}

// bottleneck
burn_in_time+decline_time+1 early() {
	sim.addSubpop("p2",0);
	migrants = sample(p1.individuals, K_bottleneck);
	p2.takeMigrants(migrants);
	
	sim.addSubpop("p3",0);
	migrants = sample(p1.individuals, K_bottleneck);
	p3.takeMigrants(migrants);
	
	p1.fitnessScaling=0.0;
}

// small population
burn_in_time+decline_time+bottleneck_time+1 early() {
	sim.addSubpop("p4",0);
	migrants = sample(p2.individuals, K_small_pop);
	p4.takeMigrants(migrants);
	p2.fitnessScaling=0.0;
}

//Kyriazis, 2024
early() {
	inds = sim.subpopulations.individuals;
	ages = inds.age;
	mortality = L[ages];
	survival = 1 - mortality;
	inds.fitnessScaling = survival;
	
	//no fitness increases due to density dependence allowed
	//use p1.individuals.tagF to keep track of population-level fitness scaling
	//as well as fitness scaling for each individual due to age
	//need quantity to divide out of cachedFitness to get unscaled absolute fitness for the population
	// burn in
	if(sim.cycle <= burn_in_time){
		p1.fitnessScaling = min(K_burn_in /(p1.individualCount * mean(survival)), 1.0);
		p1.individuals.tagF = p1.individuals.fitnessScaling*p1.fitnessScaling;
	}
	// decline
	if(sim.cycle <= burn_in_time+decline_time & sim.cycle > burn_in_time){
		p1.fitnessScaling = min(K_declien /(p1.individualCount * mean(survival)), 1.0);
		p1.individuals.tagF = p1.individuals.fitnessScaling*p1.fitnessScaling;
	}
	// bottleneck, 1914~2014
	if(sim.cycle > burn_in_time+decline_time & sim.cycle <= burn_in_time+decline_time+bottleneck_time){
		p3.fitnessScaling = min(K_bottleneck /(p3.individualCount * mean(survival)), 1.0);
		p3.individuals.tagF = p3.individuals.fitnessScaling*p3.fitnessScaling;
		
		p2.fitnessScaling = min(K_bottleneck /(p2.individualCount * mean(survival)), 1.0);
		p2.individuals.tagF = p2.individuals.fitnessScaling*p2.fitnessScaling;
	}
	//2015~future
	if(sim.cycle > burn_in_time+decline_time+bottleneck_time){
		p3.fitnessScaling = min(K_big_pop_recover /(p3.individualCount * mean(survival)), 1.0);
		p3.individuals.tagF = p3.individuals.fitnessScaling*p3.fitnessScaling;
		
		p4.fitnessScaling = min(K_small_pop_recover /(p4.individualCount * mean(survival)), 1.0);
		p4.individuals.tagF = p4.individuals.fitnessScaling*p4.fitnessScaling;
	}
}

//Robinson 2022
(burn_in_time+decline_time+bottleneck_time+1):(burn_in_time+decline_time+bottleneck_time+500) late() {

	if(p4.individuals.size() < 2){
		stats = c("NA,NA,NA,NA,NA,NA,NA,NA,NA,NA"); //cant get stats from just one individual
	}
	
	// case when p1 size is less than sample size but greater than 1
	if(p4.individuals.size() < sampleSize & p4.individuals.size() > 1){	
		stats = getStats(p4, p4.individuals.size());
	}
	//case when p1 size is greater than or equal to sample size
	if(p4.individuals.size() >= sampleSize){ 
		stats = getStats(p4, sampleSize);
	}
	
	cat(sim.cycle +","+p4.individuals.size()+ "," + stats + "\n");
	
	//end sim if pop goes extinct
	if(p4.individuals.size() < 2){
		sim.simulationFinished();
	}
}


// define function to sample a population for
// mean fitness, heterozygosity, inbreeding load (2B), mean Froh
// and avg num of mutations of different classes per individual (very str del, str del, mod del, wk del)
function (s) getStats(o pop, i sampSize)
{
	i = sample(pop.individuals, sampSize, F);
	
	m = sortBy(i.genomes.mutations, "position"); //get all mutations in sample
	m_uniq = unique(m); // get rid of redundant muts
	DAF = sapply(m_uniq, "sum(m == applyValue);"); // count number of each mut in pop
	m_uniq_polym = m_uniq[DAF != i.genomes.size()]; //remove fixed mutations
	
	//calculate mean pop heterozygosity
	meanHet = calcHeterozygosity(pop.genomes);
	
	// tally mutations of each type
	count_wkDel = mean(i.genomes.countOfMutationsOfType(m1));
	count_modDel = mean(i.genomes.countOfMutationsOfType(m2));
	count_strDel = mean(i.genomes.countOfMutationsOfType(m3));
	count_semiLet = mean(i.genomes.countOfMutationsOfType(m4));
	count_let = mean(i.genomes.countOfMutationsOfType(m5));
	
	
	// initialize vector for storing total ROH length for each individual
	ROH_length_sumPerInd_1Mb = c();
	
	for (individual in i) {
		
		indm = sortBy(individual.genomes.mutations, "position");
		indm = indm[match(indm, m_uniq_polym) >= 0];   // Check that individual mutations are not fixed 
		indm_uniq = unique(indm);
		
		genotype = sapply(indm_uniq, "sum(indm == applyValue);");
		
		//code for getting ROHs
		ID_het = (genotype == 1); //outputs T/F for genotypes if they are het or homDer
		ID_homDer = (genotype == 2);
		pos_het = indm_uniq.position[ID_het]; //outputs positions of heterozgoys genotypes
		
		startpos = c(0, pos_het); //adds 0 to beggining of vector of hets
		endpos = c(pos_het, sim.chromosome.lastPosition); //adds last position in genome to vector of hets
		pos_het_diff = endpos - startpos;
		
		//sum for ROHs > 1Mb
		ROH_startpos_1Mb = startpos[pos_het_diff > ROHcutoff]; //filter out startpos that dont correspond to ROH > 1Mb
		ROH_endpos_1Mb = endpos[pos_het_diff > ROHcutoff];
		ROH_length_1Mb = pos_het_diff[pos_het_diff > ROHcutoff]; //vector of ROHs for each individual	
		ROH_length_sum_1Mb = sum(ROH_length_1Mb);
		ROH_length_sumPerInd_1Mb = c(ROH_length_sumPerInd_1Mb, ROH_length_sum_1Mb); // add sum of ROHs for each individual to vector of ROHs for all individuals
	}
	
	//calculate 2B (inbreeding load)
	
	// minor allele frequencies
	q = i.genomes.mutationFrequenciesInGenomes(m_uniq);
	
	// get selection coefficients as positive s
	s = -m_uniq.selectionCoeff;
	
	// replace mutations with s>1.0 with 1.0 (can happen when drawing from gamma distribution)
	s[s>1.0]=1.0;
	
	// get h for each mutation
	// note that this will not work if changing h using fitness callbacks
	h=m_uniq.mutationType.dominanceCoeff;
	
	// calculate s in terms of probability of survival to age x
	age=12;
	s_gen = 1-(1-s)^age;
	
	// calculate number of diploid lethal equivalents (2B or inbreeding load)
	// equation from Morton et al 1956
	B_year = 2*(sum(q*s)-sum(q^2*s)-2*sum(q*(1-q)*s*h));
	B_gen = 2*(sum(q*s_gen)-sum(q^2*s_gen)-2*sum(q*(1-q)*s_gen*h));
	
	return(mean(pop.cachedFitness(NULL)/pop.individuals.tagF) + "," + meanHet + "," + B_year + "," + B_gen + "," + mean(ROH_length_sumPerInd_1Mb)/seqLength + "," +  count_let + "," +  count_semiLet + "," + count_strDel+ "," + count_modDel + "," + count_wkDel);
}
EOM

