#!/bin/bash

burn_in_time=${1}
out_dir=${2}


cat > ${out_dir}/panda_slim_${burn_in_time}.slim << EOM
initialize() {
	
	initializeSLiMModelType("nonWF");
	
	defineConstant("burn_in_time", ${burn_in_time});				// burn-in
	defineConstant("first_decline_time", 9000);		          // first_decline

	defineConstant("K_burn_in", 47100);					            // burn-in, carrying capacity
	defineConstant("K_decline_first", 2355);		          	// first decline, carrying capacity

	defineConstant("fitness_Scale", 1.0);				            // fitness scale
	defineConstant("sampleSize", 30); 					            // sample to stats
	defineConstant("ROHcutoff", 100000);				            // ROH to stats >1Mb
	
	defineConstant("p_repr", 0.625); 					              // annual probability of an adult female reproducing, wei,1994 
	initializeSex("A"); 									                 	// Sexual model
	initializeMutationRate(1.29e-8); 					              // zhao, 2013
	
	defineConstant("g",18530); 							                // number of genes
	defineConstant("geneLength", 1800); 				            // difine gene length
	defineConstant("seqLength", g*geneLength);		          // total gene length to caculate het
	cat("Genome length:"+seqLength+"\n");						
	
	defineConstant("L", c(0.6, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0,0.0,0.0,0.0,0.0, 0.1, 0.2, 0.3, 0.3, 0.5, 1.0)); // age-dependent increases in mortality
	
	defineConstant("h_strDel", 0.05); 
	defineConstant("h_modDel", 0.2); 
	defineConstant("h_wkDel", 0.45);
 
	// Kyriazis, 2024
	// set up discrete DFE with four mutation types coming from gamma DFE
	// augmented with recessive lethals
	// this approach for implementing an h-s relationship is much faster than using fitness callbacks in SLiM (see manual)

	//strongly deleterious mutations (s<-0.01)
	initializeMutationType("m1", h_strDel, "s", "do x=rgamma(1,-0.01314833,0.186); while (x >= -0.01); return x;");
	//moderately deleterious mutations (-0.001 > s >= -0.01)
	initializeMutationType("m2", h_modDel, "s", "do x=rgamma(1,-0.01314833,0.186); while (x < -0.01 | x >= -0.001); return x;");
	//weakly deleterious mutations (s >= -0.001)
	initializeMutationType("m3", h_wkDel, "s", "do x=rgamma(1,-0.01314833,0.186); while (x < -0.001); return x;");
	//neutral mutations
	initializeMutationType("m4", 0.5,"f",0.0);
	
	
	//ratio of different deleterious mutation types taken from vaquita DFE (sum to 100 below)
	//assume ratio of deleterious to neutral muts of 2.31:1 (Huber et al 2017) 
	//giving 100/2.31=43.3 for neutral mutations below
	initializeGenomicElementType("g1", c(m1,m2,m3,m4), c(26.2, 24.7, 49.1, 43.3));

	//number of genes on each autosome from giant panda annotations
	gene_vec=c(1298,1563,911,1580,1048,939,924,1329,609,1070,719,1402,1347,712,691,1312,410,190,212,264);
	
	gene_num=sum(gene_vec);
	
	for (i in 1:gene_num){
		initializeGenomicElement(g1, ((i-1)*geneLength)+(i-1), (i*geneLength)+(i-2) );
	}
	
	rates=NULL;
	
	//assume no recombination within genes, a rate of 1e-3 between genes, and free recombination between chroms
	for (i in 1:(size(gene_vec)-1)){
		rates=c(rates, 0, rep(c(1e-3, 0), asInteger(gene_vec[i-1]-1)), 0.5);
	}
	rates=c(rates, 0, rep(c(1e-3, 0), asInteger(gene_vec[size(gene_vec)-1]-1)));
	
	ends=NULL;
	for (i in 1:gene_num){
		ends=c(ends, (i*geneLength)+(i-2), (i*geneLength)+(i-1));
	}
	ends=ends[0:(size(ends)-2)];
	
	initializeRecombinationRate(rates, ends);

}


reproduction() {
	for(pop in sim.subpopulations){
		
		males = pop.individuals[pop.individuals.sex=='M'];
		females = pop.individuals[pop.individuals.sex=='F'];
		
		//get reproductive age males and females
		repro_females = females[females.age >= 6 & females.age <= 20 ];
		repro_males = males[males.age >= 7 & males.age <= 20 ];
		
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
	cat("gen,popSize,meanFitness,meanHet,B_gen,B_year,FROH_1Mb,avgStrDel,avgModDel,avgWkDel,avgHomDel,avgHetDel\n");
	sim.addSubpop("p1", K_burn_in);
	p1.individuals.age = rdunif(K_burn_in, min=0, max=26);
	defineConstant("simID", getSeed());
	//use sim.tag to enforce sims to run for 500 years after reaching min pop size
	sim.tag = 500;
}

(burn_in_time+first_decline_time) late() {		
	sim.outputFull("/media/desk16/tly6838/data/slim/burn_in" + simID + ".txt");
	sim.simulationFinished();
}

// print generation time and ages
100 early() {
	// catn(p1.individuals.age);
	females = p1.individuals[p1.individuals.sex=='F'];
	males = p1.individuals[p1.individuals.sex=='M'];
	// get age > 1
	After_1_females = females[ females.age >= 1];
	After_1_males = males[males.age >= 1];
	// get reproductive age males and females
	repro_females = females[ females.age >= 6 & females.age <= 20 ];
	repro_males = males[males.age >= 7 & males.age <= 20 ];
	
	// print mean lift span and generation time
	// catn(mean(c(females.age,males.age)));
	// catn(mean(c(After_1_females.age,After_1_males.age)));
	// catn(mean(c(repro_females.age,repro_males.age)));
}

//track statistics before second decline every 500 years
1:(burn_in_time+first_decline_time) late() {
	if (p1.individuals.size()<sampleSize){
		cat(sim.cycle+ ",Fail\n");
		sim.simulationFinished();	
	}
	if (sim.cycle % 500 == 0) {
		stats = getStats(p1, sampleSize);
		cat(sim.cycle+","+p1.individuals.size() + "," + stats + "\n");
	}
}


early() {
	inds = sim.subpopulations.individuals;
	ages = inds.age;
	mortality = L[ages];
	survival = 1 - mortality;
	inds.fitnessScaling = survival;
	//sim.outputMutations(sim.mutationsOfType(m3));
	
	//no fitness increases due to density dependence allowed
	//use p1.individuals.tagF to keep track of population-level fitness scaling
	//as well as fitness scaling for each individual due to age
	//need quantity to divide out of cachedFitness to get unscaled absolute fitness for the population
	// burn in
	if(sim.cycle <= burn_in_time){
		p1.fitnessScaling = min(K_burn_in /(p1.individualCount * mean(survival)), fitness_Scale);
		p1.individuals.tagF = p1.individuals.fitnessScaling*p1.fitnessScaling;
	}
	// first decline
	if(sim.cycle <= burn_in_time+first_decline_time & sim.cycle > burn_in_time){
		p1.fitnessScaling = min(K_decline_first /(p1.individualCount * mean(survival)), fitness_Scale);
		p1.individuals.tagF = p1.individuals.fitnessScaling*p1.fitnessScaling;
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

	count_strDel = mean(i.genomes.countOfMutationsOfType(m1));
	count_modDel = mean(i.genomes.countOfMutationsOfType(m2));
	count_wkDel = mean(i.genomes.countOfMutationsOfType(m3));
	
	count_hom_Del = c();
	count_het_Del = c();
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
	
	// 
	del_muts = c(individual.genomes.mutationsOfType(m1),individual.genomes.mutationsOfType(m2),individual.genomes.mutationsOfType(m3));
	het_del_num = 0;
	hom_del_num = 0;
	if (del_muts.length()>0) {
		for(m in del_muts){
			//check if mut is heterozygous
			if(individual.genomes.mutationCountsInGenomes(m)==1){
				het_del_num=het_del_num+1;
			}
			else{
				hom_del_num=hom_del_num+1;		
			}
		}
	}
	count_hom_Del=c(count_hom_Del,hom_del_num);
	count_het_Del=c(count_het_Del,het_del_num);	

	return(mean(pop.cachedFitness(NULL)/pop.individuals.tagF) + "," + meanHet + "," + B_gen + "," + B_year + "," + mean(ROH_length_sumPerInd_1Mb)/seqLength + "," + count_strDel+ "," + count_modDel + "," + count_wkDel+","+mean(count_hom_Del)+","+mean(count_het_Del));
}



EOM


