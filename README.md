# Genomic evidence for divergent co-infections of SARS-CoV-2 variants

Here is the code and data for *Hang-Yu Zhou, et al.* "Genomic evidence for divergent co-infections of SARS-CoV-2 variants." Journal (2021); doi:xxx

***

## Overview
![image](https://github.com/wuaipinglab/SARS-CoV-2_co-infection/blob/main/img/Fig2.png)

#### Hypergeometric Distribution
![image](https://github.com/wuaipinglab/SARS-CoV-2_co-infection/blob/main/img/formula.png)

where
* N is the population size (the total number of nonsynonyumous mutations that occur in 2.5 million SARS-CoV-2 consensus genomes), 
* K is the number of success states in the population (the number of feature variations of one lineage),
* n is the number of draws (the number of remaining undefined mutations of sample),
* k is the number of observed successes (the number of remaining undefined mutations that occur both in sample and lineage feature variations)

For every lineage, hypergeometric test computes the probability (p-value) of observing **k** or more remaining undefined mutations that occur both in sample and this lineage feature variations under null hypothesis that is nothing special about this lineage. If this probability is sufficiently low, we can decide to reject the null hypothesis as too unlikely - mutations in sample are highly correlated with feature variations of this lineage.

#### Statistical Significance
Lineages with lower p-value is more likely assigned to the sample.

#### Mutation Frequency Uniformity
Frequencies of mutations that occur both in sample and one lineage feature variations should have a standard deviation less than 20.

#### Mutation Concentration
Number of mutations that occur both in sample and one lineage feature variations should be greater than 6.
Number of mutations that occur both in sample and one lineage feature variations divided by number of the lineage feature variations should be greater than 0.3.

## Usage
1.Run `get_lineagesFV_and_mutationNum.py` to get lineages feature variations and number of global SARS-CoV-2 nonsynonymous mutations.

2.Run `get_candidate_lineages.py` to get candidate with-in host lineages.

3.Run `identify_co-infection_lineages.py` to get defined with-in host lineages, and distinguish co-infection from quasispecies in one virus.
