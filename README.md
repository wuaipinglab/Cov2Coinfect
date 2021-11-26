# Cov2Coinfect

Here is the code and data for *Hang-Yu Zhou, et al.* "Genomic evidence for divergent co-infections of severe acute respiratory syndrome coronavirus 2 lineages"

***

## Overview
![image](https://github.com/wuaipinglab/SARS-CoV-2_co-infection/blob/main/img/Figure2.png)

### Hypergeometric Distribution
![image](https://github.com/wuaipinglab/SARS-CoV-2_co-infection/blob/main/img/formula.png)

where
* N is the population size (the total number of nonsynonymous mutations that occur in 2.5 million SARS-CoV-2 consensus genomes), 
* K is the number of success states in the population (the number of feature variations of one lineage),
* n is the number of draws (the number of remaining undefined mutations of sample),
* k is the number of observed successes (the number of remaining undefined mutations that occur both in sample and lineage feature variations)

For each lineage, hypergeometric test computes the probability (p-value) of observed successes (mutations that occurred both in sample and lineage feature variations) under “null hypothesis”: the hypothesis that there is nothing special about the lineage. If the p-value is sufficiently low, we can reject the null hypothesis as too impossible and conclude that the sample are highly correlated with the tested lineage.

### Statistical Significance
* Lineages with lower p-value are more likely assigned to the sample.

### Mutation Frequency Uniformity
* Frequencies of mutations that occur both in sample and one lineage feature variations should have a standard deviation less than 20.

### Mutation Concentration
* Number of mutations that occur both in sample and one lineage feature variations should be greater than 6.
* Number of mutations that occur both in sample and one lineage feature variations divided by number of the lineage feature variations should be greater than 0.3.

## Usage
1. Detect variant in SRA file using `wuaipinglab/sra2variant`.

2. Run `get_lineagesFV_and_mutationNum.py` to get lineages feature variations and number of global SARS-CoV-2 nonsynonymous mutations.

3. Run `get_candidate_lineages.py` to get candidate with-in host lineages.

4. Run `identify_co-infection_lineages.py` to get defined with-in host lineages and identify co-infection cases.
