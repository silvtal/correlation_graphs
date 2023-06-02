## Description

Co-abundance correlations in the experimental passage data were evaluated following the same strategy recently employed by [Goyal et al (2022)](https://doi.org/10.7554/eLife.74987) which assesses if two populations abundances are more coupled than expected by chance. The method is grounded on the assumptions that OTU abundance fluctuations are independent of each other (i.e. no interactions) and follow a gamma distribution. First, the temporal abundance trajectories of all OTUs are calculated. Then, the Pearson correlation coefficient is calculated for each pair of OTUs based on their trajectories. The statistical significance of the correlations are calculated against the expected correlation distribution produced following a null model. 

- **For the null model, a gamma distribution is constructed for each OTU based on the experimental data.**

- **Then, random communities are generated from the gamma distributions and community compositions are renormalized by dividing each individual abundance by the communities' total sum.** 

- **The simulated communities are then arranged in passage trajectories and the correlations between each pair of OTUs obtained as per the experimental data.**

- **Finally, the P values of the experimentally observed correlations are obtained by comparing with the expected null distribution of 1000 simulated correlations.**

Interactions are deemed as present if its P value is below 0.05, and then depicted as a network using R package iGraph. In our case, we assessed the existence of interactions within the last 5 passages, when communities had been shown to remain phylogenetically stable (i.e. the PCGs had emerge and their relative abundance was comparatively stable). The analysis was conducted independently for experimental communities arising from different initial communities, and OTUs not detected in all samples of a passage trajectory were removed to avoid spurious correlations based on double zero co-abundances.
