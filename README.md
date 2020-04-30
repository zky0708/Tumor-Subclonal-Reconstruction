# Subclonal reconstruction from tumor genomic variation data using a cascade ensemble model
The repository contains the code for the winning model we developed in the [ICGC-TCGA-DREAM Somatic Mutation Calling -- Tumor Heterogeneity Challenge (SMC-Het)](https://www.synapse.org/#!Synapse:syn2813581/wiki/303137). The model is designed to accurately infer intra-tumor heterogeneity, including the number of subclones, cellular proportion and mutation assignment of each subclone, as well as their phylogenic relationships. The paper describing this Challenge has been published on [Nature Biotechnology](https://www.nature.com/articles/s41587-019-0364-z).

## Inputs and Outputs
Input data: 
1. A variant call format ([VCF](http://www.internationalgenome.org/wiki/Analysis/variant-call-format)) file that contains data on the somatic mutations found in the tumor sample. 
2. A [Battenberg](https://github.com/cancerit/cgpBattenberg) file that contains data on the copy number of each chromosome within the tumor sample.

Output data:
1. Tumor cellularity (*c*), i.e. the proportion of cells in the sample that are cancerous.
2. Number of cancerous subpopulations (*k*), which are defined by a unique set of somatic mutations.
3. Cellular proportions of each identified subpopulation (*p*).
4. Mutation assignments to subpopulations (*a*).
5. Phylogenetic tree structure among these subpopulations (*T*).


## Model description
We developed a cascade ensemble model, incorporating extensive use of the genomic attributes of tumors to detect the false positives of the mutation calls and adjust the predictions of the subclones. This model consists of four sequentially connected modules (M1 to M4). 

![image of flowchart](https://github.com/zky0708/tumor-subclonal-reconstruction/blob/master/flowchart.png)

1. The first module (M1) of the algorithm derives an estimate of the cellularity using the phased genotype of a tumor sample. 
2. The second (M2) and third (M3) modules predict the number of subclones and their proportions in the samples, based on a modified **truncated Dirichlet process**. Specifically, M2 decomposes each sample using the somatic mutations with a deterministic ratio of variant and reference alleles and the estimated cellularity from the first module. 
3. To reduce the effect of the false positives (FPs), M3 filters the mutations by the status of their genomic positions using the database of single nucleotide polymorphisms (dbSNP) and their base read quality and requires the positions in the genome of the matched normal samples to be variant-free. It then uses only the mutations that pass these filters, as well as the estimated number of subclones from the second module, to decompose the samples, and it corrects the predictions using the information of copy number variations. Among the resulting estimates of the three modules, the algorithm outputs the most conservative estimate of cellularity along with the corresponding characteristics of subclones and assignment of mutations. 
4. Finally, the last module (M4) reconstructs the evolutionary relationships of these subclones using a heuristic tree building method. 

## References
Here are some useful articles for the Dirichlet process techniques.

[[1]](https://doi.org/10.1198/016214501750332758) Ishwaran, H. and James, L. F. Gibbs sampling methods for stick-breaking priors. *J Am Stat Assoc*, 96, 453 (Mar 2001), 161-173.

[[2]](http://papers.nips.cc/paper/4108-tree-structured-stick-breaking-for-hierarchical-data) Adams, R. P., Ghahramani, Z. and Jordan, M. I. Tree-structured stick breaking for hierarchical data. *Proceedings of the 23rd International Conference on Neural Information Processing Systems - Volume 1* (Vancouver, British Columbia, Canada, 2010). 

[[3]](http://www.cs.columbia.edu/~blei/papers/BleiJordan2004.pdf) Blei, D. M. and Jordan, M. I. Variational Inference for Dirichlet Process Mixtures. Bayesian Anal, 1, 1 (2006), 121-143.
