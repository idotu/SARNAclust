# SARNAclust
This repository contains code for both RNApeakFold and SARNAclust.
RNApeakFold is an RNA secondary structure prediction tool tailored to CLIP peaks.
SARNAclust is a clustering algorithm for CLIP peaks that typically uses the
secondary structure prediction from RNApeakFold.

## RNApeakFold
RNApeakFold requires RNAfold from the Vienna Package for calculating RNA base pairing
probabilities. With these probabilities, RNApeakFold performs Nussinov-like secondary
structure prediction on a given set of RNA sequences.
As input, RNApeakFold takes a fasta file with a list of RNA sequences. Fasta comments need to 
have at least 4 fields (separated by | ) that completely identify each sequence (see test.fa).
RNApeakFold generates several .ps and .fa files that can be removed afterwards by the user.
The result is printed on the terminal.

## SARNAclust
SARNAclust requires the following python libraries:
* networkx
* [EDeN](https://github.com/fabriziocosta/EDeN) (graph kernel)
* sklearn 
* matplotlib
* forgi (bulge graph)
* Biotools

It also requires ClustalW to be installed, and the path indicated in variable *clustalw_exe*.
Moreover, it imports dpcluster, an in-house implementation [1] of Density Clustering also available here.
There are also a few fixed parameters that the user can change in the code:
* MINL = 10 - min length of a peak for it to be considered
* MAXL = 80 - max length of a peak for it to be considered
* MAXP = 2000 - maximun number of peaks to be clustered in the same iteration

SARNAclust inputs an extended fasta file (i.e., the output from RNApeakFold) where there has to be
a set of sequence/structures each composed of three lines: fasta comment (see RNApeakFold above), sequence
and secondary structure. It also requires a set of parameters:
* r - radius (see [EDeN](https://github.com/fabriziocosta/EDeN) library)
* d - distance (see [EDeN](https://github.com/fabriziocosta/EDeN) library)
* cluster_option - depends on the method selected, but the format is "(method param1 param2 ... paramN)":
  * 0 - MeanShift (sklearn)
  * 1 - DBSCAN (sklearn)  
    * similarity - maximum distance between 2 sequences in the same cluster (0:1)
    * minsamples - minimum number of samples for a cluster to be considered
  * 2 - Affinity Propagation (sklearn)                        
  * 3 - Prints similarity matrix                      
  * 4 - KMeans (sklearn)                                 
  * 5 - Spectral Clustering (sklearn)                            
  * 6 - Density Clustering (dclust in-house implementation of [1], read paper to understand parameters) 
    * Flag for adding halo as part of the cluster (0,1) 
    * lowerbound percentile
    * dc
* graph_option:
* thClus - distance threshold for merging clusters (0:1)
* iterations - number of iterations
* verbose 
* debug - plots the graph transformations

Since python can only handle a certain number of graphs loaded in memory (around 2000 depending on the machine
running the code), SARNAclust performs a number of *iterations* selecting each time *MAXP* random sequences
and then merging clusters among iterations if the average distance between sequences in 2 clusters is lower
than *thClus*.

At the end, SARNAclust returns a set of clusters along with their consensus sequence and structure.

# References

* [1] - Rodriguez A, Laio A. 2014. Machine learning. Clustering by fast search and find of density peaks. Science 344: 1492–1496.

