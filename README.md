# SARNAclust
This repository contains code for both RNApeakFold and SARNAclust.
RNApeakFold is an RNA secondary structure prediction tool tailored to CLIP peaks.
SARNAclust is a clustering algorithm for CLIP peaks that typically uses the
secondary structure prediction from RNAfold.

## SARNAclust
SARNAclust requires the following python libraries:
* networkx
* [EDeN](https://github.com/fabriziocosta/EDeN) (graph kernel)
* sklearn 
* matplotlib
* forgi (bulge graph)
* Biotools

It also requires Locarna 1.8.12 to be installed.
Moreover, it imports dpcluster, an in-house implementation [1] of Density Clustering also available here.
There are also a few fixed parameters that the user can change in the code:
* MINL = 10 - min length of a peak for it to be considered
* MAXL = 320 - max length of a peak for it to be considered
* MAXP = 1400 - maximun number of peaks to be clustered in the same iteration

SARNAclust inputs an extended fasta file (-file FileName, default GAGAandRandom.efa)) (i.e., the output from RNAfold) where there has to be
a set of sequence/structures each composed of three lines: fasta comment, sequence
and secondary structure. It also requires a set of parameters:
* -r - radius (see [EDeN](https://github.com/fabriziocosta/EDeN) library, default 2)
* -d - distance (see [EDeN](https://github.com/fabriziocosta/EDeN) library, default 3)
* -cluster_option - depends on the method selected, but the format is "(method,param1,param2, ... ,paramN)" (default DBSCAN '(1,0.6,10)'):
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
* -graph_option (default 9):
  * 0 - Complete graph with no structural information
  * 1 - graphProt minus directional                     
  * 2 - Bulge graph plus Complete Graph                 
  * 3 - Bulge graph                                     
  * 4 - Bulge graph plus sequence in HP loops           
  * 5 - Bulge graph plus sequence in IL loops           
  * 6 - Bulge graph plus sequence in external loops     
  * 7 - Bulge graph plus sequence in double stranded    
  * 8 - Bulge graph plus sequence in HP and IL          
  * 9 - Bulge graph plus sequence in single stranded    
  * 10 - Bulge graph plus sequence everywhere           
  * 11 - Same as 2 but no sequence in double stranded  
* -thClus - distance threshold for merging clusters (0:1, default 0.6)
* -iterations - number of iterations (default 5)
* -verbose (default 0)
* -alpha - alphabet (0- A,C,G,U; 1- G,R, default 0)

Since python can only handle a certain number of graphs loaded in memory (around 800 depending on the machine
running the code), SARNAclust performs a number of *iterations* selecting each time *MAXP* random sequences
and then merging clusters among iterations if the average distance between sequences in 2 clusters is lower
than *thClus*. Each sequence/structure is decomposed into substructures that are automatically detected and the
clustering is performed on those, unless graph option *0* is selected. To disallow structure decomposition just
change readFile(fn,gopt) in the main function for readFile(fn,0). This should allow the user to up *MAXP* to
around 2000.

At the end, SARNAclust returns a set of clusters along with their consensus sequence and structure.
You can try with our example file GAGAandRandom.efa and DBSCAN (default options), for instance, with the following command
line:

python SARNAclust.py

# References

* [1] - Rodriguez A, Laio A. 2014. Machine learning. Clustering by fast search and find of density peaks. Science 344: 1492–1496.

