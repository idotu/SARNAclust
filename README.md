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
