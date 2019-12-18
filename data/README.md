# Data files

We rely on the Multi-Locus Sequence Typing scheme for *Borrelia* available 
at <a href="https://pubmlst.org/borrelia/">pubmlst</a>, downloaded on 
October 11, 2019.

This scheme contains 8 genes, clpA, clpX, nifS, pepX, pyrG, recG, rplB, uvrA. 
For each of these genes, the available alleles were aligned using MAFFT to 
form a multiple sequence alignment (MSA). 
The alleles sequences for gene X are available in the files 'X.fasta' and the 
MSA in the file 'X_mafft.fasta'.

We downloaded also from <a href="https://pubmlst.org/borrelia/"pubmlst</a> a list 
of knwon strain types (ST), each being composed of a set of 8 alleles, one per
locus, and originate from North America, Asia and Europe. 
They are available in the file 'bigsdb.txt'. In this file each line represents one
strain type, defined as 8 alleles, one per locus of the MLST scheme.

To map short reads, we defined one genome per ST as follows:
We first downloaded a fully assembled Borrelia burgdorferi strain, the strain B331 
(https://www.ncbi.nlm.nih.gov/nuccore/CP017201, file CP017201.fasta). Then for each of 
the 8 loci of the MLST scheme, we aligned, using BLASTn, all alleles onto this strain 
and recorded the location of the best alignment, correspondign to the allele present 
in the strain. Then we extracted from the genome 75bp before and after each such 
alignment (called the flanking regions from now). Finally for each ST, for each locus, 
we flanked the allele present in the ST by the corresponding flanking regions and concatenated 
all these sequences into a single sequence by a segment of 30 Ns. These genomes are
available in the file 'bigsdb_flanked.fasta'.

Lasst, the file 'samples.txt' contains the <a href="https://www.ncbi.nlm.nih.gov/sra">SRA</a> 
accession numbers of the 24 isolates we analyzed.
