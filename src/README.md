# Source files

The code is in python3.

The file 'borrelia_utils.py' provides scripts to manipulate an MLST scheme and data files.

The file ILP_phasing.py implements the ILP described in the notebook '../doc/Experiments_02122019.ipynb'.
The ILP is implemented in CPLEX, using docplex.

To reproduce the results shown in the directory ../results', it can be used as follows:  
python ILP_phasing \
    STRAINS_DATA_FILE \ #File describing the known STs that can occur in a sample
    ALLELES_DATA_FILE \ #File describing the allele frequencies observed in a sample
    PHASING_DATA_FILE \ #File describing the phasing data 
    OUTPUT_PREF       \ #Prefix of the output files
    K                 \ #Maximum number of strains in the sample
    0                 \ #Unused parameter (related to Hamming distance)
    0                 \ #Unused parameter (related to Hamming distance)
    0.0               \ #Minimum abundance of strains
    1                 \ #Weight of the usage of known STs in the objective function
    -2                \ #Same for the average deviation term
    1                 \ #Same for the maximum deviation term
    0                 \ #Unused parameter (related to Hamming distance)
   
For example, one run would be  
python ./ILP_phasing.py ../results/SRR2034335_bt2_2_1_01_st_db.txt ../results/SRR2034335_bt2_2_1_01_varpos.txt ../results/SR
R2034335_bt2_2_1_01_phasing.txt ../results/SRR2034335_bt2_2_1_01_S_4_XMIN_00_phasing 4 0 0 0.0 1 -2 1 0 

The input files have the following formats:
STRAINS_DATA_FILE (SRR2034335_bt2_2_1_01_st_db.txt): each line is a ST defined in terms of the variable positions of the sample
ST<TAB>{var_pos:nucleotide}^*

ALLELES_DATA_FILE (SRR2034362_bt2_2_1_01_varpos.txt): each line respresents the content of a variable position:
var_pos<TAB>{nucleotide:frequency}^*

PHASING_DATA_FILE (SRR2034335_bt2_2_1_01_phasing.txt): each line represents one read and the variable positions it overlaps
read_id<TAB>var_pos:nucleotide

The output file is OUTPUT_PREF_sol.txt and describes the solution in a human readable format.

