# Creating a diversity input from real Borellia data

import os
import sys
import random
from Bio import Seq,SeqIO
from Bio.Alphabet import generic_dna

# ---------------------------------------------------------------------------------------
# allele_id           = string
# st_id               = string
# st_weight           = float
# locus_name          = string
# allele_name         = locus_name+_+allele_id
# position            = locus_name+_+integer
# nucleotide          = {A,C,G,T,X} 
# proportion          = floating number in [0,1]
# distance            = integer
# loci_names          = dictionary (i -> locus_name) for i from 1 to the number of loci in the MLST scheme
# st_profile          = dictionary (locus_name -> allele_id)
# sts_db              = dictionary (st_id -> st_profile)
# st_seq              = dictionary (msa_position -> nucleotide)
# sts_seq             = dictionary (st_id -> st_seq)
# sts_id              = list(st_id)
# sts_weights         = dictionary (st_id -> st_weight)
# sts_distances       = dictionary (st_id -> dictionary (st_id -> distance))
# locus_alleles_fasta = dictionary (allele_name -> SeqIO object for the fasta entry for the allele)
# loci_alleles_fasta  = dictionary (locus_name -> locus_alleles_fasta)
# positions           = list(msa_position)
# sts_proportions     = dictionary (st_id -> proportion) suming to 1
# nuc_proportions     = dictionary (nucleotide -> proportion) suming to 1
# pos_nuc_proportions = dictionary (position -> nuc_proportions)

# ---------------------------------------------------------------------------------------
# Read the profile of every strain type in a database of strain types
# Input:  sts_db_file_name: path to file containing the strain types profiles
#                           format: line starting with ST = loci names
#                           each other line is a tab separated line with alleles ids
#         loci_names
# Output: sts_db
def read_sts_db(sts_db_file_name,loci_names):
    # Loci names: from the strain types files, in the order they are in the input file
    # loci_names = {1:'clpA',2:'clpX',3:'nifS',4:'pepX',5:'pyrG',6:'recG',7:'rplB',8:'uvrA'}
    loci_nb    = len(loci_names.keys())
    # reading the strain types database
    sts_db_data = open(sts_db_file_name,'r').readlines()
    sts_db      = {}
    for data_line in sts_db_data:
        data_line_1 = data_line.rstrip().split('\t')
        if data_line_1[0] != 'ST':
            st_id      = data_line_1[0]
            st_profile = {loci_names[i]:data_line_1[i] for i in range(1,loci_nb+1)}
            sts_db[st_id] = st_profile
    return(sts_db)

# Filter from sts_db all entries of key in sts_id
def filter_sts_db(sts_db,sts_id):
    sts_db_filtered = {st_id:sts_db[st_id] for st_id in sts_bd.keys() if st_id not in sts_id}
    return(sts_db_filtered)

# ---------------------------------------------------------------------------------------
# Read the DNA+gap sequences of all alleles from the MSAs of an MLST scheme
# Input:  loci_names
#         alleles_fasta_dir: directory where are located the MSAs files for the MLST scheme alleles
#                            one file per locus of name locus_name+_msa_suffix+.+fasta_suffix
# Output: loci_alleles_fasta
def read_alleles_MSA(loci_names,alleles_fasta_dir,msa_suffix,fasta_suffix):
    loci_alleles_fasta = {}
    for locus_idx in loci_names.keys():
        locus_name = loci_names[locus_idx]
        if msa_suffix != '':
            alleles_sequences_file = alleles_fasta_dir+'/'+locus_name+'_'+msa_suffix+'.'+fasta_suffix
        else:
            alleles_sequences_file = alleles_fasta_dir+'/'+locus_name+'.'+fasta_suffix
        locus_alleles_fasta    = SeqIO.to_dict(SeqIO.parse(alleles_sequences_file,'fasta',generic_dna))
        loci_alleles_fasta[locus_name] = locus_alleles_fasta
    return(loci_alleles_fasta)

# ---------------------------------------------------------------------------------------
# Read the st_id to st_id distances from a csv file
def read_distances_csv(csv_file_name):
    sts_id,sts_distances = None,{}
    with open(csv_file_name) as f:
        for line in f:
            if sts_id is None:
                sts_id = [st_id.rstrip() for st_id in line.split(',')]
            else:
                distances = line.split(',')
                st_id = distances[0]
                sts_distances[st_id] = {sts_id[i]:int(distances[i]) for i in range(1,len(distances))}
    return(sts_distances)

# Extract k random strains at bounded distance from a source strain
def random_st_distance(sts_distances,st_id,d_max,k):
    st_id_candidates = [st_id_cand for st_id_cand,dist in sts_distances[st_id].items() if dist<=d_max and st_id_cand!=st_id]
    if k == 0:
        return(st_id_candidates)
    elif k>0 and k<=len(st_id_candidates):
        st_id_chosen = random.sample(st_id_candidates,k)
        return(st_id_chosen)
    else:
        return(None)

# ---------------------------------------------------------------------------------------
# Selects k random strain types from a satabase of strain types as read in the previous function
# Input:  sts_db
#         k
# Output: sts_db
def random_sts(sts_db,k):
    nb_sts           = len(sts_db.keys())
    sts_db_keys      = list(sts_db.keys())
    chosen_sts_index = [random.randint(0,nb_sts-1) for i in range(k)]
    chosen_sts       = {sts_db_keys[idx]: sts_db[sts_db_keys[idx]] for idx in chosen_sts_index}
    return(chosen_sts)

# Generates random floating numbers whose sum is 1.0 and each is greater than M/100
# Input:  sts_id
#         M
# Output: sts_proportions
def random_sts_proportions(M,sts_id):
    k = len(sts_id)
    l = list(range(1, 100-k*M))
    random.shuffle(l)
    l1 = sorted([l[i] for i in range(k-1)])
    l1.append(100-k*M)
    l2 = [l1[i] for i in range(k)]
    sts_proportions = {}
    for i in range(k):
        if i>0:
            l2[i]=l1[i]-l1[i-1]
        proportion             = (M+l2[i]) / 100.0
        sts_proportions[sts_id[i]] = proportion
    return(sts_proportions)

# Writes/read sts_proportions to/from a file
def write_sts_proportions(sts_proportions_file_name,sts_proportions):
    output_file = open(sts_proportions_file_name,'w')
    for st_id in sts_proportions.keys():
        proportion = sts_proportions[st_id]
        output_file.write(st_id+'\t'+str(proportion)+'\n')
    output_file.close()

def read_sts_proportions(sts_proportions_file_name):
    input_file  = open(sts_proportions_file_name,'r').readlines()
    sts_proportions = {}
    for l in input_file:
        l1 = l.rstrip().split('\t')
        st_id      = l1[0]
        proportion = float(l1[1])
        sts_proportions[st_id] = proportion
    return(sts_proportions)

# ---------------------------------------------------------------------------------------

DNA = ['A','C','G','T','X']
def read_DNA(n):
    DNA_dic = {'a':'A','A':'A','c':'C','C':'C','g':'G','G':'G','t':'T','T':'T','-':'X'}
    return(DNA_dic[n])

# Extract the DNA sequence of a strain type
# Input:  st_profile
#         loci_alleles_fasta
# Output: st_seq
def extract_st_seq(st_profile,loci_alleles_fasta):
    st_seq={}
    for locus_name in loci_alleles_fasta.keys():
        allele_id   = st_profile[locus_name]
        allele_name = locus_name+'_'+allele_id
        allele_seq  = str(loci_alleles_fasta[locus_name][allele_name].seq).rstrip()
        for i in range(len(allele_seq)):
            st_seq[locus_name+'_'+str(i)] = read_DNA(allele_seq[i])
    return(st_seq)

# Writing out the nucleotides and weights for given sts_id
# If weights are provided they are written, otherwise, they are given value 0.0
# Input:  output_file_name: format = st_id<TAB>st_weight<TAB>' '.join(list(<position:nucleotide>))
#         sts_seq
#         sts_weights: optional
def write_sts_seq_sts_weights(output_file_name,sts_seq,sts_weights=None):
    output_file = open(output_file_name,'w')
    for st_id in sts_seq.keys():
        if sts_weights==None:
            output_file.write(st_id+'\t0.0')
        else:
            output_file.write(st_id+'\t'+str(sts_weights[st_id]))
        for pos in sts_seq[st_id].keys():
            output_file.write('\t'+pos+':'+sts_seq[st_id][pos])
        output_file.write('\n')
    output_file.close()

# Reading the sequence and weight information from the same file
def read_sts_seq_sts_weights(input_file_name):
    input_file = open(input_file_name,'r').readlines()
    sts_seq     = {}
    sts_weights = {}
    for input_data in input_file:
        st_data        = input_data.rstrip().split('\t')
        st_id          = st_data[0]
        st_weight      = float(st_data[1])
        sts_seq[st_id]     = {seq_pos.split(':')[0]:seq_pos.split(':')[1] for seq_pos in st_data[2:]}
        sts_weights[st_id] = st_weight
    return((sts_seq,sts_weights))

# Separate from an sts_seq object invariant and variable positions
# Input:  sts_seq
# Output: two sts_seq object (variable and invariant positions)
def sts_seq_separate_invariant_variable_positions(sts_seq):
    sts_id            = list(sts_seq.keys())
    positions_all       = sts_seq[sts_id[0]].keys()
    positions_variable  = []
    positions_invariant = []
    for pos in positions_all:
        nucleotides = list(set([sts_seq[st_id][pos] for st_id in sts_id]))
        if len(nucleotides) > 1:            
            positions_variable.append(pos)
        else:
            positions_invariant.append(pos)
    sts_seq_variable  = {st_id:{pos:sts_seq[st_id][pos] for pos in positions_variable} for st_id in sts_id}
    sts_seq_invariant = {st_id:{pos:sts_seq[st_id][pos] for pos in positions_invariant} for st_id in sts_id}
    return((sts_seq_invariant,sts_seq_variable))

# Remove from an sts_seq sub-dictionary defined by sts_id all positions that are not in in the given positions
# and collapse the resulting identical st_id into a single one
def sts_seq_reduce_by_positions(sts_seq,positions,sts_id):
    sts_seq_reduced_wdup = {st_id:{position:sts_seq[st_id][position] for position in positions} for st_id in sts_id}
    sts_seq_reduced_aux1 = {'##'.join([k1+':'+v1 for k1,v1 in v.items()]):[] for v in sts_seq_reduced_wdup.values()}
    for st_id in sts_seq_reduced_wdup.keys():
        sts_seq_reduced_aux1['##'.join([k1+':'+v1 for k1,v1 in sts_seq_reduced_wdup[st_id].items()])].append(st_id)
    sts_seq_reduced = {'_'.join(v):{pos.split(':')[0]:pos.split(':')[1] for pos in k.split('##')} for k,v in sts_seq_reduced_aux1.items()}
    return(sts_seq_reduced)

# ---------------------------------------------------------------------------------------
# Extract the proportion of nucleotides for each position in positions and for a given list of sts_proportions
# Output:  pos_nuc_proportions
def extract_pos_nuc_proportions(sts_seq,sts_proportions):
    pos_nuc_proportions = {}
    sts_id          = list(sts_proportions.keys())
    positions       = list(sts_seq[sts_id[0]].keys())
    for position in positions:
        pos_nuc_proportions[position] = {nucleotide:0.0 for nucleotide in DNA}
        for st_id in sts_id:
            nucleotide = sts_seq[st_id][position]
            pos_nuc_proportions[position][nucleotide]+=sts_proportions[st_id]
    return(pos_nuc_proportions)

# Writing the proportion of nucleotides per position
def write_pos_nuc_proportions(output_file_name,pos_nuc_proportions):
    output_file = open(output_file_name,'w')
    for position in pos_nuc_proportions.keys():
        output_file.write(position)
        for nucleotide in pos_nuc_proportions[position].keys():
            output_file.write('\t'+nucleotide+':'+str(pos_nuc_proportions[position][nucleotide]))
        output_file.write('\n')
    output_file.close()

# Reading the proportions of nucleotides per position
def read_pos_nuc_proportions(input_file_name):
    input_file = open(input_file_name,'r').readlines()
    pos_nuc_proportions = {}
    for l in input_file:
        l1 = l.rstrip().split('\t')
        position = l1[0]
        nuc_proportions = {nuc_prop.split(':')[0]:float(nuc_prop.split(':')[1]) for nuc_prop in l1[1:]}
        pos_nuc_proportions[position] = nuc_proportions
    return(pos_nuc_proportions)

# ---------------------------------------------------------------------------------------
if __name__ == "__main__":

    TASK = sys.argv[1]

    borellia_loci_names = {1:'clpA',2:'clpX',3:'nifS',4:'pepX',5:'pyrG',6:'recG',7:'rplB',8:'uvrA'}

    
    # if TASK == 'read_dist':
    #     sts_distances = read_distances_csv(sys.argv[2])
    #     print(random_st_distance(sts_distances,'1',12,5))
    #     print(random_st_distance(sts_distances,'1',12,5))
    #     print(random_st_distance(sts_distances,'1',12,5))
    #     print(random_st_distance(sts_distances,'1',12,5))

    if TASK == 'evo_mod1':
        DATA_DIR        = sys.argv[2]
        STRAINS_DB_FILE = sys.argv[3]
        MSA_SUFFIX      = sys.argv[4]
        FASTA_SUFFIX    = sys.argv[5]
        K1              = int(sys.argv[6])
        K2              = int(sys.argv[7])
        D_MAX           = int(sys.argv[8])
        NB_NOVEL        = int(sys.argv[9])
        REPLICATES      = int(sys.argv[10])
        OUTPUT_PREF     = sys.argv[11]
    
        #SEED = int.from_bytes(os.urandom(8), byteorder="big")
        SEED=int(sys.argv[12])
        random.seed(SEED)        
        

        # Reading the strain types database and hamming distances
        borellia_sts_db = read_sts_db(DATA_DIR+'/'+STRAINS_DB_FILE,borellia_loci_names)
        borellia_sts_id = list(borellia_sts_db.keys())
        borellia_sts_distances = read_distances_csv(DATA_DIR+'/'+STRAINS_DB_FILE.replace('.txt','_'+MSA_SUFFIX+'_dist.csv'))
        
        log = open(OUTPUT_PREF+'.log','w')
        log.write('#Seed\t'+str(SEED)+'\n')
        log.write('#K1\t'+str(K1)+'\tK2'+str(K2)+'\n')
        log.write('#Replicates\t'+str(REPLICATES)+'\n')
        log.write('#Number of STs in full DB\t'+str(len(borellia_sts_id))+'\n')
        for rep in range(REPLICATES):
            log.write('##Replicate\t'+str(rep)+'\n')
            # Selecting K1 random source strains
            chosen_sts_source_db = random_sts(borellia_sts_db,K1)
            chosen_sts_source_id = list(chosen_sts_source_db.keys())
            # Selecting K2 random derived strains
            chosen_sts_derived_id_aux = []
            for j1 in chosen_sts_source_id:
                chosen_sts_derived_id_aux += [(j1,st_id) for st_id in random_st_distance(borellia_sts_distances,j1,D_MAX,0)]
            random.shuffle(chosen_sts_derived_id_aux)
            log.write('##Derived strains\t'+'\t'.join([st_id+':'+j1+':'+str(borellia_sts_distances[j1][st_id]) for (j1,st_id) in chosen_sts_derived_id_aux[0:K2]])+'\n')
            chosen_sts_derived_id = [st_id for (j1,st_id) in chosen_sts_derived_id_aux[0:K2]]
            # All chosen strains
            chosen_sts_id = chosen_sts_source_id + chosen_sts_derived_id 
            chosen_sts_db = {st_id:borellia_sts_db[st_id] for st_id in chosen_sts_id}
            # Creating random proportions for the chosen strain types
            chosen_sts_proportions = random_sts_proportions(10,chosen_sts_id)
            write_sts_proportions(OUTPUT_PREF+'_replicate_'+str(rep)+'_st_prop.txt',chosen_sts_proportions)
            log.write('##Chosen source strains\t'+'\t'.join([st_id+':'+str(chosen_sts_proportions[st_id]) for st_id in chosen_sts_source_id])+'\n')
            log.write('##Chosen derived strains\t'+'\t'.join([st_id+':'+str(chosen_sts_proportions[st_id]) for st_id in chosen_sts_derived_id])+'\n')
            # Extracting the sequence of each chosen strain type
            borellia_alleles_fasta = read_alleles_MSA(borellia_loci_names,DATA_DIR,MSA_SUFFIX,FASTA_SUFFIX)
            borellia_sts_seq       = {st_id: extract_st_seq(borellia_sts_db[st_id],borellia_alleles_fasta) for st_id in borellia_sts_id}
            chosen_sts_seq         = {st_id: extract_st_seq(chosen_sts_db[st_id],borellia_alleles_fasta) for st_id in chosen_sts_id}
            # Separating variable and invariant positions for the chosen strain types
            (chosen_sts_seq_inv,chosen_sts_seq_var) = sts_seq_separate_invariant_variable_positions(chosen_sts_seq)
            # Writing out chosen strain types sequences
            write_sts_seq_sts_weights(OUTPUT_PREF+'_replicate_'+str(rep)+'_var.txt',chosen_sts_seq_var)
            write_sts_seq_sts_weights(OUTPUT_PREF+'_replicate_'+str(rep)+'_inv.txt',chosen_sts_seq_inv)
            chosen_sts_var_pos    = chosen_sts_seq_var[chosen_sts_id[0]].keys()
            chosen_sts_inv_pos    = chosen_sts_seq_inv[chosen_sts_id[0]].keys()
            chosen_sts_nb_var_pos = len(chosen_sts_var_pos)
            chosen_sts_nb_inv_pos = len(chosen_sts_inv_pos)
            log.write('##Number of variable/invariant positions\t'+str(chosen_sts_nb_var_pos)+'\t'+str(chosen_sts_nb_inv_pos)+'\n')
            # Creating exact alleles proportions for the chosen strain types and variable columns
            chosen_pos_nuc_proportions = extract_pos_nuc_proportions(chosen_sts_seq_var,chosen_sts_proportions)
            write_pos_nuc_proportions(OUTPUT_PREF+'_replicate_'+str(rep)+'_alleles.txt',chosen_pos_nuc_proportions)
            # Filtering the strains db to keep only strains in agreement with the chosen sts variables and invariant nucleotides
            # --> Not done
            # Filtering out of the db to remove the chosen strains that are novel
            chosen_sts_novel_id = [chosen_sts_id[i] for i in random.sample(list(range(K1+K2)),K1+K2)][0:NB_NOVEL]
            borellia_sts_id_no_novel = list(set(borellia_sts_id).difference(set(chosen_sts_novel_id)))
            # Creating the strains db based on the variable positions
            borellia_var_sts_seq = sts_seq_reduce_by_positions(borellia_sts_seq,chosen_sts_var_pos,borellia_sts_id_no_novel)
            log.write('##Number of STs in reduced DB\t'+str(len(list(borellia_var_sts_seq.keys())))+'\n')
            write_sts_seq_sts_weights(OUTPUT_PREF+'_replicate_'+str(rep)+'_reduced_strains_db.txt',borellia_var_sts_seq)

    # Creating replicates of selecting K random strains and for each such simulated sample, writing the
    # file containing for each selected strain the nucleotides at variable positions as given by the alleles MSAs
    elif TASK == 'random_strains':    
        DATA_DIR        = sys.argv[2]
        STRAINS_DB_FILE = sys.argv[3]
        MSA_SUFFIX      = sys.argv[4]
        FASTA_SUFFIX    = sys.argv[5]
        K               = int(sys.argv[6])
        NB_STRAINS_2KEEP= int(sys.argv[7])
        REPLICATES      = int(sys.argv[8])
        OUTPUT_PREF     = sys.argv[9]

        #SEED = int.from_bytes(os.urandom(8), byteorder="big")
        SEED=20
        random.seed(SEED)        

        # Reading the strain types database        
        borellia_sts_db = read_sts_db(DATA_DIR+'/'+STRAINS_DB_FILE,borellia_loci_names)
        borellia_sts_id = list(borellia_sts_db.keys())

        log = open(OUTPUT_PREF+'.log','w')
        log.write('#Seed\t'+str(SEED)+'\n')
        log.write('#K\t'+str(K)+'\n')
        log.write('#Replicates\t'+str(REPLICATES)+'\n')
        log.write('#Number of STs in full DB\t'+str(len(borellia_sts_id))+'\n')
        for rep in range(REPLICATES):
            print('Replicate\t'+str(rep))
            log.write('##Replicate\t'+str(rep)+'\n')
            NB_STRAINS_2KEEP_REP = NB_STRAINS_2KEEP
            # Selecting K random strains
            chosen_sts_db = random_sts(borellia_sts_db,K)
            chosen_sts_id = list(chosen_sts_db.keys())
            # Creating random proportions for the chosen strain types
            chosen_sts_proportions = random_sts_proportions(10,chosen_sts_id)
            write_sts_proportions(OUTPUT_PREF+'_replicate_'+str(rep)+'_st_prop.txt',chosen_sts_proportions)            
            log.write('##Chosen strains\t'+'\t'.join([st_id+':'+str(prop) for st_id,prop in chosen_sts_proportions.items()])+'\n')
            # Extracting the sequence of each chosen strain type
            borellia_alleles_fasta = read_alleles_MSA(borellia_loci_names,DATA_DIR,MSA_SUFFIX,FASTA_SUFFIX)
            borellia_sts_seq       = {st_id: extract_st_seq(borellia_sts_db[st_id],borellia_alleles_fasta) for st_id in borellia_sts_id}
            chosen_sts_seq         = {st_id: extract_st_seq(chosen_sts_db[st_id],borellia_alleles_fasta) for st_id in chosen_sts_id}
            # Separating variable and invariant positions for the chosen strain types
            (chosen_sts_seq_inv,chosen_sts_seq_var) = sts_seq_separate_invariant_variable_positions(chosen_sts_seq)
            # Writing out chosen strain types sequences
            write_sts_seq_sts_weights(OUTPUT_PREF+'_replicate_'+str(rep)+'_var.txt',chosen_sts_seq_var)
            write_sts_seq_sts_weights(OUTPUT_PREF+'_replicate_'+str(rep)+'_inv.txt',chosen_sts_seq_inv)
            chosen_sts_var_pos    = chosen_sts_seq_var[chosen_sts_id[0]].keys()
            chosen_sts_inv_pos    = chosen_sts_seq_inv[chosen_sts_id[0]].keys()
            chosen_sts_nb_var_pos = len(chosen_sts_var_pos)
            chosen_sts_nb_inv_pos = len(chosen_sts_inv_pos)
            log.write('##Number of variable/invariant positions\t'+str(chosen_sts_nb_var_pos)+'\t'+str(chosen_sts_nb_inv_pos)+'\n')
            # Creating exact alleles proportions for the chosen strain types and variable columns
            chosen_pos_nuc_proportions = extract_pos_nuc_proportions(chosen_sts_seq_var,chosen_sts_proportions)
            write_pos_nuc_proportions(OUTPUT_PREF+'_replicate_'+str(rep)+'_alleles.txt',chosen_pos_nuc_proportions)
            # Filtering the strains db to keep only strains in agreement with the chosen sts variables and invariant nucleotides
            # --> Not done
            # Creating the strains db based on the variable positions
            borellia_var_sts_seq = sts_seq_reduce_by_positions(borellia_sts_seq,chosen_sts_var_pos,borellia_sts_id)
            log.write('##Number of STs in reduced DB\t'+str(len(list(borellia_var_sts_seq.keys())))+'\n')
            write_sts_seq_sts_weights(OUTPUT_PREF+'_replicate_'+str(rep)+'_reduced_strains_db.txt',borellia_var_sts_seq)
            # Reducing the strain db 1: weighting the second part of the db by 1
            borellia_var_sts_id  = list(borellia_var_sts_seq.keys())
            if len(borellia_var_sts_id)<NB_STRAINS_2KEEP_REP:
                NB_STRAINS_2KEEP_REP = len(borellia_var_sts_id)
            borellia_var_sts_weights = {st_id:1.0 for st_id in borellia_var_sts_id}
            for i in range(NB_STRAINS_2KEEP_REP):
                borellia_var_sts_weights[borellia_var_sts_id[i]] = 0.0
            write_sts_seq_sts_weights(OUTPUT_PREF+'_replicate_'+str(rep)+'_reduced_strains_db_w1.txt',borellia_var_sts_seq,borellia_var_sts_weights)
            # Reducing the strain db 2:removing the second part of the db
            for st_id in borellia_var_sts_id:
                if borellia_var_sts_weights[st_id] == 1.0:
                    del borellia_var_sts_seq[st_id]
            write_sts_seq_sts_weights(OUTPUT_PREF+'_replicate_'+str(rep)+'_reduced_strains_db_w2.txt',borellia_var_sts_seq,borellia_var_sts_weights)
            log.write('##Number of STs in reduced DB for ILP2\t'+str(len(borellia_var_sts_seq.keys()))+'\n')
