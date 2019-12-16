# Implementation of a within-host evolution model for the diversity problem
# Cedric Chauve, November 15, 2019

import os
import sys
import numpy as np
import json
import string
import borrelia_data_utils as bdu
from docplex.mp.model import Model

# ---------------------------------------------------------------------------------------

def collapse_strains(STRAINS_DATA):
    RES_AUX = {'##'.join([locus+':'+nuc for locus,nuc in profile.items()]):[] for profile in STRAINS_DATA.values()}
    for st_id,profile in STRAINS_DATA.items():
        RES_AUX['##'.join([locus+':'+nuc for locus,nuc in profile.items()])].append(st_id)
    RES = {}
    for profile,st_ids in RES_AUX.items():
        RES['_'.join(st_ids)] = {locus_nuc.split(':')[0]:locus_nuc.split(':')[1] for locus_nuc in profile.split('##')}
    return(RES)

def read_phasing(PHASING_FILE):
    PHASING_DATA = {}
    PHASING_STREAM = open(PHASING_FILE,'r').readlines()
    for read_data in PHASING_STREAM:
        read_data1 = read_data.rstrip().split('\t')
        read_id = read_data1[0]
        read_content = read_data1[1].split(' ')
        if len(read_content)>1:
            PHASING_DATA[read_id] = {r.split(':')[0]:r.split(':')[1] for r in read_content}
    return(PHASING_DATA)

def phasing_keys(PHASING_DATA,ALLELES_KEYS):
    PHASING_KEYS = []
    for r,content in PHASING_DATA.items():
        to_keep=True
        for k,a in content.items():
            if a  not in ALLELES_KEYS[k]:
                to_keep = False
        if to_keep:
            PHASING_KEYS.append(r)
    return(PHASING_KEYS)
    
# ---------------------------------------------------------------------------------------
def ct_and(model,res,cname,*vars):
    _Y = {}
    nb_vars = len(vars)
    _Y[cname] = model.binary_var(name=cname+'_AUX')
    i=1
    for var in vars:
        model.add_constraint(_Y[cname] <= var,ctname=cname+'_'+str(i))
        i+=1
    model.add_constraint(_Y[cname] >= model.sum(vars) - nb_vars + 1,ctname=cname+'_'+str(i))
    i+=1
    # Does not work as it encodes AND with 2 _Y[cname] >= var1 + var2
    #_Y[cname] = model.logical_and(var1,var2)
    return(model.add_constraint(res == _Y[cname],ctname=cname+'_'+str(i)))

def ct_or(model,res,cname,*vars):
    _Y = {}
    nb_vars = len(vars)
    _Y[cname] = model.binary_var(name=cname+'_AUX')
    i=1
    for var in vars:
        model.add_constraint(_Y[cname] >= var,ctname=cname+'_'+str(i))
        i+=1
    model.add_constraint(_Y[cname] <= model.sum(vars),ctname=cname+'_'+str(i))
    i+=1
    return(model.add_constraint(res == _Y[cname],ctname=cname+'_'+str(i)))


def diversity_lp_model3_create(STRAINS_DATA_ALL,ALLELES_DATA,PHASING_DATA,K1,K2,D_MAX,X_MIN,coeffUsage=1,coeffDev=-1,coeffDevMax=1,coeffDist=-1,coeffPhasing=1,LOCI_LIST=[0]):
    # Preprocessing input data
    if LOCI_LIST[0]==0:
        LOCI_KEYS  = list(ALLELES_DATA.keys())
    else:
        LOCI_KEYS  = [list(ALLELES_DATA.keys())[locus] for locus in LOCI_LIST]
    NB_LOCI      = len(LOCI_KEYS)
    ALLELES_KEYS = {k: [a for a in ALLELES_DATA[k].keys()] for k in LOCI_KEYS}
    STRAINS_DATA = collapse_strains({st_id:{locus:profile[locus] for locus in LOCI_KEYS} for st_id,profile in STRAINS_DATA_ALL.items()})
    STRAINS_KEYS = list(STRAINS_DATA.keys())
    PHASING_KEYS = phasing_keys(PHASING_DATA,ALLELES_KEYS)
    NB_READS     = float(len(PHASING_KEYS))

    # Setting the objective function weights
    # -- Known strains usage
    LAMBDA_USAGE = coeffUsage
    # -- Deviation from observed frequencies    
    if coeffDev == -2:  # average deviation per locus 
        LAMBDA_DEV = 1.0/(float(NB_LOCI))
    elif coeffDev == -1: # average deviation per locus normalized to belong to [0,1]
        MAX_DEV    = sum([max(ALLELES_DATA[k][a],1.0-ALLELES_DATA[k][a]) for k in LOCI_KEYS for a in ALLELES_KEYS[k]])
        LAMBDA_DEV = 1.0/(MAX_DEV) # average deviation per locus normalized to belong to [0,1]
    else:
        LAMBDA_DEV = coeffDev
    LAMBDA_MAX_DEV = coeffDevMax
    # -- Hamming distance
    if coeffDist == -2:
        LAMBDA_DIST = 1.0/(float(K2)) # average Hamming distance per derived strain
    elif coeffDist == -1:
        LAMBDA_DIST = 1.0/(float(K2*D_MAX)) # average Hamming distance per derived strain  normalized to belong to [0,1]
    else:
        LAMBDA_DIST = coeffDist

    LAMBDA_PHASING = coeffPhasing

    # Creating LP
    DIV_M3 = Model('Diversity_model_3')

    # Decision variables: strains proportions
    # X[j] = proportion of source abstract strain j
    K = K1+K2
    X_IDX  = list(range(K))
    XS_IDX = list(range(K1))   # Source abstract strains
    XD_IDX = list(range(K1,K)) # Derived abstract strains
    X      = DIV_M3.continuous_var_dict(X_IDX,lb=X_MIN,ub=1.0,name='X') 
    # The sum of the strains proportions should be at most 1
    DIV_M3.add_constraint(DIV_M3.sum(X[j] for j in X_IDX) == 1.0,ctname='X1a')
    MAX_STRAINS_USAGE = DIV_M3.continuous_var(lb=0.0,ub=1.0,name='MAX_STRAINS_USAGE')
    DIV_M3.add_constraint(MAX_STRAINS_USAGE == 1.0,ctname='OBJ_MAX_STRAIN_USAGE')
    #TOTAL_STRAINS_USAGE = DIV_M3.continuous_var(lb=0.0,ub=1.0,name='TOTAL_STRAINS_USAGE')
    #DIV_M3.add_constraint(TOTAL_STRAINS_USAGE == DIV_M3.sum(X[j] for j in X_IDX),ctname='OBJ_TOTAL_STRAIN_USAGE')

    # Forcing the proportions to be weakly increasing to reduce the search space
    for (j1,j2) in [(j1,j2) for j1 in XS_IDX for j2 in XS_IDX if j1>j2]:
        DIV_M3.add_constraint(X[j1]>=X[j2],ctname='XX1_'+str(j1)+'_'+str(j2))
    for (j1,j2) in [(j1,j2) for j1 in XD_IDX for j2 in XD_IDX if j1>j2]:
        DIV_M3.add_constraint(X[j1]>=X[j2],ctname='XX1_'+str(j1)+'_'+str(j2))

    # Decision variables: alleles content and proportions
    # Y[(k,a,j)] = proportion of allele a in position k in abstract strain j
    # Z[(k,a,j)] = 1 if allele a chosen in position k for abstract strain j
    Y_IDX = [(k,a,j) for k in LOCI_KEYS for a in ALLELES_KEYS[k] for j in X_IDX]
    Y     = DIV_M3.continuous_var_dict(Y_IDX,lb=0.0,ub=1.0,name='Y')
    Z     = DIV_M3.binary_var_dict(Y_IDX,name='Z')
    # Constraint: Ensuring one allele per locus for each abstract strain with non-zero abundance
    Z_IDX = [(k,j) for k in LOCI_KEYS for j in X_IDX]
    for (k,j) in Z_IDX:
        DIV_M3.add_constraint(DIV_M3.sum(Z[(k,a,j)] for a in ALLELES_KEYS[k]) >= X[j],ctname='Z1_'+'_'+str(k)+'_'+str(j))
        DIV_M3.add_constraint(DIV_M3.sum(Z[(k,a,j)] for a in ALLELES_KEYS[k]) <= 1.0,ctname='Z2_'+'_'+str(k)+'_'+str(j))
    # Constraints defining Y[(k,a,j)] as X[j]xZ[(k,a,j)]
    for (k,a,j) in Y_IDX:
        DIV_M3.add_constraint(Y[(k,a,j)] <= Z[(k,a,j)],ctname='YZ1_'+'_'+str(a)+'_'+str(k)+'_'+str(j))
        DIV_M3.add_constraint(Y[(k,a,j)] >= Z[(k,a,j)] - 1 + X[j],ctname='XYZ1_'+'_'+str(a)+'_'+str(k)+'_'+str(j))
        DIV_M3.add_constraint(Y[(k,a,j)] <= 1 - Z[(k,a,j)] + X[j],ctname='XYZ2_'+'_'+str(a)+'_'+str(k)+'_'+str(j))                
    # Ensuring that no two abstract strains are identical
    for (j1,j2) in [(j1,j2) for j1 in X_IDX for j2 in X_IDX if j1>j2]:
        DIV_M3.add_constraint(DIV_M3.sum(DIV_M3.sum(DIV_M3.abs(Z[k,a,j1]-Z[k,a,j2]) for a in ALLELES_KEYS[k]) for k in LOCI_KEYS) >= 1,ctname='XX2_'+str(j1)+'_'+str(j2))

    # Implied variables: allele usage per locus
    # U[(k,a)] = proportion of allele a at locus k as implied by the Y variables
    # DEV[(k,a)] = |U[(k,a)]-ALLELES_DATA[k][a]|
    U_IDX = [(k,a) for k in LOCI_KEYS for a in ALLELES_KEYS[k]]
    U     = DIV_M3.continuous_var_dict(U_IDX,lb=0.0,ub=1.0,name='U')
    DEV   = DIV_M3.continuous_var_dict(U_IDX,name='DEV')
    for (k,a) in U_IDX:
        DIV_M3.add_constraint(U[(k,a)] == DIV_M3.sum(Y[(k,a,j)] for j in X_IDX),ctname='U1_'+str(k)+'_'+str(a))
        DIV_M3.add_constraint(DEV[(k,a)] >= U[(k,a)]-ALLELES_DATA[k][a],ctname='DEV1_'+str(k)+'_'+str(a))        
        DIV_M3.add_constraint(DEV[(k,a)] >= ALLELES_DATA[k][a]-U[(k,a)],ctname='DEV2_'+str(k)+'_'+str(a))        

    # L_DEV[k] = sum of deviations on locus k
    # MAX_DEV = max deviation over a column
    K_IDX   = [k for k in LOCI_KEYS]
    L_DEV   = DIV_M3.continuous_var_dict(K_IDX,lb=0.0,ub=1.0,name='L_DEV')
    MAX_DEV = DIV_M3.continuous_var(name='MAX_DEV')
    for k in K_IDX:
        DIV_M3.add_constraint(L_DEV[k] == DIV_M3.sum(DEV[(k,a)] for a in ALLELES_KEYS[k]),ctname='L_DEV1'+str(k))
        # Note: we do not consider adding non-present nucleotides
        DIV_M3.add_constraint(L_DEV[k] <= MAX_DEV,ctname='MAX_DEV1'+str(k))

    # Implied variables: usage of known strains
    # E[(s,j)] = 0 if database strain s is not abstract strain j and X[j] otherwise
    E_IDX = [(s,j) for s in STRAINS_KEYS for j in X_IDX]
    E     = DIV_M3.continuous_var_dict(E_IDX,lb=0.0,ub=1.0,name='E')
    for (s,j) in E_IDX:
        for k in LOCI_KEYS:
            # Constraint ensuring that E is zero if some position can not match its nucleotide
            a = STRAINS_DATA[s][k] # a = allele at locus k for database strain s
            if a in ALLELES_KEYS[k]: # a is also the allele of locus k in abstract strain j: j could be s
                DIV_M3.add_constraint(E[(s,j)] <= Y[(k,a,j)],ctname='EY_'+str(s)+'_'+str(k)+'_'+str(a)+'_'+str(j))
            else: # abstract strain j can not be database strain s
                DIV_M3.add_constraint(E[(s,j)] == 0.0,ctname='EY1_'+str(s)+'_'+str(k)+'_'+str(a)+'_'+str(j))
            # E[(s,j)] = 0 if abstract strain j has nucleotide a and a is not the nucleotide at position k for database strain s
            for a in ALLELES_KEYS[k]:
                if a != STRAINS_DATA[s][k]:
                    DIV_M3.add_constraint(E[(s,j)] <= 1-Z[(k,a,j)],ctname='EY2_'+str(s)+'_'+str(k)+'_'+str(a)+'_'+str(j))

    # Implied variables: Hamming distance
    if K2>0:
        D_IDX = [(j1,j2) for j1 in XS_IDX for j2 in XD_IDX]
        PD    = DIV_M3.integer_var_dict(D_IDX,lb=0,ub=NB_LOCI,name='PD')
        D_MIN = DIV_M3.integer_var_dict(XD_IDX,lb=0,ub=D_MAX,name='D')
        # Constraints defining the pairwise Hamming distance
        for (j1,j2) in D_IDX:
            DIV_M3.add_constraint(PD[(j1,j2)] == 0.5* DIV_M3.sum(DIV_M3.abs(Z[(k,a,j1)]-Z[(k,a,j2)]) for (k,a) in U_IDX),ctname='D_'+str(j1)+'_'+str(j2))
            # Constraint defining the minimum Hamming distance
        for j2 in XD_IDX:
            DIV_M3.add_constraint(D_MIN[j2] == DIV_M3.min(PD[(j1,j2)] for j1 in XS_IDX),ctname='DMIN_'+str(j2))

    # Implied variables: reads whose phasing agrees with the chosen strains
    # PH1[(r,j)]=1 if read r agrees with abstract strain j
    PH1_IDX = [(r,j) for r in PHASING_KEYS for j in X_IDX]
    PH1     = DIV_M3.binary_var_dict(PH1_IDX,name='PH1')
    PH1_CT  = {}
    PH2_CT  = {}
    for (r,j) in PH1_IDX:
        PH1_CT[(r,j)] = ct_and(DIV_M3,PH1[(r,j)],'PH1C1_'+str(r)+'_'+str(j),*[Z[(k,PHASING_DATA[r][k],j)] for k in PHASING_DATA[r].keys()])
    # PH2[r] = 1 ifread r agrees with at least one abstract strain
    PH2 = DIV_M3.binary_var_dict(PHASING_KEYS,name='PH2')
    for r in PHASING_KEYS:
        PH2_CT[r] = ct_or(DIV_M3,PH2[r],'PH2C1_'+str(r),*[PH1[(r,j)] for j in X_IDX])
    # PHASING_TOTAL = number of reads in agreement with at least one strain
    PHASING_TOTAL = DIV_M3.integer_var(name='PHASING_TOTAL')
    DIV_M3.add_constraint(PHASING_TOTAL == DIV_M3.sum(PH2[r] for r in PHASING_KEYS),ctname='OBJ_PHASING_TOTAL1')
        

    # Objective function: deviation from input alleles proportions
    # DEV_TOTAL = sum of all DEV variables
    DEV_TOTAL = DIV_M3.continuous_var(name='OBJ_DT')
    DIV_M3.add_constraint(DEV_TOTAL == DIV_M3.sum(DEV[(k,a)] for (k,a) in U_IDX),ctname='OBJ_DEV_TOTAL1')
    
    # Objective function: usage of known strains
    KNOWN_STRAINS_USAGE = DIV_M3.continuous_var(lb=0.0,ub=1.0,name='OBJ_KSU')
    DIV_M3.add_constraint(KNOWN_STRAINS_USAGE == DIV_M3.sum(E[(s,j)] for (s,j) in E_IDX),ctname='OBJ_KNOWN_STRAIN_USAGE1')
    
    # Objective function: total Hamming distance
    TOTAL_HAMMING_D = DIV_M3.integer_var(lb=0,ub=K2*D_MAX,name='OBJ_THD')
    if K2>0:
        DIV_M3.add_constraint(TOTAL_HAMMING_D == DIV_M3.sum(D_MIN[j2] for j2 in XD_IDX),ctname='OBJ_TOTAL_HAMMING_D1')
    else:
        DIV_M3.add_constraint(TOTAL_HAMMING_D == 0,ctname='OBJ_TOTAL_HAMMING_D1')

    # Implied variables used for output of objective function
    OUT_DT  = DIV_M3.continuous_var(name='OUT_DT')
    OUT_MD  = DIV_M3.continuous_var(name='OUT_MD')
    DIV_M3.add_constraint(OUT_DT == DEV_TOTAL*LAMBDA_DEV,ctname='OBJ_DEV_TOTAL2')    
    DIV_M3.add_constraint(OUT_MD == MAX_DEV,ctname='OBJ_MAX_DEV2')    
    OUT_KSU = DIV_M3.continuous_var(name='OUT_KSU')
    DIV_M3.add_constraint(OUT_KSU == LAMBDA_USAGE*(MAX_STRAINS_USAGE-KNOWN_STRAINS_USAGE),ctname='OBJ_KNOWN_STRAIN_USAGE2')
    #DIV_M3.add_constraint(OUT_KSU == LAMBDA_USAGE*(TOTAL_STRAINS_USAGE-KNOWN_STRAINS_USAGE),ctname='OBJ_KNOWN_STRAIN_USAGE2')
    OUT_HD  = DIV_M3.continuous_var(name='OUT_HD')
    DIV_M3.add_constraint(OUT_HD == TOTAL_HAMMING_D*LAMBDA_DIST,ctname='OBJ_TOTAL_HAMMING_D2')
    OUT_PH  = DIV_M3.continuous_var(name='OUT_PH')
    DIV_M3.add_constraint(OUT_PH == LAMBDA_PHASING*(1.0-PHASING_TOTAL/NB_READS),ctname='OBJ_PHASING_TOTAL2')

    # Setting the objective function
    DIV_M3.set_objective('min',OUT_KSU+OUT_DT+OUT_MD+OUT_HD+OUT_PH)
    #LAMBDA_USAGE*(MAX_STRAINS_USAGE-KNOWN_STRAINS_USAGE)+(DEV_TOTAL*LAMBDA_DEV)+(MAX_DEV*LAMBDA_MAX_DEV)+(TOTAL_HAMMING_D*LAMBDA_DIST)+LAMBDA_PHASING*(1.0-PHASING_TOTAL/NB_READS))

    return(DIV_M3)

# ---------------------------------------------------------------------------------------
# Managing output

def read_json(json_file_name):
    with open(json_file_name) as json_data:
        RESULTS_ILP = {data['name']:data['value'] for data in json.load(json_data)['CPLEXSolution']['variables']}
    return(RESULTS_ILP)

def write_dictionary(dict):
    output='\n'.join(['_'.join(k)+':'+str(dict[k]) for k in dict.keys()])
    return(output+'\n')

def X_key(k):
    k1 = k.split('_')
    return((k1[0],'_'.join(k1[1:])))

def U_key(k):
    k1 = k.split('_')
    return((k1[0],k1[1],k1[2],k1[3]))

def E_key(k):
    k1 = k.split('_')
    l1 = len(k1)
    return((k1[0],'_'.join(k1[1:l1-1]),k1[l1-1]))

def Z_key(k):
    k1 = k.split('_')
    return((k1[0],k1[1],k1[2],k1[3],k1[4]))

def D_key(k):
    k1 = k.split('_')
    return((k1[0],k1[1]))

def PH_key(k):
    k1 = k.split('_')
    return((k1[0],k1[1],k1[2]))

def O_key(k):
    k1 = k.split('_')
    return((k1[0],k1[1]))

def read_json_3x(JSON_SOL_FILE_NAME,prec):
    # -- Reading the json file
    RESULTS = read_json(JSON_SOL_FILE_NAME)
    # -- Recording variables values in dictionaries
    # ---- Strain proportion variables
    X_VAL   = {}
    Z_VAL   = {}
    U_VAL   = {}
    E_VAL   = {}
    PH_VAL  = {}
    DEV_VAL = {}
    OBJ_VAL = {} 
    for k,v in RESULTS.items():
        val  = float(v)
        if k[0:2]=='X_' and val>=prec:            
            X_VAL[X_key(k)]=val
        elif k[0:2]=='U_' and val>=prec:            
            U_VAL[U_key(k)]=val
        elif k[0:2]=='Z_' and val>=prec:            
            Z_VAL[Z_key(k)]=val
        elif k[0:2]=='E_' and val>=prec:            
            E_VAL[E_key(k)]=val
        elif k[0:4]=='DEV_' and val>=prec:            
            DEV_VAL[U_key(k)]=val
        elif k[0:4]=='PH1_' and val>=prec:            
            PH_VAL[PH_key(k)]=val
        elif k[0:4]=='OUT_' and val>=prec:            
            OBJ_VAL[O_key(k)]=val
    return((X_VAL,U_VAL,Z_VAL,E_VAL,DEV_VAL,PH_VAL,OBJ_VAL))

def json2sol_3x(pref,prec):
    OUTPUT_FILE = open(pref+'_sol_aux.txt','w')    
    OUTPUT      = ''
    ILP_JSON = pref+'_sol.json'
    if os.path.isfile(ILP_JSON):
        (X_VAL,U_VAL,Z_VAL,E_VAL,DEV_VAL,PH_VAL,OBJ_VAL) = read_json_3x(ILP_JSON,prec)
        for DICT in  (X_VAL,U_VAL,Z_VAL,E_VAL,DEV_VAL,PH_VAL,OBJ_VAL):
            if len(list(DICT.keys()))>0:
                OUTPUT+=write_dictionary(DICT)
    else:
        X_VAL = None
        Y_VAL = None
        E_VAL = None
        Z_VAL = None
        DEV_VAL = None
        PH_VAL  = None
        OBJ_VAL = None
    OUTPUT_FILE.write(OUTPUT)
    OUTPUT_FILE.close()

def txtsol2summary_3(pref):
    LOCI = ['clpA','clpX','nifS','pepX','pyrG','recG','rplB','uvrA']

    X   = {}
    E   = {}
    SEQ = {}
    POS = {locus:[] for locus in LOCI}
    DEV = {locus:{} for locus in LOCI}
    OUT = {'DT':0.0,'MD':0.0,'KSU':0.0,'PH':0.0}
    
    # Reading the results file
    INPUT = open(pref+'_sol_aux.txt','r').readlines()
    for l in INPUT:
        l1 = l.rstrip().split(':')
        var = l1[0].split('_')
        val = float(l1[1])
        
        if var[0]=='X':
            st = var[1]
            X[st]   = val
            SEQ[st] = {locus:{} for locus in LOCI}
            E[st]   = 'Novel'
        elif var[0]=='Z':
            locus,pos,nuc,st = var[1],int(var[2]),var[3],var[4]
            SEQ[st][locus][pos]=nuc
            if st=='0':
                POS[locus].append(pos)
                DEV[locus][pos] = 0.0
        elif var[0]=='E':
            lg1 = len(var)
            st1,st2 = var[1],var[lg1-1]
            E[st2] = st1
        elif var[0]=='DEV':
            locus,pos = var[1],int(var[2])
            DEV[locus][pos]+=val
        elif var[0]=='OUT':
            OUT[var[1]] = val

    for locus in LOCI:
        POS[locus].sort()

    # Printing the output
    OUTPUT = open(pref+'_sol.txt','w')
    pos_idx=1
    OUTPUT.write('#\t')
    for locus in LOCI:
        for pos in POS[locus]:
            OUTPUT.write(str(pos_idx)+':'+locus+'_'+str(pos)+' ')
            pos_idx+=1
    OUTPUT.write('\n')
    OUTPUT.write(' '+'\t'+'    '+'\t'+'     '+'\t')
    pos_idx=1
    for locus in LOCI:
        for pos in POS[locus]:
            OUTPUT.write(str(round(pos_idx,2)).ljust(3))
            pos_idx+=1
    OUTPUT.write('\n')
    for st in X.keys():
        abundance = str(round(X[st],2))
        existing  = E[st]        
        OUTPUT.write(st+'\t'+abundance+'\t'+existing+'\t')
        for locus in LOCI:
            for pos in POS[locus]:
                OUTPUT.write(SEQ[st][locus][pos]+'  ')
        OUTPUT.write('\n')
    OUTPUT.write('Dev'+'\t\t\t')
    for locus in LOCI:
        for pos in POS[locus]:
            OUTPUT.write(str(round(DEV[locus][pos],2)).replace('0.','').ljust(3))
    OUTPUT.write('\n')
    OUTPUT.write('Objective function: '+str(OUT['KSU']+OUT['DT']+OUT['MD']+OUT['PH'])+'\n')
    OUTPUT.write('Strain usage penalty: '+str(OUT['KSU'])+'\n')
    OUTPUT.write('Normalized cumulated deviation penalty: '+str(OUT['DT'])+'\n')
    OUTPUT.write('Maximum deviation penalty: '+str(OUT['MD'])+'\n')
    OUTPUT.write('Ratio of phased reads: '+str(OUT['PH'])+'\n')
    OUTPUT.close()



# ---------------------------------------------------------------------------------------
if __name__ == "__main__":
    STRAINS_DATA_FILE = sys.argv[1]
    ALLELES_DATA_FILE = sys.argv[2]
    PHASING_DATA_FILE = sys.argv[3]
    OUTPUT_PREF       = sys.argv[4]
    K1                = int(sys.argv[5])
    K2                = int(sys.argv[6])
    D_MAX             = int(sys.argv[7])
    X_MIN             = float(sys.argv[8])
    LAMBDA_USAGE      = int(sys.argv[9])
    LAMBDA_DEV        = int(sys.argv[10])
    LAMBDA_DEV_MAX    = int(sys.argv[11])
    LAMBDA_DIST       = int(sys.argv[12])
    PREC = 0.00001

    ALLELES_DATA                  = bdu.read_pos_nuc_proportions(ALLELES_DATA_FILE)
    (STRAINS_SEQ,STRAINS_WEIGHTS) = bdu.read_sts_seq_sts_weights(STRAINS_DATA_FILE)
    STRAINS_DATA                  = {ST_ID:STRAINS_SEQ[ST_ID] for ST_ID in list(STRAINS_SEQ.keys())}
    PHASING_DATA                  = read_phasing(PHASING_DATA_FILE)
    
    M3 = diversity_lp_model3_create(STRAINS_DATA,ALLELES_DATA,PHASING_DATA,K1,K2,D_MAX,X_MIN,LAMBDA_USAGE,LAMBDA_DEV,LAMBDA_DEV_MAX,LAMBDA_DIST)
    sys.stdout = open(OUTPUT_PREF+'.log','w')
    M3.print_information()
    sys.stdout.flush()
    M3.export_as_lp(OUTPUT_PREF+'.lp')
    M3.set_time_limit(24*60*60)
    M3.parameters.mip.tolerances.mipgap = 0.02
    M3.parameters.mip.strategy.lbheur = 1
    M3.parameters.emphasis.mip = 0
    SOL3 = M3.solve(log_output=True)
    sys.stdout.flush()
    if SOL3!=None:
        print('Solution found')
        SOL3.export(OUTPUT_PREF+'_sol.json',format='json')
        #(_X,_Z,_U,_E,_DEV,_O) = read_json_3x(OUTPUT_PREF+'_sol.json',PREC)
        json2sol_3x(OUTPUT_PREF,PREC)
        txtsol2summary_3(OUTPUT_PREF)
    else:
        print('No solution found')
    sys.stdout.flush()
