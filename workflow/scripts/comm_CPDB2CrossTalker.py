# source: https://github.com/CostaLab/CrossTalkeR/blob/master/CellPhoneDB%20Tutorial.md

import sys
import pandas as pd
import csv

### INTERNAL VARIABLES

## INPUT
input_file1 = sys.argv[1]
input_file2 = sys.argv[2]
## OUTPUT
output_file = sys.argv[3]

#input_file1 = "out/comm/Groups/young/AA/combined_significant_means.tsv"
#input_file2 = "out/comm/Groups/young/AA/idx2cluster.tsv"
#output_file = "/tmp/test.csv"

def correct_lr(data):
    '''
    Invert the RL to LR and R1R2 to r2>r1
    '''
    def swap(a,b): return b,a
    data = data.to_dict('index')
    for k,v in data.items():
        if v['isReceptor_fst'] and v['isReceptor_scn']:
            v['isReceptor_fst'],v['isReceptor_scn'] = swap(v['isReceptor_fst'],v['isReceptor_scn'])
            v['Ligand'],v['Receptor'] = swap(v['Ligand'],v['Receptor'])
            v['Ligand.Cluster'],v['Receptor.Cluster'] = swap(v['Ligand.Cluster'],v['Receptor.Cluster'])
        elif v['isReceptor_fst'] and not v['isReceptor_scn']:
            v['isReceptor_fst'],v['isReceptor_scn'] = swap(v['isReceptor_fst'],v['isReceptor_scn'])
            v['Ligand'],v['Receptor'] = swap(v['Ligand'],v['Receptor'])
            v['Ligand.Cluster'],v['Receptor.Cluster'] = swap(v['Ligand.Cluster'],v['Receptor.Cluster'])
    res_df = pd.DataFrame.from_dict(data,orient='index')
    return (res_df)
def cpdb2df(data,clsmapping):
    data = data.fillna(0)
    df_data = {}
    df_data['Ligand'] = []
    df_data['Receptor'] = []
    df_data['Ligand.Cluster'] = []
    df_data['Receptor.Cluster'] = []
    df_data['isReceptor_fst'] = []
    df_data['isReceptor_scn'] = []
    df_data['MeanLR'] = []
    for i in range(data.shape[0]):
        pair = list(data['interacting_pair'])[i].split('_')
        for j in range(data.iloc[:,12:].shape[1]):
            c_pair = list(data.columns)[j+12].split('|')
            if float(data.iloc[i,j+12]) !=0.0:
                df_data['Ligand'].append(pair[0])
                df_data['Receptor'].append(pair[1])
                df_data['Ligand.Cluster'].append(clsmapping[c_pair[0]])
                df_data['Receptor.Cluster'].append(clsmapping[c_pair[1]])
                df_data['isReceptor_fst'].append(list(data['receptor_a'])[i])
                df_data['isReceptor_scn'].append(list(data['receptor_b'])[i])
                df_data['MeanLR'].append(data.iloc[i,j+12])
    data_final = pd.DataFrame.from_dict(df_data)
    return(data_final)

def cpdb2df_nocls(data):
    '''
   		When the cluster name is used on CPDB
    '''
    data = data.fillna(0)
    df_data = {}
    df_data['Ligand'] = []
    df_data['Receptor'] = []
    df_data['Ligand.Cluster'] = []
    df_data['Receptor.Cluster'] = []
    df_data['isReceptor_fst'] = []
    df_data['isReceptor_scn'] = []
    df_data['MeanLR'] = []
    for i in range(data.shape[0]):
        pair = list(data['interacting_pair'])[i].split('_')
        for j in range(data.iloc[:,12:].shape[1]):
            c_pair = list(data.columns)[j+12].split('|')
            if float(data.iloc[i,j+12]) !=0.0:
                df_data['Ligand'].append(pair[0])
                df_data['Receptor'].append(pair[1])
                df_data['Ligand.Cluster'].append(c_pair[0])
                df_data['Receptor.Cluster'].append(c_pair[1])
                df_data['isReceptor_fst'].append(list(data['receptor_a'])[i])
                df_data['isReceptor_scn'].append(list(data['receptor_b'])[i])
                df_data['MeanLR'].append(data.iloc[i,j+12])
    data_final = pd.DataFrame.from_dict(df_data)
    return(data_final)

## Load data
s1 = pd.read_csv(input_file1,sep='\t')
#dict with the mapping
with open(input_file2) as f:
    next(f) # skip the header
    reader = csv.reader(f, skipinitialspace=True, delimiter="\t")
    num_to_clust = dict(reader)
   
## Convert
s1_filtered = cpdb2df(s1,num_to_clust)
s1_filtered = correct_lr(s1_filtered)

## Save
s1_filtered.to_csv(output_file)
