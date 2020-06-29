# Editor : AJ

######### STRING API - Python implementation #########

# Document and resource used : 
#
#        -STRING API ref.
#        -https://string-db.org/cgi/help.pl?subpage=api%23getting-the-string-network-interactions

############
# purpose : analysis of protein-protein interactions of genes by utilizing STRING database
# parameters : csv files (containing EnsemblGeneID of the differentially expressed genes )
# return value : 
#                  - Tabels of mapped EnsemblGene identifiers to the STRING identifiers
#                  (containing queried EnsemblGeneIDs, STRING identifier, NCBI taxon identifier,
#                  species name, common protein name and protein's annotation),
#                  - STRING interaction network for identified proteins (Containing stringId_A (for protein A),
#                  stringId_B for (protein B), common protein name (protein A), common protein name (protein B),
#                  NCBI taxon identifier, combined score, gene neighborhood score, gene fusion score
#                  phylogenetic profile score, coexpression score, experimental score, database score
#                  and , textmining score)
# details :
#
############

# import Libraries

from pathlib import Path
import pandas as pd
import requests 

# Creating related directories

Path("../../Results/String/").mkdir(parents=True, exist_ok=True)

'''
Defining number of conditions (treatments), and intial
STRING parameters
'''

treat_count = 3

string_api_url = "https://string-db.org/api"
#output_format = "tsv-no-header"
output_format = "tsv"
id_method = "get_string_ids"
net_method = "network"

NCBI_species_Id=7227

# Reading the data files

gene_dfl = []
gene_id = []
symbol =[]

for i in range(1,treat_count+1):
    gene_df = pd.read_csv("../../Results/tabels/final_df_" + str(i) + "_.csv")
    gene_dfl.append(gene_df)
    gene_id.append(gene_dfl[i-1].EnsemblGeneID)
    symbol.append(gene_dfl[i-1].symbol)

#del gene_df

## Set parameters for STRING api identifier Mapping 

lid_params = []

for i in range(0,treat_count):
    id_params = {
        "identifiers" : "\r".join(gene_id[i]), # your protein list
        "species" : NCBI_species_Id, # species NCBI identifier 
        "limit" : 1, # only one (best) identifier per input protein
        "echo_query" : 1, # see your input identifiers in the output
        #"caller_identity" : "app_name" # your app name
    }
    lid_params.append(id_params)

## Define function to read and parse the results from identifier mapping api

def string_res(result):
    result_df = pd.DataFrame()
    i = 0
    for line in result.text.strip().split("\n"):
        l = line.split("\t")
        if i == 0 :
            columns = l
        else :
            result_df = result_df.append([l], ignore_index=True)
        i = 1
    result_df.columns = columns
    return result_df

## Construct URL for STRING identifier mapping api request

id_request_url = "/".join([string_api_url, output_format, id_method])

## Call STRING

lid_results = []
lnet_results = []
stringId = []

for i in range(0,treat_count):
    id_results = requests.post(id_request_url, data=lid_params[i])
    lid_results.append(id_results)
    stringId.append(string_res(lid_results[i])['stringId'])
    
## Set parameters for STRING's interaction network's api

lnet_params = []

for i in range(0,treat_count):
    net_params = {
        "identifiers" : "%0d".join(stringId[i]), # your protein list
        "species" : NCBI_species_Id, # species NCBI identifier 
        #"echo_query" : 1, # see your input identifiers in the output
        #"caller_identity" : "app_name" # your app name
    }
    lnet_params.append(net_params)
    
## Set parameters for STRING's interaction network's api

net_request_url = "/".join([string_api_url, output_format, net_method])

for i in range(0,treat_count):
    net_results = requests.post(net_request_url, data=lnet_params[i])
    lnet_results.append(net_results)
    
for i in range(0,treat_count):
    string_res(lid_results[i]).to_csv('../../Results/String/string_id_res_' + str(i+1) +
                                      '.tsv',sep='\t',index=False)
    string_res(lnet_results[i]).to_csv('../../Results/String/string_net_res_' + str(i+1) +
                                       '.tsv',sep='\t',index=False)
    
#END