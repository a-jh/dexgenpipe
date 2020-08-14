# Editor : AJ

######### STRING API - R implementation #########

# Document and resource used : 
#
#        -STRING API ref.
#        -https://string-db.org/cgi/help.pl?subpage=api%23getting-the-string-network-interactions
#

############
# purpose : analysis of protein-protein interactions of genes by utilizing STRING database
#           and adding String proteomics information to previous genomic analysis information
# parameters : csv files (containing EnsemblGeneID of the differentially expressed genes ),
#              number of conditions, species, ontology domain, STRING Api's input data attributes
# return value : 
#                  - Tabels of proteogenomic informations created by adding STRING proteomics information
#                  to previous genomic analysis information
#                  - STRING interaction network for identified proteins (Containing stringId_A (for protein A),
#                  stringId_B for (protein B), common protein name (protein A), common protein name (protein B),
#                  NCBI taxon identifier, combined score, gene neighborhood score, gene fusion score
#                  phylogenetic profile score, coexpression score, experimental score, database score
#                  and , textmining score)
# details :
#           - Species, ontology domain, addrress of the directory which contains genomic information 
#           should be provided in the section : 
#           ##  Setting initil Analysis' info  ##
#           - Data should be provided as a .csv file containing ensembl gene IDs
#
############

# Loading needed Libraries

libs = c("httr", "stringi", "stringr", "caret", "data.table")

for (i in libs){
  if( !is.element(i, .packages(all.available = TRUE)) ) {
    install.packages(i)
  }
  library(i,character.only = TRUE)
}

# Creating related directories

dir.create(file.path("../../Results/","Proteomics_interaction_tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("../../Results/","Proteogenomic_info"), recursive = TRUE, showWarnings = FALSE)

# Creatnig a 'named list' containing mapping between NCBI species Ids and
# supported species to be used for downloading proteomics data from STRING 

NCBI_id_list <- list(9606,10088,10114,4932,7227,3701,7955,10224,27592,9031,9615,9823,
                     9544,83333,8353,7164,9598,36329,83334,32)

list_names <- c("human","mouse","rat","yeast","fly","arabidopsis","zebrafish","worm","bovinae",
                "chicken","canine","pig","rhesus","ecoli_k12","xenopus","anopheles","chimp",
                "malaria","ecoli_sakai","myxococcus")

names(NCBI_id_list) <- list_names

##       Setting initil info       ##
#####################################

# define the preferred species' name from the names in list_names 
species = "fly"

# define the preferred ontology domain
ont_domain = "molecular_function"   #  possible values : "cellular_component","molecular_function",
                                                                                                                #" biological_process"

##   Setting the input data info   ##
#####################################

# Defining number of conditions (treatments) except the base condition,
# and intial STRING parameters
treat_count = 3

##  Setting STRING input data info ##
#####################################

string_api_url = "https://string-db.org/api"
#output_format = "tsv-no-header"
output_format = "tsv"
id_method = "get_string_ids"
net_method = "network"

# required_score and additional nodes are not used if default is "true"
default = "false"      # values : "true", "false"
required_score = 900 # threshold of significance to include a interaction,
                      # a number between 0 and 1000 (default depends on the network)
add_nodes = 0         # adds a number of proteins with to the network 
                      # based on their confidence score

NCBI_species_Id=NCBI_id_list[[species]]

# Load data from file

gene_dfl = list()

for (idx in 1:treat_count){
    gene_dfl[[idx]] <- read.table(file = paste("../../Results/Enrichment_tabels/", species,
                                               "/final_df_", ont_domain,"_", idx, ".csv", sep = "")
                                  ,sep = ",", header = TRUE,  stringsAsFactors = FALSE)
}

## Set parameters for STRING api identifier Mapping 

id_params = list()
lid_params = list()

for (idx in 1:treat_count){
    id_params$identifiers = paste(gene_dfl[[idx]]$EnsemblGeneID, collapse = "\r")
    id_params$species = NCBI_species_Id  # species NCBI identifier 
    id_params$limit = 1 # only one (best) identifier per input protein
    id_params$echo_query = 1 # see your input identifiers in the output
    #id_params$caller_identity" : "app_name" # your app name
    
    lid_params[[idx]] <- id_params
}

## Construct URL for STRING identifier mapping api request

id_request_url = paste(string_api_url,output_format,id_method, sep = "/")

## Call STRING

lid_results = list()

for (idx in 1:treat_count){
    id_results <- POST(id_request_url, body = lid_params[[idx]])
    lid_results[[idx]] <- id_results
    content(lid_results[[idx]], as = "text", encoding = "UTF-8") # automatically parses JSON
}

## Define function to read and parse the tsv results from STRING api

string_res <- function(results){
    
    string_result = list()
    check = 1
    
    lines = str_split(results, "\n")
    len = lengths(lines)-1
    
    lines = unlist(lines, use.names = FALSE)
    
    for (idx in 1:len){
        #string_result[[idx]] <- unlist(str_split(lines[[idx]], "\t"), use.names = FALSE)
        if (check == 1) {
            col_name = unlist(str_split(lines[[idx]], "\t"), use.names = FALSE)
        } else {
            data_row = unlist(str_split(lines[[idx]], "\t"), use.names = FALSE)
            string_result[[idx]] <- data_row
        }
        check = 0
    }
    df = as.data.frame(t(as.data.table(string_result)))
    colnames(df) <- col_name
    return(df)
}

# parsing the results from identifier mapping api and creating a dataframe

lid_parsed_results = list()

for (idx in 1:treat_count){
    lid_parsed_results[[idx]] = string_res(content(lid_results[[idx]], as = "text", encoding = "UTF-8"))
}

# Creating a new data frame by adding String proteomics information 
# to previous genomic analysis information

full_df = list()
for (idx in 1:treat_count){
    full_df[[idx]] = merge(gene_dfl[[idx]], lid_parsed_results[[idx]], by.x = "ID", by.y = "queryItem")
    full_df[[idx]]$queryIndex <- NULL
}

## Set parameters for STRING's interaction network's api

net_params = list()
lnet_params = list()

for (idx in 1:treat_count){
    
    if (default == "false"){
        net_params$required_score = required_score
        net_params$add_nodes = add_nodes        
    }
    
    net_params$identifiers = paste(lid_parsed_results[[idx]]$stringId, collapse = "%0d")
    net_params$species = NCBI_species_Id  # species NCBI identifier 
    #net_params$caller_identity" : "app_name" # your app name

    lnet_params[[idx]] <- net_params
}

## Construct URL for STRING interaction network's api request

net_request_url = paste(string_api_url,output_format,net_method, sep = "/")

## Call STRING

lnet_results = list()

for (idx in 1:treat_count){
    net_results <- POST(net_request_url, body = lnet_params[[idx]])
    lnet_results[[idx]] <- net_results
    content(lnet_results[[idx]], as = "text", encoding = "UTF-8") # automatically parses JSON
}

# parsing the results from interaction network's api and creating a dataframe
col_name = c("stringId_A", "stringId_B", "preferredName_A", "preferredName_B", "ncbiTaxonId",
            "combined_score", "neighborhood_score", "fusion_score", "phylogenetic_profile_score",
            "coexpression_score", "experimental_score", "database_score", "textmining_score")

lnet_parsed_results = list()

for (idx in 1:treat_count){
    lnet_parsed_results[[idx]] = string_res(content(lnet_results[[idx]], as = "text", encoding = "UTF-8"))
    colnames(lnet_parsed_results[[idx]]) <- col_name
}

# Saving the STRING api query results

for (idx in 1:treat_count){
    ending = paste(species, ont_domain, sep = "_")
    write.csv(full_df[[idx]], paste("../../Results/Proteogenomic_info/", ending,
                                           "_proteogenomic_info_", idx, ".csv", sep = ""), row.names = FALSE)
    
    write.csv(lnet_parsed_results[[idx]], paste("../../Results/Proteomics_interaction_tables/",
                                    ending, "_Proteomics_interaction_table_", idx, ".csv", sep = ""), row.names = FALSE)
}

### END

'''
sessionInfo()

R version 3.6.3
Platform: x86_64-conda_cos6-linux-gnu (64-bit)

Matrix products: default
BLAS/LAPACK: ../anaconda3/envs/Bioconda/lib/libopenblasp-r0.3.10.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] data.table_1.12.2 caret_6.0-83      ggplot2_3.1.1     lattice_0.20-38  
[5] stringr_1.4.0     stringi_1.4.6     httr_1.4.0       

loaded via a namespace (and not attached):
 [1] pbdZMQ_0.3-3       tidyselect_0.2.5   repr_0.19.2        purrr_0.3.2       
 [5] reshape2_1.4.3     splines_3.6.3      colorspace_1.4-1   generics_0.0.2    
 [9] stats4_3.6.3       htmltools_0.3.6    base64enc_0.1-3    survival_2.44-1.1 
[13] prodlim_2018.04.18 rlang_0.3.4        ModelMetrics_1.2.2 pillar_1.3.1      
[17] glue_1.3.1         withr_2.1.2        uuid_0.1-2         foreach_1.4.4     
[21] plyr_1.8.4         lava_1.6.5         timeDate_3043.102  munsell_0.5.0     
[25] gtable_0.3.0       recipes_0.1.5      codetools_0.2-16   evaluate_0.13     
[29] curl_3.3           class_7.3-15       IRdisplay_0.7.0    Rcpp_1.0.1        
[33] scales_1.0.0       IRkernel_0.8.15    ipred_0.9-8        jsonlite_1.6      
[37] digest_0.6.18      dplyr_0.8.0.1      grid_3.6.3         tools_3.6.3       
[41] magrittr_1.5       lazyeval_0.2.2     tibble_2.1.1       crayon_1.3.4      
[45] pkgconfig_2.0.2    MASS_7.3-51.3      Matrix_1.2-17      lubridate_1.7.4   
[49] gower_0.2.0        assertthat_0.2.1   iterators_1.0.10   R6_2.4.0          
[53] rpart_4.1-15       nnet_7.3-12        nlme_3.1-139       compiler_3.6.3    
'''