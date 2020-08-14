# Editor : AJ

######### CyREST - R implementation #########

# Document and resource used : 
#
#        -RCy3: Network biology using Cytoscape from within R ref.
#         https://f1000research.com/articles/8-1774/v1
#        -RCy3 Cytoscape functions ref.
#         https://cytoscape.org/RCy3/reference/index.html
#        -cytoscape-automation ref.
#         https://github.com/cytoscape/cytoscape-automation/blob/master/for-scripters/R/notebooks/Cancer-networks-and-data.Rmd
#        -RCy3: Functions to Access and Control Cytoscape ref.
#         https://rdrr.io/bioc/RCy3/
#        -https://cytoscape.org/cytoscape-tutorials/presentations/bioc2018_Rcy3_intro.html#/
#        -bioconductor RCy3 Documentation ref.
#         https://bioconductor.org/packages/release/bioc/html/RCy3.html
#        -clusterMaker2: Creating and Visualizing Cytoscape Clusters ref.
#         https://www.rbvi.ucsf.edu/cytoscape/clusterMaker2/clusterMaker2.shtml#mcl
#        -AutoAnnotate Cytoscape App 1.3 ref.
#         https://autoannotate.readthedocs.io/en/latest/index.html
#

############
# purpose : Creating analysis network from proteogenomic and protein-protein interactions' information 
#           obtained from differential expression, gene ontology, pathway analysis and STRING database
#           queries
# parameters : csv files (containing proteogenomic and proteomics interactions' information ),
#              number of conditions, species, ontology domain, network creation attributes,
#              visual styles attributes, clustring attributes, auto-annotation attributes,
#              network's layput attributes
# return value : 
#                - Summary of network and clustering evaluation
#                -
#                
# details :
#        - Cytoscape must be running whenever using RCy3 is used
#        - Initial setting such as Species, ontology domain, addrress of the directory which contains proteogenomic information,
#        visual styles attributes, mcl clustring attributes, and annotation attributes should be provided in the section : 
#        ##  Setting initil Analysis' info  ##
#        - Helpful Cytoscape commands :
#                commandsAPI(), commandsHelp("app name"), getLayoutNames(),
#                commandsRun('app name'), commandsHelp("app name subcommand")
#
############

############################################
#            Network analysis :
############################################ 

# Loading needed Libraries

libs = c("RCy3", "RColorBrewer")

for (i in libs){
  if( !is.element(i, .packages(all.available = TRUE)) ) {
    install.packages(i)
  }
  library(i,character.only = TRUE)
}

# Install needed Cytoscape apps

# App list
cyto_app_l <- c("clustermaker2", "autoannotate", "wordcloud", "Legend Creator")

# check if the apps are already installed
for(idx in 1:length(cyto_app_l)){
    
    check <- grep(commandsGET(paste("apps status app=\"", cyto_app_l[idx],"\"", sep="")),
                    pattern = "status: Installed")
    
    # install app if it is not installed
    if(length(check) == 0){
        print(cyto_app_l[idx])
        installApp(cyto_app_l[idx])
    }
}

# Creating related directories

dir.create(file.path("../../Results/","Cy_networks"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("../../Results/","Cy_manual_networks"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("../../Results/","Cy_network_images"), recursive = TRUE, showWarnings = FALSE)

# Creatnig a 'named list' containing mapping between NCBI species Ids and
# supported species to be used for downloading proteomics data from STRING 

NCBI_id_list <- list(9606,10088,10114,4932,7227,3701,7955,10224,27592,9031,9615,9823,
                     9544,83333,8353,7164,9598,36329,83334,32)

list_names <- c("human","mouse","rat","yeast","fly","arabidopsis","zebrafish","worm","bovinae",
                "chicken","canine","pig","rhesus","ecoli_k12","xenopus","anopheles","chimp",
                "malaria","ecoli_sakai","myxococcus")

names(NCBI_id_list) <- list_names

##         Setting initil info          ##
##########################################

# define the preferred species' name from the names in list_names 
species = "fly"

# define the preferred ontology domain
ont_domain = "molecular_function"           #  possible values : "cellular_component","molecular_function",
                                                                                                                        #" biological_process"

##     Setting the input data info      ##
##########################################

# Defining number of conditions (treatments) except the base condition,
treat_count = 3

##     Network creation attributes      ##
##########################################
id_name = "preferredName"               # Name of nodes' input table column being used as
                                        # nodes' input table id
interaction_source = "preferredName_A"  # Name of edges(interaction)' input table column being used as
                                        # edges' table source (Example of interaction input table can be 
                                        # the result of an query from STRING database)
interaction_target = "preferredName_B"  # Name of edges(interaction)' input table column being used as
                                        # edges' table target (Example of interaction input table can be 
                                        # the result of an query from STRING database)
interaction_weight = "combined_score"   # Name of edges(interaction)' input table column being used as
                                        # edges' table weight (Example of interaction input table can be 
                                        # the result of an query from STRING database)

## Setting the visual styles attributes ##
##########################################
BackgroundColor = '#888888'          # Hexadecimal8-bit per channel 
NodeShape = "ellipse"                # Use getArrowShapes() to see the list of shapes
NodeSize = 60                        # Numeric
NodeDefaultColor = "#AAAAAA"         # Hexadecimal8-bit per channel
EdgeDefaultLineWidth = 2

size_col = "F"                       # Numeric column - Name of nodes' table column being used 
                                     # to map the size of the nodes
color_col = "score"                  # Numeric column - Name of nodes' table column being used
                                     # to map the color of the nodes
NodeLabel = "preferredName"          # Name of nodes' table column being used to get the node names from

##   setting mcl clustring attributes   ##
##########################################
edgeCutOff = 0.4                     # Allows the user to set a custoff for edge weights
                                     # in MC clustering very dense networks
edgeWeighter = "None"                # Conversions that might be applied to the edge weights before clustering.
                                     # list values : 
                                     #      None: Don't do any conversion
                                     #      1/value: Use the inverse of the value.
                                     #      LOG(value): Take the log of the value.
                                     #      -LOG(value): Take the negative log of the value.
                                     #      SCPS: (Spectral Clustering of Protein Sequences) algorithm
inflation_parameter = 2.5            # The Granularity Parameter is also known as the inflation 
                                     # parameter or inflation value.This is the power used to 
                                     # inflate the matrix. Reasonable values range from ~1.8 to about 2.5

##   setting auto annotate attributes   ##
##########################################

clusterIdColumn = "__mclCluster"     # Allows to define clusters using node attributes for AutoAnnotate 
                                     # This allows clustering algorithms provided by other Cytoscape Apps 
                                     # or by external scripts to be used with AutoAnnotate.
createSingletonClusters = "false"    # If "true" a cluster will be created for each un-clustered node.
maxWords = 3                         # Number of words to be extrcted by annotation algorithm and used for annotation
useClusterMaker = "false"            # Alows AutoAnnotate to use other cytoscape apps for clustering the nodes

labelColumn = "Go_Term"              # Select a Node Column that will be used to calculate the cluster labels.

## setting network's layput attributes  ##
##########################################
layout = "force-directed" # other possible lay outs : attribute-circle, stacked-node-layout,
                          # degree-circle, circular , attributes-layout, kamada-kawai,
                          # force-directed, cose, grid, hierarchical, fruchterman-rheingold,
                          # isom

defaultSpringCoefficient = 0.00008
defaultSpringLength = 100

##  setting the results' output format  ##
##########################################

image_type = "SVG"           # PNG (default), JPEG, PDF, SVG, PS (PostScript).
network_type = "SIF"         # File type. SIF (default), CX, cyjs, graphML, NNF, xGMML. 

# Testing Cytoscape connection

base.url = "http://localhost:1234"
cytoscapeVersionInfo ()
cytoscapePing ()

## Load data from file

ending = paste(species, ont_domain, sep = "_")

# Load Proteogenomic information tabels 

Proteogenomic_dfl = list()

for (idx in 1:treat_count){
    Proteogenomic_dfl[[idx]] <- read.table(file = paste("../../Results/Proteogenomic_info/",
                                                      ending, "_proteogenomic_info_", idx, ".csv", sep = ""),
                                  ,sep = ",", header = TRUE,  stringsAsFactors = FALSE)
    Proteogenomic_dfl[[idx]]$score <- Proteogenomic_dfl[[idx]]$F * Proteogenomic_dfl[[idx]]$regulation
}

# Load Proteomics interaction tables' data

interaction_df = list()
node_id = list()
target_id = list()

for (idx in 1:treat_count){
    interaction_df[[idx]] <- read.table(file = paste("../../Results/Proteomics_interaction_tables/",
                                                 ending, "_Proteomics_interaction_table_", idx, ".csv", sep = ""),
                                  ,sep = ",", header = TRUE,  stringsAsFactors = FALSE)
}

# Creating intial proteomics interaction networks using imported tables' data

network_id = list()

for (idx in 1:treat_count){
    nodes <- data.frame(id=Proteogenomic_dfl[[idx]][,paste(id_name)],
                        Proteogenomic_dfl[[idx]],
                        stringsAsFactors=FALSE)
    edges <- data.frame(source=interaction_df[[idx]][,paste(interaction_source)],
                        target=interaction_df[[idx]][,paste(interaction_target)],
                        interaction_df[[idx]],
                        #interaction = interaction_list # optional to provide a list interaction types
                        weight=interaction_df[[idx]][,paste(interaction_weight)], # numeric
                        stringsAsFactors=FALSE)
    
    network_id[[idx]] <- createNetworkFromDataFrames(nodes, edges, 
                                                      title = paste("Proteogenomic_net_", idx, sep = ""),
                                                      collection = "String_DataFrames")
    #layoutNetwork('force-directed')
}

## Definning the visualiziation expression and setting the network Style

# Define new visual styles and style's default values

style.name = "new_style"
createVisualStyle(style.name)
setBackgroundColorDefault(BackgroundColor, style.name)
setNodeShapeDefault(NodeShape, style.name) 
setNodeSizeDefault(NodeSize, style.name)
setNodeColorDefault(NodeDefaultColor, style.name)
setEdgeLineWidthDefault(EdgeDefaultLineWidth, style.name)

for (idx in 1:treat_count){
    
    # Choosing the network
    setCurrentNetwork(network = getNetworkName(suid=as.numeric(network_id[[idx]])))
    
    #Setting new visual styles
    setVisualStyle(style.name)
    
    # Get the min and max differentiation score from node table and defining the range of values
    # of visualiziation size node expression data mapping.
    
    size = getTableColumns("node", size_col)
    size.values = c(min(size[,1],na.rm=TRUE),mean(size[,1],na.rm=TRUE),max(size[,1],na.rm=TRUE))
    
    # Get the min and max protein interaction score from edge table and defining the range of values
    # of visualiziation the edge line width expression data mapping .
    
    edge_weight = getTableColumns('edge', "combined_score" )
    min.edge_weight = min(edge_weight[,1],na.rm=TRUE)
    max.edge_weight = max(edge_weight[,1],na.rm=TRUE)
    edge_weight.values = c(min.edge_weight,
                           (min.edge_weight+max.edge_weight)*0.5,
                           max.edge_weight)

    # Get the min and max differentiation score from node table and defining the range of values
    # of visualiziation color node expression data mapping .
    
    node_col = getTableColumns('node', color_col )
    
    min.node_col = min(node_col[,1],na.rm=TRUE)
    max.node_col = max(node_col[,1],na.rm=TRUE)
    color.data.values = c(min.node_col,0,max.node_col)
    
    # picking colors with RColorBrewer package
    # display.brewer.all(length(data.values), colorblindFriendly=TRUE, type="div") # div,qual,seq,all
    node.colors <- c(rev(brewer.pal(length(color.data.values), "RdBu")))
    
    # constructing CyREST data objects for style mappings 
    # Mapping node size
    setNodeSizeMapping(size_col , size.values, 
                       sizes = c(30, 70, 150), mapping.type = "c", 
                       style.name = style.name)
    # Mapping node labels
    setNodeLabelMapping(NodeLabel, style.name)
    # Mapping node colors
    setNodeColorMapping(color_col, color.data.values, node.colors,
                        style.name=style.name)
    # Mapping edge line width
    setEdgeLineWidthMapping("combined_score", edge_weight.values, widths = c(1,3,6),
                            mapping.type = "c", style.name=style.name)
}

## Clustring with mcl Algorithm

new_node_table = list()
new_edge_table = list()

# setting mcl clustring attributes

adjustLoops = "true"
undirectedEdges = "true"
selectedOnly = "false"
createGroups = "true"
showUI = "true"

for (idx in 1:treat_count){
        
    # Choosing the network
    network = paste("Proteogenomic_net_", idx, sep = "")
    setCurrentNetwork(network = getNetworkName(suid=as.numeric(network_id[[idx]])))
    
    # creating clustering command
    clustermaker_url <- paste("cluster mcl edgeCutOff=", edgeCutOff, "network=", network, 
                              "edgeWeighter=", edgeWeighter, "adjustLoops=", adjustLoops, 
                              "inflation_parameter=", inflation_parameter, "edgeWeighter=",
                              edgeWeighter, "showUI=", showUI, "undirectedEdges=",
                              undirectedEdges)
    
    # running clustring command
    # commandsRun(clustermaker_url)
    commandsGET(clustermaker_url)
    
    #get the clustering results
    network = paste("Proteogenomic_net_", idx, "--clustered", sep = "")
    new_node_table[[idx]] <- getTableColumns(table= "node",network = as.numeric(network_id[[idx]]))
    new_edge_table[[idx]] <- getTableColumns(table= "edge",network = as.numeric(network_id[[idx]]))
}

# Applying network layput

if (TRUE){
    
    for (idx in 1:treat_count){
        
        # Choosing the network
        network = paste("Proteogenomic_net_", idx, "--clustered", sep = "")
        setCurrentNetwork(network = network)
        # Creating layout command
        
        # constructing CyREST data object for style mappings 
        if(layout == "force-directed"){
            
            layout <- paste(layout, "defaultSpringCoefficient=",
                              defaultSpringCoefficient,
                              "defaultSpringLength=",
                              defaultSpringLength)
            layoutNetwork(layout)
            }
        else {
            layoutNetwork(layout)
        }
    }
}

# Preforming auto-annotation

for (idx in 1:treat_count){
    
    # Choosing the network
    network = paste("Proteogenomic_net_", idx, "--clustered", sep = "")
    setCurrentNetwork(network = network)
    
    # Add list of words to be ignored during annotation 
    
    wordcloud_ignore1 <- paste("wordcloud ignore add", "value=", "null")
    wordcloud_ignore2 <- paste("wordcloud ignore add", "value=", "na")

    commandsGET(wordcloud_ignore1)
    commandsGET(wordcloud_ignore2)

    # Run the AutoAnnotate command
    annotate <- paste("autoannotate annotate-clusterBoosted", "network=", network,
                      "clusterIdColumn=", clusterIdColumn, "labelColumn=", labelColumn,
                      "maxWords=", maxWords, "useClusterMaker=", useClusterMaker)
    commandsGET(annotate)
    
    annotate <- paste("autoannotate layout", "network=", network,
                      "layout=", "grid")
    commandsGET(annotate)
}

# General Network and Cluster analysis and evaluation informtion

base_net = list()
clustered_net = list()

analyze_cmd <- paste("analyzer analyze directed= false")

for (idx in 1:treat_count){
    
    network = paste("Proteogenomic_net_", idx, sep = "")
    setCurrentNetwork(network = network)
    
    base_net[[idx]] <- commandsGET(analyze_cmd)
    
    network = paste("Proteogenomic_net_", idx, "--clustered", sep = "")
    setCurrentNetwork(network = network)
    
    clustered_net[[idx]] <- commandsGET(analyze_cmd)
}

# Saving the results

ending = paste(species, ont_domain, "-annotated_by", labelColumn, sep = "_")
AddrnName = paste("../../Results/Cy_networks/",
                      "Network_session_", ending, sep = "")
saveSession( filename = AddrnName)

for (idx in 1:treat_count){
    network = paste("Proteogenomic_net_", idx, "--clustered", sep = "")
    AddrnName = paste("../../Results/Cy_network_images/",
                      "Network_images_", ending, idx, sep = "")
    exportImage(
        filename = AddrnName,
        type = image_type,       # PNG (default), JPEG, PDF, SVG, PS (PostScript).
        resolution = 600,   # (numeric) The resolution of the exported image, in DPI. 
                            # Valid only for bitmap formats. The possible values are:
                            # 72 (default), 100, 150, 300, 600.
        network = network
        )
    
    AddrnName = paste("../../Results/Cy_networks/",
                      "Network_", ending, idx, sep = "")
    exportNetwork(
        filename = AddrnName, # Full path or path relavtive to current working directory
                              # in addition to the name of the file. 
        type = network_type,         # File type. SIF (default), CX, cyjs, graphML, NNF, xGMML.
        network = network     # (optional) Name or SUID of a network or view. Default is the "current" network active in Cytoscape.
    )
}

### END

'''
sessionInfo()

R version 3.6.3 (2020-02-29)
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
[1] RColorBrewer_1.1-2 RCy3_2.6.0        

loaded via a namespace (and not attached):
 [1] igraph_1.2.4.1      graph_1.64.0        Rcpp_1.0.1         
 [4] magrittr_1.5        BiocGenerics_0.32.0 uuid_0.1-2         
 [7] R6_2.4.0            httr_1.4.0          tools_3.6.3        
[10] parallel_3.6.3      R.oo_1.22.0         htmltools_0.3.6    
[13] digest_0.6.18       crayon_1.3.4        RJSONIO_1.3-1.1    
[16] IRdisplay_0.7.0     repr_0.19.2         base64enc_0.1-3    
[19] R.utils_2.8.0       curl_3.3            IRkernel_0.8.15    
[22] evaluate_0.13       pbdZMQ_0.3-3        compiler_3.6.3     
[25] R.methodsS3_1.7.1   stats4_3.6.3        XML_3.98-1.19      
[28] jsonlite_1.6        pkgconfig_2.0.2
'''