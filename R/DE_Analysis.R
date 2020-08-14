# Editor : AJ

######### EdgeR - R implementation #########

# Document and resource used : 
#
#        -Adrian Alexa and Jorg Rahnenfuhrer (2019). topGO: Enrichment Analysis for Gene Ontology
#        -Adrian Alexa, Jorg Rahnenfuhrer (2020). Gene set enrichment analysis with topGO
#        -Marc Carlson (2019). GO.db: A set of annotation maps describing the entire Gene Ontology.
#        -Y. Chen, D. McCarthy, M. Ritchie, M. Robinson, and G. Smyth, “edgeR: 
#        differential analysis of sequence read count data  User’s Guide,”
#        -Defined colMap, ref. https://rpubs.com/aemoore62/TopGo_colMap_Func_Troubleshoot
#        -Gene Ontology or KEGG Pathway Analysis ref.
#         https://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/html/goana.html
#        -Table of Top GO Terms or Top KEGG Pathways ref.
#         https://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/html/topGO.html
#        -AnnotationDb-class: AnnotationDb objects and their progeny, methods etc. ref.
#         https://rdrr.io/bioc/AnnotationDbi/man/AnnotationDb-class.html
#        -Charity Law, et al. (2018). RNA-seq analysis is easy as 1-2-3 with limma, Glimma and edgeR ref.
#         https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
#
############
# purpose :      Differential expression, gene ontology and pathway analysis 
# parameters :   xlsx file (containing raw counts of the genes), preferred species, ontology domain,
#                 number of conditions and condition replications, index of the base condition
# return value : Differential expression analysis tabels (with and without ontologies and
#                pathway informations), Gene ontology analysis tabels, Pathway analysis
#                tabels, log-fold change against log-counts per million plots, dispersion
#                estimate plot, quasi-likelihood (QL) F-test plot
# details :
#           - In input data the first column should provide the gene ids and it should be 
#           ensembl gene id
#           - Number of replications for each condition should be the same
#           - Species, ontology domain, addrress of the directory which contains the genes' read 
#           counts, number of conditions (treatments), number of condition replications, and
#           index of the base condition should be provided in the section : 
#           ##  Setting initil Analysis' info  ##

############

############################################
#      Preparing the environment :
############################################

# Loading needed Libraries

libs = c("BiocManager", "statmod", "caret", "readxl", "ggplot2", "RColorBrewer", "stringr")
Bioc_libs = c("topGO","Rgraphviz","limma","edgeR","GO.db","AnnotationDbi")

for (i in libs){
  if( !is.element(i, .packages(all.available = TRUE)) ) {
    install.packages(i)
  }
  library(i,character.only = TRUE)
}

for (i in Bioc_libs){
  if( !is.element(i, .packages(all.available = TRUE)) ) {
    BiocManager::install(i)
  }
  library(i,character.only = TRUE)
}

##  Setting initil Analysis' info  ##
#####################################

# define the preferred species' name from the names in list_names 
species = "fly"

# define the preferred ontology domain
ont_domain = "cellular_component"           #  possible values : "cellular_component","molecular_function",
                                                                                                                        #" biological_process"

##   Setting the input data info   ##
#####################################

# Address of the excel file containing the genes' read counts
data_set_addr = "../../Table S3.xlsx" 

# Index or name sheet of the which contain the which contains the genes' read counts 
# if the spreadsheet contains only one sheet it should be 'NULL'
sheet_num = 5 

# Defining number of conditions (treatments) except the base condition,
# experiment's replications, and index of reference experiment in data set

treat_count = 3       # number of conditions
replic_count = 3      # number of experiment's replications
ref_idx = 4           # index of reference experiment in data set

#   Expression analysis thresholds  #
#####################################
p.value = 0.05       # P-value threshold
lfc = 2              # Fold change threshold

# Creatnig a 'named list' containing annotation mapping for supported species
# to be used for downloading prefered database from prefered species
#
# Supported species : human, mouse, rat, yeast, fly, arabidopsis, zebrafish, worm,
# bovine, chicken, canine, pig, rhesus, E coli strain K12 (ecoli_k12), xenopus, anopheles,
# anopheles, anopheles, chimp, malaria, E coli strain Sakai (ecoli_sakai) and Myxococcus xanthus DK
# (myxococcus)
#
# ###
#
# Note : KEGG database contains specific lists for different species of organisms. The genomic 
#        databases used in here may refer to a different species of the desired organisms
#

species_db_list <- list("org.Hs.eg.db","org.Mm.eg.db","org.Rn.eg.db","org.Sc.sgd.db","org.Dm.eg.db",
                        "org.At.tair.db","org.Dr.eg.db","org.Ce.eg.db","org.Bt.eg.db","org.Gg.eg.db",
                        "org.Cf.eg.db","org.Ss.eg.db","org.Mmu.eg.db","org.EcK12.eg.db","org.Xl.eg.db",
                        "org.Ag.eg.db","org.Pt.eg.db","org.Pf.plasmo.db","org.EcSakai.eg.db","org.Mxanthus.db")

kegg_db_list <- list("hsa","mmu","rno","sce","dme","ath","dre","sko","bta","gga",
                     "cfa","ssc","mcc","ecok","xla","aga","ptr","pfa","ecs","mxa")

list_names <- c("human","mouse","rat","yeast","fly","arabidopsis","zebrafish","worm","bovinae","chicken",
                "canine","pig","rhesus","ecoli_k12","xenopus","anopheles","chimp","malaria","ecoli_sakai",
                "myxococcus")

names(kegg_db_list) <- list_names
names(species_db_list) <- list_names

# Creatnig a 'named list' for mapping ontology domain names with TopGO ontology identifier

topgo_id <- list("CC","MF","BP")

names(topgo_id) <- c("cellular_component","molecular_function","biological_process")

# Installing the species' required database

if (!requireNamespace(species_db_list[[species]], quietly = TRUE))
    BiocManager::install(name)

library(species_db_list[[species]],character.only = TRUE)

# Creating related directories

dir.create(file.path("../../Results/","Data_process"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("../../Results/","Enrichment_tabels", species), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("../../Results/","GO_tabels", species, ont_domain), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("../../Results/","GO_figures", species, ont_domain), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("../../Results/","GO_enrich_terms", species, ont_domain), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("../../Results/","Pathway_enrich", species), recursive = TRUE, showWarnings = FALSE)

# Importing genes' read counts 

raw_df = read_excel(data_set_addr, sheet_num)

############################################
#    Differential expression analysis :
############################################

# Creating genes counts matrix, list of gene names

options(java.parameters = "-Xmx10000m")
genes=raw_df[,1]                            # list of gene names
raw_counts = as.matrix(raw_df[,2:(1+((treat_count+1)*replic_count))])       # genes counts matrix
row.names(raw_counts) <- unlist(genes)
remove(raw_df)

# Creating groups for edgeR DGEList

y = list()
group = c()
for (idx in 1:(treat_count+1)){
    for (idx2 in 1:replic_count){
    group <-c(group,idx)
    }
}

# Creating edgeR DGEList
y <- DGEList(counts=raw_counts, group=group, genes = genes)
names(y$genes)[1]<-paste("EnsemblGeneID")

remove(raw_counts)
remove(genes)

# Setting reference group as reference_idx

y$samples$group <- relevel(y$samples$group, ref=ref_idx)
levels(y$samples$group)

# Defining the design matrix 
design <- model.matrix(~group, data=y$samples)
design

# Filtering
keep <- filterByExpr(y,design)
y <- y[keep, , keep.lib.sizes=FALSE]

# Normalization
y <- calcNormFactors(y)
#y$samples

# Initial Data exploration

jpeg(paste("../../Results/Data_process/initial_data_exploration.jpg"), width = 800, height = 600)
plotMDS(y)
dev.off()
#plotMDS(y)

# Estimating dispersions
y <- estimateDisp(y, design)
#print(y$common.dispersion)

# Plotting dispersion estimates

jpeg(paste("../../Results/Data_process/dispersion_estimate.jpg"), width = 800, height = 600)
plotBCV(y)
dev.off()
#plotBCV(y)

## Starting Generalized Linear Models (glm) analysis ##
#######################################################

# Quasi-likelihood (QL) F-test

#fit <- glmFit(y, design, robust=TRUE) 

fit <- glmQLFit(y, design, robust=TRUE) 

# Quasi-likelihood (QL) F-test visualization

jpeg(paste("../../Results/Data_process/QL_F_test.jpg"), width = 800, height = 600)
plotQLDisp(fit)
dev.off()
#plotQLDisp(fit)

# Differential Expression test

qlf = list()
reg_res = list()    # regulation results
FT_res = list()     # FTest results

for (idx in 1:treat_count){
    qlf[[idx]] <- glmQLFTest(fit, coef=(idx+1))
    reg_res[[idx]] = as.data.frame(decideTests(qlf[[idx]], adjust.method = "BH", p.value = p.value, lfc = lfc))
    colnames(reg_res[[idx]])[1] <- "regulation"
    sum = summary(decideTests(qlf[[idx]], adjust.method = "BH", p.value = p.value, lfc = lfc))
    print(sum)
    write.csv(sum[c(1:3),], paste("../../Results/Enrichment_tabels/", species,
                                  "/DE_test_sum_", idx, ".csv", sep = ""), row.names = TRUE)
    FT_res[[idx]] = as.data.frame(topTags(qlf[[idx]], n=nrow(reg_res[[idx]])))
    remove(sum)
}
remove(fit)

# Plotting fold changes

for (idx in 1:treat_count){
    jpeg(paste("../../Results/Data_process/fold_change_", idx, ".jpg", sep = ""), width = 800, height = 600)
    plotMD(qlf[[idx]],main = "")
    abline(h=c(-1, 1), col="blue")
    dev.off()
}

# Creating differential expression output tables from
# regulation and FTest results

final_df = list()

for (idx in 1:treat_count){
    FT_res[[idx]] <- cbind(ID = rownames(FT_res[[idx]]), FT_res[[idx]])
    reg_res[[idx]] <- cbind(ID = rownames(reg_res[[idx]]), reg_res[[idx]])
    merge_df = merge(reg_res[[idx]], FT_res[[idx]], by = 'ID')
    final_df[[idx]] = subset(merge_df, merge_df$regulation!=0 )
    remove(merge_df)
}
remove(FT_res)
remove(reg_res)

############################################
#         Gene annotation :
############################################

# Adding Gene annotations to the data set

for (idx in 1:treat_count){
    keys = dplyr::pull(final_df[[idx]]["EnsemblGeneID"])
    multiVals <- function(x) paste(x,collapse=";")
    Symbol <- mapIds(get(species_db_list[[species]]), keys=keys, keytype="ENSEMBL", column="SYMBOL", multiVals=multiVals)
    EntrezID <- mapIds(get(species_db_list[[species]]), keys=keys, keytype="ENSEMBL", column="ENTREZID", multiVals=multiVals)
    final_df[[idx]]["symbol"] <- data.frame(Symbol=Symbol, stringsAsFactors=FALSE)
    final_df[[idx]]["EntrezID"] <- data.frame(EntrezID=EntrezID, stringsAsFactors=FALSE)
}

# Saving the raw differential expression analysis Results

for (idx in 1:treat_count){
    write.csv(final_df[[idx]], paste("../../Results/Enrichment_tabels/", species,
                                     "/final_df_raw_", idx, ".csv", sep = ""), row.names = FALSE)
}

############################################
# Gene ontology (GO) and pathway analysis :
############################################

########################################
# Gene ontology (GO) enrichment analysis

# Data preparation for Gene ontology (GO) enrichment analysis

go_list = list()

for (idx in 1:treat_count){
    go_list[[idx]] = final_df[[idx]][,"PValue"]
    names(go_list[[idx]]) = final_df[[idx]][,"symbol"]
}

# Defining function that returns TRUE/FALSE for p-values<0.05
test_func <- function(p_val){ return(p_val < 0.01)}

# GO to Symbol mappings Using the defined species' GO annotations database
go2sym <- annFUN.org(topgo_id[ont_domain], mapping = species_db_list[[species]], ID = "symbol")

# Performing GO enrichment analysis with both classic and
# eliminating genes approach using Kolmogorov–Smirnov-like
# statistic and aggregating the results in allRes table

GOdata = list()
go_res_cl_ks = list()
go_res_elim_ks = list()
allRes = list()

for (idx in 1:treat_count){
    GOdata[[idx]] <- new("topGOdata", ontology = paste(topgo_id[ont_domain]),
                         allGenes = go_list[[idx]], geneSelectionFun=test_func,
                         GO2genes=go2sym, nodeSize = 10, annot=annFUN.GO2genes)
    go_res_cl_ks[[idx]] <- runTest(GOdata[[idx]], algorithm ="classic", statistic = "ks")
    go_res_elim_ks[[idx]] <- runTest(GOdata[[idx]], algorithm="elim", statistic="ks")
    x <- geneData(go_res_elim_ks[[idx]])
    allRes[[idx]] = GenTable(GOdata[[idx]], classicKS = go_res_cl_ks[[idx]], elimKS = go_res_elim_ks[[idx]],
                      orderBy = "elimKS", ranksOf = "classicKS", topNodes = x['SigTerms'])
}

# Saving the GO enrichment analysis Results

for (idx in 1:treat_count){
    write.csv(allRes[[idx]], paste("../../Results/GO_tabels/", species, "/", ont_domain, 
                                   "/allRes_", idx, ".csv", sep = ""), row.names = TRUE)
    #print(GOdata[[idx]])
    #print(go_res_cl_ks[[idx]])
    #print(go_res_elim_ks[[idx]])
}

########################################
# Visualizing Results of GO enrichment analysis 

# Defining functions for visualization of GO enrichment analysis Results

colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}

# Visualization function
visualize <- function(go_res_cl_ks,go_res_elim_ks,GOdata) {
    
    pValue.classic <- score(go_res_cl_ks)
    pValue.elim <- score(go_res_elim_ks)[names(pValue.classic)]
    gstat <- termStat(GOdata, names(pValue.classic))
    gSize <- gstat$Annotated / max(gstat$Annotated) * 4
    gCol <- colMap(gstat$Significant)
    
    jpeg(paste("../../Results/GO_figures/", species, "/", ont_domain,"/elim_vs_classic_methods_differences_",
               idx, ".jpg", sep = ""), width = 800, height = 600)
    plot(pValue.classic, pValue.elim, xlab = "p-value classic", ylab = "p-value elim",
         pch = 19, cex = gSize, col = gCol)
    dev.off()
    
    #showSigOfNodes(GOdata, score(go_res_cl_ks), firstSigNodes = 10, useInfo = 'all')
    printGraph(GOdata, go_res_cl_ks, firstSigNodes = 10,
               fn.prefix = paste("../../Results/GO_figures/", species, "/", ont_domain, "/classic_method_GO_graph_",
                                 idx, sep = ""),useInfo = "all", pdfSW = TRUE)
    
    #showSigOfNodes(GOdata, score(go_res_elim_ks), firstSigNodes = 10, useInfo ='all')
    printGraph(GOdata, go_res_elim_ks, firstSigNodes = 10,
               fn.prefix = paste("../../Results/GO_figures/", species, "/", ont_domain,"/elim_method_GO_graph_",
                                 idx, sep = ""), useInfo = "all", pdfSW = TRUE)
}

# Visualizing GO enrichment analysis Results

for (idx in 1:treat_count){
    visualize(go_res_cl_ks[[idx]],go_res_elim_ks[[idx]],GOdata[[idx]])
}

########################################
# Adding GO Terms to the final data sets

# Creating a map of symbols to GO IDs
sym2go = inverseList(go2sym)

# Creating a map of symbols to GO terms

sym2term = list()

for (idx in 1:length(sym2go)){
    sym2term[names(sym2go[idx])] <- list(paste(Term(sym2go[[idx]])))
}

# Creating a matrix from sym2go map
sym2term_mx = as.matrix(sym2term)

# Defining the function to add GO Terms to the final data sets

term_add <- function(sym2term_mx,df) {
    sym2term_df = list()
    sym2term_df$a <- names(as.data.frame(sym2term_mx)[,"V1"])
    sym2term_df$V1 <- as.data.frame(sym2term_mx)
    sym2term_df = as.data.frame(sym2term_df)
    names(sym2term_df) = c("symbol","Go_Term")
    merge_df = merge(df, sym2term_df, by = 'symbol', all.x= TRUE)
    return(merge_df)
}

# Adding GO Terms to the final data sets using term_add function

for (idx in 1:treat_count){
    final_df[[idx]] = term_add(sym2term_mx,final_df[[idx]])
}

########################################
# Extraction of significantly enrichmented GO terms

# Defining functions for extraction of significantly enrichmented
# GO terms and saving the significantly enrichmented GO terms

go_enrich_func <- function(result,GOdata,idx,sign) {
    go_enrich <- GenTable(GOdata, KS=result, orderBy="KS", topNodes=20)
    #go_enrich <- go_enrich[go_enrich$KS<0.05,]
    go_enrich <- go_enrich[,c("GO.ID","Go_Term","KS")]
    go_enrich$KS <- as.numeric(go_enrich$KS)
    write.csv(go_enrich, paste("../../Results/GO_enrich_terms/", species, "/", ont_domain,
                               "/go_enrich_", sign, idx, ".csv", sep = ""), row.names = FALSE)
    #print(go_enrich)
}

# Extracting the significantly enrichmented GO terms and saving
# the significantly enrichmented GO terms using the go_enrich_func

for (idx in 1:treat_count){
    go_enrich_func(go_res_cl_ks[[idx]],GOdata[[idx]],idx,sign="classic")
    go_enrich_func(go_res_elim_ks[[idx]],GOdata[[idx]],idx,sign="elim")
}

############################################
# Pathway analysis with online KEGG database
############################################

## Preparing data for pathway analysis with kegg
# getting equivalent symbols and entrez-ids

#Symbol = list()
EntrezID = list()

for (idx in 1:treat_count){
    keys = dplyr::pull(qlf[[idx]]$genes)
    multiVals <- function(x) paste(x,collapse=";")
    EntrezID[[idx]] <- mapIds(get(species_db_list[[species]]), keys=keys,
                              keytype="ENSEMBL", column="ENTREZID", multiVals=multiVals)
    #Symbol[[idx]] <- mapIds(get(species_db_list[[species]]), keys=keys,
    #                         keytype="ENSEMBL", column="SYMBOL", multiVals=multiVals)
    #print(head(EntrezID[[idx]]))
}

## Preparing data for pathway analysis with kegg
# changing index in qlf from from ensembl-id to entrez-id

for (idx in 1:treat_count){    
    row.names(qlf[[idx]]$coefficients) <- EntrezID[[idx]]
    row.names(qlf[[idx]]$fitted.values) <- EntrezID[[idx]]
    row.names(qlf[[idx]]$unshrunk.coefficients) <- EntrezID[[idx]]
    names(qlf[[idx]]$deviance) <- EntrezID[[idx]]
    names(qlf[[idx]]$df.prior) <- EntrezID[[idx]]
    names(qlf[[idx]]$var.post) <- EntrezID[[idx]]
    qlf[[idx]]$genes <- EntrezID[[idx]]
    names(qlf[[idx]]$df.total) <- EntrezID[[idx]]
    #print(qlf[[idx]])
}

## Preparing data for pathway analysis with kegg
# Changing index in qlf$table from from ensembl-id to entrez-id
# and removing duplicate entrez-ids

for (idx in 1:treat_count){
    qlf[[idx]]$table$EntrezID <- EntrezID[[idx]]
    qlf[[idx]]$table <- qlf[[idx]]$table[!duplicated(qlf[[idx]]$table[,"EntrezID"]),]
    row.names(qlf[[idx]]$table) <- qlf[[idx]]$table$EntrezID
    qlf[[idx]]$table$EntrezID <- NULL
}

# performing the pathway analysis with online kegg database :
keg = list()

for (idx in 1:treat_count){
    keg[[idx]] <- kegga(qlf[[idx]], restrict.universe = FALSE, species.KEGG=kegg_db_list[[species]],
                        gene.pathway = NULL, pathway.names = NULL, convert=TRUE)
}

# Saving the significantly enrichmented pathways

path_number = 20 # number of significantly differentiated pathways to be 
                 # to be saved (use "Inf" to save all the pathways)

for (idx in 1:treat_count){
    up_path <- topKEGG(keg[[idx]], number=path_number, sort = "up")
    write.csv(up_path, paste("../../Results/Pathway_enrich/", species,
                               "/up_path_", idx, ".csv", sep = ""), row.names = FALSE)
    
    down_path <- topKEGG(keg[[idx]], number=path_number, sort = "down")
    write.csv(down_path, paste("../../Results/Pathway_enrich/", species,
                               "/down_path_", idx, ".csv", sep = ""), row.names = FALSE)
}

########################################
# Adding pathways to the final data sets

# Creating the map between keggid of each condition
# to pathway as a data frame

lkeggid2path = list()

for (idx in 1:treat_count){
    keggid2path = list()
    keggid2path$kegg_id <- rownames(keg[[idx]])
    keggid2path$kegg_id <- str_sub(keggid2path$kegg_id, start = -5, end = -1)
    keggid2path$pathway <- keg[[idx]]$Pathway
    keggid2path = as.data.frame(keggid2path)
    lkeggid2path[[idx]] <- keggid2path
}

# Extracting the database containing entrez_ids
# and KEGG_ids from the species database

species_db_str <- paste(species_db_list[[species]]) 
str_sub(species_db_str, -3) <- ""
path_db <- get(paste(species_db_str,'PATH',sep=""))

# Get the entrez_ids that are mapped to a KEGG_ids
mapped_genes <- mappedkeys(path_db)

# Creating a map of entrez_ids to KEGG_ids
entrez2keggid <- as.list(path_db[mapped_genes])

# inversing the entrez_ids to KEGG_ids map and Creating 
# the map between kegg_ids to entrez_ids as a data frame

keggid2entrez_df <- as.data.frame(as.matrix(inverseList(entrez2keggid)))
keggid2entrez_df$kegg_id <- rownames(keggid2entrez_df)

# Creating the map between entrez_ids of each condition to pathways
# as a matrix by merging lkeggid2path with keggid2entrez_df

entrez2path_mx = list()

for (idx in 1:treat_count){
    merge_df = merge(lkeggid2path[[idx]], keggid2entrez_df, by = 'kegg_id', all.x= TRUE)
    names(merge_df$V1) <- merge_df$pathway
    entrez2path_mx[[idx]] <- as.matrix(inverseList(merge_df$V1))
}

# Defining the function to add pathways to the final data sets
# using entrez2path_mx matrix

path_add <- function(entrez2path_mx,df) {
    entrez2kegg_df = list()
    entrez2kegg_df$a <- names(as.data.frame(entrez2path_mx)[,"V1"])
    entrez2kegg_df$V1 <- as.data.frame(entrez2path_mx)
    entrez2kegg_df = as.data.frame(entrez2kegg_df)
    names(entrez2kegg_df) = c("EntrezID","pathways")
    merge_df = merge(final_df[[idx]], entrez2kegg_df, by = 'EntrezID', all.x= TRUE)
    merge_df$'NA' <- merge_df$'NA.1' <- NULL
    return(merge_df)
}

# Adding pathways to the final data sets using path_add function

for (idx in 1:treat_count){
    final_df[[idx]] = path_add(entrez2path_mx[[idx]],final_df[[idx]])
}

# Saving the differential expression analysis Results

for (idx in 1:treat_count){
    final_df[[idx]]$Go_Term <- vapply(final_df[[idx]]$Go_Term, paste, collapse = ", ", character(1L))
    final_df[[idx]]$pathways <- vapply(final_df[[idx]]$pathways, paste, collapse = ", ", character(1L))
    write.csv(final_df[[idx]], paste("../../Results/Enrichment_tabels/", species,
                                     "/final_df_", ont_domain,"_", idx, ".csv", sep = ""), row.names = FALSE)
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
 [1] grid      stats4    parallel  stats     graphics  grDevices utils    
 [8] datasets  methods   base     

other attached packages:
 [1] org.Dm.eg.db_3.10.0  edgeR_3.28.0         limma_3.42.0        
 [4] Rgraphviz_2.30.0     topGO_2.37.0         SparseM_1.77        
 [7] GO.db_3.10.0         AnnotationDbi_1.48.0 IRanges_2.20.0      
[10] S4Vectors_0.24.0     Biobase_2.46.0       graph_1.64.0        
[13] BiocGenerics_0.32.0  stringr_1.4.0        RColorBrewer_1.1-2  
[16] readxl_1.3.1         caret_6.0-83         ggplot2_3.1.1       
[19] lattice_0.20-38      statmod_1.4.34       BiocManager_1.30.10 

loaded via a namespace (and not attached):
 [1] bit64_0.9-7        jsonlite_1.6       splines_3.6.3      foreach_1.4.4     
 [5] prodlim_2018.04.18 assertthat_0.2.1   blob_1.1.1         cellranger_1.1.0  
 [9] ipred_0.9-8        pillar_1.3.1       RSQLite_2.1.1      glue_1.3.1        
[13] uuid_0.1-2         digest_0.6.18      colorspace_1.4-1   recipes_0.1.5     
[17] htmltools_0.3.6    Matrix_1.2-17      plyr_1.8.4         timeDate_3043.102 
[21] pkgconfig_2.0.2    purrr_0.3.2        scales_1.0.0       gower_0.2.0       
[25] lava_1.6.5         tibble_2.1.1       generics_0.0.2     withr_2.1.2       
[29] repr_0.19.2        nnet_7.3-12        lazyeval_0.2.2     survival_2.44-1.1 
[33] magrittr_1.5       crayon_1.3.4       memoise_1.1.0      evaluate_0.13     
[37] nlme_3.1-139       MASS_7.3-51.3      class_7.3-15       tools_3.6.3       
[41] data.table_1.12.2  matrixStats_0.54.0 locfit_1.5-9.4     munsell_0.5.0     
[45] compiler_3.6.3     rlang_0.3.4        pbdZMQ_0.3-3       iterators_1.0.10  
[49] IRkernel_0.8.15    base64enc_0.1-3    gtable_0.3.0       ModelMetrics_1.2.2
[53] codetools_0.2-16   DBI_1.0.0          reshape2_1.4.3     R6_2.4.0          
[57] lubridate_1.7.4    dplyr_0.8.0.1      bit_1.1-14         stringi_1.4.6     
[61] IRdisplay_0.7.0    Rcpp_1.0.1         rpart_4.1-15       tidyselect_0.2.5  
'''