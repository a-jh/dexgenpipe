# Editor : AJ

######### EdgeR - R implementation #########

# Document and resource used : 
#
#        -Adrian Alexa and Jorg Rahnenfuhrer (2019). topGO: Enrichment Analysis for Gene Ontology
#        -Marc Carlson (2019). GO.db: A set of annotation maps describing the entire Gene Ontology.
#        -Y. Chen, D. McCarthy, M. Ritchie, M. Robinson, and G. Smyth, “edgeR: 
#        differential analysis of sequence read count data  User’s Guide,”
#        -Defined colMap, ref. https://rpubs.com/aemoore62/TopGo_colMap_Func_Troubleshoot

############
# purpose : gene ontology and differential expression analysis 
# parameters : xlsx file (containing raw counts of the genes)
# return value : differential expression analysis tabels, Gene ontology analysis tabels
#                  log-fold change against log-counts per million plots, dispersion
#                  estimate plot, quasi-likelihood (QL) F-test plot
# details :
#
############

# differential expression analysis :

# Loading needed Libraries

library(topGO)
library(statmod)
library(Rgraphviz)

library(limma)
library(edgeR)

library(readxl)

library(GO.db)
library(org.Dm.eg.db)
library(org.Hs.eg.db)

# Creating related directories

dir.create("../../Results/Figures/", showWarnings = FALSE)
dir.create("../../Results/Enrichment_tabels/", showWarnings = FALSE)
dir.create("../../Results/GO_tabels/", showWarnings = FALSE)
dir.create("../../Results/GO_figures/", showWarnings = FALSE)
dir.create("../../Results/GO_enrich_terms/", showWarnings = FALSE)

# Defining number of conditions (treatments), experiment's replications,
# and index of reference experiment in data set

treat_count = 3
replic_count = 3
ref_idx = 4

# Importing genes read counts 

raw_df = read_excel("../../Table S3.xlsx", 5)

# creating genes counts matrix, list of gene names

options(java.parameters = "-Xmx10000m")
genes=raw_df[,1]                            # list of gene names
raw_counts = data.matrix(raw_df[,2:13])     # genes counts matrix
remove(raw_df)

# Creating groups edgeR DGEList

y = list()
group = c()
for (idx in 1:(treat_count+1)){
    for (idx2 in 1:replic_count){
    group <-c(group,idx)
    }
}

# Creating edgeR DGEList
y <- DGEList(counts=raw_counts, group=group, genes=genes)
remove(raw_counts)
remove(genes)

# Setting reference group as reference_idx

y$samples$group <- relevel(y$samples$group, ref=ref_idx)
levels(y$samples$group)

# defining the design matrix 
design <- model.matrix(~group, data=y$samples)
design

# Filtering
keep <- filterByExpr(y,design)
y <- y[keep, , keep.lib.sizes=FALSE]

# Normalization
y <- calcNormFactors(y)
#y$samples

# Initial Data exploration

jpeg(paste("../../Results/Figures/initial_data_exploration.jpg"), width = 600, height = 600)
plotMDS(y)
dev.off()
#plotMDS(y)

# Estimating dispersions
y <- estimateDisp(y, design)
#print(y$common.dispersion)

# Plotting dispersion estimates

jpeg(paste("../../Results/Figures/dispersion_estimate.jpg"), width = 600, height = 600)
plotBCV(y)
dev.off()
#plotBCV(y)

# Starting Generalized Linear Models (glm) analysis

# quasi-likelihood (QL) F-test

fit <- glmQLFit(y, design, robust=TRUE) 

# quasi-likelihood (QL) F-test visualization

jpeg(paste("../../Results/Figures/QL_F_test.jpg"), width = 600, height = 600)
plotQLDisp(fit)
dev.off()
#plotQLDisp(fit)

# Differential Expression test

qlf = list()
reg_res = list()    # regulation results
FT_res = list()     # FTest results

for (idx in 1:treat_count){
    qlf[[idx]] <- glmQLFTest(fit, coef=(idx+1))
    reg_res[[idx]] = as.data.frame(decideTests(qlf[[idx]], adjust.method = "BH", p.value = 0.05, lfc = 2))
    colnames(reg_res[[idx]])[1] <- "regulation"
    sum = summary(decideTests(qlf[[idx]], adjust.method = "BH", p.value = 0.05, lfc = 2))
    print(sum)
    write.csv(sum[c(1:3),], paste("../../Results/Enrichment_tabels/DE_test_sum", idx, ".csv", sep = "_"), row.names = TRUE)
    FT_res[[idx]] = as.data.frame(topTags(qlf[[idx]], n=1000000))
    remove(sum)
}
remove(fit)

# plotting fold changes

for (idx in 1:treat_count){
    jpeg(paste("../../Results/Figures/fold change", idx, ".jpg", sep = "_"), width = 600, height = 600)
    plotMD(qlf[[idx]],main = "")
    abline(h=c(-1, 1), col="blue")
    dev.off()
}

# Creating differential expression output tables from
# regulation anf FTest results

final_df = list()

for (idx in 1:treat_count){
    FT_res[[idx]] <- cbind(ID = rownames(FT_res[[idx]]), FT_res[[idx]])
    reg_res[[idx]] <- cbind(ID = rownames(reg_res[[idx]]), reg_res[[idx]])
    merge_df = merge(reg_res[[idx]], FT_res[[idx]], by = 'ID')
    colnames(merge_df)[3] <- "EnsemblGeneID"
    final_df[[idx]] = subset(merge_df, merge_df$regulation!=0 )
    remove(merge_df)
}
remove(FT_res)
remove(reg_res)

# Gene annotation :

# Gene annotation

for (idx in 1:treat_count){
    keys = dplyr::pull(final_df[[idx]]["EnsemblGeneID"])
    multiVals <- function(x) paste(x,collapse=";")
    Symbol <- mapIds(org.Dm.eg.db, keys=keys, keytype="ENSEMBL", column="SYMBOL", multiVals=multiVals)
    final_df[[idx]]["symbol"] <- data.frame(Symbol=Symbol, stringsAsFactors=FALSE)
}

# Saving the differential expression analysis Results

for (idx in 1:treat_count){
    write.csv(final_df[[idx]], paste("../../Results/Enrichment_tabels/final_df", idx, ".csv", sep = "_"), row.names = FALSE)
}

# Gene ontology (GO) and pathway analysis :

# data preparation for Gene ontology (GO) enrichment analysis

go_list = list()

for (idx in 1:treat_count){
    go_list[[idx]] = final_df[[idx]][,"PValue"]
    names(go_list[[idx]]) = final_df[[idx]][,"symbol"]
}

# function that returns TRUE/FALSE for p-values<0.05
test_func <- function(p_val){ return(p_val < 0.01)}

# GO to Symbol mappings Using the org.Dm.eg.db annotations
Dm <- annFUN.org("BP", mapping = "org.Dm.eg.db", ID = "symbol")

# Performing GO enrichment analysis with both classic and
# eliminating genes approach using Kolmogorov–Smirnov-like
# statistic and aggregating the results in allRes table

GOdata = list()
go_res_cl_ks = list()
go_res_elim_ks = list()
allRes = list()

for (idx in 1:treat_count){
    GOdata[[idx]] <- new("topGOdata", ontology = "BP", allGenes = go_list[[idx]], geneSelectionFun=test_func,
                  GO2genes=Dm, nodeSize = 10, annot=annFUN.GO2genes)
    go_res_cl_ks[[idx]] <- runTest(GOdata[[idx]], algorithm ="classic", statistic = "ks")
    go_res_elim_ks[[idx]] <- runTest(GOdata[[idx]], algorithm="elim", statistic="ks")
    allRes[[idx]] = GenTable(GOdata[[idx]], classicKS = go_res_cl_ks[[idx]], elimKS = go_res_elim_ks[[idx]],
                      orderBy = "elimKS", ranksOf = "classicKS", topNodes = 150)
}

# Saving the GO enrichment analysis Results

for (idx in 1:treat_count){
    print(GOdata[[idx]])
    print(go_res_cl_ks[[idx]])
    print(go_res_elim_ks[[idx]])
    write.csv(allRes[[idx]], paste("../../Results/GO_tabels/allRes", idx, ".csv", sep = "_"), row.names = FALSE)
}

# defining functions for visualization of GO enrichment analysis Results

colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}

# visualization function
visualize <- function(go_res_cl_ks,go_res_elim_ks,GOdata) {
    
    pValue.classic <- score(go_res_cl_ks)
    pValue.elim <- score(go_res_elim_ks)[names(pValue.classic)]
    gstat <- termStat(GOdata, names(pValue.classic))
    gSize <- gstat$Annotated / max(gstat$Annotated) * 4
    gCol <- colMap(gstat$Significant)
    
    jpeg(paste("../../Results/GO_figures/elim_vs_classic_methods_differences", idx, ".jpg", sep = "_"), width = 600, height = 600)
    plot(pValue.classic, pValue.elim, xlab = "p-value classic", ylab = "p-value elim",
         pch = 19, cex = gSize, col = gCol)
    dev.off()
    
    
    #showSigOfNodes(GOdata, score(go_res_cl_ks), firstSigNodes = 10, useInfo = 'all')
    printGraph(GOdata, go_res_cl_ks, firstSigNodes = 10,
               fn.prefix = paste("../../Results/GO_figures/classic_method_GO_graph",
                                 idx, sep = "_"),
               useInfo = "all", pdfSW = TRUE)
    
    #showSigOfNodes(GOdata, score(go_res_elim_ks), firstSigNodes = 10, useInfo ='all')
    printGraph(GOdata, go_res_elim_ks, firstSigNodes = 10,
               fn.prefix = paste("../../Results/GO_figures/elim_method_GO_graph",
                                 idx, sep = "_"),
               useInfo = "all", pdfSW = TRUE)
}

# visualizing GO enrichment analysis Results

for (idx in 1:treat_count){
    visualize(go_res_cl_ks[[idx]],go_res_elim_ks[[idx]],GOdata[[idx]])
}

# defining functions for extraction of significantly
# enrichmented GO terms

go_enrich_func <- function(result,GOdata,idx,sign) {
    go_enrich <- GenTable(GOdata, KS=result, orderBy="KS", topNodes=20)
    #go_enrich <- go_enrich[go_enrich$KS<0.05,]
    go_enrich <- go_enrich[,c("GO.ID","Term","KS")]
    go_enrich$KS <- as.numeric(go_enrich$KS)
    write.csv(go_enrich, paste("../../Results/GO_enrich_terms/go_enrich", sign, idx, ".csv", sep = "_"), row.names = FALSE)
    #print(go_enrich)
}

# saving the significantly enrichmented GO terms

for (idx in 1:treat_count){
    go_enrich_func(go_res_cl_ks[[idx]],GOdata[[idx]],idx,sign="classic")
    go_enrich_func(go_res_elim_ks[[idx]],GOdata[[idx]],idx,sign="elim")
}

#END