library(Matrix)
library(scran)
library(Rtsne)
library(BiocParallel)



ncores = 16
mcparam = SnowParam(workers = ncores)
register(mcparam)



counts = readMM("./raw_counts.mtx")



genes = read.table("./genes.tsv", stringsAsFactors = F)
meta = read.table("./meta.tab", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "$")
  
rownames(counts) = genes[,1] #ensembl
colnames(counts) = meta$cell

sce = SingleCellExperiment(assays = list("counts" = counts))

list_cells = read.csv('list_cells.csv',header = FALSE)

sce2 = sce[,list_cells$V1]



sce_g = readRDS('./gata1_sce.Rds')

rownames(sce2) = rownames(sce_g)[1:29452]


sce_g2 = sce_g[1:29452,]

mat = cbind(counts(sce2),counts(sce_g2))

big_sce <- SingleCellExperiment(assays = list("counts" = mat))

clusts = as.numeric(quickCluster(big_sce, method = "igraph", min.size = 100, BPPARAM = mcparam))

min.clust = min(table(clusts))/2
new_sizes = c(floor(min.clust/3), floor(min.clust/2), floor(min.clust))
big_sce = computeSumFactors(big_sce, clusters = clusts, sizes = new_sizes, max.cluster.size = 3000)


write.table(sizeFactors(big_sce), quote = F, col.names = F, row.names = F, file = "./sizefactors_big.tab")


saveRDS(file="unnorm_big_sce.rds", big_sce)

sfs = read.table("./sizefactors_big.tab", stringsAsFactors = F)[,1]

sizeFactors(big_sce) = sfs


big_sce <- scater::logNormCounts(big_sce,size_factors=sfs)

saveRDS(file="norm_big_sce.rds", big_sce)


ychr_e = read.table("./ygenes.tab", stringsAsFactors = F)[,1]
genes_meta = read.table("./genes.tsv", stringsAsFactors = F)
rownames(genes_meta) = genes_meta$V1
ychr = genes_meta[ychr_e,]$V2

p.val_th = 0.05

print(big_sce)

big_sce   <- big_sce[!rownames(big_sce) %in% c('Xist', ychr),]

 
print(big_sce)

dec <- modelGeneVar(big_sce)
hvg <- getTopHVGs(dec, fdr.threshold = p.val_th)

write.csv(hvg, file = "./hvg_noy.csv")


