# This code creates a lognormalised object with the chimera unspliced counts

library(Matrix)
library(scran)
library(Rtsne)
library(BiocParallel)
ncores = 16
mcparam = SnowParam(workers = ncores)
register(mcparam)


sce_s = readRDS('./sce_spliced.rds')

counts(sce_s) = assay(sce_s)

clusts = as.numeric(quickCluster(sce_s, method = "igraph", min.size = 100, BPPARAM = mcparam))

min.clust = min(table(clusts))/2
new_sizes = c(floor(min.clust/3), floor(min.clust/2), floor(min.clust))
sce_s = computeSumFactors(sce_s, clusters = clusts, sizes = new_sizes, max.cluster.size = 3000)

write.table(sizeFactors(sce_s), quote = F, col.names = F, row.names = F, file = "./sizefactors_sp.tab")

sfs = sizeFactors(sce_s)
sce_s <- scater::logNormCounts(sce_s,size_factors=sfs)

ychr_e = read.table("./ygenes.tab", stringsAsFactors = F)[,1]
genes_meta = read.table("./genes.tsv", stringsAsFactors = F)
rownames(genes_meta) = genes_meta$V1
ychr = genes_meta[ychr_e,]$V2

sce_s <- sce_s[!rownames(sce_s) %in% c('Xist', ychr),]

dec <- modelGeneVar(sce_s)
hvg <- getTopHVGs(dec, n = 5000)

write.csv(hvg, file = "./hvg_s.csv")

sce_s = sce_s[hvg,]

saveRDS(file="sce_s_ann.rds", sce_s)




