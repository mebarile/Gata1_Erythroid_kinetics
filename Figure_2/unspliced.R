# This code creates a lognormalised object with the chimera unspliced counts

library(Matrix)
library(scran)
library(Rtsne)
library(BiocParallel)
ncores = 16
mcparam = SnowParam(workers = ncores)
register(mcparam)


sce_u = readRDS('sce_unspliced.rds')

counts(sce_u) = assay(sce_u)

clusts = as.numeric(quickCluster(sce_u, method = "igraph", min.size = 100, BPPARAM = mcparam))

min.clust = min(table(clusts))/2
new_sizes = c(floor(min.clust/3), floor(min.clust/2), floor(min.clust))
sce_u = computeSumFactors(sce_u, clusters = clusts, sizes = new_sizes, max.cluster.size = 3000)

sfs = sizeFactors(sce_u)

sce_u <- scater::logNormCounts(sce_u,size_factors=sfs)

ychr_e = read.table("./ygenes.tab", stringsAsFactors = F)[,1]
genes_meta = read.table("./genes.tsv", stringsAsFactors = F)
rownames(genes_meta) = genes_meta$V1
ychr = genes_meta[ychr_e,]$V2


sce_u <- sce_u[!rownames(sce_u) %in% c('Xist', ychr),]


dec <- modelGeneVar(sce_u)
hvg <- getTopHVGs(dec, n = 5000)

write.csv(hvg, file = "./hvg_u.csv")

sce_u = sce_u[hvg,]

print(sce_u)

saveRDS(file="sce_u_ann.rds", sce_u)




