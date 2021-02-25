# This code corrects the mofa actors for batch effects

library(Matrix)
library(scran)
library(Rtsne)
library(BiocParallel)
require(irlba)
BPPARAM = SerialParam()
library(batchelor)

mata = read.csv('./files/factors_df_all.csv',row.names = 1, header= TRUE)
mata = t(mata[,3:10])


meta_a = read.table("/meta.tab", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "$")
meta_a$comb = paste(meta_a$stage,meta_a$sample,sep = '_')


list_cells = read.csv('list_cells.csv',header = FALSE)

rownames(meta_a) =meta_a$cell
meta_a = meta_a[list_cells$V1,]


met65_1 = mata[,meta_a$comb == 'E6.5_1']
met65_18 = mata[,meta_a$comb == 'E6.5_18']
met65_5 = mata[,meta_a$comb == 'E6.5_5']


met675_7 = mata[,meta_a$comb == 'E6.75_7']

met7_10 = mata[,meta_a$comb == 'E7.0_10']
met7_14 = mata[,meta_a$comb == 'E7.0_14']
met7_15 = mata[,meta_a$comb == 'E7.0_15']
met7_30 = mata[,meta_a$comb == 'E7.0_30']
met7_31 = mata[,meta_a$comb == 'E7.0_31']
met7_32 = mata[,meta_a$comb == 'E7.0_32']

met725_23 = mata[,meta_a$comb == 'E7.25_23']
met725_26 = mata[,meta_a$comb == 'E7.25_26']
met725_27 = mata[,meta_a$comb == 'E7.25_27']

met75_19 = mata[,meta_a$comb == 'E7.5_19']
met75_2 = mata[,meta_a$comb == 'E7.5_2']
met75_20 = mata[,meta_a$comb == 'E7.5_20']
met75_3 = mata[,meta_a$comb == 'E7.5_3']
met75_4 = mata[,meta_a$comb == 'E7.5_4']
met75_6 = mata[,meta_a$comb == 'E7.5_6']


met775_12 = mata[,meta_a$comb == 'E7.75_12']
met775_13 = mata[,meta_a$comb == 'E7.75_13']
met775_8 = mata[,meta_a$comb == 'E7.75_8']
met775_9 = mata[,meta_a$comb == 'E7.75_9']

met8_16 = mata[,meta_a$comb == 'E8.0_16']
met8_33 = mata[,meta_a$comb == 'E8.0_33']
met8_34 = mata[,meta_a$comb == 'E8.0_34']
met8_35 = mata[,meta_a$comb == 'E8.0_35']

met825_24 = mata[,meta_a$comb == 'E8.25_24']
met825_25 = mata[,meta_a$comb == 'E8.25_25']
met825_28 = mata[,meta_a$comb == 'E8.25_28']



met85_17 = mata[,meta_a$comb == 'E8.5_17']
met85_29 = mata[,meta_a$comb == 'E8.5_29']
met85_36 = mata[,meta_a$comb == 'E8.5_36']
met85_37 = mata[,meta_a$comb == 'E8.5_37']


atlas_pca = t(cbind(met65_1,met65_18,met65_5,met675_7,met7_10,met7_14,met7_15,met7_30,met7_31,met7_32 ,
                        met725_23,met725_26,met725_27,met75_19,met75_2,met75_20,met75_3,met75_4,met75_6,
                        met775_12,met775_13,met775_8,met775_9,met8_16,met8_33,met8_34,met8_35,met825_24,
                        met825_25,met825_28,met85_17,met85_29,met85_36,met85_37))


meta_a = meta_a[rownames(atlas_pca),]

order_df = meta_a[!duplicated(meta_a$sample), c("stage", "sample")]
order_df$ncells = sapply(order_df$sample, function(x) sum(meta_a$sample == x))
order_df$stage = factor(order_df$stage,
                          levels = rev(c("E8.5",
                                         "E8.25",
                                         "E8.0",
                                         "E7.75",
                                         "E7.5",
                                         "E7.25",
                                         "E7.0",
                                         "E6.75",
                                         "E6.5")))
order_df = order_df[order(order_df$stage, order_df$ncells, decreasing = TRUE),]
order_df$stage = as.character(order_df$stage)


set.seed(42)


timepoints = meta_a$stage
samples = meta_a$sample
timepoint_order = order_df$stage
sample_order = order_df$sample
  


#create nested list
  pc_list = lapply(unique(timepoints), function(tp){
    sub_pc = atlas_pca[timepoints == tp, , drop = FALSE]
    sub_samp = samples[timepoints == tp]
    list = lapply(unique(sub_samp), function(samp){
     sub_pc[sub_samp == samp, , drop = FALSE]
    })
    names(list) = unique(sub_samp)
    return(list)
  })

  names(pc_list) = unique(timepoints)
print('dd')
  #arrange to match timepoint order
  pc_list = pc_list[order(match(names(pc_list), timepoint_order))]
  pc_list = lapply(pc_list, function(x){
    x[order(match(names(x), sample_order))]
  })


# corr atlas
# corr tps

corr8_5 = reducedMNN(pc_list[[1]][[1]],pc_list[[1]][[2]],pc_list[[1]][[3]],pc_list[[1]][[4]],BPPARAM = BPPARAM)$corrected

corr8_25 = reducedMNN(pc_list[[2]][[1]],pc_list[[2]][[2]],pc_list[[2]][[3]],BPPARAM = BPPARAM)$corrected

corr8_0 = reducedMNN(pc_list[[3]][[1]],pc_list[[3]][[2]],pc_list[[3]][[3]],pc_list[[3]][[4]],BPPARAM = BPPARAM)$corrected
corr7_75 = reducedMNN(pc_list[[4]][[1]],pc_list[[4]][[2]],pc_list[[4]][[3]],pc_list[[4]][[4]],BPPARAM = BPPARAM)$corrected
corr7_5 = reducedMNN(pc_list[[5]][[1]],pc_list[[5]][[2]],pc_list[[5]][[3]],pc_list[[5]][[4]],pc_list[[5]][[5]],pc_list[[5]][[6]],BPPARAM = BPPARAM)$corrected
print('hh')
corr7_25 = reducedMNN(pc_list[[6]][[1]],pc_list[[6]][[2]],pc_list[[6]][[3]],BPPARAM = BPPARAM)$corrected
corr7_0 = reducedMNN(pc_list[[7]][[1]],pc_list[[7]][[2]],pc_list[[7]][[3]],pc_list[[7]][[4]],pc_list[[7]][[5]],pc_list[[7]][[6]],BPPARAM = BPPARAM)$corrected
corr6_75 = pc_list[[8]][[1]]
print('ff')

corr6_5 = reducedMNN(pc_list[[9]][[1]],pc_list[[9]][[2]],pc_list[[9]][[3]],BPPARAM = BPPARAM)$corrected


# corrct atlas
correct_atlas = reducedMNN(corr8_5,corr8_25,corr8_0,corr7_75,corr7_5,corr7_25,corr7_0,corr6_75,corr6_5 ,BPPARAM = BPPARAM)$corrected

write.csv(correct_atlas,'./corrected_pca_all_atlas_reducedMNN.csv')




