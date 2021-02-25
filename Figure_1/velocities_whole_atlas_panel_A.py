######################### This code plots Figure 1 A and B
#########################
import numpy as np
import scvelo as scv
import scanpy as sc
import matplotlib.pyplot as plt

celltype_colours = {
  "Allantois" : "#532C8A",
  "Anterior Primitive Streak" : "#c19f70",
  "Blood progenitors 1" : "#f9decf",
  "Blood progenitors 2" : "#c9a997",
  "Cardiomyocytes" :  "#B51D8D",#
  "Caudal epiblast" : "#9e6762",
  "Caudal Mesoderm" : "#3F84AA",
  "Def. endoderm" : "#F397C0",#
  "Nascent mesoderm" :  "#C594BF",#
  "Mixed mesoderm" :  "#DFCDE4",#
  
  "Endothelium" :  "#eda450",#
  "Epiblast" :  "#635547",#
  "Erythroid1" :  "#C72228",#
  "Erythroid2" :  "#EF4E22",#
  "Erythroid3" : "#f77b59",
  "ExE ectoderm" :  "#989898",#
  
  "ExE endoderm" : "#7F6874",#
  "ExE mesoderm" :  "#8870ad",#
  
  
  "Rostral neurectoderm" :  "#65A83E",#
  "Forebrain/Midbrain/Hindbrain" : "#647a4f",
  
  
  "Gut" :  "#EF5A9D",
  "Haematoendothelial progenitors" :  "#FBBE92",#
  "Caudal neurectoderm": "#354E23",
  
  
  "Intermediate mesoderm" :  "#139992",#
  "Neural crest": "#C3C388",
  
  "NMP" :  "#8EC792",#
  "Notochord" :  "#0F4A9C",#
  "Paraxial mesoderm" :  "#8DB5CE",#
  "Parietal endoderm" :  "#1A1A1A",
  "PGC" :  "#FACB12",
  
  "Pharyngeal mesoderm" :  "#C9EBFB",#
  "Primitive Streak" :  "#DABE99",#
  "Mesenchyme" : "#ed8f84",
  "Somitic mesoderm" :  "#005579",#
  "Spinal cord" :  "#CDE088",#
  "Surface ectoderm" : "#BBDCA8",#
  
  
  "Visceral endoderm" : "#F6BFCB",#
  "Mes1": "#c4a6b2",#
  "Mes2":"#ca728c",#
  
  "Cardiomyocytes" :  "#B51D8D",
}
  
stage_colours = { "E6.5" : "#D53E4F",
                  "E6.75" : "#F46D43",
                  "E7.0" : "#FDAE61",
                  "E7.25" : "#FEE08B",
                  "E7.5" : "#FFFFBF",
                  "E7.75" : "#E6F598",
                  "E8.0" : "#ABDDA4",
                  "E8.25" : "#66C2A5",
                  "E8.5" : "#3288BD",
}

adata = sc.read('adata_umap_pca.h5',cache=True)

adata_proc = adata.copy()
scv.pp.filter_and_normalize(adata_proc, n_top_genes=3000, min_shared_counts=20)

scv.pp.moments(adata_proc, n_pcs=50, n_neighbors=30)

scv.tl.recover_dynamics(adata_proc)

scv.tl.velocity(adata_proc, mode="dynamical")

scv.tl.velocity_graph(adata_proc)
#scv.tl.velocity_graph(adata_proc, approx=True) 

#adata_proc.write('adata_dynamical_top3000.h5')

## Plots for manuscript
#adata_proc = sc.read('adata_dynamical_top3000.h5',cache=True)

colPalette_celltype = [celltype_colours[i] for i in sorted(np.unique(adata_proc.obs['celltype']))]
adata_proc.uns['celltype_colors'] = colPalette_celltype

colPalette_stage = [stage_colours[i] for i in sorted(np.unique(adata_proc.obs['stage']))]
adata_proc.uns['stage_colors'] = colPalette_stage

fig, axs = plt.subplots(1,1, figsize = (20,20))

scv.pl.velocity_embedding_stream(adata_proc, basis="umap", color="celltype", size=50, ax = axs)

fig, axs = plt.subplots(1,1, figsize = (20,20))

scv.pl.velocity_embedding_stream(adata_proc, basis="umap", color="stage", size=50, ax = axs)

scv.tl.velocity_pseudotime(adata_proc)

scv.pl.scatter(adata_proc, color='velocity_pseudotime', cmap='gnuplot')
