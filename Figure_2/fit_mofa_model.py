# This code takes as input the prepared MOFA objects for spliced and unspliced and outputs the fitted model

from mofapy2.run.entry_point import entry_point
ent = entry_point()
import scanpy as sc

adata_u = sc.read('./input_mofa_uc_all.h5') % unspliced object
adata_s = sc.read('./input_mofa_sc_all.h5') % spliced object

G = 1
M = 2

data_mat = [[None for g in range(G)] for m in range(M)]

data_mat[0][0] = adata_s.X.todense()
data_mat[1][0] = adata_u.X.todense()


ent.set_data_options(
    scale_groups = False,
    scale_views = True
)

ent.set_data_matrix(data_mat, views_names = ['spliced','unspliced'])

ent.set_model_options(
    factors = 10,
    spikeslab_weights = True,
    ard_factors = True,
    ard_weights = True
)

ent.set_train_options(
    iter = 1000, 
    convergence_mode = "fast", 
    startELBO = 1, 
    freqELBO = 1, 
    dropR2 = 0.001, 
    gpu_mode = True, 
    verbose = False, 
    seed = 1
)

ent.build()

ent.run()

# Save the output
outfile = "./mofa_out_cos_all.hdf5"
ent.save(outfile)

