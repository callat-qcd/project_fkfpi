import os
import itertools
import pandas as pd
import numpy as np
import sys

# Get collection of fits to prune
if len(sys.argv) > 1:
    collection_name = sys.argv[1]
else:
    collection_name = input('Name for fit collection: ')

list_of_models = [] # models included in average

# Now we prune unused models
# First prune results
project_path = os.path.normpath(os.path.join(os.path.realpath(__file__), os.pardir, os.pardir))
filepath = project_path +'/results/'+ collection_name +/'results.csv'
df_fit = pd.read_csv(filepath, index_col=0, header=0)

# Gets all models not in list_of_models
prune_list = list(set(df_fit['name']) - set(list_of_models))
indices = [np.argwhere(df_fit['name'].values == model_to_prune).item() for model_to_prune in prune_list]
df_fit = df_fit.drop(indices)

# Reset indices, save results
df_fit = df_fit.reset_index(drop=True)
df_fit.to_csv(filepath)

####
# Next prune pickles
for file in os.listdir(project_path +'/results/'+ collection_name +'/pickles/'):
    if file.endswith('.p'):
        if not file.split('.')[0] in list_of_models:
            os.remove(project_path +'/results/'+ collection_name +'/pickles/'+ file)

###
# Finally prune priors
filepath = project_path +'/results/'+ collection_name +' creat/prior.csv'
df_prior = pd.read_csv(filepath, index_col=0)

#    # Gets all models not in list_of_models
prune_list = list(set(df_prior['name']) - set(list_of_models))
indices = [np.argwhere(df_prior['name'].values == model_to_prune).item() for model_to_prune in prune_list]
df_prior = df_prior.drop(indices)

df_prior = df_prior.reset_index(drop=True)
df_prior.to_csv(filepath)
