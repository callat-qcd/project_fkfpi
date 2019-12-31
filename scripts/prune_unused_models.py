import os
import itertools
import pandas as pd
import numpy as np

# Get all proper model names
p_dict = {
    'order' : {
        'fit' : 'nnlo', # 'nlo', 'nnlo', or 'nnnlo'
    },
}

choices = {
    'fit' : ['nnlo'],
    'fit_type' : ['ma', 'ma-ratio', 'xpt', 'xpt-ratio'],
    'F2' : ['FKFK', 'FKFpi', 'FpiFpi'],
    'vol' : [0, 10],

    # semi-nnlo corrections
    'include_alpha_s' : [False, True],
    'semi-nnlo_corrections' : ['nnlo-full_bijnens', 'nnlo-full', 'nnlo-logSq_ct', 'nnlo-ct'],

    # nnnlo corrections
    'include_latt_n3lo' : [False, True],
    'output_name' : 'optimized_priors',
}

list_of_models = [] # models included in average
# Get all enumerations of these choices
for j, choice in enumerate(dict(zip(choices, x)) for x in itertools.product(*choices.values())):

    p_dict = {}
    p_dict['order'] = {}

    p_dict['order']['fit'] = choice['fit']
    p_dict['fit_type'] = choice['fit_type']
    p_dict['F2'] = choice['F2']

    p_dict['order']['vol'] = choice['vol']
    p_dict['order']['include_alpha_s'] = choice['include_alpha_s']
    p_dict['order']['include_latt_n3lo'] = choice['include_latt_n3lo']

    if choice['semi-nnlo_corrections'] == 'nnlo-full_bijnens':
        p_dict['order']['include_log'] = True
        p_dict['order']['include_log2'] = True
        p_dict['order']['include_sunset'] = True
        p_dict['use_bijnens_central_value'] = True

    elif choice['semi-nnlo_corrections'] == 'nnlo-full':
        p_dict['order']['include_log'] = True
        p_dict['order']['include_log2'] = True
        p_dict['order']['include_sunset'] = True
        p_dict['use_bijnens_central_value'] = False

    elif choice['semi-nnlo_corrections'] == 'nnlo-logSq_ct':
        p_dict['order']['include_log'] = False
        p_dict['order']['include_log2'] = True
        p_dict['order']['include_sunset'] = True
        p_dict['use_bijnens_central_value'] = False

    elif choice['semi-nnlo_corrections'] == 'nnlo-ct':
        p_dict['order']['include_log'] = False
        p_dict['order']['include_log2'] = False
        p_dict['order']['include_sunset'] = False
        p_dict['use_bijnens_central_value'] = False


    name = p_dict['fit_type'] +'_'+ p_dict['F2']+'_'+p_dict['order']['fit']
    if p_dict['order']['include_log']:
        name = name + '_log'
    if p_dict['order']['include_log2']:
        name = name + '_logSq'
    if p_dict['order']['include_sunset']:
        name = name + '_sunset'
    if p_dict['order']['include_alpha_s']:
        name = name + '_alphaS'
    if p_dict['order']['include_latt_n3lo']:
        name = name + '_a4'
    if p_dict['order']['vol'] > 6:
        name = name + '_FV'
    if p_dict['use_bijnens_central_value']:
        name = name + '_bijnens'

    list_of_models.append(name)


# Now we prune unused models
# First prune ./results/ folder
project_path = os.path.normpath(os.path.join(os.path.realpath(__file__), os.pardir, os.pardir))
filepath = project_path +'/results/'+ output_name +'.csv'
df_fit = pd.read_csv(filepath, index_col=0, header=0)

# Gets all models not in list_of_models
prune_list = list(set(df_fit['name']) - set(list_of_models))
indices = [np.argwhere(df_fit['name'].values == model_to_prune).item() for model_to_prune in prune_list]
df_fit = df_fit.drop(indices)

# Reset indices, save results
df_fit = df_fit.reset_index(drop=True)
df_fit.to_csv(filepath)

####
# Next prune ./pickles/ folder
for file in os.listdir(project_path +'/pickles/'+ output_name):
    if file.endswith('.p'):
        if not file.split('.')[0] in list_of_models:
            os.remove(project_path+'/pickles/'+ output_name +'/'+ file)

###
# Finally prune ./priors/ folder
for fit_type in ['ma', 'ma-ratio', 'xpt', 'xpt-ratio']:
    filepath = project_path +'/priors/'+ fit_type +'.csv'
    df_prior = pd.read_csv(filepath, index_col=0)

    # Gets all models not in list_of_models
    prune_list = list(set(df_prior['name']) - set(list_of_models))
    indices = [np.argwhere(df_prior['name'].values == model_to_prune).item() for model_to_prune in prune_list]
    df_prior = df_prior.drop(indices)

    df_prior = df_prior.reset_index(drop=True)
    df_prior.to_csv(filepath)
