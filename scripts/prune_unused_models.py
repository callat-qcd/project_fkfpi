import os
import itertools
import pandas as pd

# Get all proper model names
p_dict = {
    'order' : {
        'fit' : 'nnlo', # 'nlo', 'nnlo', or 'nnnlo'
    },
}

choices = {
    'fit_type' : ['ma', 'ma-ratio', 'xpt', 'xpt-ratio'],
    'F2' : ['FKFK', 'FKFpi', 'FpiFpi'],
    'vol' : [0, 10],

    # semi-nnlo corrections
    'include_alpha_s' : [False, True],
    'semi-nnlo_corrections' : ['nnlo-full_bijnens', 'nnlo-full', 'nnlo-logSq_ct', 'nnlo-ct'],

    # nnnlo corrections
    'include_latt_n3lo' : [False, True],
}

list_of_models = []
# Get all enumerations of these choices
for j, choice in enumerate(dict(zip(choices, x)) for x in itertools.product(*choices.values())):

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

# Get all model names in csv file
#project_path = os.path.normpath(os.path.join(os.path.realpath(__file__), os.pardir, os.pardir))
project_path = './../'
filepath = project_path +'/results/fit_results.csv'
df_fit = pd.read_csv(filepath, header=0)

# Gets all models not in list_of_models
list(set(df_fit['name']) - set(list_of_models))
