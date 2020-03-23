import itertools
import sys

sys.path.append("../")
import fitter.data_loader as dl

### Edit these lines
collection_name = 'test'
choices = {
    'fit_type' : ['ma', 'ma-ratio', 'xpt', 'xpt-ratio'],
    'order' : ['nnlo'],
    'F2' : ['FKFK', 'FKFpi', 'FpiFpi'],

    # Vol corrections
    'include_FV' : [True],

    # semi-nnlo corrections
    'include_alpha_s' : [False, True],
    'semi-nnlo_corrections' : ['nnlo-ct', 'nnlo-logSq_ct', 'nnlo-full', 'nnlo-full_bijnens'],

    # nnnlo corrections
    'include_latt_n3lo' : [False, True],
}

### Don't edit
model_list = []
for j, choice in enumerate(dict(zip(choices, x)) for x in itertools.product(*choices.values())):

    # Get model info from choices
    model_info = {}
    model_info['fit_type'] = choice['fit_type']
    model_info['order'] = choice['order']
    model_info['F2'] = choice['F2']

    model_info['include_FV'] = choice['include_FV']
    model_info['include_alpha_s'] = choice['include_alpha_s']
    model_info['include_latt_n3lo'] = choice['include_latt_n3lo']

    if choice['semi-nnlo_corrections'] == 'nnlo-full_bijnens':
        model_info['include_log'] = True
        model_info['include_log2'] = True
        model_info['include_sunset'] = True
        model_info['use_bijnens_central_value'] = True

    elif choice['semi-nnlo_corrections'] == 'nnlo-full':
        model_info['include_log'] = True
        model_info['include_log2'] = True
        model_info['include_sunset'] = True
        model_info['use_bijnens_central_value'] = False

    elif choice['semi-nnlo_corrections'] == 'nnlo-logSq_ct':
        model_info['include_log'] = False
        model_info['include_log2'] = True
        model_info['include_sunset'] = True
        model_info['use_bijnens_central_value'] = False

    elif choice['semi-nnlo_corrections'] == 'nnlo-ct':
        model_info['include_log'] = False
        model_info['include_log2'] = False
        model_info['include_sunset'] = False
        model_info['use_bijnens_central_value'] = False


    # Get model name from info
    name = model_info['fit_type'] +'_'+ model_info['F2'] +'_'+ model_info['order']
    if model_info['include_log']:
        name = name + '_log'
    if model_info['include_log2']:
        name = name + '_logSq'
    if model_info['include_sunset']:
        name = name + '_sunset'
    if model_info['include_alpha_s']:
        name = name + '_alphaS'
    if model_info['include_latt_n3lo']:
        name = name + '_a4'
    if model_info['include_FV']:
        name = name + '_FV'
    if model_info['use_bijnens_central_value']:
        name = name + '_bijnens'

    model_list.append(name)

data_loader = dl.data_loader()
data_loader.save_model_names(collection_name, model_list)
