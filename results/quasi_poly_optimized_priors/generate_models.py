import itertools
import sys

sys.path.append("../../")
import fitter.data_loader as dl
data_loader = dl.data_loader()

### Edit these lines
choices = {
    'fit_type' : ['poly'],
    'order' : ['nlo', 'n2lo', 'n3lo'],
    'latt_ct' : ['nlo', 'n2lo', 'n3lo', 'n4lo'],
    'F2' : ['FKFK', 'FKFpi', 'FpiFpi'],

    # Vol corrections
    'include_FV' : [False],

    # semi-nnlo corrections
    'include_alpha_s' : [False],
    'semi-nnlo_corrections' : ['nnlo-ct'],
}

### Don't edit

if len(sys.argv) > 1:
    collection_name = sys.argv[1]
else:
    collection_name = input('Name for fit collection: ')

model_list = []
for j, choice in enumerate(dict(zip(choices, x)) for x in itertools.product(*choices.values())):

    # Get model info from choices
    model_info = {}
    model_info['fit_type'] = choice['fit_type']
    model_info['order'] = choice['order']
    model_info['latt_ct'] = choice['latt_ct']
    model_info['F2'] = choice['F2']

    model_info['include_FV'] = choice['include_FV']
    model_info['include_alpha_s'] = choice['include_alpha_s']

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
    name = data_loader.get_model_name_from_model_info(model_info)

    model_list.append(name)

data_loader.save_model_names(collection_name, model_list)