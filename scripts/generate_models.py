import itertools
import sys

sys.path.append("../../")
import fitter.data_loader as dl

data_loader = dl.data_loader()

### Edit these lines
choices = []

choices.append({
    'fit_type' : ['ma', 'ma-ratio', 'xpt', 'xpt-ratio'],
    'order' : ['n2lo'],
    'latt_ct' : ['n2lo', 'n3lo'],
    'F2' : ['FKFK', 'FKFpi', 'FpiFpi'],

    # Vol corrections
    'include_FV' : [True],

    # semi-n2lo corrections
    'include_alpha_s' : [False, True],
    'semi-n2lo_corrections' : ['n2lo-ct', 'n2lo-full_bijnens'],
})

choices.append({
    'fit_type' : ['ma', 'ma-ratio', 'xpt', 'xpt-ratio'],
    'order' : ['n3lo'],
    'latt_ct' : ['n3lo'],
    'F2' : ['FKFK', 'FKFpi', 'FpiFpi'],

    # Vol corrections
    'include_FV' : [True],

    # semi-n2lo corrections
    'include_alpha_s' : [False, True],
    'semi-n2lo_corrections' : ['n2lo-full_bijnens'],
})

choices.append({
    'fit_type' : ['ma', 'ma-ratio', 'xpt', 'xpt-ratio'],
    'order' : ['n3lo'],
    'latt_ct' : ['n4lo'],
    'F2' : ['FKFK', 'FKFpi', 'FpiFpi'],

    # Vol corrections
    'include_FV' : [True],

    # semi-n2lo corrections
    'include_alpha_s' : [False],
    'semi-n2lo_corrections' : ['n2lo-full_bijnens'],
})

### Don't edit

if len(sys.argv) > 1:
    collection_name = sys.argv[1]
else:
    collection_name = input('Name for fit collection: ')

model_list = []
for k in range(len(choices)):
    for j, choice in enumerate(dict(zip(choices[k], x)) for x in itertools.product(*choices[k].values())):

        # Get model info from choices
        model_info = {}
        model_info['fit_type'] = choice['fit_type']
        model_info['order'] = choice['order']
        model_info['F2'] = choice['F2']

        model_info['include_FV'] = choice['include_FV']
        model_info['include_alpha_s'] = choice['include_alpha_s']
        model_info['latt_ct'] = choice['latt_ct']

        if choice['semi-n2lo_corrections'] == 'n2lo-full_bijnens':
            model_info['include_log'] = True
            model_info['include_log2'] = True
            model_info['include_sunset'] = True
            model_info['use_bijnens_central_value'] = True

        elif choice['semi-n2lo_corrections'] == 'n2lo-full':
            model_info['include_log'] = True
            model_info['include_log2'] = True
            model_info['include_sunset'] = True
            model_info['use_bijnens_central_value'] = False

        elif choice['semi-n2lo_corrections'] == 'n2lo-logSq_ct':
            model_info['include_log'] = False
            model_info['include_log2'] = True
            model_info['include_sunset'] = True
            model_info['use_bijnens_central_value'] = False

        elif choice['semi-n2lo_corrections'] == 'n2lo-ct':
            model_info['include_log'] = False
            model_info['include_log2'] = False
            model_info['include_sunset'] = False
            model_info['use_bijnens_central_value'] = False


        # Get model name from info
        name = data_loader.get_model_name_from_model_info(model_info)

        model_list.append(name)


data_loader.save_model_names(collection_name, model_list)
