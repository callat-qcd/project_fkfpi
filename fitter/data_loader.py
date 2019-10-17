import pandas as pd
import numpy as np
import gvar as gv
import sys
import time
import os
import h5py
from matplotlib.backends.backend_pdf import PdfPages
from collections import OrderedDict

class data_loader(object):

    def __init__(self):
        self.project_path = os.path.normpath(os.path.join(os.path.realpath(__file__), os.pardir, os.pardir))
        self.file_h5 = os.path.normpath(self.project_path+'/FK_Fpi_data.h5')

    def get_ensembles(self):
        with h5py.File(self.file_h5, "r") as f:
            ensembles = f.keys()
        return sorted(ensembles)

    def get_fit_data(self):
        data = {}
        with h5py.File(self.file_h5, "r") as f:
            for ensemble in f.keys():
                data[ensemble] = {}
                for key in f[ensemble].keys():
                    dset = f[ensemble][key]
                    shape = dset.shape
                    if shape != ():
                        data[ensemble][key] = dset[:]
                    else:
                        data[ensemble][key] = dset[()]
        return data
    # whose: 'mine', 'others'
    def get_fit_info(self, whose=None):
        filepath = None
        if whose=='mine' or whose is None:
            filepath = os.path.normpath(self.project_path + '/results/fit_results.csv')
        elif whose=='others':
            filepath = os.path.normpath(self.project_path + '/results/other_results.csv')

        if not os.path.isfile(filepath):
            return None


        df_fit = pd.read_csv(filepath, header=0)
        cols = df_fit.columns.values
        #cols = ['fit', 'logGBF', 'chi2/df', 'Q', 'vol corr', 'latt corr']
        #cols = np.intersect1d(cols, df_fit.columns.values)

        fit_types = df_fit['name'].values
        output_dict = OrderedDict()
        for name in sorted(fit_types):
            index = np.argwhere(df_fit['name'].values == name)

            output_dict[name]= {key: np.asscalar(df_fit[key].values[index]) for key in cols}

        return output_dict

    def get_prior(self, fit_type, F2):
        filepath = self.project_path + '/priors/'+fit_type+'_'+F2+'.csv'
        if not os.path.isfile(filepath):
            return None
        df_read = pd.read_csv(filepath, index_col=0)
        return gv.gvar({key : df_read.to_dict("index")[key]['0']
                for key in df_read.to_dict("index").keys()})

    def get_variable_names(self):
        names = []
        with h5py.File(self.file_h5, "r") as f:
            for ensemble in f.keys():
                for key in f[ensemble].keys():
                    names.append(key)
        return sorted(np.unique([names]))

    def save_fit_info(self, fit_info):
        print "Saving..."

        if not os.path.exists(self.project_path + '/results/'):
            os.makedirs(self.project_path + '/results/')
        filepath = os.path.normpath(self.project_path + '/results/fit_results.csv')

        # get fit info
        cols = np.array(['name', 'fit', 'delta_su2', 'logGBF', 'chi2/df', 'Q', 'vol corr', 'latt corr'])
        lecs_cols = ['L_4', 'L_5', # nlo terms
                     'A_a', 'A_k', 'A_p', 'A_loga'] # nnlo terms
                     #'A_aa', 'A_ak', 'A_ap', # nnnlo terms
                     #'A_kk', 'A_kp', 'A_pp'] # more nnnlo terms

        # append LEC results
        cols = np.concatenate((cols,lecs_cols), axis=0)
        #diff = sorted(list(set(cols).symmetric_difference(fit_info.keys()))) # Gets set difference (ie, LEC keys)
        #cols = np.concatenate((cols, diff), axis=0) # Append LEC keys to columns
        #print "\n\n---"
        #print cols
        #print"---"

        # fit_info keys not in cols -> create key in fit_info
        for key in cols:
            if key not in fit_info.keys():
                #print "fit info key not in col: ", key
                fit_info[key] = np.nan

        # Add result to file if file exists
        if os.path.isfile(filepath):
            df_best_fits = pd.read_csv(filepath, index_col=0).to_dict()
            #print df_best_fits.keys()
            #print sorted(list(set(cols).symmetric_difference(df_best_fits.keys())))


            # Get keys in df but not in cols
            #diff = sorted(list(set(cols).symmetric_difference(df_best_fits.keys())))
            #cols = np.concatenate((cols, diff), axis=0)
            #print "\n\n---"
            #print cols
            #print"---"

            output_dict = {}
            for key in df_best_fits.keys():
                output_dict[key] = np.array(df_best_fits[key].values(), dtype="object")


            # df keys not in cols -> create keys in df
            for key in cols:
                if key not in output_dict.keys():
                    print "df key not in col: ", key
                    output_dict[key] = np.repeat(np.nan, len(output_dict['name']))



            if fit_info['name'] in output_dict['name']:
                index = np.asscalar(np.argwhere(output_dict['name'] == fit_info['name']))
                for key in fit_info.keys():
                    output_dict[key][index] = fit_info[key]
            else:
                for key in output_dict.keys():
                    output_dict[key] = np.append(output_dict[key], fit_info[key])

            df = pd.DataFrame.from_dict(output_dict)
            df = df[cols]
            df.sort_values('name')
            df.to_csv(filepath)

            #print
            #print df

        # Create new file if it doesn't exist
        else:
            output_dict = {key : [fit_info[key]] for key in fit_info.keys()}
            df = pd.DataFrame.from_dict(output_dict)
            df = df[cols] # rearrange in logical order
            df.to_csv(filepath)

        print "Done."
        return None

    def save_plots(self, figs=None, output_filename=None):

        if not os.path.exists(os.path.normpath(self.project_path+'/tmp/')):
            os.makedirs(os.path.normpath(self.project_path+'/tmp/'))

        if output_filename is None:
            output_file = os.path.normpath(self.project_path+'/tmp/temp.pdf')
        else:
            output_file = os.path.normpath(self.project_path+'/tmp/'+output_filename+'.pdf')

        if not os.path.exists(os.path.dirname(output_file)):
            os.makedirs(os.path.dirname(output_file))

        output_pdf = PdfPages(output_file)
        try:
            for fig in figs:
                output_pdf.savefig(fig)
        except TypeError: # save figs directly if figs is just a figure, not an array
            output_pdf.savefig(figs)

        output_pdf.close()
        print "Done."

        return None

    def save_prior(self, prior, fit_type, F2):
        filepath = self.project_path + '/priors/'+fit_type+'_'+F2+'.csv'

        current_prior = self.get_prior(fit_type, F2)

        #for key in out_prior.keys():
        #    out_prior[key] = [str(prior[out_prior])]

        out_prior = {}
        if current_prior is not None:
            for key in current_prior.keys():
                out_prior[key] = [str(current_prior[key])]

        for key in prior.keys():
            out_prior[key] = [str(prior[key])]

        df = pd.DataFrame.from_dict(out_prior).T
        df.to_csv(filepath)

        return None
