import os
import argparse
import tables as h5
import pandas as pd
import numpy as np
import gvar as gv
import matplotlib as mpl
import matplotlib.pyplot as plt
# now module for Madras-Sokal autocorr time
from emcee import autocorr

# Figure formatting for paper
fig_width = 6.75 # in inches, 2x as wide as APS column
gr        = 1.618034333 # golden ratio
fig_size  = (fig_width, fig_width / gr)
fig_size2 = (fig_width, fig_width * 1.6)
plt_axes  = [0.14,0.14,0.85,0.85]
axes1     = [0.13,0.13,0.69,0.85]
axes2     = [0.821,0.13,0.12,0.85]
fs_text   = 18 # font size of text
fs_leg    = 16 # legend font size
mrk_size  = '5' # marker size
tick_size = 16 # tick size
lw        = 1 # line width
# saving as rc params
mpl.rcParams['figure.figsize'] = fig_size
mpl.rcParams['lines.linewidth'] = lw
mpl.rcParams['lines.markersize'] = mrk_size
mpl.rcParams['font.size'] = fs_text
mpl.rcParams['legend.fontsize'] = fs_leg
mpl.rcParams['axes.labelsize'] = fs_text
mpl.rcParams['xtick.labelsize'] = tick_size
mpl.rcParams['ytick.labelsize'] = tick_size
mpl.rcParams['backend'] = 'Agg'

# streams and cfg offsets for the 3 ensembles
streams = {
    'a15m135XL_s':['b','c','d','e'],
    'a09m135_s'  :['a','b'],
    'a06m310L_s' :['b','c'],
}
## offset in units of MDTU
offset = {
    'a15m135XL_s':[0, 2700, 5400, 8100],
    'a09m135_s'  :[0, 3422],
    'a06m310L_s' :[0, 3200],
}

pbp_reps = 5  # number of stochastic point sources for PBP

labels = {'a15m135XL_s':'a15m135XL', 'a09m135_s':'a09m135', 'a06m310L_s':'a06m310L'}
lattices = {'a15m135XL':'a15', 'a09m135':'a09', 'a06m310L':'a06'}
colors = {'a15':'#ec5d57', 'a12':'#70bf41', 'a09':'#51a7f9', 'a06':'#00FFFF'}

# we have a total of 6 observables
## The first 4 are measured every MDTU
## The last 2 are measured every TRAJ
obs = ['deltaS','pbp_c','pbp_s','pbp_l','acc_rej','plaq']
# we use MDTU and TRAJ as indices for observables
index = ['mdt','traj']

######
parser = argparse.ArgumentParser(description='Extract HMC observable for an ensemble')
parser.add_argument('--ens', type=str, default='a06m135XL_s', help='pick an ensemble [%(default)s]')
parser.add_argument('--autocorr', default=True, action='store_const', const=False, help='Compute autocorr times? [%(default)s]')

if __name__ == '__main__':
    ### Setup
    args = parser.parse_args()
    ensemble = args.ens
    acorr = args.autocorr

    print('Extract HMC observable for {}'.format(ensemble))
    print('This ensemble has {} streams: {}'.format(len(streams[ensemble]),streams[ensemble]))
    prex = labels[ensemble]
    filename_obs = 'hmc_obs_'+ensemble+'.h5'
    print('Getting data from {}'.format(filename_obs))

    if ! os.path.isfile(filename_obs):
        print('File {} not found'.format(filename_obs))
        return
    ### Get data
    dataset_mdt = dict()
    dataset_trj = dict()
    dset_mdt = dict()
    dset_trj = dict()
    ## Open file
    f = h5.open_file(filename_obs)
    for s in streams[ensemble]:
        node_name = "/{}_{}/{}".format(prex,s,index[0])
        index_mdt = f.get_node(node_name).read()
        node_name = "/{}_{}/{}".format(prex,s,index[1])
        index_trj = f.get_node(node_name).read()
        for o in obs:
            node_name = "/{}_{}/{}".format(prex,s,o)
            data = f.get_node(node_name).read()
            if o in ['pbp_c','pbp_s','pbp_l']:
                data = np.reshape(data,(-1,pbp_reps*2))
                dataset_mdt[o] = data.mean(axis=1)
            if o in ['deltaS']:
                dataset_mdt[o] = data
            if o in ['acc_rej']:
                dataset_trj[o] = data[:,1]
            if o in ['plaq']:
                dataset_trj[o] = data.mean(axis=1)/3.

        dset_mdt[s] = pd.DataFrame(data=dataset_mdt,index=index_mdt)
        dset_trj[s] = pd.DataFrame(data=dataset_trj,index=index_trj)
        assert(dset_mdt[s].shape[1] == 4)  # we have 4 observables for each MDTU
        assert(dset_trj[s].shape[1] == 2)  # we have 2 observable for each TRAJ
    # Close file
    f.close()

    ###

    ### Raw analysis of autocorrelations of pbp
    for s in dset_mdt.keys():
        for pbp in ['pbp_l','pbp_s','pbp_c']:
            series = dset_mdt[s][pbp].values
            tau = autocorr.integrated_time(series,c=5)
            print('{}/{}: {}'.format(s,pbp,tau))
    ### Raw analysis of autocorrelations of plaquette
    for s in dset_trj.keys():
        series = dset_trj[s]['plaq'].values
        tau = autocorr.integrated_time(series,c=5)
        print('{}/{}: {}'.format(s,'plaq',tau))
    #### plotting
    ax = plt.gca()
    for s in dset_mdt.keys():
        dset_mdt[s].pbp_l.plot(ax=ax)
    ### Save data into txt files for unew analysis
