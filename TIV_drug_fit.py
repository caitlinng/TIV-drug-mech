'''
Generate simulation data
'''


# Module import
import TIV_drug_funcs as TIV
import test_fit as test
import pandas as pd
import numpy as np
import scipy as sci
from scipy import stats as st
import matplotlib.pyplot as plt
import datetime
import math
import MCMC_chain as MCMC
import TIV_plot

import TIV_drug_sim_data as sim_data

import random
random.seed(42)

import time
start_time = time.time()

'''
MCMC Set-up
'''

# Set random seed for reproducibility
#np.random.seed(42)

# Save chain ('y') or not (anything else)
save = 'y'

# Generate simulation data
sim_V_data = sim_data.fitted_data.y[2]
drug = 'plc'

# Import viral load data (will try fitting first with plc_n1)
import TIV_Rubi_data
real_V_data = TIV_Rubi_data.plc_n1_raw
drug = 'plc'

# Number of iterates
n_iterates = 1000

# Prior functions (what is our prior belief about beta, gamma)

def prior_param_belief_uni(min, max):
    return st.uniform(loc=min, scale=max-min)

def prior_param_belief_normal(mean, var):
    return st.norm(loc=mean, scale=var)


# Setting parameters
param_names = ['beta_dot', 'beta', 'deltaI', 'pV', 'deltaV', 'V0']

# Input own parameters (e.g. from previously run chain)
init_param = sim_data.final_fit_param

# Proposal widths (1 for each parameter)
w = [1e-3, 1e-3, 1.5, 1, 1, 2]

# Setting prior belief as that parameters are within a reasonable given range
beta_dot_prior_fun = prior_param_belief_uni(1e-12, 3)
beta_prior_fun = prior_param_belief_uni(1e-12, 3) #(0, 0.5) # (e-12, e-4)  #(0, 0.5)
deltaI_prior_fun = prior_param_belief_uni(0, 1e+6)
pV_prior_fun = prior_param_belief_uni(0, 1e+6) #(0, 24) #(-6, 6) #(150, 250)
deltaV_prior_fun = prior_param_belief_uni(0, 1e+6)
V0_prior_fun = prior_param_belief_uni(0, 1e+4)

prior_funcs = [beta_dot_prior_fun, beta_prior_fun, deltaI_prior_fun, pV_prior_fun, deltaV_prior_fun, V0_prior_fun]

# Choose to fit TIV or TLIV model
model_ll = TIV.TIV_drug_ll
model = TIV.TIV_drug_model

# Start simulation
V_data = sim_V_data
chain_vals = MCMC.run_chain(model_ll, param_names, init_param, V_data, drug, len(V_data), n_iterates, w, prior_funcs, save='y', calibrate_truth=False)

chain = chain_vals[0]
best_MCMC_param = chain_vals[1]
last_MCMC_param = chain_vals[2]

# Produce graph showing how well simulation-generated parameters fits true data
test.test_fit(model=model, drug='plc', true_data=V_data, fitted_param=best_MCMC_param, max_time=len(V_data))

# Plot MCMC chain (log likelihood, parameter chains and histograms)
TIV_plot.plot_chain(chain)

print("--- %s seconds ---" % (time.time() - start_time))