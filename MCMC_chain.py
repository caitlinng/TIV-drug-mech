"""
MCMC chain
"""

import scipy as sci
import numpy as np
import datetime
import math
import pandas as pd

'''
tune_w()
Checks acceptance rate for this parameter (need to input param-specific values and modifies w accordingly
    If too low (<0.3) multiply parameter's w by 0.75
    If too high (>0.5), multiply w by 1.5
Returns w
'''


def tune_w(accepts, current_iter, w, param_name):
    print('w was ' + str(w))

    acceptance_rate = accepts / current_iter

    if acceptance_rate < 0.3:
        w *= 0.75

    if acceptance_rate > 0.5:
        w *= 1.5

    print('Tuned w is now ' + str(w))
    print('Acceptance rate for ' + param_name + ' is ' + str(acceptance_rate))

    return w


'''
Runs MCMC chain
'''


def run_chain(model_ll, param_names, init_param, V_data, drug, max_time, n_iterates, w, prior_funcs, save, calibrate_truth):  # NOTE: Need to input TIV_funcs.model_ll
    now = datetime.datetime.now()  # Record time chain started
    calibrate = calibrate_truth  # Calibrate or not

    # Calculate the log likelihood of the initial guess
    init_ll = model_ll(V_data, drug, init_param, max_time)

    # And log likelihood of all subsequent guesses
    param = init_param.copy()
    ll = init_ll.copy()

    # Establish data storage space for chain
    # Where first column is ll and 1: are the n-parameters [ll, beta, gamma]
    chain = np.zeros((n_iterates, len(param) + 1))
    chain[0, 0] = ll
    chain[0, 1:] = param

    accepts = np.zeros(len(init_param))

    # Run MCMC
    for i in range(n_iterates):

        # Print status every 100 iterations
        if i % 100 == 0:
            print('Iteration ' + str(i) + ' of ' + str(n_iterates))

        # Gibbs loop over number of parameters (j = 0 is beta, j = 1 is gamma)
        for j in range(len(param)):

            # For every 100 iterations, calibrate w
            if calibrate is True:
                if i % 100 == 0 and i > 0:
                    w[j] = tune_w(accepts=accepts[j], current_iter=i, w=w[j], param_name=param_names[j])
                    print('w = ' + str(w))
                    accepts[j] = 0  # Reset acceptance rate

            # Propose a parameter value within prev. set widths
            prop_param = param.copy()

            # Take a random step size from a uniform distribution (that has width of w)
            prop_param[j] = prop_param[j] - (sci.stats.uniform.rvs(loc=-w[j] / 2, scale=w[j], size=1))
            prop_param[j] = np.ndarray.item(prop_param[j])  # Converts param value from single element array into a scalar

            # Deal with invalid proposals by leaving ll, param unchanged
            if prop_param[j] <= 0:  # Invalid, so try next parameter proposal
                prop_ll = -1 * np.inf

            else:
                # Calculate LL of proposal
                prop_ll = model_ll(V_data, drug, prop_param, max_time)  # TODO: Testing with TLIV fitting

            # Decide on accept/reject
            prior_fun = prior_funcs[j]  # Grab correct prior function

            # Likelihood ratio st.norm.rvs(1, 1)
            r = np.exp(prop_ll - ll) * prior_fun.pdf(prop_param[j]) / prior_fun.pdf(param[j])

            if math.isnan(r):
                alpha = 0
                print('r was a NaN')

            else:
                # Is likelihood ratio less than or equal to one
                alpha = min(1, r)

                print('r = ' + str(r))

            # Random number between 0 to 1
            # So will have weighted chance of possibly accepting depending on how likely the new parameter is
            test = np.random.uniform(0, 1)
            # Maybe accept
            if test < alpha:
                ll = prop_ll.copy()
                param = prop_param.copy()
                accepts[j] += 1

                print('Accept new parameters ' + str(param))


            # "Else" reject, though nothing to write
            else:
                print('Reject parameters ' + str(prop_param))

            # Store iterate
            chain[i, 0] = ll
            chain[i, 1:] = param

            # Update chain.txt file with new iterate as chain is generated (discards first half though)
            '''
            if save == 'y' and i > int(n_iterates/2):
                f = open('chain.txt', 'a')
                f.write('\n' + str(chain[i]))
                f.close()
            '''

    if calibrate is True:
        print('w = ' + str(w))
        for i in range(len(accepts)):
            acceptance_rate = accepts[i]/ n_iterates
            print('Acceptance rate for ' + param_names[i] is str(acceptance_rate))

    # Modifying data to make
    # easier to retrieve values
    chain_half = pd.DataFrame(chain[int(n_iterates / 2):, :], columns=['ll'] + param_names)
    MCMC_param = chain[-1, 1:]  # TODO: Look at last estimated ll

    best_ll_index = chain_half[['ll']].idxmax()
    best_ll_row = chain_half.iloc[best_ll_index, :]

    if save == 'y':
        f = open('chain_info.txt', 'a')
        f.write('\n')
        f.write('\n' + now.strftime("%Y-%m-%d %H:%M"))
        f.write('\n' + 'init_param = ' + str(init_param))
        f.write('\n' + 'w = ' + str(w))
        f.write('\n' + 'n_iterates = ' + str(n_iterates))
        f.write('\n' + str(chain_half))
        f.write('\n' + 'Last tested parameters were: ' + str(MCMC_param))
        f.write('\n' + 'best_ll row = ' + str(best_ll_row))
        f.close()

    return [chain_half, best_ll_row[1:], MCMC_param]
