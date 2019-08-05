"""
Show MCMC graphs
"""

import matplotlib.pyplot as plt
import pandas as pd


def plot_chain(chain):
    # Creating a list of dataframe columns
    columns = list(chain)

    # Plotting ll
    chain.plot(kind= 'line', y = 'll')
    plt.ylabel('Log Likelihood')
    plt.xlabel('Iterate')
    plt.savefig('ll_chain.png')  # Save as .png

    # Plotting each parameter
    for i in columns[1:]:
        param = i

        # Plotting chain
        chain.plot(kind='line', y=param)
        plt.ylabel('Estimate for ' + str(param))
        plt.xlabel('Iterate')
        plt.savefig(str(param) + '_chain' + '.png')  # Save as .png

        # Plotting histogram
        chain[[param]].plot(kind='hist')
        plt.savefig(str(param) + '_freq' + '.png')  # Save as .png

    return plt.show()
