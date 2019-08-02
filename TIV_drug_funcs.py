'''
Returns viral titre measurements (8 days) of solved TIV model variations
'''

from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

def TIV_drug_model(drug, param, max_time):  # where drug = [type, epsilon], param = [

    # Write differential equations describing TIV drug model (right-hand side)
    def rhs(t, y):
        T, I, V = y
        return [g * T * (1 - (T + I) / T0) - (beta * V * T),
                beta * V * T - (deltaI * I),
                pV * I - (deltaV * V) - (beta * V * T)]

    # Input parameters
    g = param[0]
    beta = param[1]
    deltaI = param[2]
    pV = param[3]
    deltaV = param[4]

    # If drug present, adjust respective parameter that drug acts on by multiplying with drug scaling factor
    if drug[0] == 'rep':
        pV = pV *  (1 - drug[1])

    if drug[0] == 'ent':
        beta = beta * (1 - drug[1])

    # Else leave parameters unchanged (e.g. placebo)

    # Initial conditions
    T0 = 4e+8  # Fixing this parameter
    I0 = 0
    V0 = param[5]

    y_init = [T0, I0, V0]

    # Create an array of measurement time points (Day 1, 2, 3, ... 8)
    measurement_times = np.arange(0, max_time, 1)

    # Solve TIV
    sol = solve_ivp(rhs, t_span=(0, max_time), y0=y_init, method='RK45', t_eval=measurement_times)
    return sol  # NOTE: Returns raw (not log) values

def PK_model(drug, param, max_time):
    rhs = epsilon_max * D
    return drug  # [type, epsilon]
# Generate placebo data

# Parameters expressed as log10 (e.g. deltaV = 1e+5)
pV = 12.6  # This is technically production rate of infectious virions though not all V are infectious
beta = 5e-7
V0 = 1e+4
T0 = 7e+7
g = 0.8
deltaV = 4
deltaI = 2

param = [g, beta, deltaI, pV, deltaV, V0]

plc_data = TIV_drug_model(drug='plc', param=param, max_time=8)

# Generate entry-blocking drug data
#ent_drug = [ent, #epsilon]

# Generate replication-blocking drug data
#rep_drug = [rep, #epsilon]

