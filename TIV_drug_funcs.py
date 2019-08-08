'''
Returns viral titre measurements (8 days) of solved TIV model variations
'''

from scipy.integrate import solve_ivp
import numpy as np
from scipy import stats as st

def TIV_drug_model(Drug, param, max_time):  # where drug = [type, epsilon], param = [

    # Write differential equations describing TIV drug model (right-hand side)
    def rhs(t, y):
        T, I, V = y
        return [g * T * (1 - (T + I) / T0) - (beta_dot * V * T),
                beta_dot * V * T - (deltaI * I),
                pV * I - (deltaV * V) - (beta * V * T)]

    # Input parameters
    g = 0.8
    beta_dot = param[0]
    beta = param[1]
    deltaI = param[2]
    pV = param[3]
    deltaV = param[4]

    # If drug present, adjust respective parameter that drug acts on by multiplying with drug scaling factor
    if Drug.type == 'rep':
        pV = pV * (1 - Drug.epsilon)

    if Drug.type == 'ent':
        beta_dot = beta_dot * (1 - Drug.epsilon)

    # Initial conditions
    T0 = 4e+8  # Fixing this parameter 7e+7
    I0 = 0
    V0 = param[5]

    y_init = [T0, I0, V0]

    # Create an array of measurement time points (Day 1, 2, 3, ... 8)
    measurement_times = np.arange(0, max_time, 1)

    # Solve TIV
    sol = solve_ivp(rhs, t_span=(0, max_time), y0=y_init, method='BDF', t_eval=measurement_times)
    return sol  # NOTE: Returns raw (not log) values


def TIV_Drug_model(Drug, param, max_time):  # where drug = [type, epsilon], param = [

    # Write differential equations describing TIV drug model (right-hand side)
    def rhs(t, y):
        T, I, V, D = y
        return [g * T * (1 - (T + I) / T0) - (beta_dot * V * T),
                beta_dot * V * T - (deltaI * I),
                pV * (1 - epsilon) * I - (deltaV * V) - (beta * V * T)]

    # Input parameters
    g = 0.8
    beta_dot = param[0]
    beta = param[1]
    deltaI = param[2]
    pV = param[3]
    deltaV = param[4]

    # If drug present, adjust respective parameter that drug acts on by multiplying with drug scaling factor
    '''
    if Drug.type == 'rep':
        pV = pV * (1 - Drug.epsilon)

    if Drug.type == 'ent':
        beta_dot = beta_dot * (1 - Drug.epsilon)
    '''
    # Else leave parameters unchanged (e.g. placebo)

    # Initial conditions
    T0 = 4e+8  # Fixing this parameter 7e+7
    I0 = 0
    V0 = param[5]

    y_init = [T0, I0, V0]

    # Create an array of measurement time points (Day 1, 2, 3, ... 8)
    measurement_times = np.arange(0, max_time, 1)

    # Solve TIV
    sol = solve_ivp(rhs, t_span=(0, max_time), y0=y_init, method='BDF', t_eval=measurement_times)
    return sol  # NOTE: Returns raw (not log) values

def TLIV_drug_model(drug, param, max_time):  # where drug = [type, epsilon], param = [

    # Write differential equations describing TIV drug model (right-hand side)
    def rhs(t, y):
        T, L, I, V = y
        return [g * T * (1 - (T + I) / T0) - (beta * V * T),
                beta * V * T - (gamma * L),
                (gamma * L) - (deltaI * I),
                pV * I - (deltaV * V) - (beta * V * T)]

    # Input parameters
    g = 0.8
    beta = param[0]
    deltaI = param[1]
    pV = param[2]
    deltaV = param[3]
    gamma = param[4]

    # If drug present, adjust respective parameter that drug acts on by multiplying with drug scaling factor
    if drug[0] == 'rep':
        pV = pV * (1 - drug[1])

    if drug[0] == 'ent':
        beta = beta * (1 - drug[1])

    # Else leave parameters unchanged (e.g. placebo)

    # Initial conditions
    T0 = 4e+8  # Fixing this parameter
    L0 = 0
    I0 = 0
    V0 = param[5]

    y_init = [T0, L0, I0, V0]

    # Create an array of measurement time points (Day 1, 2, 3, ... 8)
    measurement_times = np.arange(0, max_time, 1)

    # Solve TIV
    sol = solve_ivp(rhs, t_span=(0, max_time), y0=y_init, method='RK45', t_eval=measurement_times)
    return sol  # NOTE: Returns raw (not log) values



def TIV_drug_ll(V_data, drug, param, max_time):  # Where I_data = I (infected individuals) as retrieved from data

    ll = 0

    sol = TIV_drug_model(drug, param, max_time)
    if sol.success == False:  #Reject parameter if solve_ivp failed
        ll = -1 * np.inf
        print('solve_ivp failed for param ' + str(param) + ' so rejecting value')
        return ll

    V_model = np.log(sol.y[2])  # Obtain model values for I, given new parameters
    exp_sd = 1.5  # sd of experimental measurements of viral titre

    for k in range(len(V_data)):
        new_ll = st.norm.logpdf(V_data[k], loc=V_model[k], scale=exp_sd)  # norm.logpdf(i, loc=mu, scale=sd)

        ll = ll + new_ll

    return ll


def TLIV_drug_ll(V_data, drug, param, max_time):  # Where I_data = I (infected individuals) as retrieved from data

    ll = 0

    sol = TLIV_drug_model(drug, param, max_time)
    if sol.success == False:  #Reject parameter if solve_ivp failed
        ll = -1 * np.inf
        print('solve_ivp failed for param ' + str(param) + ' so rejecting value')
        return ll

    V_model = sol.y[2]  # Obtain model values for I, given new parameters
    exp_sd = 1.5  # sd of experimental measurements of viral titre

    for k in range(len(V_data)):
        new_ll = st.norm.logpdf(np.log(V_data[k]), loc=np.log(V_model[k]), scale=exp_sd)  # norm.logpdf(i, loc=mu, scale=sd)

        ll = ll + new_ll

    return ll
