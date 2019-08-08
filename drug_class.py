"""
Create class for drug
"""

from scipy.integrate import solve_ivp
import numpy as np

'''
Pharmacokinetics model (currently keeping dose constant though)
'''

def PK_model(drug, max_time):
    '''
    def D_model(param, max_time):
        # Write differential equations describing TIV drug model (right-hand side)
        def rhs(t, y):
            D = y
            return w * ka * D_admin * e ** (ka * (t - tadmin)) - ke * D

        # Input parameters
        w = 4.63
        ka = 11.04
        ke = 2.64

        # Initial conditions
        D_admin = 75 # mg
        y_init = D_admin

        # Create an array of measurement time points (Day 1, 2, 3, ... 8)
        measurement_times = np.arange(0, max_time, 1)

        # Solve TIV
        sol = solve_ivp(rhs, t_span=(0, max_time), y0=y_init, method='RK45', t_eval=measurement_times)
        return sol  # NOTE: Returns raw (not log) values
    '''

    epsilon_max = param[0]  #0.98
    EC50 = param[1]  # 30
    #D = D_model(param, max_time)
    D = param[2]

    epsilon = (epsilon_max * D) / (D + EC50)

    return epsilon  # [type, epsilon]


class Drug:
    # Instance attributes
    def __init__(self, type, EC50, epsilon_max, dose_admin, ka, ke, omega, t_admin):
        self.type = type
        self.EC50 = EC50
        self.epsilon_max = epsilon_max
        self.dose_admin = dose_admin
        self.ka = ka
        self.ke = ke
        self.omega = omega
        self.t_admin = t_admin


    def solve_epsilon(self, time):
        def solve_dose(self, time):
            # Solve d for given time point
            return self.omega * self.ka * self.dose_admin * e ** (self.ka * (time - self.tadmin)) - (self.ke * D)

        D = solve_dose(time)
        return (self.epsilon_max * D) / (D + self.EC50)

OST = Drug(type='ent_drug', EC50=36.1, epsilon_max=0.98, dose_admin=75, ka=11.04, ke=2.64, omega=4.63, t_admin=)

class EntDrug(Drug):
    pass


class RepDrug(Drug):
    pass


# Generate OST (rep drug) parameters
epsilon_max = 0.98  # From (Cao et al, 2017)
EC50 = 36.1  # IC50 against IBV = 36.1 nm
D = 75e+6  # 75 mg (though may need to do titration)

OST_param = [epsilon_max, EC50, D]

OST_epsilon = PK_model(OST_param, max_time=8)


class Param:
    pass