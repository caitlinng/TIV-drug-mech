"""
Create class for drug
"""
import numpy as np
import math
from scipy.integrate import solve_ivp


class Drug:
    # Instance attributes
    def __init__(self, type, EC50, epsilon_max, dose_admin, t_admin, ka, ke, omega):
        self.type = type
        self.EC50 = EC50
        self.epsilon_max = epsilon_max
        self.dose_admin = dose_admin
        self.t_admin = t_admin
        self.ka = ka
        self.ke = ke
        self.omega = omega

    def solve_TIV(self, param, max_time):  # where drug = [type, epsilon], param = [

        # Input parameters
        g = 0.8
        beta_dot = param[0]
        beta = param[1]
        deltaI = param[2]
        pV = param[3]
        deltaV = param[4]

        # Initial conditions
        D0 = self.dose_admin
        T0 = 4e+8  # Fixing this parameter 7e+7
        I0 = 0
        V0 = param[5]

        y_init = [D0, T0, I0, V0]

        # If drug present, adjust respective parameter that drug acts on by multiplying with drug scaling factor
        if self.type == 'rep':
            # Write differential equations describing TIV drug model (right-hand side)
            def rhs(t, y):
                D, T, I, V = y
                return [self.omega * self.ka * self.dose_admin * math.e ** (self.ka * (t - self.t_admin)) - self.ke * D,
                        g * T * (1 - (T + I) / T0) - (beta_dot * V * T),
                        beta_dot * V * T - (deltaI * I),
                        pV * (1 - (self.epsilon_max * D) / (D + self.EC50)) * I - (deltaV * V) - (beta * V * T)]

        # Create an array of measurement time points (Day 1, 2, 3, ... 8)
        measurement_times = np.arange(0, max_time, 1)

        # Solve TIV
        sol = solve_ivp(rhs, t_span=(0, max_time), y0=y_init, method='BDF', t_eval=measurement_times)
        return sol  # NOTE: Returns raw (not log) values


OST = Drug(type='rep', EC50=36.1, epsilon_max=0.98, dose_admin=75, t_admin=0, ka=11.04, ke=2.64, omega=4.63)


sol = OST.solve_TIV(param=[0.000217, 0.000751, 3.3, 75, 35, 6], max_time=8)

# Plot

class EntDrug(Drug):
    pass

'''
class RepDrug(Drug):
    def solve_TIV(param, max_time):  # where drug = [type, epsilon], param = [

        # Input parameters
        g = 0.8
        beta_dot = param[0]
        beta = param[1]
        deltaI = param[2]
        pV = param[3]
        deltaV = param[4]

        # Initial conditions
        D0 = dose_admin
        T0 = 4e+8  # Fixing this parameter 7e+7
        I0 = 0
        V0 = param[5]

        y_init = [D0, T0, I0, V0]

        # If drug present, adjust respective parameter that drug acts on by multiplying with drug scaling factor
        if Drug.type == 'rep':
            # Write differential equations describing TIV drug model (right-hand side)
            def rhs(t, y):
                D, T, I, V = y
                return [omega * ka * dose_admin * math.e ** (ka * (t - t_admin)) - ke * D,
                        g * T * (1 - (T + I) / T0) - (beta_dot * V * T),
                        beta_dot * V * T - (deltaI * I),
                        pV * (1 - (epsilon_max * D) / (D + EC50)) * I - (deltaV * V) - (beta * V * T)]

        # Create an array of measurement time points (Day 1, 2, 3, ... 8)
        measurement_times = np.arange(0, max_time, 1)

        # Solve TIV
        sol = solve_ivp(rhs, t_span=(0, max_time), y0=y_init, method='BDF', t_eval=measurement_times)
        return sol  # NOTE: Returns raw (not log) values


class Param:
    def __init__(self, name, w, prior, guess):
        self.name = name
        self.w = w
        self.prior = prior
        self.guess = guess
'''