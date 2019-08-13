"""
Create class for drug
"""
import numpy as np
import math
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


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


    def drug_factor(self, D):
        return 1 - ((self.epsilon_max * D) / (D + self.EC50))


#    def t_admin(self, t):


    def solve_TIV(self, param, max_time):
        # Input parameters
        g = 0.8
        beta_dot = param[0]
        beta = param[1]
        deltaI = param[2]
        pV = param[3]
        deltaV = param[4]

        # Initial conditions
        G0 = self.dose_admin
        D0 = 0
        T0 = 4e+8  # Fixing this parameter 7e+7
        I0 = 0
        V0 = param[5]

        y_init = [G0, D0, T0, I0, V0]

        # If drug present, adjust respective parameter that drug acts on by multiplying with drug scaling factor
        #if self.type == 'rep':
            # Write differential equations describing TIV drug model (right-hand side)
        def rhs(t, y):
            G, D, T, I, V = y  # Add new equation
            return [-self.ka * G,
                    (self.omega * self.ka * G) - (self.ke * D),
                    g * T * (1 - (T + I) / T0) - (beta_dot * V * T),
                    beta_dot * V * T - (deltaI * I),
                    pV * self.drug_factor(D) * I - (deltaV * V) - (beta * V * T)]


        # Create an array of measurement time points (Every hour until the next point of administration)
        measurement_times = np.arange(0, self.t_admin, 1)

        # Solve TIV

        # Solve ODEs for initial conditions until the next time point for drug administration, add solutions to all_sol
        t = 0  # time tracker
        all_sol = solve_ivp(rhs, t_span=(0, self.t_admin), y0=y_init, method='BDF', t_eval=measurement_times)
        t += self.t_admin

        while t < max_time:
            # Update initial conditions to the conditions at the end of the previously solved ODEs
            G0 = all_sol.y[0][-1].copy() + self.dose_admin  # G (drug administered) is the amount of drug leftover + new amount of drug administered
            D0 = all_sol.y[1, -1].copy()
            T0 = all_sol.y[2, -1].copy()
            I0 = all_sol.y[3, -1].copy()
            V0 = all_sol.y[4, -1].copy()

            y_init = [G0, D0, T0, I0, V0]

            # Solve again for the next __ hours until the drug is administered again
            sol = solve_ivp(rhs, t_span=(0, self.t_admin), y0=y_init, method='BDF', t_eval=measurement_times)

            # Add solution to previous chains of solutions
            all_sol.t = np.append(all_sol.t, np.arange(all_sol.t[-1] + 1, all_sol.t[-1] + self.t_admin + 1))

            all_sol.y = np.hstack((all_sol.y, sol.y))

            # Update time
            t += self.t_admin

        print (all_sol)
        return all_sol  # NOTE: Returns raw (not log) values

    # Plot pharmacokinetics
    def plot_PK(self, param, max_time):
        sol = self.solve_TIV(param, max_time)

        fig1, ax1 = plt.subplots()
        ax1.plot(sol.t, sol.y[0], 'o', color='blue')  # Plot guess fitted parameter solution

        ax1.set_title('PK')
        ax1.set_xlabel('Hours')
        ax1.set_ylabel('D (ng/mL)')

        #    plt.savefig(str(MCMC_param) + '_Vcurve' + '.png')  # Save as .png

        return plt.show()

    # Step-wise dose_admin function

    # For loop, Iteration every 12 hours, y_init = last conditions and D = leftover (d_admin * math.e ** (-k * t (i.e. 12)) + new dose admin (75 mg)

    # Plot viral load curve given param
    def plot_TIV(self, param, max_time):
        sol = self.solve_TIV(param, max_time)

        fig1, ax1 = plt.subplots()
        ax1.plot(sol.t, sol.y[4], '-', color='blue')  # Plot guess fitted parameter solution

        ax1.set_title('Viral Load')
        ax1.set_xlabel("Hours")
        ax1.set_ylabel("Viral load (log TCID50)")
        ax1.set_yscale('log')

        #    plt.savefig(str(MCMC_param) + '_Vcurve' + '.png')  # Save as .png

        return plt.show()


# Defining drug and drug parameters
OST = Drug(type='rep', EC50=36.1, epsilon_max=0.98, dose_admin=75, t_admin=12, ka=0.46, ke=0.11, omega=4.63)

# Plot pharmacokinetics of OST
sol = OST.plot_PK(param=[0.000217, 0.000751, 3.3, 75, 35, 6], max_time=8*12)

# Plot viral load curve of OST-treated ferret
OST.plot_TIV(param=[0.000217, 0.000751, 3.3, 75, 35, 6], max_time=8*12)


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