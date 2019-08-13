'''
Generate simulation data:
Measurement of (log?) viral titre over 8 days
For placebo-treated and 2 drug-treated (entry-blocking and replication-blocking drug)
'''

import TIV_drug_funcs as TIV
from test_fit import test_fit


# Generate placebo data
final_fit_param = [0.000217, 0.000751, 3.3, 75, 35, 6]

# Parameters expressed as log10 (e.g. deltaV = 1e+5)
T0 = 7e+7
beta_dot = 3e-8
beta = 5e-7
deltaI = 2
pV = 150 #12.6  # This is technically production rate of infectious virions though not all V are infectious
deltaV = 5
V0 = 1e+2

plc_param = [beta_dot, beta, deltaI, pV, deltaV, V0]

# Generate placebo data with Cao's parameters
#plc_data = TIV.TIV_drug_model(Drug='plc', param=plc_param, max_time=8)
#test_fit(TIV.TIV_drug_model, 'plc', plc_data.y[2], plc_param, len(plc_data.y[2]))


# Generate placebo data with fitted parameters from real data
#fitted_data = TIV.TIV_drug_model(Drug='plc', param=final_fit_param, max_time=8)
#test_fit(TIV.TIV_drug_model, 'plc', fitted_data.y[2], fitted_param, len(plc_data.y[2]))

'''
# Generate placebo data with TLIV model
gamma = 4
plc_TLIV_param = plc_param + [gamma]
print(plc_TLIV_param)
TLIV_plc_data = TIV.TLIV_drug_model(drug='plc', param=plc_TLIV_param, max_time=8)
test_fit(TIV.TLIV_drug_model, 'plc', TLIV_plc_data.y[2], plc_TLIV_param, len(TLIV_plc_data))
'''


# Generate entry-blocking drug data with fitted parameters
#ent_drug = ['ent', OST_epsilon]
#ent_drug_data = TIV.TIV_drug_model(drug=ent_drug, param=fitted_param, max_time=8)

#test_fit(TIV.TIV_drug_model, ent_drug, ent_drug_data.y[2], fitted_param, len(ent_drug_data.y[2]))


# Generate replication-blocking drug data with fitted parameters
#rep_drug = ['rep', OST_epsilon]
#rep_drug_data = TIV.TIV_drug_model(drug=rep_drug, param=fitted_param, max_time=8)

#test_fit(TIV.TIV_drug_model, rep_drug, rep_drug_data.y[2], fitted_param, len(rep_drug_data.y[2]))


# Testing real data
import TIV_Rubi_data
V_data = TIV_Rubi_data.plc_n1_raw

test_param = [0.000217, 0.000751, 3.3, 75, 35, 6]

#test_fit(model=TIV.TIV_drug_model, drug='plc', true_data=V_data, fitted_param=test_param, max_time=len(V_data))

