'''
Generate simulation data:
Measurement of (log?) viral titre over 8 days
For placebo-treated and 2 drug-treated (entry-blocking and replication-blocking drug)
'''

import TIV_drug_funcs as TIV
from test_fit import test_fit


# Generate placebo data

# Parameters expressed as log10 (e.g. deltaV = 1e+5)
T0 = 7e+7
g = 0.8
beta = 5e-7
deltaI = 2
pV = 12.6  # This is technically production rate of infectious virions though not all V are infectious
deltaV = 4
V0 = 1e+4

plc_param = [g, beta, deltaI, pV, deltaV, V0]

plc_data = TIV.TIV_drug_model(drug='plc', param=plc_param, max_time=8)

#test_fit(TIV.TIV_drug_model, 'plc', [plc_data.t, plc_data.y[2]], plc_param, len(plc_data.y[2]))


# Generate OST (rep drug) parameters
epsilon_max = 0.98  # From (Cao et al, 2017)
EC50 = 36.1  # IC50 against IBV = 36.1 nm
D = 75e+3  # 75 mg (though may need to do titration)

OST_param = [epsilon_max, EC50, D]

OST_epsilon = TIV.PK_model(OST_param, max_time=8)



# Generate entry-blocking drug data
ent_drug = ['ent', OST_epsilon]
ent_drug_data = TIV.TIV_drug_model(drug=ent_drug, param=plc_param, max_time=8)

#test_fit(TIV.TIV_drug_model, ent_drug, [ent_drug_data.t, ent_drug_data.y[2]], plc_param, len(ent_drug_data.y[2]))


# Generate replication-blocking drug data
rep_drug = ['rep', OST_epsilon]
rep_drug_data = TIV.TIV_drug_model(drug=rep_drug, param=plc_param, max_time=8)

#test_fit(TIV.TIV_drug_model, rep_drug, [rep_drug_data.t, rep_drug_data.y[2]], plc_param, len(rep_drug_data.y[2]))
