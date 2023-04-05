import numpy as np
import pandas as pd

def mle(model, obs, obs_e):
    model = model.astype('float')
    lk = 1.0/((2*np.pi*obs_e**2.0)**(0.5)) * np.exp(0.0 - (model - obs)**2.0/(2.0*obs_e**2.0))
    return lk

def get_age(obs, model):
    obs.insert(len(obs.columns), 'gyro_age', value = -9999)
    obs.insert(len(obs.columns), 'egyro_age_m', value = -9999)
    obs.insert(len(obs.columns), 'egyro_age_p', value = -9999)
    for index, row in obs.iterrows():
        bp_rp0, e_bp_rp0, prot, e_prot =  row['bprp0'], row['e_bprp0'],row['prot'], row['e_prot']
        starset = []
        sigma = 3
        starset = model[((np.abs(model['bp_rp0'] - bp_rp0) <= sigma*e_bp_rp0) &
                    (np.abs(model['Prot'] - prot) <= sigma*e_prot))]
        starset.index = range(len(starset))
        starset['lk_prot'] = mle(starset['Prot'], obs = prot, obs_e = e_prot)
        starset['lk_bprp0'] = mle(starset['bp_rp0'], obs = bp_rp0, obs_e = e_bp_rp0)
        starset['lk_total'] = starset['lk_bprp0']*starset['lk_prot']
        try:
            starset = starset.sample(frac = 0.01, weights = 'lk_total', random_state = 1000)
            low, center, high = np.percentile(starset['age'], [16,50,84])
            obs.loc[index, ['gyro_age', 'egyro_age_m', 'egyro_age_p']] =center, center - low,  high - center
        except:
            continue
    return obs
