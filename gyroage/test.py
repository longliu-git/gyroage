
import numpy as np
import pandas as pd
import age_grids
import gyro_age.get_gyro_age as get_gyro_age

obs = pd.DataFrame({'bprp0':0.826, 'e_bprp0':0.1, 'prot': 25, 'e_prot': 0.1},index=[0])
model = age_grids.get_age_grids(age_range=np.linspace(700, 6000, 100))
result = get_gyro_age.get_age(obs,model)
print(result)
