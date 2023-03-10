import pandas as pd
import xarray as xr
import numpy as np

def get_datasources():
    datasources = {
        'nacl': {
            'o': '/Users/ianbillinge/Documents/yiplab/projects/tsse/oxyanions/nacl/NaCl_25C_equil/sodium chloride dipa equil - export_org_25C.csv',
            'a': '/Users/ianbillinge/Documents/yiplab/projects/tsse/oxyanions/nacl/NaCl_25C_equil/sodium chloride dipa equil - export_aq_25C.csv'},
        'nano3': {
            'o': '/Users/ianbillinge/Documents/yiplab/projects/tsse/oxyanions/nano3/NaNO3_25C_equil/sodium nitrate dipa equil - export_org_25C.csv',
            'a': '/Users/ianbillinge/Documents/yiplab/projects/tsse/oxyanions/nano3/NaNO3_25C_equil/sodium nitrate dipa equil - export_aq_25C.csv'},
        'naac': {
            'o': '/Users/ianbillinge/Documents/yiplab/projects/tsse/oxyanions/naac/sodium acetate dipa equil - export_org_25C.csv',
            'a': '/Users/ianbillinge/Documents/yiplab/projects/tsse/oxyanions/naac/sodium acetate dipa equil - export_aq_25C.csv'},
        'na2so4': {
            'o': '/Users/ianbillinge/Documents/yiplab/projects/tsse/oxyanions/na2so4/SS_25C_equil/sodium sulfate dipa equil - export_org_25C.csv',
            'a': '/Users/ianbillinge/Documents/yiplab/projects/tsse/oxyanions/na2so4/SS_25C_equil/sodium sulfate dipa equil - export_aq_25C.csv'},
        'na2co3': {
            'o': '/Users/ianbillinge/Documents/yiplab/projects/tsse/oxyanions/na2co3/sodium carbonate dipa equil - export_org_25C.csv',
            'a': '/Users/ianbillinge/Documents/yiplab/projects/tsse/oxyanions/na2co3/sodium carbonate dipa equil - export_aq_25C.csv'},
        'na3po4': {
            'o': '/Users/ianbillinge/Documents/yiplab/projects/tsse/oxyanions/na3po4/sodium phosphate dipa equil - export_org_25C.csv',
            'a': '/Users/ianbillinge/Documents/yiplab/projects/tsse/oxyanions/na3po4/sodium phosphate dipa equil - export_aq_25C.csv'}
    }
    return datasources

def get_na_params():
    params = {
        'mw': 22.989,
        'pz': 0.258,
        'H_sol': -416,
        'G_sol': -375,
        'Hf':-240,
        'Gf':-261
    }
    return params

def get_params(normalized):
    params = {  # [nacl, nano3, naac, na2so4, na2co3, na3po4]
        'label': ['NaCl', r'NaNO$_{3}$', r'NaCH$_{3}$CHOO', r'Na$_{2}$SO$_{4}$', r'Na$_{2}$CO$_{3}$',
                  r'Na$_{3}$PO$_{4}$'],
        'mw': [58.44, 84.99, 83.04, 142.04, 105.99, 163.94],
        'mw_i': [35.45, 62.004, 2*12.011 + 2*15.999 + 3*1.007, 96.06, 60.008, 94.9714],
        'charge': [1, 1, 1, 2, 2, 3],
        'tchem_r_i': [0.175, 0.179, 0.162, 0.258, 0.178, 0.238],
        'Hsol': [3.88, 20.50, -17.32,
                 # CRC handbook, 5-108 [kJ/mol, 25 C]. First three values (Nacl, nano3, naac are verified.
                 -2.343,
                 # [kJ/mol] #Pitzer 1938. He reported -568 cal/mol assuming he meant cal/mol and not kcal/mol. I then converted to kJ/mol.
                 -26.61,  # probably na2co3 based on wikipedia
                 None  # haven't been able to find a value for Na3PO4.
                 ],
        'Hf_i': [-167.2, -207.4, -486.0, -909.3, -677.1, -1277.4],  # CRC handbook, 5-65 [kJ/mol]
        'S_in': [56.5, 146.4, 86.6, 20.1, -56.9, -220.5],  # CRC handbook, 5-65 [J/molK]
        'Gf_i': [-131.2, -111.3, -369.3, -744.5, -527.8, -1018.7],  # CRC handbook, 5-65 [kJ/mol]
        'Cp_i': [-136.4, -86.6, -6.3, -293.0, None, None],  # CRC Handbook 5-65 [J/molK]
        'r2water': [0.318, 0.345, 0.35, 0.381, None, 0.306],
        'w_sat': [0.2645, 0.477, 124. / (124. + 100.), 0.2194, 0.235, 0.126],
        # CRC handbook 5-179. NaaC value comes from Wikipedia and is interpolation.
        'jones_dole': [-0.007, -0.046, 0.250, 0.208, 0.294, 0.590],
        'pz_an_chemicalize': [2.11, 3.68, 5.03, 6.28, 3.88, 5.56],
        # polarizability of anion Calculated using Chemicalize.com. [Å^3]
        'pz_chemicalize': [2.02, 3.68, 4.91, 6.04, 3.64, 5.56],# polarizability of whole salt Calculated using Chemicalize.com. [Å^3]
        'pz_i' : [3.42, 4.13, 5.50, 5.47, 4.54, 5.99], #from Ion properties, Markus, [A^3]
        'H_sol_i': [-367, -312, -425, -909, -677, -1018], #from Ion properties [kJ/mol]
        'G_sol_i':[-347, -306, -373, -1090, -479, -2773], #from Ion properties [kJ/mol]
    }
    extra_params = {
        'nacl': {
            'Hsol': 3.88,
            'a25': 3.532,
            'aP': 2.54,
            'Hsol_Parker': 929,  # cal/mol. From Parker, 1965
            'adhn': 0
        },
        'nano3': {
            'a25': 4.313,
            'aP': 3.011,
            'Hsol_Parker': 4900,  # cal/mol #Parker et al 1965? NIST
        },
        'naac': {
            'Hsol': -3.973,  # NOTE: kcal/mole Harriet Snell et al.
            'Hsol_Parker': -4140,  # cal/mol Parker
        },
        'na2so4': {
            'a25': 4.907,
            'aP': 4.280,
        },
        'na2co3': {
        },
        'na3po4': {
            'adhn': 5.1}
    }
    xsat = []
    for wsat, mw in zip(params['w_sat'], params['mw']):
        xsat.append(wsat/mw / (wsat/mw + (1. - wsat)/18.015))
    params['x_sat'] = xsat

    norm_params = {}
    for key, val in params.items():
        try:
            std = np.std(np.array(val))
            # print(key)
            # print('std is {:1.3f}'.format(std))
            norm_params[key] = [v / std for v in val]
        except TypeError as te:
            print(te)

    params = get_difference(params) #note: the difference parameters will not work if normalized is True.

    if normalized == True:
        return norm_params
    else:
        return params

def get_difference(params):
    na_params = get_na_params()
    params['d_pz'] = [pz - na_params['pz'] for pz in params['pz_i']]
    params['d_pz_scaled'] = [pz - z * na_params['pz'] for pz, z in zip(params['pz_i'], params['charge'])]
    params['H_sol'] = [z * na_params['H_sol'] + hsol for z, hsol in zip(params['charge'], params['H_sol_i'])]
    params['G_sol'] = [z * na_params['G_sol'] + Gsol for z, Gsol in zip(params['charge'], params['G_sol_i'])]

 


    return params

def load_data_params_to_ds(normalized = False):
    datasources = get_datasources()
    params = get_params(normalized)
    ccolumns = {'Sample': 'sample',
                'w_{water}': 'w_w',
                '+-': 'dw_w',
                'w_{dipa}': 'w_a',
                '+-.1': 'dw_a',
                'w_{salt}': 'w_s',
                '+-.2': 'dw_s',
                'x_{dipa}': 'x_a',
                'x_{salt}': 'x_s',
                'x_{water}': 'x_w',
                'Phase': 'phase'}
    d = ['salt', 'T', 'sample', 'phase']
    first = 0
    for slt, salt_vals in datasources.items():
        for phase, filepath in salt_vals.items():
            df = pd.read_csv(filepath).rename(columns=ccolumns).dropna()
            df['salt'] = slt
            df = df.set_index(d)

            if phase == 'o':
                ds_o = df.to_xarray()
            elif phase == 'a':
                ds_a = df.to_xarray()

        ds_salt = xr.concat([ds_o, ds_a], dim='phase')

        if first == 0:
            ds = ds_salt
            first += 1
        else:
            ds = xr.concat([ds, ds_salt], dim='salt')

    ds = ds.drop_duplicates('salt')
    ds = ds.drop_duplicates(dim='phase')
    ds = ds.drop_duplicates(dim='sample').drop_duplicates(dim='T')
    for param, values in params.items():
        ds = ds.assign_coords({param: ('salt', values)})

    if normalized == True:
        ds['x_s'] = (ds['x_s'] - ds['x_s'].mean())/ds['x_s'].std()
    return ds

def main():
    ds = load_data_params_to_ds(normalized = False)
    print(ds.sel({'T': 25, 'phase': 'o'})['x_s'].dropna(dim = 'sample', how = 'all'))
    ds = load_data_params_to_ds(normalized = True)
    print(ds.sel({'T': 25, 'phase': 'o'})['x_s'].dropna(dim = 'sample', how = 'all'))
    return
#
if __name__ == '__main__':
    main()