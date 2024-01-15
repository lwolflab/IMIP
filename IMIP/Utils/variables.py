import os,subprocess,re,json
from IMIP.Utils.util import *
# directory of some important input files
dir_path = os.path.dirname(os.path.realpath(__file__))
config_path = os.path.join(dir_path, 'Config.json')
tmol_config_path = os.path.join(dir_path, 'Turbomole/Config.json')
ovito_config_path= os.path.join(dir_path, 'Visualizer/Config.json')
xtb_config_path = os.path.join(dir_path, 'xTB/Config.json')
den_path = os.path.join(dir_path, 'Multiwfn/den.txt')
multiwfn_path = os.path.join(dir_path, 'Multiwfn/multiwfn.txt')
molden_path = os.path.join(dir_path, 'Turbomole/molden.inp')
ridft_path = os.path.join(dir_path, 'Turbomole/ridft.yaml')
radii_path = os.path.join(dir_path, 'Vdw/radii.csv')
cosmo_path = os.path.join(dir_path, 'Turbomole/cosmo')
element_path = os.path.join(dir_path, 'Visualizer/elements.csv')
element_dict_path = os.path.join(dir_path, 'Visualizer/elements_dict.csv')
color_code_path = os.path.join(dir_path, 'Visualizer/colormap_spectral.png')
working_directory, path_1, path_2, path_3 = get_paths()
xtb_dir = subprocess.run(['which', 'xtb'], stdout=subprocess.PIPE).stdout.decode().strip()
xtbiff_dir = subprocess.run(['which','xtbiff'],stdout=subprocess.PIPE).stdout.decode().strip()
multiwfn_dir = subprocess.run(['which', 'Multiwfn'], stdout=subprocess.PIPE).stdout.decode().strip()




# lists of EDA energy terms
tmol_eda_terms =   ['Total Interaction energy', 'Electrostatic Interaction',
                  'Nuc---Nuc', '1-electron', '2-electron',
                  'Exchange-Repulsion', 'Exchange Int.', 'Repulsion',
                  'Orbital Relaxation', 'Correlation Interaction',
                  'Dispersion Interaction','Corr_Disp','Steric','NBO_charge']

tmol_header =      ["Grid_Index", "Tot", "Electro", "Nuc_Nuc", "1e", "2e", "Exc_Rep", "Exc",
                    "Rep", "Orb_Relax", "Corr", "Disp","Corr_Disp",'Steric','NBO_charge']

xtb_header  =      ['Grid_Index','E_Pauli', 'E_disp_ATM', 'E_disp_2B', 'E_disp_total', 'E_ES_atom',
                    'E_ES_LMO','E_ES_total', 'E_induction', 'E_CT', 'Eint_total,gas']

ovito_grid_header = {'turbomole':['Particle Type', 'Position.X', 'Position.Y', 'Position.Z',
                    'x_norm', 'y_norm', 'z_norm','Tot', 'Electro', 'Nuc_Nuc',
                    '1e', '2e', 'Ex-Rep','Ex', 'Rep', 'Orb', 'Corr', 'Disp', 'Disp_Corr',
                    'Steric', 'NBO'],
                     'xtb':['Particle Type', 'Position.X', 'Position.Y', 'Position.Z',
                    "x_norm", "y_norm", "z_norm",'Pauli', 'E_disp_ATM', 'E_disp_2B', 'Disp', 'E_ES_atom',
                    'E_ES_LMO','Electro', 'Induction', 'Orb', 'Tot']}

ovito_title_labels= {'turbomole':{'Tot': '<sub>&Delta;</sub>E<sub>int<\sub>',
                    'Electro': '<sub>&Delta;</sub>E<sub>es<\sub>',
                    'Ex-Rep': '<sub>&Delta;</sub>E<sub>ex_rep<\sub>',
                    'Orb': '<sub>&Delta;</sub>E<sub>orb<\sub>',
                    'Corr': '<sub>&Delta;</sub>E<sub>corr<\sub>',
                    'Disp': '<sub>&Delta;</sub>E<sub>disp<\sub>',
                    'Disp_Corr': '<sub>&Delta;</sub>E<sub>disp_corr<\sub>',
                    'Steric': '<sub>&Delta;</sub>E<sub>steric<\sub>',
                    'NBO': 'NBO'},
                     'xtb':{'Tot': '<sub>&Delta;</sub>E<sub>int<\sub>',
                    'Electro': '<sub>&Delta;</sub>E<sub>es<\sub>',
                    'Disp': '<sub>&Delta;</sub>E<sub>disp<\sub>',
                    'Induction': '<sub>&Delta;</sub>E<sub>ind<\sub>',
                    'Orb': '<sub>&Delta;</sub>E<sub>CT<\sub>',
                    'Pauli': '<sub>&Delta;</sub>E<sub>Pauli<\sub>',
                    'E_disp_ATM': '<sub>&Delta;</sub>E<sub>dispATM<\sub>',
                    'E_disp_2B': '<sub>&Delta;</sub>E<sub>disp2B<\sub>',
                    'E_ES_atom': '<sub>&Delta;</sub>E<sub>ESatom<\sub>',
                    'E_ES_LMO': '<sub>&Delta;</sub>E<sub>ESlmo<\sub>'}}