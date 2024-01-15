import os,subprocess,re,json
from dataclasses import dataclass
import concurrent.futures
import pandas as pd
import numpy as np
from numpy import linspace
from IMIP.Utils.util import *
from IMIP.Utils.variables import *
from IMIP.grid import Grid
from ase import Atoms
from ase.calculators.turbomole import Turbomole
from sklearn.cluster import KMeans
from xtb.ase.calculator import XTB
from turbomoleio.input.define import DefineRunner
from turbomoleio.input.utils import get_define_template
from turbomoleio.core.datagroups import DataGroups
from tabulate import tabulate
from scipy.spatial.distance import cdist
from scipy.spatial import distance
from tqdm import tqdm

@dataclass
class TurbomoleRunner(Grid):
    """
    Author: Wentong Zhou

    A class to run Turbomole LMO-EDA calulations and extract the grid data from the output file.
    For detailed information about Turbomole, please refer to Turbomole documentation.
    """
    functional: str = 'r2scan'
    basis: str = 'def2-SVP'
    radsize: int = 5
    damp: list = None
    shift: float = 0.3
    marij: bool = False
    lastdiag: bool = False
    solvent: str = 'None'
    restart: int = 0

    @classmethod
    def load_file(cls,tmol_configs=None):
        with open(tmol_config_path, 'r') as file:
            data = json.load(file)
        if tmol_configs:
            for k, v in tmol_configs.items():
                data[k] = v
        with open('configs', 'w') as file:
            json.dump(data, file, indent=4)
        with open('configs', 'r') as file:
            return cls(**json.load(file))


    def __post_init__(self):
        [os.makedirs(dir, exist_ok=True) for dir in ['1', '2', '3']]
        self.eda_energy_df = pd.DataFrame()
        if self.restart != 0:
            self.eda_energy_df = pd.read_csv(f"energy_values_{self.molecule.split('.')[0]}.out",
                                        skiprows=2,sep='\s+',header=None)
            self.restart = len(self.eda_energy_df)
            self.eda_energy_df = self.eda_energy_df.iloc[:,1:]
            self.eda_energy_df.columns = tmol_eda_terms
            self.extractor()
            self.gridpoints_filtered = pd.read_csv(f'{self.molecule.split(".")[0]}_filtered.xyz',
                                                   skiprows=2, sep='\s+', header=None)
        else:
            super().__post_init__()
        self.tmol_mol()
        with tqdm(total=(len(self.gridpoints_filtered)),initial= self.restart + 1,desc='Grid Initiated',unit="gridpoint") as pbar:
           for i, gridpoint in enumerate(self.gridpoints_filtered.iloc[self.restart:].to_numpy()):
                pbar.set_description(f"grid {i + 1} of {len(self.gridpoints_filtered.iloc[self.restart:].to_numpy())} ")
                self.tmol_probe(i + self.restart)
                self.tmol_eda(i + self.restart)
                df = self.tmol_parse(i + self.restart)
                self.tmol_export(df)
                pbar.update(1)
        self.tmol_xyz()



    def tmol_mol(self):
        """
        This is the function to run single point calculation on the molecule

        """
        os.chdir(path_1)
        os.remove('control') if os.path.exists('control') else None
        os.system(f'x2t ../{self.molecule} > coord')
        dp_a = get_define_template(ridft_path)
        dp_a['basis'] = self.basis
        dp_a['functional'] = self.functional
        dp_a['charge'] = self.chrg
        dr_a = DefineRunner(parameters=dp_a)
        dr_a.run_full()
        a_c = DataGroups.from_file('control')
        self.options = {"radsize": f" radsize  {self.radsize}"}
        a_c.mdgo('dft', options=self.options)
        a_c.cdg('ricore', '5000')
        a_c.adg('disp4', '')
        a_c.adg('pop', 'nbo')
        a_c.to_file('control')
        if self.solvent != 'None':
            os.system(f'cosmoprep < {cosmo_path} > cosout 2> out1')
            a_c = DataGroups.from_file('control')
            option = {"epsilon": "epsilon = Infinity ion",
                      "use_contcav": "use_contcav"}
            a_c.mdgo('cosmo', options=option)
            a_c.adg('dcosmo_rs', f'file={self.solvent}_25.pot')
            a_c.to_file('control')
        os.system('ridft > ridft.out 2> out1')
        os.chdir(working_directory)


    def tmol_probe(self,num):
        """
        This is the function to run single point calculation on the probe

        num: the index of the probe

        """
        os.chdir(path_2)
        os.remove('control') if os.path.exists('control') else None
        os.system(f'cp {working_directory}/{self.molecule.split(".")[0]}_probe/probe_{num}.xyz coord.xyz')
        os.system(f'x2t coord.xyz > coord')
        dp_b = get_define_template(ridft_path)
        dp_b['charge'] = self.probe_chrg
        dp_b['basis'] = self.basis
        dp_b['functional'] = self.functional
        dr_b = DefineRunner(parameters=dp_b, timeout=2)
        dr_b.run_full()
        b_c = DataGroups.from_file('control')
        b_c.mdgo('dft', options=self.options)
        b_c.cdg('ricore', '5000')
        b_c.adg('disp4', '')
        b_c.adg('pop', 'nbo')
        b_c.to_file('control')
        if self.solvent != 'None':
            os.system(f'cosmoprep < {cosmo_path} > cosout 2> out1')
            b_c = DataGroups.from_file('control')
            option = {"epsilon": "epsilon = Infinity ion",
                      "use_contcav": "use_contcav"}
            b_c.mdgo('cosmo', options=option)
            b_c.adg('dcosmo_rs', f'file={self.solvent}_25.pot')
            b_c.to_file('control')
        os.system('ridft > ridft.out 2> out1')
        os.chdir(working_directory)


    def tmol_eda(self,num):
        """
        This is the function to run LMO-EDA calculation on the molecule and the probe

        num: the index of the probe

        """
        os.chdir(path_3)
        os.remove('control') if os.path.exists('control') else None
        os.system(f'cp {working_directory}/{self.molecule.split(".")[0]}_probe/mol_probe_{num}.xyz coord.xyz')
        os.system(f'x2t coord.xyz > coord')
        dp_c = get_define_template(ridft_path)
        dp_c['charge'] = self.chrg + self.probe_chrg
        dp_c['basis'] = self.basis
        dp_c['functional'] = self.functional
        dr_c = DefineRunner(parameters=dp_c, timeout=2)
        dr_c.run_full()
        c_c = DataGroups.from_file('control')
        c_c.mdgo('dft', options=self.options)
        c_c.cdg('ricore', '5000')
        c_c.adg('disp4', '')
        c_c.cdg('scfdamp', f'   start={self.damp[0]}  step={self.damp[1]}  min={self.damp[2]}\n')
        c_c.cdg('scforbitalshift', f'noautomatic    closedshell={self.shift}\n')
        c_c.adg('subsystems', f'\n molecule#1 file={path_1}/control \n molecule#2 file={path_2}/control')
        if self.lastdiag == True:
            c_c.adg('lastdiag', '')
        if self.marij == True:
            c_c.adg('marij', '')
        c_c.adg('pop', 'nbo')
        c_c.to_file('control')
        if self.solvent != 'None':
            os.system(f'cosmoprep < {cosmo_path} > cosout 2> out1')
            c_c = DataGroups.from_file('control')
            option = {"epsilon": "epsilon = Infinity ion",
                      "use_contcav": "use_contcav"}
            c_c.mdgo('cosmo', options=option)
            c_c.adg('dcosmo_rs', f'file={self.solvent}_25.pot')
            c_c.to_file('control')
        os.system('promowa > out 2> out1')
        os.system('ridft > ridft.out 2> out1')
        os.chdir(working_directory)

    def tmol_parse(self,num) -> pd.DataFrame:
        """
        This is the function to extract the EDA energy terms from the output file

        :param num: the index of the probe
        :return: pd.DataFrame of values of EDA energy terms
        """
        os.chdir(path_3)
        df = extract_energy_values('ridft.out', 'dscf_problem')
        df['Corr_Disp'] = float(df['Correlation Interaction'].iloc[0]) + float(df['Dispersion Interaction'].iloc[0])
        df['Steric'] = float(df['Exchange-Repulsion'].iloc[0]) + float(df['Correlation Interaction'].iloc[0]) + float(
            df['Dispersion Interaction'].iloc[0])
        with open('ridft.out', 'r') as f:
            lines = f.readlines()
            for i,line in enumerate(lines):
                if 'Summary of Natural Population Analysis:' in line:
                    start_index = i
                    break
                else:
                    start_index = -1
            if start_index != -1:
                probe = pd.read_csv(f'{dir_path}/Probes/{self.probe}.xyz',
                                    skiprows=2,sep='\s+',header=None)
                nbo_text = lines[start_index + 5:start_index + 5 + len(self.mol_coords) + len(probe)]
                atoms = [float(line.split()[2]) for line in nbo_text]
                nbo_value = np.array(atoms[-len(probe):]).sum()
            else:
                nbo_value = self.probe_chrg
        df['NBO_charge'] = nbo_value - self.probe_chrg
        df["grid_index"] = num
        os.chdir(working_directory)
        if float(df['Steric'].iloc[0]) + float(df['Electrostatic Interaction'].iloc[0]) + float(
                df['Orbital Relaxation'].iloc[0]) == 0:
            os.makedirs("failed", exist_ok=True)
            os.system(f'cp -rf 3 failed/probe_3_{num}')
            os.system(f'cp -rf 2 failed/probe_2_{num}')
        df.set_index("grid_index", inplace=True)
        # convert hartree to kcal/mol convert df.iloc[:,-14:-1] to float
        df.iloc[:, -14:-1] = df.iloc[:, -14:-1].astype(float)
        hartree2kcal = 627.509
        df.iloc[:, -14:-1] = df.iloc[:, -14:-1] * hartree2kcal
        return df

    def tmol_export(self,df_single:pd.DataFrame):
        """
        This is the function to export the EDA energy terms to an output file

        """
        self.eda_energy_df = pd.concat([self.eda_energy_df, df_single], axis=0)
        with open(f"energy_values_{self.molecule.split('.')[0]}.out", "w") as f:
                f.write(tabulate(self.eda_energy_df,headers=tmol_header,tablefmt="simple", floatfmt=".10f",
                                 colalign=("center",)*len(tmol_header)))

    def tmol_xyz(self):
        """
        This is the function to export the grid points and energy values to an extended xyz file

        """
        df2 = pd.read_csv(f"energy_values_{self.molecule.split('.')[0]}.out",
                          skiprows=2, sep='\s+', header=None)
        df2 = df2.drop([0], axis=1)
        new = np.concatenate((self.gridpoints_filtered, df2), axis=1)
        with open(f'{self.molecule.split(".")[0]}_turbomole.xyz', 'w') as f:
            f.write(str(len(df2)))
            f.write('\n\n')
            np.savetxt(f, new, fmt='%s')
        os.system('rm -rf 1 2 3')
















