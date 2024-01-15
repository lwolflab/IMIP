import os,subprocess,re,json
from dataclasses import dataclass
import concurrent.futures
import pandas as pd
import numpy as np
from numpy import linspace
from IMIP.Utils.util import *
from IMIP.Utils.variables import *
from IMIP.grid import Grid
from IMIP.vis_ovito import Ovito
from tqdm import tqdm
from tabulate import tabulate



@dataclass
class xTBRunner(Grid):
    gfn: int = 2
    '''
    Author: Wentong Zhou

    A class to run xTB-IFF calulations and extract the grid data from the output file.
    For detailed information about xTB abd xTB-IFF, please refer to xTB documentation.
    '''

    @classmethod
    def load_file(cls, xtb_configs=None):
        with open(xtb_config_path, 'r') as file:
            data = json.load(file)
        if xtb_configs:
            for k, v in xtb_configs.items():
                data[k] = v
        with open('configs', 'w') as file:
            json.dump(data, file, indent=4)
        with open('configs', 'r') as file:
            return cls(**json.load(file))

    def __post_init__(self):
        [os.makedirs(dir, exist_ok=True) for dir in ['1', '2', '3']]
        self.eda_energy_df = pd.DataFrame()
        super().__post_init__()
        self.xtb_mol()
        with tqdm(total=(len(self.gridpoints_filtered)), initial=1, desc='Grid Initiated',
                  unit="gridpoint") as pbar:
            for i, gridpoint in enumerate(self.gridpoints_filtered.to_numpy()):
                pbar.set_description(f"grid {i + 1} of {len(self.gridpoints_filtered.to_numpy())} ")
                self.xtb_probe(i)
                df = self.xtb_eda(i)
                self.xtb_export(df)
                pbar.update(1)
        self.xtb_xyz()


    def xtb_mol(self):
        os.chdir(path_1)
        os.system(f'cp -f {working_directory}/{self.molecule} .')
        subprocess.run(f'{xtb_dir} {self.molecule} --chrg {self.chrg} '
                       f'--gfn {self.gfn} --uhf {(self.multiplicity -1)//2} --lmo', shell=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        os.rename('xtblmoinfo', '1')
        os.system(f'cp -f 1 {path_3}')
        os.chdir(working_directory)

    def xtb_probe(self,num):
        os.chdir(path_2)
        os.system(f'cp -f {working_directory}/{self.molecule.split(".")[0]}_probe/probe_{num}.xyz gridpoint.xyz')
        subprocess.run(f'xtb gridpoint.xyz --chrg {self.probe_chrg} --gfn {self.gfn} --uhf {(self.multiplicity -1)//2} --lmo', shell=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        if 'xtblmoinfo' in os.listdir():
            os.rename('xtblmoinfo', '2')
            os.system(f'cp -f 2 {path_3}')
        else:
            pass
        os.chdir(working_directory)



    def xtb_eda(self,num):
        os.chdir(path_3)
        with open('int.txt', 'w') as f:
            a = subprocess.run([xtbiff_dir, '1', '2', '-sp'], stdout=f, stderr=subprocess.PIPE)
            # read stderr to check if the calculation is normal
            if 'severe' not in a.stderr.decode():
                with open('int.txt', 'r') as f:
                    lines = f.readlines()
                    if '==' in lines[-1]:
                        gridpoint_xTB = np.array([line.split(':')[1] for line in lines[-11:-1]]).astype(float)
                    else:
                        gridpoint_xTB = np.zeros(10)
            else:
                gridpoint_xTB = np.zeros(10)
        df = pd.DataFrame(gridpoint_xTB.reshape(1, -1))
        df['Grid_Index'] = num
        df.set_index("Grid_Index", inplace=True)
        os.chdir(working_directory)
        return df


    def xtb_export(self, df_single:pd.DataFrame):
        self.eda_energy_df = pd.concat([self.eda_energy_df, df_single], axis=0)
        with open(f"energy_values_xtb_{self.molecule.split('.')[0]}.out", "w") as f:
            f.write(tabulate(self.eda_energy_df, headers=xtb_header, tablefmt="simple", floatfmt=".10f",
                             colalign=("center",) * len(xtb_header)))

    def xtb_xyz(self):
        df2 = pd.read_csv(f"energy_values_xtb_{self.molecule.split('.')[0]}.out",
                          skiprows=2, sep='\s+', header=None)
        df2 = df2.drop([0], axis=1)
        new = np.concatenate((self.gridpoints_filtered, df2), axis=1)
        with open(f'{self.molecule.split(".")[0]}_xtb.xyz', 'w') as f:
            f.write(str(len(df2)))
            f.write('\n\n')
            np.savetxt(f, new, fmt='%s')
        os.system('rm -rf 1 2 3')