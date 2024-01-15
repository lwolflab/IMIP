import numpy as np
import pandas as pd
import EDA_analyzer.analyzer as eda
from EDA_analyzer.util import *
import shutil
import os
import glob

class frag(eda.EDA_analyzer):
    def __init__(self,supra,frag_A=["3"],run='supra',comm=False,probe='Li',boundary=10,grid_spacing=0.02,grids_source='turbomole',isovalue=0.00006,range=1.03,redirection=False):
        self.supra = supra
        self.line_numbers_or_ranges = frag_A
        self.path = os.getcwd()
        self.split_xyz_files(self.line_numbers_or_ranges, path=self.path)
        os.chdir(f'{self.path}/{supra.split(".")[0]}')
        test = super().__init__(molecule=supra,probe='Li',grids_source=grids_source,grid_spacing=grid_spacing,isovalue=isovalue,range=range,redirection=redirection)
        os.system(f'cp -f {supra.split(".")[0]}_filtered.xyz {supra.split(".")[0]}_grids.xyz ')
        os.system(f'cp -f {supra.split(".")[0]}_grids.xyz ../{supra.split(".")[0]}_A/{supra.split(".")[0]}_A_grids.xyz')
        os.system(f'cp -f {supra.split(".")[0]}_grids.xyz ../{supra.split(".")[0]}_B/{supra.split(".")[0]}_B_grids.xyz')
        os.chdir('../')
        if os.path.exists(f'{run}.py'):
            print('eda script found')
            for dir in [f'{supra.split(".")[0]}',f'{supra.split(".")[0]}_A',f'{supra.split(".")[0]}_B']:
                shutil.copy(f'{run}.py',dir)
                shutil.copy(f'{run}.sh',dir)
        else:
            pass
        if comm == True:
            command = f"find . -mindepth 1 -maxdepth 1 -type d \( ! -name . \) -exec bash -c \"cd '{{}}' && sbatch {run}.sh\" \;"
            subprocess.run(command, shell=True, check=True)
        else:
            pass


    def split_xyz_files(self, line_numbers_or_ranges, path=None):
        for file_path in glob.glob(os.path.join(path, '*.xyz')):
            base_name = os.path.splitext(os.path.basename(file_path))[0]
            fragment_A = os.path.join(path, f'{base_name}_A')
            fragment_B = os.path.join(path, f'{base_name}_B')
            Supra = os.path.join(path, f'{base_name}')

            if not os.path.exists(fragment_A):
                os.makedirs(fragment_A)
            if not os.path.exists(fragment_B):
                os.makedirs(fragment_B)
            if not os.path.exists(Supra):
                os.makedirs(Supra)

            shutil.copy(file_path, Supra)

            with open(file_path, 'r') as file:
                lines = file.readlines()
                copied_lines = []
                not_copied_lines = []
                for i, line in enumerate(lines):
                    if any(start <= i + 1 <= end for start, end in
                           [map(int, x.split('-')) for x in line_numbers_or_ranges if '-' in x]):
                        copied_lines.append(line)
                    elif str(i + 1) in line_numbers_or_ranges:
                        copied_lines.append(line)
                    else:
                        not_copied_lines.append(line)

                new_file_path = os.path.join(fragment_A, base_name + '_A.xyz')
                with open(new_file_path, 'w') as new_file:
                    new_file.write(f"{len(copied_lines)}\n")
                    new_file.write("\n")
                    new_file.writelines(copied_lines)

                not_copied_file_path = os.path.join(fragment_B, base_name + '_B.xyz')
                with open(not_copied_file_path, 'w') as not_copied_file:
                    not_copied_file.write(f"{len(not_copied_lines) - 2}\n")
                    not_copied_file.write("\n")
                    not_copied_file.writelines(not_copied_lines[2:])


