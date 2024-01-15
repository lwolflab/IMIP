import os,subprocess,re,json
from dataclasses import dataclass
import concurrent.futures
import pandas as pd
import numpy as np
from numpy import linspace
from IMIP.Utils.util import *
from IMIP.Utils.variables import *
from ase import Atoms
from ase.calculators.turbomole import Turbomole
from sklearn.cluster import KMeans
from xtb.ase.calculator import XTB
from turbomoleio.input.define import DefineRunner
from turbomoleio.input.utils import get_define_template
from turbomoleio.core.datagroups import DataGroups
from scipy.spatial.distance import cdist
from scipy.spatial import distance
from tqdm import tqdm


@dataclass
class Grid:
    '''
    Author: Wentong Zhou

    This class is used to generate a grid for a given molecule.

    The methodology was described in the paper:

    Amin Kiani, Wentong Zhou, and Lawrence M. Wolf "Using molecular interaction
    potential maps derived from energy decomposition analysis to interpret
    electronic structure and reactivity. " 2023 In Manuscript.

    '''
    # MolConfig attributes
    molecule: str
    chrg: int
    multiplicity: int

    # GridConfig attributes
    probe: str
    probe_chrg: int
    type: str
    spacing: float
    isovalue: float
    range_ve: float
    rad_val: float
    num_gridpoints: int

    # FilterConfig attributes
    filtration: float
    tangent_radius: float = 0.5
    tangent_max_nn: int = 40
    redirection: bool = True
    range_redir: float = 0.0001

    # ScanConfig attributes
    tight: bool = False


    @classmethod
    def load_file(cls, configs=None):
        """Update and save a Grid instance."""
        with open(config_path, 'r') as file:
            data = json.load(file)
        if configs:
            for k, v in configs.items():
                data[k] = v
        with open('configs', 'w') as file:
            json.dump(data, file, indent=4)
        with open('configs', 'r') as file:
            # save all arguments to series of variables
            return cls(**json.load(file))

    def __post_init__(self):
        self.extractor()
        types = {
            'cubic': self.cubic,
            'vdw': self.vdw,
            'molden': self.molden,
            'turbomole': self.turbomole,
            'xtb': self.xTB,
        }
        type = types.get(self.type) if self.type in types else self.readgrids
        if type is not None:
            type()
        if self.type not in types:
            pass
        elif self.redirection == False or self.type == 'vdw':
            self.tangent()
            self.filter()
            self.position()
        else:
            self.filter()
            self.redirect()
            self.position()
        try:
            self.xyz(self.gridpoints_filtered, name='filtered')
            print(15 * '-' + 'Grid generation finished' + 15 * '-')
            print(f'Number of gridpoints: {len(self.gridpoints_filtered)}')
        except:
            print(15 * '-' + 'Molecule loaded' + 15 * '-')


    def extractor(self):
        self.mol_coords = []
        f = open(self.molecule, 'r')
        lines = f.readlines()
        del lines[0:2]
        for line in lines:
            self.mol_coords.append(line.split())
        self.mol_coords = pd.DataFrame(self.mol_coords)
        self.mol_coords.columns = ['atom_name', 'X', 'Y', 'Z']
        self.mol_coords.dropna(axis=0, inplace=True)
        self.origin = [pd.to_numeric(self.mol_coords.iloc[:, i]).mean() for i in [1, 2, 3]]
        coords = self.mol_coords.iloc[:, 1:4].to_numpy(dtype=float)
        dist = [distance.euclidean(self.origin, coord) for coord in coords]
        self.grid_length = np.array(dist).max() * 1.414 + 8


    def cubic(self):
        length = self.grid_length
        spacing = self.spacing
        x = linspace(self.origin[0] - length / 2, self.origin[0] + length / 2, int(length / spacing) + 1)
        y = linspace(self.origin[1] - length / 2, self.origin[1] + length / 2, int(length / spacing) + 1)
        z = linspace(self.origin[2] - length / 2, self.origin[2] + length / 2, int(length / spacing) + 1)
        self.X, self.Y, self.Z = np.meshgrid(x, y, z)
        grid = np.stack((self.X, self.Y, self.Z), axis=-1)
        self.gridpoints_coords = pd.DataFrame(np.reshape(grid, (-1, 3)))
        self.gridpoints_coords.insert(0, 'atom_name', self.probe)
        self.gridpoints_coords.columns = ['atom_name', 'X', 'Y', 'Z']


    def vdw(self):
        elements = self.mol_coords['atom_name']
        coords = self.mol_coords.iloc[:, 1:4].to_numpy(dtype=float)
        points = []
        for i, element in enumerate(elements):
            inc = np.pi * (3 - np.sqrt(5))
            off = 2 / float(self.num_gridpoints)
            for k in range(self.num_gridpoints):
                y = k * off - 1 + (off / 2)
                r = np.sqrt(1 - y ** 2)
                phi = k * inc
                x = np.cos(phi) * r
                z = np.sin(phi) * r
                points.append([x, y, z])
        points = np.array(points)
        points_per_center = len(points) // len(elements)
        radii = pd.read_csv(radii_path, header = None, index_col = 0).to_dict()[1]
        radius_list = []
        for i, element in enumerate(elements):
            radius_list.append((radii[element] + radii[self.probe]) * self.rad_val)
            points[i * points_per_center:(i + 1) * points_per_center] *= (radii[element] + radii[self.probe]) * self.rad_val
            points[i * points_per_center:(i + 1) * points_per_center] += coords[i]
        distances = distance.cdist(points, coords, metric='euclidean')
        # Filter out points inside the radius of any center, with tolerance
        dist_arr = np.reshape(distances, (-1, len(elements)))
        mask = np.any(dist_arr - np.array(radius_list)[None, :] < -0.001, axis=1)
        self.gridpoints_coords = pd.DataFrame(points[~mask])
        self.gridpoints_coords.insert(0, 'atom_name', self.probe)
        self.gridpoints_coords.columns = ['atom_name', 'X', 'Y', 'Z']

    def process_origin(self, origin, name):
        grids_test = volumn_element(origin, length=self.spacing * 22, probe=self.probe,
                                    spacing=self.spacing, name=name)
        grids_test = grids_test[
            (grids_test['ele_density'] > self.isovalue) & (grids_test['ele_density'] < self.isovalue * self.range_ve)]
        return grids_test

    def rediection_origin(self, origin, name):
        count = 0
        e_diff = [1]
        while np.min(e_diff) > self.isovalue * self.range_redir:
            grids_test = volumn_element(origin, length=0.05 + count, probe=self.probe,
                                        spacing=0.005 * (1 - count), name=name)
            e_diff = np.abs(grids_test['ele_density'].to_numpy() - self.isovalue)
            min_diff_idx = np.argmin(e_diff)
            count += 0.02
        grid_min = grids_test.iloc[min_diff_idx:min_diff_idx + 1, :]
        return grid_min

    def xyz(self, gridpoints, name='grids'):
        """
        Export the gridpoints coordinates to a .xyz file
        gridpoints: pandas dataframe containing the coordinates of gridpoints
        """
        with open(f'{self.molecule.split(".")[0]}_{name}.xyz', 'w') as f:
            f.write(str(len(gridpoints)))
            f.write('\n\n')
            np.savetxt(f, gridpoints.to_numpy(), fmt='%s')

    def run_multiwfn(self):
        """
        Run Multiwfn to get the initial grid around the molecule

        """

        print(15 * '-' + 'wavefunction file found' + 15 * '-')
        molden_file = f'{self.molecule.split(".")[0]}.molden'

        # modify the molden file if the molecule contains ECP
        molden_modify(molden_file, num=len(self.mol_coords))

        # run Multiwfn to get the density
        with open(den_path, 'r') as f:
            lines = f.readlines()
        num_gridpoints = int(self.grid_length / (self.spacing * 13)) + 1
        # convert to Bohr because Multiwfn only accepts Bohr
        boundary = self.grid_length * 1.889725988 / 2

        # write the commands to Multiwfn input file
        lines.insert(3, f"{self.origin[0]} {self.origin[1]} {self.origin[2]}\n")
        lines.insert(4, f"{num_gridpoints} {num_gridpoints} {num_gridpoints}\n")
        lines.insert(5, f"{boundary} {boundary} {boundary}\n")
        with open(f'{self.molecule.split(".")[0]}_den.txt', 'w') as f:
            f.writelines(lines)
        den_cmd = f'{multiwfn_dir} {molden_file} < {self.molecule.split(".")[0]}_den.txt'
        process = subprocess.Popen(den_cmd,shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        stdout, stderr = process.communicate()

        # get the gridpoints coordinates and density
        df = pd.DataFrame(np.loadtxt('output.txt'))
        # the commented codes are for generating original cubic grid coordicates
        # not used in real calculations
        # self.origin_cubic = pd.DataFrame(np.loadtxt('output.txt'))
        # self.origin_cubic.insert(0, 'atom_name', self.probe)
        # self.origin_cubic.columns = ['atom_name', 'X', 'Y', 'Z', 'ele_density']
        # self.xyz(self.origin_cubic, name='init_cubic')
        self.origin_den = df[(df[3] > self.isovalue * 0.2) & (df[3] < self.isovalue * 1.8)]
        self.origin_den = self.origin_den.reset_index(drop=True)
        self.origin_den.insert(0, 'atom_name', self.probe)
        self.origin_den.columns = ['atom_name', 'X', 'Y', 'Z', 'ele_density']
        self.xyz(self.origin_den, name='init_den')
        self.grid_seeds = self.origin_den.iloc[:, 1:4].to_numpy()
        os.system(f'rm output.txt {self.molecule.split(".")[0]}_den.txt')


    def molden(self):
        """
        Generate the gridpoints using molden file
        """
        self.run_multiwfn()
        os.system(f'cp -f {self.molecule.split(".")[0]}.molden 1.molden')
        with concurrent.futures.ProcessPoolExecutor() as executor:
            # Use a list comprehension to create a list of Future objects
            futures = [executor.submit(self.process_origin, origin, i) for i, origin in
                       enumerate(self.grid_seeds)]
            results = []
            for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures),
                               desc='Gridpoints generation started', unit="cubes"):
                results.append(future.result())
        self.gridpoints_coords_den = pd.concat(results, axis=0, ignore_index=True)
        self.gridpoints_coords_den = self.gridpoints_coords_den.drop_duplicates()
        self.gridpoints_coords = self.gridpoints_coords_den.iloc[:, 0:4]


    def turbomole(self):
        """
        Generate the gridpoints using Turbomole single point calculation
        """
        os.mkdir('turbomole_temp')
        os.system(f'cp -f {self.molecule} turbomole_temp')
        os.chdir('turbomole_temp')
        os.system(f'x2t ../{self.molecule} > coord')
        dp_a = get_define_template(ridft_path)
        dp_a['charge'] = self.chrg
        dr_a = DefineRunner(parameters=dp_a)
        dr_a.run_full()
        a_c = DataGroups.from_file('control')
        a_c.cdg('ricore', '0')
        a_c.adg('disp4', '')
        a_c.to_file('control')
        os.system('ridft > ridft.out 2> out1')
        os.system(f'tm2molden norm < {molden_path} > out 2> out2')
        os.system(f'mv molden_std.input ../{self.molecule.split(".")[0]}.molden')
        os.chdir('..')
        os.system('rm -rf turbomole_temp')
        self.molden()


    def xTB(self):
        """
        Generate the gridpoints using xTB single point calculation

        """
        os.mkdir('xtb_temp')
        os.system(f'cp -f {self.molecule} xtb_temp')
        os.chdir('xtb_temp')
        subprocess.run([f'{xtb_dir}', self.molecule, '--molden', '--gfn 2', f'--chrg {self.chrg}'],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        os.system(f'cp -f molden.input ../{self.molecule.split(".")[0]}.molden')
        os.chdir('..')
        os.system('rm -rf xtb_temp')
        self.molden()

    def readgrids(self):
        """
        Read the gridpoints coordinates from a .xyz file

        """
        try:
            self.gridpoints_coords = pd.read_csv(f'{self.molecule.split(".")[0]}_{self.type}.xyz',
                                                                sep='\s+', skiprows=2, header=None)
            # check column number of self.gridpoints_coordinate
            if self.gridpoints_coords.shape[1] == 4:
                self.gridpoints_filtered = self.gridpoints_coords.iloc[:, 0:4]
                self.gridpoints_filtered.columns = ['atom_name', 'X', 'Y', 'Z']
            else:
                self.gridpoints_filtered = self.gridpoints_coords.iloc[:, 0:7]
                self.gridpoints_filtered.columns = ['atom_name', 'X', 'Y', 'Z', 'xn', 'yn', 'zn']
        except:
            pass

    @timer
    def filter(self):
        cluster = len(self.gridpoints_coords) // 10000
        cluster = 1 if cluster == 0 else cluster
        self.gridpoints_filtered = MIF_filter(clusters=cluster, threshold=self.filtration,
                                              gridpoints=self.gridpoints_coords,
                                              probe=self.probe)

    @timer
    def tangent(self):
        self.gridpoints_coords = tangent_coordinate(self, radius=self.tangent_radius, max_nn=self.tangent_max_nn)


    def redirect(self):
        self.origins_redircted = self.gridpoints_filtered.iloc[:,1:4].to_numpy()
        with concurrent.futures.ProcessPoolExecutor() as executor:
            # Use a list comprehension to create a list of Future objects
            futures = [executor.submit(self.rediection_origin,origin,i)
                       for i, origin in enumerate(self.origins_redircted)]
            results = []
            for future in tqdm(concurrent.futures.as_completed(futures),
                               total=len(futures),
                               desc='Redirection started',
                               unit="cubes"):
                results.append(future.result())
        self.gridpoints_redirected = pd.concat(results, axis=0, ignore_index=True)
        self.xyz(self.gridpoints_redirected, name='redirection')

        self.gridpoints_coords = pd.concat([self.gridpoints_redirected.iloc[:, 0:4],
                                           self.gridpoints_coords_den.iloc[:, 0:4]],
                                          axis=0, ignore_index=True)
        self.tangent()
        self.gridpoints_filtered = self.gridpoints_coords.iloc[0:len(self.gridpoints_redirected),:]
        self.gridpoints_filtered = self.gridpoints_filtered.drop_duplicates()

    def rotate_scan(self,stepsize:int,start_deg:int,rot_deg:int,normal,probe_coord,num):
        energy = []
        self.errors_opt = []
        for index in range(int(start_deg / stepsize), int(rot_deg / stepsize)):
            index_1 = stepsize * index
            axis = normal
            rotated, rotated_probe = rotate_coordinates(probe_coord,axis,index_1)
            new_molecule = pd.concat([self.mol_coords, rotated_probe], ignore_index=True)
            positions = np.vstack((self.mol_coords.iloc[:, 1:4].to_numpy().astype(float), rotated))
            atoms = new_molecule.iloc[:, 0].to_numpy()
            # convert the molecule to ase Atoms object
            numbers = np.array([self.element_dict[atom] for atom in atoms])
            atoms = Atoms(numbers=numbers, positions=positions)
            # call xTB python api to do rigid body scan single point calculations
            initial_charges = atoms.get_initial_charges()
            initial_charges[-1] = self.chrg + self.probe_chrg
            atoms.set_initial_charges(initial_charges)
            initial_uhf = atoms.get_initial_magnetic_moments()
            initial_uhf[-1] = (self.multiplicity - 1) // 2
            atoms.set_initial_magnetic_moments(initial_uhf)
            atoms.calc = XTB(method='GFN2-xTB', max_iterations=512)
            # check if calc report errors
            try:
                single_energy = atoms.get_potential_energy()
            except Exception as e:
                # export e as outputs
                with open('errors_scan.txt', 'a') as f:
                    f.write(f'the error was at {num} gridpoints and {index_1} degree\n')
                self.errors_opt.append([num, index_1])
                single_energy = 1000
            if len(energy) < 2:
                try:
                    energy.append(single_energy)
                except UnboundLocalError:
                    energy.append(1000)
            elif len(energy) == 2:
                if energy[-1] == energy[-2]:
                    break
                else:
                    energy.append(single_energy)
                    continue
            else:
                if single_energy > energy[-1] and energy[-1] < energy[-2]:
                    break
                else:
                    energy.append(single_energy)
                    continue
        try:
            # find the minimum energy and the corresponding degree
            min_energy_index = start_deg + np.argmin(energy) * stepsize
        except ValueError:
            min_energy_index = 0
        rotated, rotated_probe = rotate_coordinates(probe_coord, axis, min_energy_index)
        return energy, rotated, rotated_probe



    def position(self):
        '''
        Position the probe at each gridpoint and export the coordinates to series of .xyz files
        :return:
        '''
        print(15 * '-' + 'Grid positioning started' + 15 * '-')
        with tqdm(total=(len(self.gridpoints_filtered)), desc='Scan Initiated', unit="gridpoint") as pbar:
            for i, gridpoint in enumerate(self.gridpoints_filtered.to_numpy()):
                pbar.set_description(f"grid {i + 1} of {len(self.gridpoints_filtered)}")
                normal = self.gridpoints_filtered.iloc[:, 4:7].to_numpy()
                b = self.gridpoints_filtered.iloc[:, 1:4].to_numpy()
                self.element_dict = pd.read_csv(element_dict_path, index_col=0).to_dict()['0']
                if self.probe in self.element_dict.keys():
                    atom_probe = pd.DataFrame([f'{self.probe}', 0, 0, 0]).to_numpy()
                    xyz_generator("1",atom_probe.T,name=f'{dir_path}/Probes/{self.probe}.xyz')
                x, y, z, new_trans, new_probe = probe_gen_(normal=normal[i], b=b[i], probe=self.probe)
                if self.probe in set(['ch3+','ben','Sben','Tben','c2h2','c2h4','bf3','bh3']):

                    energy, rotated, rotated_probe = self.rotate_scan(10,0,90,
                                                                      normal[i],new_probe,i)
                    energy_new_range = (len(energy) - 1) * 10
                    energy, rotated, rotated_probe = self.rotate_scan(1,energy_new_range - 10,
                                                                      energy_new_range + 10,
                                                                      normal[i],new_probe,i)
                    if self.tight == True:
                        energy_new_range = (len(energy) - 1) * stepsize + energy_new_range - 10
                        energy, rotated, rotated_probe = self.rotate_scan(0.1,energy_new_range - 1,
                                                                          energy_new_range + 1,
                                                                          normal[i],new_probe,i)
                else:
                    rotated_probe = new_probe
                # create a new folder to store all positioned proble coordinates
                if not os.path.exists(f'{working_directory}/{self.molecule.split(".")[0]}_probe'):
                    os.mkdir(f'{working_directory}/{self.molecule.split(".")[0]}_probe')
                os.chdir(f'{working_directory}/{self.molecule.split(".")[0]}_probe')
                xyz_generator(f'{len(rotated_probe)}', rotated_probe.to_numpy(), name=f'probe_{i}.xyz')
                combined_coordinates = np.vstack((self.mol_coords.to_numpy(), rotated_probe.to_numpy()))
                xyz_generator(f'{len(self.mol_coords) + len(rotated_probe)}', combined_coordinates,
                              name=f'mol_probe_{i}.xyz')
                pbar.update(1)
        os.chdir(working_directory)
        print(15 * '-' + 'Grid positioning finished' + 15 * '-')


















