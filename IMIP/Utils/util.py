import os,re,time
import pandas as pd
import numpy as np
from numpy import linspace
import matplotlib.pyplot as plt
import open3d as o3d
from sklearn.cluster import KMeans
from scipy.linalg import norm
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation as R
def timer(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"{func.__name__} elapsed time: {elapsed_time} seconds")
        return result
    return wrapper

# directory of config.json
dir_path = os.path.dirname(os.path.realpath(__file__))
multiwfn_path = os.path.join(dir_path, 'Multiwfn/multiwfn.txt')


def get_paths_(name,args=('1','2','3')):
    os.makedirs(os.path.join(os.getcwd(), name),exist_ok=True)
    dir = [os.getcwd(),os.path.join(os.getcwd(),name)]
    for arg in args:
        os.makedirs(os.path.join(os.getcwd(),name,arg),exist_ok=True)
        dir.append(os.path.join(os.getcwd(),name,arg))
    return tuple(dir)


def get_paths():
    working_directory = os.getcwd()
    path_1 = working_directory + '/1'
    path_2 = working_directory + '/2'
    path_3 = working_directory + '/3'
    return working_directory, path_1, path_2, path_3



def output_parser(file_name: str, error_file: str, output_types: list) -> pd.DataFrame:
        output_values = {output_type: 0 for output_type in output_types}
        if os.path.isfile(error_file):
            df = pd.DataFrame(output_values, index=["0"])
            return df
        else:
            with open(file_name, 'r') as f:
                lines = f.read()
                output_types_re = "|".join(output_types)
                matches = re.finditer(output_types_re, lines)
                for match in matches:
                    output_type = match.group()
                    match = re.search(r'-?\d+\.\d+', lines[match.start():])
                    if match:
                        output_values[output_type] = match.group()
            df = pd.DataFrame(output_values, index=["0"])
            return df


def xyz_generator(num:str,contents:np.array,name:str='coord.xyz'):
    with open(name, 'w') as f:
        f.write(num)
        f.write('\n\n')
        np.savetxt(f, contents , fmt='%s')


def read_xyz(file_name):
    atoms = []
    coordinates = []
    with open(file_name, "r") as file:
        lines = file.readlines()
        for line in lines[2:]:
            line_split = line.split()
            atoms.append(line_split[0])
            coordinates.append([float(x) for x in line_split[1:]])
    return atoms, np.array(coordinates)





@timer
def MIF_filter(clusters=1,threshold=0.2,gridpoints=None,probe='F'):
    gridpoints = gridpoints.iloc[:,1:].to_numpy()
    kmeans = KMeans(clusters,n_init=10)
    kmeans.fit(gridpoints[:,0:3])
    labels = kmeans.labels_
    filtered = []
    for i in range(0,clusters,1):
        new = gridpoints[labels==i]
        distances = cdist(gridpoints[:,0:3][labels==i],gridpoints[:,0:3][labels==i])
        mask = np.ones(len(gridpoints[labels==i]), dtype=bool)
        for j in range(len(gridpoints[labels==i])):
            if mask[j]:
                mask[(j+1):][distances[j, (j+1):] < threshold] = False
        filtered.append(gridpoints[labels==i][mask])
    filtered = np.concatenate(filtered,axis=0)
    # check if filtered have three columns
    if filtered.shape[1] == 3:
        filtered = pd.DataFrame(filtered)
        filtered.insert(0,'atom_name', probe)
        filtered.columns = ['atom_name','X','Y','Z']
    elif filtered.shape[1] == 4:
        filtered = pd.DataFrame(filtered)
        filtered.insert(0, 'atom_name', probe)
        filtered.columns = ['atom_name','X','Y','Z','ele_density']
    elif filtered.shape[1] == 6:
        filtered = pd.DataFrame(filtered)
        filtered.insert(0, 'atom_name', probe)
        filtered.columns = ['atom_name','X','Y','Z','xn','yn','zn']
    return filtered


def extract_energy_values(file_name: str, error_file: str) -> pd.DataFrame:
    energy_types = ['Total Interaction energy', 'Electrostatic Interaction', 'Nuc---Nuc', '1-electron', '2-electron',
                    'Exchange-Repulsion', 'Exchange Int.', 'Repulsion', 'Orbital Relaxation', 'Correlation Interaction',
                    'Dispersion Interaction']
    energy_values = {energy_type: 0 for energy_type in energy_types}

    if os.path.isfile(error_file):
        df = pd.DataFrame(energy_values, index=["0"])
        return df
    else:
        with open(file_name, 'r',encoding='utf-8') as f:
            lines = f.read()
            energy_types_re = "|".join(energy_types)
            matches = re.finditer(f"{energy_types_re}(?=.*-?\d+\.\d+)", lines)
            for match in matches:
                energy_type = match.group()
                energy_values[energy_type] = re.search(r'-?\d+\.\d+', lines[match.end():]).group()

        df = pd.DataFrame(energy_values, index=["0"])
        return df


def probe_transform(A,B):
    BA = np.subtract(A, B)
    new_z = BA / norm(BA)
    random_vector = np.random.rand(3)
    proj = random_vector - np.dot(random_vector, new_z) * new_z
    new_x = proj / norm(proj)
    new_y = np.cross(new_x, new_z)
    return new_x,new_y,new_z

def probe_transform_(normal):
    new_z = normal/norm(normal)
    random_vector = np.random.rand(3)
    proj = random_vector - np.dot(random_vector, new_z) * new_z
    new_x = proj / norm(proj)
    new_y = np.cross(new_x, new_z)
    return new_x,new_y,new_z


def probe_gen(a,b,probe='ch3+'):
    module_dir = os.path.dirname(__file__)
    file_path = os.path.join(module_dir, f'{probe}.xyz')
    df = pd.read_csv(f'{file_path}', sep='\s+', header=None, skiprows=2)
    atom = df.iloc[:, 0]
    df = df.drop(columns=[0])
    x, y, z = probe_transform(a, b)
    new = df.to_numpy()
    new_trans = np.zeros(new.shape)
    for id in range(len(new)):
        R = np.column_stack((x, y, z))
        new_trans[id] = R @ (new[id].T) + b
    new_probe = pd.DataFrame(new_trans)
    new_probe.columns = ['X', 'Y', 'Z']
    new_probe.insert(0, 0, atom)
    return x,y,z,new_trans,new_probe


def probe_gen_(normal,b,probe='ch3+'):
    module_dir = os.path.dirname(__file__)
    file_path = os.path.join(module_dir, f'Probes/{probe}.xyz')
    df = pd.read_csv(f'{file_path}', sep='\s+', header=None, skiprows=2)
    atom = df.iloc[:, 0]
    df = df.drop(columns=[0])
    x, y, z = probe_transform_(normal)
    new = df.to_numpy()
    new_trans = np.zeros(new.shape)
    for id in range(len(new)):
        R = np.column_stack((x, y, z))
        new_trans[id] = R @ (new[id].T) + b
    new_probe = pd.DataFrame(new_trans)
    new_probe.columns = ['X', 'Y', 'Z']
    new_probe.insert(0, 0, atom)
    return x,y,z,new_trans,new_probe


def rotate_coordinates(probe:str, axis, angle):
    atoms = probe.iloc[:, 0]
    coords = probe.iloc[:, 1:4].to_numpy().astype(float)
    centroid = np.mean(coords, axis=0)
    translated_coords = coords - centroid
    rotation = R.from_rotvec(np.radians(angle) * axis)
    rotated_coords = rotation.apply(translated_coords)
    final_coords = rotated_coords + centroid
    final_probe = pd.DataFrame(final_coords)
    final_probe.columns = ['X', 'Y', 'Z']
    final_probe.insert(0, 0, atoms)
    final_probe.columns = ['atom_name', 'X', 'Y', 'Z']
    return final_coords,final_probe

def tangent_coordinate(test,radius=0.1,max_nn=40):
    new_group = test.gridpoints_coords.iloc[:,1:4].to_numpy()
    point_cloud = o3d.geometry.PointCloud()
    point_cloud.points = o3d.utility.Vector3dVector(new_group)
    search_param=o3d.geometry.KDTreeSearchParamHybrid(radius=radius, max_nn=max_nn)
    point_cloud.estimate_normals(search_param=search_param,fast_normal_computation=False)
    point_cloud.orient_normals_towards_camera_location(camera_location=test.origin)
    normals = pd.DataFrame(np.asarray(point_cloud.normals))
    new_coords = pd.concat([test.gridpoints_coords.iloc[:,0:4],normals],axis=1)
    new_coords.columns = ['atom_name', 'X', 'Y', 'Z','xn','yn','zn']
    return new_coords

def probe_visualizer(test,index):
    if len(test.gridpoints_filtered.columns) == 7:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(*test.gridpoints_filtered.iloc[:, 1:4].to_numpy().T, s=1, alpha=0.1, c='g')
        for i in index:
            normal = test.gridpoints_filtered.iloc[i,4:7].to_numpy()
            coordinate = test.gridpoints_filtered.iloc[i,1:4].to_numpy().reshape(1,-1)
            ax.scatter(*coordinate.T,c='r',s=10)
            ax.quiver(*test.gridpoints_filtered.iloc[i,1:4],*-normal,normalize=True,color='b')
# num is the length of the coordinates (atom list)
def molden_modify(test,num):
    with open(test, 'r') as f:
        lines_ = f.readlines()
        atom_line = lines_.index('[Atoms] AU\n')+1
        atom_lists = lines_[atom_line:atom_line+num]
    with open(test, 'w') as f:
        for i, atom in enumerate(atom_lists):
            atom_info_ = atom.split(' ')
            atom_info = [i for i in atom_info_ if i]
            core_potential = float(atom_info[-4])
            ECP = 0
            if 36 < core_potential < 54:
                ECP = 28
            elif 53 < core_potential < 58:
                ECP = 46
            elif 72 < core_potential < 86:
                ECP = 60
            lines_[atom_line+i] = atom.replace(atom_info[-4], f'{int(atom_info[-4]) - ECP}')
        f.writelines(lines_)

# parse the nocv output file from ORCA
def parse_nocv(test,list=5):
    with open(test, 'r') as f:
        lines_ = f.readlines()
        start_line = lines_.index('                     NOCV/ETS analysis                     \n')
    del_energy = float(lines_[start_line-4].split(' ')[-1])
    nocv_list = lines_[start_line+5:start_line+5+list]
    split_nocv_list = [i.split(' ') for i in nocv_list]
    eigen_list = [float(i[4]) for i in split_nocv_list]
    energy_list = [float(i[-2]) for i in split_nocv_list]
    energy_list.append(del_energy)
    return eigen_list,energy_list

def mol_to_xyz(mol, filename):
    atom_info = [f"{mol.atom_symbol(i)} {' '.join(map(str, mol.atom_coord(i)))}" for i in range(mol.natm)]
    with open(filename, 'w') as f:
        f.write(f"{mol.natm}\n\n")
        f.write('\n'.join(atom_info))


def volumn_element(origin,length,spacing,probe='Li',name=1):
    x = linspace(origin[0] - length/2, origin[0] + length/2, int(length/spacing) + 1)
    y = linspace(origin[1] - length/2, origin[1] + length/2, int(length/spacing) + 1)
    z = linspace(origin[2] - length/2, origin[2] + length/2, int(length/spacing) + 1)
    X, Y, Z = np.meshgrid(x, y, z)
    grid_ = np.stack((X, Y, Z), axis=-1)
    grid_ = np.reshape(grid_, (-1, 3))
    # round all the coordinates to 3 decimal places
    grid_ = np.round(grid_, 3)
    # convert grid from angstrom to bohr
    grid = grid_ / 0.529177
    gridpoints_coordinate = pd.DataFrame(grid_)
    gridpoints_coordinate.insert(0,'atom_name',probe)
    gridpoints_coordinate.columns = ['atom_name','X','Y','Z']
    with open(f'volumn_{name}.txt', 'w') as f:
        f.write(str(len(grid)) + '\n')
        np.savetxt(f,grid,fmt='%s',delimiter=' ')
    with open(multiwfn_path,'r') as f:
        lines = f.readlines()
        lines[3] = f'volumn_{name}.txt\n'
        lines[4] = f'grids_volumn_{name}.txt\n'
    with open(f'multiwfn_{name}.txt','w') as f:
        f.writelines(lines)
    os.system(f'Multiwfn 1.molden < multiwfn_{name}.txt > /dev/null')
    den = pd.read_csv(f'grids_volumn_{name}.txt',sep='\s+',header=None,skiprows=1)
    den = den.iloc[:,3]
    # insert den into gridpoints_coordinate as the last column
    gridpoints_coordinate_den = gridpoints_coordinate.copy()
    gridpoints_coordinate_den.insert(4,'ele_density',den)
    gridpoints_coordinate_den.columns = ['atom_name','X','Y','Z','ele_density']
    # remove txt files
    os.system(f'rm -rf volumn_{name}.txt grids_volumn_{name}.txt multiwfn_{name}.txt')
    return gridpoints_coordinate_den





    