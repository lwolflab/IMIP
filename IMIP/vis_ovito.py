import os,subprocess,re,json
from dataclasses import dataclass
import pandas as pd
import numpy as np
from IMIP.Utils.util import *
from IMIP.Utils.variables import *
from IMIP.grid import Grid
os.environ['OVITO_GUI_MODE'] = '1'
from ovito.io import *
from ovito.modifiers import *
from ovito.data import *
from ovito.pipeline import *
from ovito.vis import *


@dataclass
class Ovito(Grid):
    renderer: OSPRayRenderer = OSPRayRenderer(refinement_iterations = 5,max_ray_recursion = 2,
                                              direct_light_intensity = 4.0,
                                              direct_light_angular_diameter = 0.009,
                                              ambient_brightness = 0.1,
                                              material_shininess = 5.0,
                                              material_specular_brightness = 0.0)
    vp: Viewport = Viewport(type=Viewport.Type.Top, fov=5,
                             camera_pos=(0, 0, 0), camera_dir=(0, 0, 1))
    grid_type: str = 'turbomole'
    grid_term: str = 'Tot'
    grid_surface: bool = False
    grid_r_scale: float = 1
    grid_l_scale: float = 1
    '''
    Author: Wentong Zhou

    This is a class for Ovito visualization of IMIP surfaces and molecules
    '''

    @classmethod
    def load_file(cls, ovito_configs=None):
        """Update and save a Grid instance."""
        with open(ovito_config_path, 'r') as file:
            data = json.load(file)
        if ovito_configs:
            for k, v in ovito_configs.items():
                data[k] = v
        with open('configs', 'w') as file:
            json.dump(data, file, indent=4)
        with open('configs', 'r') as file:
            # save all arguments to series of variables
            return cls(**json.load(file))


    def __post_init__(self):
        super().__post_init__()
        self.vis_mol()
        # check if mol_turbomole.xyz exists
        if (os.path.exists(f'{self.molecule.split(".")[0]}_turbomole.xyz') or
        os.path.exists(f'{self.molecule.split(".")[0]}_xtb.xyz')):
            self.vis_grid()




    def vis_mol(self):
        self.ovito_mol = import_file(self.molecule,columns=['Particle Type', 'Position.X', 'Position.Y', 'Position.Z'])
        def create_atoms(frame: int, data: DataCollection):
            atom_type = pd.read_csv(element_path)
            total_atoms = np.unique(self.mol_coords['atom_name'].to_numpy())
            for atom in total_atoms:
                data.particles_.particle_types_.type_by_name_(atom).radius = \
                    atom_type.loc[atom_type['atom_name'] == atom]['size'].to_numpy()[0]
                data.particles_.particle_types_.type_by_name_(atom).color = tuple(
                    [atom_type.loc[atom_type['atom_name'] == atom][i].to_numpy()[0] for i in ['R', 'G', 'B']])
        self.ovito_mol.modifiers.append(create_atoms)
        data = self.ovito_mol.compute()
        data.particles.vis.radius = 0.2
        data.cell.vis.enabled = False
        del data
        self.vis_bonds()
        self.vp.camera_pos = self.origin



    def vis_bonds(self):
        '''
        This is a function to create bonds between atoms in the molecule for Ovito visualization

        '''
        mod = CreateBondsModifier()
        mod.vis.width = 0.2
        mod.mode = CreateBondsModifier.Mode.VdWRadius
        self.ovito_mol.modifiers.append(mod)


    def vis_grid(self):
        '''
        This is a function to visualize the IMIP surfaces of the current molecule

        '''
        self.grid_eda = pd.read_csv(f'{self.molecule.split(".")[0]}_{self.grid_type}.xyz',
                                    header=None,sep='\s+',names=ovito_grid_header[self.grid_type],skiprows=2)
        median = self.grid_eda[self.grid_term].median()
        std = self.grid_eda[self.grid_term].std()
        self.start_value = median - std * self.grid_l_scale
        self.end_value = median + std * self.grid_r_scale
        self.ovito_grid = import_file(f'{self.molecule.split(".")[0]}_{self.grid_type}.xyz',
                                      columns=ovito_grid_header[self.grid_type])
        def create_grid(frame: int, data: DataCollection):
            radius = 0.2 if self.grid_surface else 0.1
            data.particles_.particle_types_.type_by_id_(1).radius = radius
        self.ovito_grid.modifiers.append(create_grid)
        data2 = self.ovito_grid.compute()
        if self.grid_surface:
            data2.particles.vis.enabled = False
        data2.cell.vis.enabled = False
        del data2
        self.vis_color(self.start_value,self.end_value)

        if self.grid_surface:
            self.vis_color(self.start_value, self.end_value)
            self.vis_surf()



    def vis_color(self,start_value,end_value):
        self.start_value = start_value
        self.end_value = end_value
        self.ovito_grid.modifiers.append(ColorCodingModifier(
            property=self.grid_term,
            gradient=ColorCodingModifier.Image(color_code_path),
            start_value=start_value,
            end_value=end_value))



    def vis_surf(self):
        mod = ConstructSurfaceModifier()
        mod.method = ConstructSurfaceModifier.Method.GaussianDensity
        mod.transfer_properties = True
        mod.grid_resolution = 400
        mod.isolevel = 0.1
        mod.vis.surface_color = (1.0, 0.04, 0.4)
        mod.vis.show_cap = False
        self.ovito_grid.modifiers.append(mod)

    def vis_legend(self):
        self.legend_spacing = (self.end_value - self.start_value) / 4
        self.legend = ColorLegendOverlay(offset_y = 0.1,
                                legend_size = 0.5,
                                aspect_ratio = 11.0,
                                font_size = 0.12,
                                label_size = 0.5,
                                format_string = '%.1f',
                                title = ovito_title_labels[self.grid_type][self.grid_term],
                                border_enabled = True,
                                font = 'Ubuntu,11,-1,5,700,0,0,0,0,0,0,0,0,0,0,1',
                                ticks_enabled = True,
                                ticks_spacing = self.legend_spacing)
        self.legend.modifier = self.ovito_grid.modifiers[1]
        self.vp.overlays.append(self.legend)


    def render(self):
        '''
        This is a function to render the current molecule
        '''
        self.ovito_mol.add_to_scene()
        self.vp.render_image((800,800),0,f'{self.molecule.split(".")[0]}.png',
                             renderer=self.renderer,alpha=True,crop=True)
        self.ovito_mol.remove_from_scene()


    def render_grid(self):
        '''
        This is a function to render the IMIP surfaces of the current molecule

        '''
        self.vis_legend()
        self.ovito_grid.add_to_scene()
        self.vp.render_image((800,800),0,
                             f'{self.molecule.split(".")[0]}_{self.grid_type}_{self.grid_term}.png',
                             renderer=self.renderer,alpha=True,crop=True)
        self.vp.overlays.remove(self.vp.overlays[0])
        self.ovito_grid.remove_from_scene()




# def batch_mol_render(camera_dir:tuple=(0,0,1),fov:int=6,zoom_all:bool=False):
#     '''
#     This is a function to render in batch all molecules in the current directory.
#     Make sure all molecules in the current folder are molecules, not grids!!!
#
#
#     :param camera_dir: camera direction
#     :param fov: deeper fov, more zoomed in
#     :param zoom_all: all molecules are zoomed in
#
#     '''
#     for mol in glob.glob('*.xyz'):
#         df_mol = pd.read_csv(mol,sep='\s+',header=None,skiprows=2)
#         if df_mol.shape[1] == 4:
#             ovi_configs = {"molecule": f"{mol}"}
#             mol = Ovito.load_file(ovito_configs=ovi_configs)
#             mol.vp.camera_dir = camera_dir
#             mol.vp.fov = fov
#             if zoom_all:
#                 mol.vp.zoom_all()
#             mol.render()
#
# if __name__ == '__main__':
#     parser = argparse.ArgumentParser(description='Molecules renderer')
#     parser.add_argument('-c', '--camera_dir', type=float, nargs=3, help='camera direction')
#     parser.add_argument('-f', '--fov', type=int, help='fov')
#     parser.add_argument('-z', '--zoom_all', action='store_true', help='zoom all')
#     args = parser.parse_args()
#     batch_mol_render(camera_dir=args.camera_dir,fov=args.fov,zoom_all=args.zoom_all)
