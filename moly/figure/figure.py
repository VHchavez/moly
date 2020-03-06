"""

Creates main figure 

"""
import numpy as np
import qcelemental as qcel
import plotly.graph_objects as go

from ..molecule.shapes import get_sphere
from ..molecule.shapes import get_single_cylinder
from ..molecule.shapes import rotation_matrix
#from ..molecule.molecule_factory import molecule_factory

from ..layers.bonds import get_bond_mesh
from ..layers.geometry import get_sphere_mesh
from ..layers.blob import get_blob, cube_to_molecule, get_cubes, get_cubes_traces, get_buttons
from .layouts import get_layout


class Figure():
    def __init__(self, surface="matte", figsize=None, **kwargs):

        self.fig = go.Figure()
        self.molecules = []
        self.molecule_labels = []
        self.surface = surface
        self.resolution = figsize
        self.max_range = np.zeros(3)
        self.min_range = np.zeros(3)

    def show(self):
        self.fig.show()

    def add_molecule(self, name, molecule):
        bonds = get_connectivity(molecule)
        add_bonds(molecule, bonds, self.fig, self.surface)
        add_atoms(molecule, self.fig, self.surface)
        self.molecules.append(molecule)
        self.molecule_labels.append(molecule)
        self.assert_range(molecule.geometry)

        self.fig.update_layout(get_layout(molecule.geometry, self.resolution, self.max_range, self.min_range))

    def add_cubes(self, directory, iso=0.03):
        cubes, details = get_cubes(directory)
        geometry, symbols, atomic_numbers, spacing, origin = cube_to_molecule(details[0]["name"]+".cube")
        bonds = qcel.molutil.guess_connectivity(symbols, geometry)
        add_bonds(geometry, symbols, bonds, self.fig, self.surface)
        add_atoms(geometry, atomic_numbers, symbols, self.fig, self.surface)

        geometry_traces = len(self.fig.data)

        traces = get_cubes_traces(cubes, spacing, origin, iso)
        for vol in traces:
            self.fig.add_traces(vol)

        
        button_list = get_buttons(details, geometry_traces)

        self.fig.update_layout(updatemenus=[dict(showactive=False,
                                                buttons=button_list,
                                                font={"family": "Helvetica",
                                                      "size" : 18},
                                                borderwidth=0
            ),
        ])

        self.fig.update_layout(get_layout(molecule.geometry, self.resolution, self.max_range * 1.5, self.min_range * 1.5))

    def add_blob(self, index=0, iso=0.01, color="Portland", opacity=0.2):
        volume = get_blob(self.molecules[index], iso, opacity, color)
        self.fig.add_trace(volume)

    def add_layer(self, trace):
        self.fig.add_trace(trace)

    def assert_range(self, geometry):
        
        mol_max_range = [max(geometry[:,axis]) for axis in range(3)]
        mol_min_range = [min(geometry[:,axis]) for axis in range(3)]

        self.max_range = [max(value) for value in zip(mol_max_range, self.max_range)]
        self.min_range = [min(value) for value in zip(mol_min_range, self.min_range)]



def add_bonds(geometry, symbols, bonds, figure, surface):

    for idx1, idx2 in bonds:

        vec1 = geometry[idx1]
        vec2 = geometry[idx2]
        length = np.linalg.norm(vec2-vec1)
        R = rotation_matrix(np.array([0,0,1]), vec2 - vec1)

        if symbols[idx1] == symbols[idx2]:

            cyl = get_single_cylinder()
            cyl[:,2] *= length
            cyl = R.dot(cyl.T).T
            cyl += vec1

            mesh = get_bond_mesh(cyl, idx1, symbols, surface)
            figure.add_trace(mesh)

        if symbols[idx1] != symbols[idx2]:

            cyl = get_single_cylinder()
            cyl[:,2] *= length / 2
            cyl = R.dot(cyl.T).T
            cyl_1 = cyl + vec1
            cyl_2 = cyl + (vec1+vec2)/2

            mesh = get_bond_mesh(cyl_1, idx1, symbols, surface)
            figure.add_trace(mesh)
            mesh = get_bond_mesh(cyl_2, idx2, symbols, surface)
            figure.add_trace(mesh)
            

def add_atoms(geometry, atomic_numbers, symbols, figure, surface):
    sphere = np.array(get_sphere())
    for atom, xyz in enumerate(geometry):
        mesh = get_sphere_mesh(sphere * (atomic_numbers[atom]/30 + 0.6), symbols[atom], xyz, surface)
        figure.add_trace(mesh)

def get_connectivity(molecule):

    mol_dict = molecule.dict()
    symbols = molecule.symbols
    geometry = molecule.geometry
    
    if not "connectivity" in mol_dict:
        return qcel.molutil.guess_connectivity(symbols, geometry)
    
    elif mol_dict["connectivity"] is None:
        return  qcel.molutil.guess_connectivity(symbols, geometry)

    elif "connectivity" in mol_dict:
        return self.dict["connectivity"]
