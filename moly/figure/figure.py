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
from ..layers.blob import get_blob
from .layouts import get_layout


class Figure():
    def __init__(self, surface="matte", figsize=None, **kwargs):

        self.fig = go.Figure()
        self.molecules = []
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
        self.assert_range(molecule.geometry)

        self.fig.update_layout(get_layout(molecule.geometry, self.resolution, self.max_range, self.min_range))
        self.fig.show()

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



def add_bonds(molecule, bonds, figure, surface):

    for idx1, idx2 in bonds:

        vec1 = molecule.geometry[idx1]
        vec2 = molecule.geometry[idx2]
        length = np.linalg.norm(vec2-vec1)
        R = rotation_matrix(np.array([0,0,1]), vec2 - vec1)

        if molecule.symbols[idx1] == molecule.symbols[idx2]:

            cyl = get_single_cylinder()
            cyl[:,2] *= length
            cyl = R.dot(cyl.T).T
            cyl += vec1

            mesh = get_bond_mesh(cyl, idx1, molecule.symbols, surface)
            figure.add_trace(mesh)

        if molecule.symbols[idx1] != molecule.symbols[idx2]:

            cyl = get_single_cylinder()
            cyl[:,2] *= length / 2
            cyl = R.dot(cyl.T).T
            cyl_1 = cyl + vec1
            cyl_2 = cyl + (vec1+vec2)/2

            mesh = get_bond_mesh(cyl_1, idx1, molecule.symbols, surface)
            figure.add_trace(mesh)
            mesh = get_bond_mesh(cyl_2, idx2, molecule.symbols, surface)
            figure.add_trace(mesh)
            

def add_atoms(molecule, figure, surface):
    sphere = np.array(get_sphere())
    #sphere = np.array(sphere)
    for atom, xyz in enumerate(molecule.geometry):
        mesh = get_sphere_mesh(sphere * (molecule.atomic_numbers[atom]/30 + 0.6), molecule.symbols[atom], xyz, surface)
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