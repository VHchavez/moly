"""

Creates main figure 

"""
import numpy as np
import plotly.graph_objects as go

from ..molecule.shapes import get_sphere
from ..molecule.shapes import get_single_cylinder
from ..molecule.shapes import rotation_matrix
from ..molecule.molecule_factory import molecule_factory

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
        self.max_range = None
        self.min_range = None

    def show(self):
        self.fig.show()

    def add_molecule(self, label, **kwargs):
        molecule = molecule_factory(label, **kwargs)
        add_bonds(molecule, self.fig, self.surface)
        add_atoms(molecule, self.fig, self.surface)
        self.molecules.append(molecule)
        self.fig.update_layout(get_layout(molecule.geometry, self.resolution))

    def add_blob(self, index=0, iso=0.01, color="Portland", opacity=0.2):
        volume = get_blob(self.molecules[index], iso, opacity, color)
        self.fig.add_trace(volume)

    def add_layer(self, trace):
        self.fig.add_trace(trace)


def add_bonds(molecule, figure, surface):

    for idx1, idx2 in molecule.bonds:

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
