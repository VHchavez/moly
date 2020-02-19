"""

Creates main figure 

"""
#import sys
#sys.path.append("..")

import plotly.graph_objects as go
#from .molecule.adapters import molecule_factory
from ..layers.bonds import get_bond_mesh
from ..layers.geometry import get_sphere_mesh
from ..layers.blob import get_blob
from .layouts import get_layout

def plot(molecule, surface="matte"):

    fig = go.Figure()

    #Add bonds
    for bond in range(len(molecule.bonds)):
        mesh_a, mesh_b = get_bond_mesh(molecule.cilinders, molecule.bonds, bond, molecule.symbols, surface)
        fig.add_trace(mesh_a)
        fig.add_trace(mesh_b)

    #Add atoms
    for atom, xyz in enumerate(molecule.geometry):
        mesh = get_sphere_mesh(molecule.spheres[atom], molecule.symbols[atom], xyz, surface)
        fig.add_trace(mesh)

    
    # #Add Volume
    # if kwargs["iso"] != None:
    #     print("Yas queen")

    #Layout
    fig.update_layout(get_layout(molecule.geometry))

#    fig.show()
    return fig

class Figure():

    def __init__(self, surface="matte", resolution=None, **kwargs):

        #self.molecule = molecule_factory(label, **kwargs)
        self.fig        = go.Figure()
        self.molecules  = []
        self.surface    = surface
        self.resolution = resolution
        self.max_range  = None
        self.min_range  = None
        #self.add_geometry(label, **kwargs)
        #self.add_bonds()
        #self.add_atoms()

    def show(self):
        self.fig.show()


    def add_bonds(self):
        add_bonds(self.molecule, self.fig, self.surface)
        
        # for bond in range(len(self.molecule.bonds)):
        #     mesh_a, mesh_b = get_bond_mesh(self.molecule.cilinders, self.molecule.bonds, bond, self.molecule.symbols, self.surface)
        #     self.fig.add_trace(mesh_a)
        #     self.fig.add_trace(mesh_b)

    def add_atoms(self):
        add_atoms(self.molecule, self.fig, self.surface)
        # for atom, xyz in enumerate(self.molecule.geometry):
        #     mesh = get_sphere_mesh(self.molecule.spheres[atom], self.molecule.symbols[atom], xyz, self.surface)
        #     self.fig.add_trace(mesh)

    def add_blob(self, index=0,
                       iso=0.001, 
                       color="Portland_r",
                       opacity=0.2):
        volume = get_blob(self.molecules[index], iso, opacity, color)
        self.fig.add_trace(volume)

    def add_molecule(self, molecule):
        #molecule = molecule_factory(label, **kwargs)
        add_bonds(molecule, self.fig, self.surface)
        add_atoms(molecule, self.fig, self.surface)
        self.molecules.append(molecule)
        self.fig.update_layout(get_layout(molecule.geometry, self.resolution))

    def add_layer(trace):
        self.fig.add_trace(trace)


def add_bonds(molecule, figure, surface):
    for bond in range(len(molecule.bonds)):
        mesh_a, mesh_b = get_bond_mesh(molecule.cilinders, molecule.bonds, bond, molecule.symbols, surface)
        figure.add_trace(mesh_a)
        figure.add_trace(mesh_b)

def add_atoms(molecule, figure, surface):
    for atom, xyz in enumerate(molecule.geometry):
        mesh = get_sphere_mesh(molecule.spheres[atom], molecule.symbols[atom], xyz, surface)
        figure.add_trace(mesh)











