"""

Creates main figure 

"""
import numpy as np
import qcelemental as qcel
import plotly.graph_objects as go

from ..layers.bonds import get_bonds
from ..layers.geometry import get_atoms
from ..layers.measurements import get_angle, get_line
from ..layers.cube import get_cubes, cube_to_molecule, get_cubes, get_cube_trace
from .layouts import get_layout, get_range
from .widgets import get_buttons, get_buttons_wfn

from ..advanced import cubeprop


class Figure():
    def __init__(self, surface="matte", figsize=None, **kwargs):

        self.fig = go.Figure()
        self.molecules = {}
        self.geometries = []
        self.surface = surface
        self.resolution = figsize
        self.min_range = 0.0
        self.max_range = 0.0

    def show(self):
        self.fig.show()

    def assert_range(self, geometry):
        self.min_range = self.min_range if self.min_range < np.min(geometry) else np.min(geometry)
        self.max_range = self.max_range if self.max_range > np.max(geometry) else np.max(geometry)

        range_layout = get_range(self.min_range, self.max_range)
        self.fig.update_layout(range_layout)

    def get_connectivity(self, molecule):

        mol_dict = molecule.dict()
        symbols = molecule.symbols
        geometry = molecule.geometry
        
        if not "connectivity" in mol_dict:
            return qcel.molutil.guess_connectivity(symbols,geometry)
        
        elif mol_dict["connectivity"] is None:
            return qcel.molutil.guess_connectivity(symbols,geometry)

        elif "connectivity" in mol_dict:
            return self.dict["connectivity"]

    #Basic Traces

    def add_molecule(self, name, molecule, style="ball_and_stick"):

        self.molecules[name] = molecule

        bonds = self.get_connectivity(molecule)
        bond_list = get_bonds(molecule.geometry, molecule.symbols, bonds, style, self.surface)
        atom_list = get_atoms(molecule.geometry, molecule.atomic_numbers, molecule.symbols, style, self.surface)

        #Add traces
        for bond in bond_list:
            self.fig.add_trace(bond)

        for atom in atom_list:
            self.fig.add_trace(atom)

        #Update layout
        self.fig.update_layout(get_layout(self.resolution))
        self.assert_range(molecule.geometry)

    def add_cube(self, file, iso=0.01, plot_geometry=True, 
                 colorscale="portland", opacity=0.2, style="ball_and_stick", 
                 iso_steps=3):

        geometry, symbols, atomic_numbers, spacing, origin, cube = cube_to_molecule(file)
    
        if plot_geometry is True:
            bonds = qcel.molutil.guess_connectivity(symbols, geometry)
            bond_list = get_bonds(geometry, symbols, bonds, style, self.surface)
            atom_list = get_atoms(geometry, atomic_numbers, symbols, style, self.surface)

            #Add traces
            for bond in bond_list:
                self.fig.add_trace(bond)

            for atom in atom_list:
                self.fig.add_trace(atom)

        #Single value of iso
        if type(iso) == float: 
            trace, min_range, max_range = get_cube_trace(cube, spacing, origin, iso, colorscale, opacity) 
            #Add traces
            self.fig.add_trace(trace)

        #Update layout
        self.fig.update_layout(get_layout(self.resolution))
        self.assert_range([min_range, max_range])

    def add_measurement(self, mol_label, m, 
                        line_width=20,
                        line_color='grey'):

        molecule = self.molecules[mol_label]
        measurement = molecule.measure(m)

        if len(m) == 2:
            line = get_line(str(round(measurement, 2)) ,molecule.geometry[m[0]], molecule.geometry[m[1]], 
                    line_width,
                    line_color)

            self.fig.add_trace(line)

            # markers = go.Scatter3d(x=[molecule.geometry[m[0]][0],molecule.geometry[m[0]][0]],
            #                        y=[molecule.geometry[m[0]][1],molecule.geometry[m[0]][1]],  
            #                        z=[molecule.geometry[m[0]][2],molecule.geometry[m[0]][2]])

            # self.fig.add_trace(markers)

            print(f"Distance between {molecule.symbols[m[0]]} and {molecule.symbols[m[1]]} is {measurement}")

        elif len(m) == 3:
            angle = get_angle(str(round(measurement,2)), 
                    molecule.geometry[m[0]], molecule.geometry[m[1]], molecule.geometry[m[2]], 
                    line_color)

            self.fig.add_trace(angle)

        elif len(m) == 4:
            print("Unable to add dihedral")

    def add_cubes(self, directory=".", iso=0.03, style="ball_and_stick", colorscale="portland", opacity=0.3):
        cubes, details = get_cubes(directory)
        geometry, symbols, atomic_numbers, spacing, origin, _ = cube_to_molecule(details[0]["name"]+".cube")
        bonds = qcel.molutil.guess_connectivity(symbols, geometry)


        bond_list = get_bonds(geometry, symbols, bonds, style, self.surface)
        atom_list = get_atoms(geometry, atomic_numbers, symbols, style, self.surface)

        for bond in bond_list:
            self.fig.add_trace(bond)
        for atom in atom_list:
            self.fig.add_trace(atom)

        geometry_traces = len(self.fig.data)

        cube_list = []
        min_list = []
        max_list = []
        for cube in cubes:
            trace, min_range, max_range = get_cube_trace(cube, spacing, origin, iso, colorscale,opacity, visible=False)
            cube_list.append(trace)
            min_list.append(min_range)
            max_list.append(max_range)
        
        for cube in cube_list:
            self.fig.add_traces(cube)

        
        button_list = get_buttons(details, geometry_traces, directory)

        self.fig.update_layout(updatemenus=[dict(showactive=True,
                                                buttons=button_list,
                                                font={"family": "Helvetica",
                                                      "size" : 18},
                                                borderwidth=0
            ),
        ])

        #Update layout
        self.fig.update_layout(get_layout(self.resolution))
        self.assert_range([min(min_list), max(max_list)])


    #Psi4 Traces

    def add_density(self, name, wfn, iso=0.03, colorscale="portland", opacity=0.3, geometry=True, 
                    spacing=[0.2, 0.2, 0.2], overage=[4.0, 4.0, 4.0]):

        molecule = qcel.models.Molecule.from_data(wfn.basisset().molecule().save_string_xyz())
        self.add_molecule("name"+"_geometry", molecule)

        D = spacing
        L = overage

        O, N =  cubeprop.build_grid(wfn, L, D) 
        block, points, nxyz, npoints =  cubeprop.populate_grid(wfn, O, N, D)
        volume = cubeprop.compute_density(O, N, D, npoints, points, nxyz, block, wfn.Da())
        trace, min_range, max_range = get_cube_trace(volume, spacing, O, iso, colorscale, opacity)

        self.fig.add_trace(trace)
        self.fig.update_layout(get_layout(self.resolution))
        self.assert_range([min_range, max_range])

   








