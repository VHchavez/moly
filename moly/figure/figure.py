"""

Creates main figure 

"""
import numpy as np
import qcelemental as qcel
import psi4
import plotly.graph_objects as go

from ..layers.bonds import get_bond_mesh, get_bonds
from ..layers.geometry import get_sphere_mesh, get_atoms
from ..layers.measurements import get_angle, get_line
from ..layers.cube import get_cubes, cube_to_molecule, get_cubes, get_cube_trace, get_surface, get_cubes_surfaces, cube_to_array

from .layouts import get_layout, get_range
from .widgets import get_slider, get_buttons, get_buttons_wfn

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

        #Slider with multiple iso values
        elif type(iso) == list:
            traces = []
            geometry_traces = len(self.fig.data)
            iso_range = np.linspace(iso[0], iso[1], iso_steps)
            for i, iso_i in enumerate(iso_range):
                if i == 0:
                    trace, min_range, max_range = get_cube_trace(cube, spacing, origin, iso_i, colorscale, opacity, visible=True)
                elif i != 0:
                    trace, min_range, max_range = get_cube_trace(cube, spacing, origin, iso_i, colorscale, opacity, visible=False)
                self.fig.add_trace(trace)
            slider = get_slider(iso_range, iso_steps, geometry_traces)
            self.fig.update_layout(sliders=slider)

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

    def add_orbitals(self, wfn, orbitals, iso, quality=0.2, overage=2.0, style="ball_and_stick"):
        import psi4

        if max(orbitals) > wfn.basisset().nbf():
            raise ValueError("The maximum MO is greater than the avaliable number of basis functions")


        molecule = qcel.models.Molecule.from_data(wfn.molecule().save_string_xyz())
        bonds = self.get_connectivity(molecule)
        add_bonds(molecule.geometry, molecule.symbols, bonds, self.fig, style, self.surface)
        add_atoms(molecule.geometry, molecule.atomic_numbers, molecule.symbols, self.fig, style, self.surface)


        L = [4.0, 4.0, 4.0]
        D = [quality, quality, quality]

        O, N =  cubeprop.build_grid(wfn, L, D) 
        block, points, nxyz, npoints =  cubeprop.populate_grid(wfn, O, N, D)
        trace_list_a = []
        trace_list_b = []

        # for orb_i in orbitals:

        #     C = wfn.Ca()
        #     C_req = C.np[:,orbitals[0]-1:orbitals[-1]]
        #     C_req = psi4.core.Matrix.from_array(C_req)

        #     print(C_req)
        #     #matrix = np.einsum('p,q->pq', wfn.Ca().np[:,orb_i], wfn.Ca().np[:,orb_i])
        #     orb_i_np = cubeprop.compute_orbitals(O, N, D, npoints, points, nxyz, block, C_req, orbitals, iso)

        # return orb_i_np

        #for orb_i in orbitals:
        #    orbital_a = np.zeros_like(wfn.Ca().np)
        #    orbital_a[:,orb_i] = wfn.Ca().np[:,orb_i]    
        #    orbital_b = np.zeros_like(wfn.Cb().np)
        #    orbital_b[:,orb_i] = wfn.Cb().np[:,orb_i]
        #    #orbital_a = np.einsum('p,q->pq', wfn.Ca().np[:,orb_i], wfn.Ca().np[:,orb_i])
        #    #orbital_b = np.einsum('p,q->pq', wfn.Cb().np[:,orb_i], wfn.Cb().np[:,orb_i])
        #    orbital_a = psi4.core.Matrix.from_array(orbital_a)
        #    orbital_b = psi4.core.Matrix.from_array(orbital_b)
        #    orbital_mesh_a = cubeprop.compute_density(O, N, D, npoints, points, nxyz, block, orbital_a, iso)
        #    orbital_mesh_b = cubeprop.compute_density(O, N, D, npoints, points, nxyz, block, orbital_b, iso)
        #    trace_a = get_surface(orbital_mesh_a, quality, O)
        #    trace_b = get_surface(orbital_mesh_b, quality, O)
        #    trace_list_a.append(trace_a)
        #    trace_list_b.append(trace_b)

        if len(orbitals) == 1:
           self.fig.add_trace(trace_list_a[0])   
           self.fig.add_trace(trace_list_b[0])  

        if len(orbitals) > 1:
            geometry_traces = len(self.fig.data)

            orbital_list = [str(orb) for orb in orbitals] + [str(orb) for orb in orbitals]
            trace_list = trace_list_a + trace_list_b

            button_list = get_buttons_wfn(orbital_list, geometry_traces)
            self.fig.update_layout(updatemenus=[dict(showactive=True,
                                                     buttons=button_list,
                                                     font={"family": "Helvetica", "size" : 18},
                                                     borderwidth=0
            ),
            ])

            for trace in trace_list:
                self.fig.add_trace(trace)

        self.assert_range(molecule.geometry)
        self.fig.update_layout(get_layout(molecule.geometry, self.resolution, self.min_range * overage, self.max_range * overage))

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
        trace, min_range, max_range = get_cubes_traces(volume, spacing, O, iso, colorscale, opacity)

        self.fig.add_trace(trace)
        self.fig.update_layout(get_layout(self.resolution))
        self.assert_range([min_range, max_range])

   








