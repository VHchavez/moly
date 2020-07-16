import numpy as np
import plotly.graph_objects as go

from ..figure.colors import colors
from ..figure.layouts import *
from ..molecule.shapes import get_sphere


def get_sphere_mesh(sphere, sym, xyz, surface):

    lightning = surface_materials[surface]

    mesh = go.Mesh3d({
            'x': sphere[0] + xyz[0]  , 
            'y': sphere[1] + xyz[1]  , 
            'z': sphere[2] + xyz[2]  , 
            'alphahull': 0,
            'color'    : colors[sym][0],
            'flatshading' : False,
            "cmin"     :-7,# atrick to get a nice plot (z.min()=-3.31909)
            "lighting" : lightning,
            "lightposition" : {"x":100,
                               "y":200,
                               "z":0}     
    })

    return mesh


def get_atoms(geometry, atomic_numbers, symbols, style, surface):
    """
    Obtain sphere mesh for each atom based on the style

    Parameters
    ----------
    geometry: numpy arrays
        Geometry of structure
    atomic_numbers: List
        List of atomic numbers
    symbols: List
        List of atomic symbols
    style: str
        Style of structure: {ball_and_stick, tubes, spacefiling, wireframe}
    surface: dict
        Specification of style for surface

    Returns
    -------
    trace_list = List
        List with plotly graph objects of atomic species
    """
    trace_list = []
    sphere = np.array(get_sphere())

    for atom, xyz in enumerate(geometry):
        if style == "ball_and_stick":
            reshaped_sphere = sphere * (atomic_numbers[atom]/30 + 0.6)
        elif style == "tubes":
            reshaped_sphere = sphere * 0.3
        elif style is "spacefilling":
            #Uses van der waals radii
            reshaped_sphere = sphere * (colors[symbols[atom]][2])
        elif style is "wireframe":
            reshaped_sphere = sphere * 0.06
        else:
            raise ValueError("Only avaliable styles are \"ball_and_stick\", \"tubes\", \"spacefilling\" and \"wireframe\" ")
        mesh = get_sphere_mesh(reshaped_sphere,symbols[atom], xyz, surface)
        trace_list.append(mesh)

    return trace_list






