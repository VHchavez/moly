import plotly.graph_objects as go

from ..figure.colors import *
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
    trace_list = []
    sphere = np.array(get_sphere())

    for atom, xyz in enumerate(geometry):
        if style is "ball_and_stick":
            reshaped_sphere = sphere * (atomic_numbers[atom]/30 + 0.6)
        elif style is "tubes":
            reshaped_sphere = sphere * 0.3
        elif style is "spacefilling":
            reshaped_sphere = sphere * (atomic_numbers[atom]/20 + 1.5)
        elif style is "wireframe":
            reshaped_sphere = sphere * 0.06
        else:
            raise ValueError("Only avaliable styles are \"ball_and_stick\", \"tubes\", \"spacefilling\" and \"wireframe\" ")
        mesh = get_sphere_mesh(reshaped_sphere,symbols[atom], xyz, surface)
        trace_list.append(mesh)

    return trace_list






