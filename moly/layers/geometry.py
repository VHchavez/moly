import plotly.graph_objects as go

from ..figure.colors import *
from ..figure.layouts import *


def get_sphere_mesh(sphere,sym, xyz, surface):

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
