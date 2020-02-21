import sys
sys.path.append("..")

import plotly.graph_objects as go
from ..figure.colors import *
from ..figure.layouts import *

def get_bond_mesh(cilinder,bond,symbols, surface):

	lighting = surface_materials[surface]

	mesh = go.Mesh3d({  'x':cilinder[:,0].flatten(),
						'y':cilinder[:,1].flatten(),
						'z':cilinder[:,2].flatten(),
						'color': colors[symbols[bond][0]][0],
						'alphahull' :0,
						'flatshading' : True,
						"cmin"     :-7,
						'lighting' : lighting,
						'lightposition' : {"x":100,
										   "y":200,
										   "z":0}})

	return mesh
