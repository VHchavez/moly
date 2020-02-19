import sys
sys.path.append("..")

import plotly.graph_objects as go
from ..figure.colors import *
from ..figure.layouts import *

ambient   = 0.65
diffuse   = 0.35
specular  = 0

def get_bond_mesh(cilinders,bonds,bond,symbols, surface):

	lighting = surface_materials[surface]

	mesh_i = go.Mesh3d({'x':cilinders[0][bond][:,0].flatten(),
	'y':cilinders[0][bond][:,1].flatten(),
	'z':cilinders[0][bond][:,2].flatten(),
	'color': colors[symbols[bonds[bond][0]]][0],
	'alphahull' :0,
	'flatshading' : True,
	"cmin"     :-7,
	'lighting' : lighting,
	'lightposition' : {"x":100,
	"y":200,
	"z":0}})

	mesh_j = go.Mesh3d({'x':cilinders[1][bond][:,0].flatten(),
	'y':cilinders[1][bond][:,1].flatten(),
	'z':cilinders[1][bond][:,2].flatten(),
	'color': colors[symbols[bonds[bond][1]]][0],
	'alphahull' :0,
	'flatshading' : True,
	'cmin'     :-7,
	'lighting' : lighting,
	'lightposition' : {"x":100,
	"y":200,
	"z":0}
	})

	return mesh_i, mesh_j
