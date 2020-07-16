import numpy as np
import plotly.graph_objects as go
from ..molecule.shapes import rotation_matrix
from ..molecule.shapes import get_single_cylinder
from ..figure.colors import colors
from ..figure.layouts import surface_materials


def get_bond_mesh(cilinder, bond, symbols, surface):

    lighting = surface_materials[surface]
    mesh = go.Mesh3d({
        'x': cilinder[:, 0].flatten(),
        'y': cilinder[:, 1].flatten(),
        'z': cilinder[:, 2].flatten(),
        'color': colors[symbols[bond]][0],
        'alphahull': 0,
        'flatshading': True,
        "cmin": -7,
        'lighting': lighting,
        'lightposition': {
            "x": 100,
            "y": 200,
            "z": 0
        }
    })

    return mesh


def get_bonds(geometry, symbols, bonds, style, surface):

    trace_list = []
    for idx1, idx2 in bonds:

        vec1 = geometry[idx1]
        vec2 = geometry[idx2]
        length = np.linalg.norm(vec2 - vec1)
        R = rotation_matrix(np.array([0, 0, 1]), vec2 - vec1)

        if style is "ball_and_stick" or style == "tubes":
            r = 0.3
        elif style is "wireframe":
            r = 0.06
        elif style is "spacefilling":
            return []
        else:
            raise ValueError(
                "Only avaliable styles are \"ball_and_stick\", \"tubes\", \"spacefilling\" and \"wireframe\" ")

        if symbols[idx1] == symbols[idx2]:

            cyl = get_single_cylinder(radius=r)
            cyl[:, 2] *= length
            cyl = R.dot(cyl.T).T
            cyl += vec1

            mesh = get_bond_mesh(cyl, idx1, symbols, surface)
            trace_list.append(mesh)

        if symbols[idx1] != symbols[idx2]:

            cyl = get_single_cylinder(radius=r)
            cyl[:, 2] *= length / 2
            cyl = R.dot(cyl.T).T
            cyl_1 = cyl + vec1
            cyl_2 = cyl + (vec1 + vec2) / 2

            mesh = get_bond_mesh(cyl_1, idx1, symbols, surface)
            trace_list.append(mesh)
            mesh = get_bond_mesh(cyl_2, idx2, symbols, surface)
            trace_list.append(mesh)

    return trace_list
