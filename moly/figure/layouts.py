import numpy as np


def get_axis(geometry, axis, max_range, min_range, overage=1.2):

    axis = {
        "autorange": False,
        "range": (overage * min_range, overage * max_range),
        "showgrid": False,
        "zeroline": False,
        "showline": False,
        "title": "",
        "ticks": '',
        "showticklabels": False,
        "showbackground": False,
        "showspikes": False
    }

    return axis


def get_layout(geometry, figsize, max_range, min_range):
    layout = {
        "scene_aspectmode": "manual",
        "scene_aspectratio": {
            "x": 1,
            "y": 1,
            "z": 1
        },
        "scene_xaxis_showticklabels": False,
        "scene_yaxis_showticklabels": False,
        "scene_zaxis_showticklabels": False,
        "dragmode": "orbit",
        "template": "plotly_white",
        "showlegend": False,
        "hovermode": False,
        "scene": {
            "xaxis": get_axis(geometry, 0, max_range[0], min_range[0]),
            "yaxis": get_axis(geometry, 1, max_range[1], min_range[1]),
            "zaxis": get_axis(geometry, 2, max_range[2], min_range[2])
        }
    }

    if figsize is not None:
        layout["height"] = figsize[0]
        layout["width"] = figsize[1]

    return layout


#Materials
surface_materials = {
    "matte": {
        "ambient": 0.60,
        "diffuse": 0.35,
        "fresnel": 0.05,
        "specular": 0.03,
        "roughness": 0.05,
        "facenormalsepsilon": 1e-15,
        "vertexnormalsepsilon": 1e-15
    },
    "shiny": {
        "ambient": 0.18,
        "diffuse": 0.85,
        "fresnel": 0.10,
        "specular": 0.70,
        "roughness": 0.05,
        "facenormalsepsilon": 1e-15,
        "vertexnormalsepsilon": 1e-15
    }
}
