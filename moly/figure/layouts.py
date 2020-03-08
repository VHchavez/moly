import numpy as np

def get_axis(min_range, max_range):

    axis = {
        "autorange": False,
        "range": (min_range, max_range),
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


def get_layout(geometry, figsize, min_range, max_range, overage=1.5):
    axis = get_axis(min_range * overage, max_range * overage)
    layout = {
        "scene_aspectmode": "manual",
        "scene_aspectratio": {
            "x": 1,
            "y": 1,
            "z": 1
        },
        "scene_xaxis_showticklabels": True,
        "scene_yaxis_showticklabels": True,
        "scene_zaxis_showticklabels": True,
        "dragmode": "orbit",
        "template": "plotly_white",
        "showlegend": False,
        "hovermode": False,
        "scene": {
            "xaxis": axis,
            "yaxis": axis,
            "zaxis": axis,
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
        "ambient": 0.3,
        "diffuse": 0.85,
        "fresnel": 0.10,
        "specular": 0.70,
        "roughness": 0.05   ,
        "facenormalsepsilon": 1e-15,
        "vertexnormalsepsilon": 1e-15
    }
}
