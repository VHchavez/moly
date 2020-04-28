"""

Generates layouts and surface materials

"""

def get_range(min_range, max_range, overage=2.5):

    axis = {
        "autorange": False,
        "range": (min_range * overage, max_range * overage),
    }

    layout = {"scene" : {"xaxis":axis, 
                         "yaxis":axis,
                         "zaxis":axis}, 
                         
             "scene_camera" : {"up":{"x":0, "y":0, "z":1},
                               "center":{"x":0, "y":0, "z":0},
                               "eye":{"x":2.2/max_range, "y":2.2/max_range, "z":2.2/max_range}}}

    return layout


def get_layout(figsize=None):

    axis = {
        "showgrid": False,
        "zeroline": False,
        "showline": False,
        "title": "",
        "ticks": '',
        "showticklabels": False,
        "showbackground": False,
        "showspikes": False
    }

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
        "scene" : {"xaxis": axis ,
                   "yaxis": axis, 
                   "zaxis": axis,  
                    },
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
        "roughness": 0.05,
        "facenormalsepsilon": 1e-15,
        "vertexnormalsepsilon": 1e-15
    },
     "orbs": {
        "ambient": 0.3,
        "diffuse": 0.6,
        "fresnel": 0.10,
        "specular": 0.70,
        "roughness": 0.9,
        "facenormalsepsilon": 1e-15,
        "vertexnormalsepsilon": 1e-15
    },
    "glass" : { 
        "ambient": 0.3,
        "diffuse": 0.1,
        "fresnel": 0.45,
        "specular": 0.70,
        "roughness": 0.1,
        "facenormalsepsilon": 1e-15,
        "vertexnormalsepsilon": 1e-15}
}
