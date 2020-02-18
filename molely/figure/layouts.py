
import numpy as np

def get_axis(geometry, axis):

    max_range = np.abs(geometry).max()

    axis = {            "autorange" : False,
                        "range" : (-1.6*max_range, 1.6*max_range),
                        "showgrid": False,
                        "zeroline": False,
                        "showline": False,
                        "title": "",
                        "ticks": '',
                        "showticklabels": False,
                        "showbackground": False,  
                        "showspikes" : False
    }

    return axis


def get_layout(geometry, resolution):
    layout = {
            "scene_aspectmode" : "manual",
            "scene_aspectratio" : {"x":1, "y":1, "z":1},
            "scene_xaxis_showticklabels"    :  False,
            "scene_yaxis_showticklabels"    :  False,
            "scene_zaxis_showticklabels"    :  False,
            "dragmode"                      :  "orbit",
            "template"                      :  "plotly_white",
            "showlegend":False,
            "hovermode":False,
                "scene":{
                    "xaxis": get_axis(geometry,0),
                    "yaxis": get_axis(geometry,1),
                    "zaxis": get_axis(geometry,2)
                        }}

    if not resolution == None:
        layout["height"] = resolution[0]
        layout["width"]  = resolution[1]

    return layout

                


#Materials
surface_materials = {"matte":{"ambient"            : 0.60,
                            "diffuse"              : 0.35,
                            "fresnel"              : 0.05,
                            "specular"             : 0.03,
                            "roughness"            : 0.05,
                            "facenormalsepsilon"   : 1e-15,
                            "vertexnormalsepsilon" : 1e-15},

                    "shiny":{"ambient"             : 0.18,
                            "diffuse"              : 0.85,
                            "fresnel"              : 0.10,
                            "specular"             : 0.70,
                            "roughness"            : 0.05,
                            "facenormalsepsilon"   : 1e-15,
                            "vertexnormalsepsilon" : 1e-15}}