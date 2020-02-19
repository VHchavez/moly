import plotly.graph_objects as go
import numpy as np

def get_blob(molecule, iso, opacity, color):

    x, y, z = np.mgrid[:molecule.blob.shape[0], :molecule.blob.shape[1], :molecule.blob.shape[2]]
    x_r = x * molecule.spacing[0] + molecule.origin[0]
    y_r = y * molecule.spacing[1] + molecule.origin[1]
    z_r = z * molecule.spacing[2] + molecule.origin[2]


    mesh = go.Isosurface(x = x_r.flatten(),
                         y = y_r.flatten(), 
                         z = z_r.flatten(), 
                        value = molecule.blob.flatten(),
                        surface_count = 2,
                        colorscale = color,
                        showscale=False,
                        isomin=-1 * iso,
                        isomax= 1 * iso,
                        opacity = opacity)

    return mesh

