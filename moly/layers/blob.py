import plotly.graph_objects as go
import numpy as np
import qcelemental as qcel
import glob
import os

from ..figure.layouts import surface_materials

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



def get_cubes(folder):
    cube_list = [f for f in glob.glob(folder+"/*.cube")]
    cubes = []
    meta = []
    for cube_file in cube_list:
        cube_np, details = cube_to_array(cube_file)
        details.update({"name":cube_file[2:-5]})
        cubes.append(cube_np)
        meta.append(details)
        
    return cubes, meta


def get_cubes_traces(cubes, spacing, origin, iso=0.03):
    x,y,z = np.mgrid[:cubes[0].shape[0], :cubes[0].shape[1], :cubes[0].shape[2]]

    x_r = x * spacing[0] + origin[0]
    y_r = y * spacing[1] + origin[1]
    z_r = z * spacing[2] + origin[2]


    traces = []
    for i, cube in enumerate(cubes):
        value = cube.flatten()
        trace = go.Isosurface(x = x_r.flatten(),
                              y = y_r.flatten(), 
                              z = z_r.flatten(), 
                              value = cube.flatten(),
                              surface_count = 2,
                              colorscale = "Portland",
                              visible = True if i == 0 else False,
                              showscale = False,
                              isomin= -1 * iso,
                              isomax= 1 * iso, 
                              flatshading = True,
                              lighting = surface_materials["matte"], 
                              caps=dict(x_show=False, y_show=False, z_show=False),
                              opacity=0.4)
        
        traces.append(trace)
        
    return traces


def get_buttons(meta, geo_traces):
    buttons =  []
    for cube_i in meta:
        button = dict(label=cube_i["name"],
                         method="update",
                         args=[{"visible": [True for traces in range(geo_traces)] + [True if cube_i['name'] == cube_j['name'] else False for cube_j in meta]},
                               {"title": "",
                                "annotations": []}])
        buttons.append(button)
    return buttons


def cube_to_molecule(cube_file):

    _ , meta = cube_to_array(cube_file)
    origin = meta["origin"]
    atoms = meta["geometry"]
    spacing = [meta["xvec"][0], meta["yvec"][1], meta["zvec"][2]]

    geometry = []
    for atom in atoms:
        geometry.append(atom[1][1:])
    geometry = np.array(geometry)

    symbols = []
    atomic_numbers = []
    for atom in atoms:
        symbols.append(qcel.periodictable.to_symbol(atom[0]))
        atomic_numbers.append(atom[0])
    symbols = np.array(symbols)
    atomic_numbers = np.array(atomic_numbers)
    
    return geometry, symbols, atomic_numbers, spacing, origin


def cube_to_array(file):
    """
    Read cube file into numpy array
    Parameters
    ----------
    fname: filename of cube file
    Returns
    --------
    (data: np.array, metadata: dict)
    """
    cube_details = {}
    with open(file, 'r') as cube:
        cube.readline()
        cube.readline()  # ignore comments
        natm, cube_details['origin'] = _getline(cube)
        nx, cube_details['xvec'] = _getline(cube)
        ny, cube_details['yvec'] = _getline(cube)
        nz, cube_details['zvec'] = _getline(cube)
        cube_details['geometry'] = [_getline(cube) for i in range(natm)]
        data = np.zeros((nx * ny * nz))
        idx = 0
        for line in cube:
            for val in line.strip().split():
                data[idx] = float(val)
                idx += 1
    data = np.reshape(data, (nx, ny, nz))
    cube.close()

    return data, cube_details

def _getline(cube):
    """
    Read a line from cube file where first field is an int
    and the remaining fields are floats.
    Parameters
    ----------
    cube: file object of the cube file
    Returns
    -------
    (int, list<float>)
    """
    l = cube.readline().strip().split()
    return int(l[0]), list(map(float, l[1:]))


