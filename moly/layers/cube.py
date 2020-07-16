import plotly.graph_objects as go
import numpy as np
import qcelemental as qcel
import glob


from ..figure.layouts import surface_materials

def get_volume(cube, spacing, origin, iso, opacity, color):

    x, y, z = np.mgrid[:cube.shape[0], :cube.shape[1], :cube.shape[2]]
    x_r = x * spacing[0] + origin[0]
    y_r = y * spacing[1] + origin[1]
    z_r = z * spacing[2] + origin[2]

    mesh = go.Isosurface(x = x_r.flatten(),
                         y = y_r.flatten(), 
                         z = z_r.flatten(), 
                        value = molecule.cube.flatten(),
                        surface_count = 2,
                        colorscale = color,
                        showscale=False,
                        isomin=-1 * iso,
                        isomax= 1 * iso,
                        opacity = opacity)


    max_range_x = np.max(x_r)
    max_range_y = np.max(y_r)
    max_range_z = np.max(z_r)
    min_range_x = np.min(x_r)
    min_range_y = np.min(y_r)
    min_range_z = np.min(z_r)
    max_range = max(max_range_x, max_range_y, max_range_z)
    min_range = min(min_range_x, min_range_y, min_range_z)

    return mesh, min_range, max_range

def get_surface(grid, spacing, origin):

    x = np.array(grid[0])
    y = np.array(grid[1])
    z = np.array(grid[2])


    x = x * spacing + origin[0]
    y = y * spacing + origin[1]
    z = z * spacing + origin[2]

    mesh = go.Mesh3d({
            'x': x, 
            'y': y, 
            'z': z, 
            'alphahull': 0,
            'color'    : 'turquoise',
            'opacity' : 0.20,
            'visible' : False,
            'flatshading' : False,
            "lighting" : surface_materials["glass"],
            "lightposition" : {"x":100,
                                "y":200,
                                    "z":0}
    })
    return mesh

def get_cubes(folder):
    cube_list = [f for f in glob.glob(folder+"/*.cube")]
    if not cube_list:
        raise ValueError("Directory does not contain cube files") 
    cubes = []
    meta = []
    for cube_file in cube_list:
        cube_np, details = cube_to_array(cube_file)
        details.update({"name":cube_file[:-5]})
        cubes.append(cube_np)
        meta.append(details)
        
    return cubes, meta

def get_cubes_surfaces(cubes, spacing, origin, iso, colorscale, opacity, alpha, atol):
    """
    Warning: Surfaces are only functional for self-contained blobs. Multiple blobs that 
    Exists in space will be merged into a massive undesired blob. 
    """

    cubes_surfaces = []
    traces = []

    if iso > 0.0:
        color = 'mediumturquoise'
    else:
        color = 'red'

    for cube in cubes:
        it = np.nditer(cube, flags=['multi_index'])
        x, y, z = [], [], []
        while not it.finished:
            if np.isclose(it[0],iso,atol=atol):
                x.append(it.multi_index[0])
                y.append(it.multi_index[1])
                z.append(it.multi_index[2])
            it.iternext()
        cubes_surfaces.append([x,y,z])
    
    for surface in cubes_surfaces:

        x = np.array(surface[0])
        y = np.array(surface[1])
        z = np.array(surface[2])
        x = x * spacing[0] + origin[0]
        y = y * spacing[1] + origin[1]
        z = z * spacing[2] + origin[2]

        mesh = go.Mesh3d({
                'x': x, 
                'y': y, 
                'z': z, 
                'alphahull': alpha,
                'color'    : color,
                'opacity' : opacity,
                'visible' : True,
                'flatshading' : False,
                "lighting" : surface_materials["glass"],
                "lightposition" : {"x":100,
                                    "y":200,
                                    "z":0}
        })

        traces.append(mesh)

    max_range_x = np.max(x)
    max_range_y = np.max(y)
    max_range_z = np.max(z)
    min_range_x = np.min(x)
    min_range_y = np.min(y)
    min_range_z = np.min(z)
    max_range = max(max_range_x, max_range_y, max_range_z)
    min_range = min(min_range_x, min_range_y, min_range_z)

    return traces[0], min_range, max_range

def get_cube_trace(cube, spacing, origin, iso, colorscale, opacity, visible=True):
    x,y,z = np.mgrid[:cube.shape[0], :cube.shape[1], :cube.shape[2]]

    x_r = x * spacing[0] + origin[0]
    y_r = y * spacing[1] + origin[1]
    z_r = z * spacing[2] + origin[2]

    #value = cube.flatten()
    trace = go.Isosurface(  x = x_r.flatten(),
                            y = y_r.flatten(), 
                            z = z_r.flatten(), 
                            value = cube.flatten(),
                            surface_count = 2,
                            colorscale = colorscale,
                            visible = visible,
                            showscale = False,
                            isomin= -1 * iso,
                            isomax= 1 * iso, 
                            flatshading = False,
                            lighting = surface_materials["matte"], 
                            caps=dict(x_show=False, y_show=False, z_show=False),
                            opacity=opacity)
        
    max_range_x = np.max(x_r)
    max_range_y = np.max(y_r)
    max_range_z = np.max(z_r)
    min_range_x = np.min(x_r)
    min_range_y = np.min(y_r)
    min_range_z = np.min(z_r)
    max_range = max(max_range_x, max_range_y, max_range_z)
    min_range = min(min_range_x, min_range_y, min_range_z)

    return trace, min_range, max_range

def cube_to_molecule(cube_file):

    cube , meta = cube_to_array(cube_file)
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
    
    return geometry, symbols, atomic_numbers, spacing, origin, cube

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


