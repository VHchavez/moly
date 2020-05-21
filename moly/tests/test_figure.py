import pytest
import moly

@pytest.fixture()
def he_dimer():
    he2 = """
    0 1
    He 0.0 0.0 0.0
    He 0.0 0.0 6.0
    units bohr
    symmetry c1
    """
    return he2

def test_init_figure():
    fig = moly.Figure()

def test_add_molecule_from_data(he_dimer):
    fig = moly.Figure()
    mol = moly.Molecule.from_data(he_dimer)
    fig.add_molecule("he2", mol)

def test_add_molecule_from_data_wireframe(he_dimer):
    fig = moly.Figure()
    mol = moly.Molecule.from_data(he_dimer)
    fig.add_molecule("he2", mol, style="wireframe")

def test_add_molecule_from_data_spacefilling(he_dimer):
    fig = moly.Figure()
    mol = moly.Molecule.from_data(he_dimer)
    fig.add_molecule("he2", mol,  style="spacefilling")

def test_add_molecule_from_data_tubes(he_dimer):
    fig = moly.Figure()
    mol = moly.Molecule.from_data(he_dimer)
    fig.add_molecule("he2", mol,  style="tubes")

# def test_add_molecule_from_file():
#     fig = moly.Figure()
#     mol = moly.Molecule.from_file("water.xyz")

def test_figsize():
    fig = moly.Figure(figsize=(300,300))

# def test_add_cube():
#     fig = moly.Figure()
#     fig.add_cube("Da.cube", iso=0.03, colorscale="rdbu", opacity=0.2)

# def test_add_cube_slider():
#     fig = moly.Figure()
#     fig.add_cube("/mnt/c/Users/victo/Dropbox/PHD/moly/moly/moly/tests/Da.cube", iso=[0.01, 0.1], colorscale="rdbu", opacity=0.2)

# def test_add_cubes():
#     fig = moly.Figure()
#     fig.add_cubes(directory='test_files', colorscale="portland", iso=0.1)
#     fig.show()




