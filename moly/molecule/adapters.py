
from .molecule_adapter import *
from .molecule_factory import *
from .shapes import *

import qcelemental as qcel
import numpy as np

class QCMolecule(MoleculeAdapter):

    def __init__(self, molecule):
        self.molecule = molecule
        self.dict     = molecule.dict() 
        self.symbols  = self.get_symbols()
        self.atomic_numbers = self.dict["atomic_numbers"]
        self.geometry = self.get_geometry()
        self.bonds    = self.get_connectivity()
        self.cilinders = get_bonds_cilinders(self.geometry, self.bonds)
        self.spheres   = get_atoms_spheres(self.geometry, self.atomic_numbers)

    def get_geometry(self):
        return self.molecule.geometry
        
    def get_symbols(self):
        return self.dict["symbols"]

    def get_connectivity(self):

        if not "connectivity" in self.dict:
            return qcel.molutil.guess_connectivity(self.symbols, self.geometry)
        
        elif self.dict["connectivity"] is None:
            return  qcel.molutil.guess_connectivity(self.symbols, self.geometry)

        elif "connectivity" in self.dict:
            return self.dict["connectivity"]

class Psi4Molecule(MoleculeAdapter):

    def __init__(self, molecule):
        self.molecule = qcel.models.Molecule.from_data(molecule.save_string_xyz())
        self.dict     = self.molecule.dict()
        self.symbols  = self.get_symbols()
        self.atomic_numbers = self.dict["atomic_numbers"]
        self.geometry = self.get_geometry()
        self.bonds    = self.get_connectivity()
        self.cilinders = get_bonds_cilinders(self.geometry, self.bonds)
        self.spheres   = get_atoms_spheres(self.geometry, self.atomic_numbers)


    def get_geometry(self):
        return self.molecule.geometry

    def get_symbols(self):
        return self.dict["symbols"]

    def get_connectivity(self):
        return qcel.molutil.guess_connectivity(self.symbols, self.geometry)


class xyzMolecule(MoleculeAdapter):

    def __init__(self, file):  
        self.molecule = qcel.models.Molecule.from_file(file)
        self.dict     = self.molecule.dict()
        self.geometry = self.get_geometry()
        self.symbols  = self.get_symbols()
        self.bonds    = self.get_connectivity()
        self.atomic_numbers = self.dict["atomic_numbers"]
        self.cilinders = get_bonds_cilinders(self.geometry, self.bonds)
        self.spheres   = get_atoms_spheres(self.geometry, self.atomic_numbers)
        
    def get_geometry(self):
        return self.molecule.geometry

    def get_symbols(self):
        return self.dict["symbols"]

    def get_connectivity(self):
        bonds = qcel.molutil.guess_connectivity(self.symbols, self.geometry)
        return bonds

class cubeMolecule(MoleculeAdapter):

    def __init__(self, file):
        self.file = file
        self.blob = None
        self.details = None
        self.cube_to_array()
        self.origin = self.details["origin"]
        self.spacing = [self.details['xvec'][0], self.details['yvec'][1], self.details['zvec'][2]]
        self.geometry = self.get_geometry()
        self.symbols, self.atomic_numbers = self.get_symbols()
        self.bonds = self.get_connectivity()
        self.cilinders = get_bonds_cilinders(self.geometry, self.bonds)
        self.spheres   = get_atoms_spheres(self.geometry, self.atomic_numbers, cube=False)

    def get_geometry(self):
        geometry = []
        for atom in self.details['atoms']:
            # atom[1][1] = (atom[1][1] - self.origin[0]) / self.spacing[0]
            # atom[1][2] = (atom[1][2] - self.origin[1]) / self.spacing[1]
            # atom[1][3] = (atom[1][3] - self.origin[2]) / self.spacing[2]
#            atom[1][1] = (atom[1][1]) / self.spacing[0]
#            atom[1][2] = (atom[1][2]) / self.spacing[1]
#            atom[1][3] = (atom[1][3]) / self.spacing[2]
            geometry.append(atom[1][1:])
        geometry = np.array(geometry)
        return geometry

    def get_symbols(self):
        symbols = []
        atomic_numbers = []
        for atom in self.details['atoms']:
            symbols.append(qcel.periodictable.to_symbol(atom[0]))
            atomic_numbers.append(atom[0])
        symbols = np.array(symbols)
        atomic_numbers = np.array(atomic_numbers)

        return symbols, atomic_numbers

        return

    def get_connectivity(self):
        return qcel.molutil.guess_connectivity(self.symbols, self.geometry)

    def cube_to_array(self):
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
        with open(self.file, 'r') as cube:
            cube.readline()
            cube.readline()  # ignore comments
            natm, cube_details['origin'] = self._getline(cube)
            nx, cube_details['xvec'] = self._getline(cube)
            ny, cube_details['yvec'] = self._getline(cube)
            nz, cube_details['zvec'] = self._getline(cube)
            cube_details['atoms'] = [self._getline(cube) for i in range(natm)]
            data = np.zeros((nx * ny * nz))
            idx = 0
            for line in cube:
                for val in line.strip().split():
                    data[idx] = float(val)
                    idx += 1
        data = np.reshape(data, (nx, ny, nz))
        cube.close()

        self.blob = data
        self.details = cube_details

        return data, cube_details

    def _getline(self, cube):
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




register_provenance('QC', QCMolecule)
register_provenance('Psi4', Psi4Molecule)
register_provenance('xyz', xyzMolecule)
register_provenance("Cube", cubeMolecule)
    

    

        
