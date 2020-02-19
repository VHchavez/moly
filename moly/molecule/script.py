from adapters import molecule_factory
from adapters import QCMolecule
import qcelemental as qcel

#mol = qcel.models.Molecule.from_data("pubchem:caffeine")
Molecule = molecule_factory("xyz", file='h2o.xyz')

print(Molecule.geometry)
print(Molecule.symbols)
print(Molecule.bonds)
