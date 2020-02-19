from .molecule_adapter import MoleculeAdapter

_provenance = {}

def register_provenance(provenance, provenance_class):

    if not issubclass(provenance_class, MoleculeAdapter):
        raise TypeError(f"{provenance_class} is not a Molecule Adapter")

    _provenance[provenance] = provenance_class


def molecule_factory(provenance, **kwargs):

    molecule_provenance = _provenance[provenance](**kwargs)

    return molecule_provenance


