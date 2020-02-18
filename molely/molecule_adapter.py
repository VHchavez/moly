from abc import ABC, abstractmethod

class MoleculeAdapter(ABC):

    @abstractmethod
    def get_geometry(self):
        pass

    @abstractmethod
    def get_symbols(self):
        pass

    @abstractmethod
    def get_connectivity(self):
        pass

    