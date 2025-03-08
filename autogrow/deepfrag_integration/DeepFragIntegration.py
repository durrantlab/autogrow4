from abc import ABC, abstractmethod
import numpy as np


class DeepFragBase(ABC):

    @abstractmethod
    def get_prediction_for_parent_receptor(self, parent_mol, receptor, branching_point):
        pass

    @abstractmethod
    def get_fingerprints_for_fragment(self, fragment):
        pass


class DeepFragFake(DeepFragBase):

    def get_prediction_for_parent_receptor(self, parent_mol, receptor, branching_point):
        return np.random.rand(2048)

    def get_fingerprints_for_fragment(self, fragment):
        return np.random.randint(2, size=2048)
