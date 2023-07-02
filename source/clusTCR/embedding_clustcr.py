from profile.properties import OPTIMAL
from tools import make_profiles
import pandas as pd
from typing import Union

class Embeddingclustcr:

    def __init__(self, max_sequence_size=None, properties: list = OPTIMAL, n_cpus: Union[str, int] = 1):
        self.properties = properties
        self.max_sequence_size = max_sequence_size
        self.n_cpus = n_cpus
        self.profiles = None

    def embed(self, data: pd.Series):

        matrix = make_profiles(data, self.properties, self.max_sequence_size, self.n_cpus)

        return matrix