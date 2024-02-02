# use to inspect the pickled data

from descriptor_generator.descriptor_generator import Molecule_Aggregate
import json
import pickle

from os import path

# helper functions to help cache the data of calculated molecules
def cache_molecule(molecule:Molecule_Aggregate, filename: str)->None:
    """
    This function save the molecule as cache to continue from incomplete calculation
    """

    with open(filename, 'wb') as f:

        pickle.dump(molecule, f)

def load_molecule(filename:str)->Molecule_Aggregate:

    with open(filename, 'rb') as f:

        molecule = pickle.load(f)

        return molecule

with open("config.json") as f:
    configs = json.load(f)


    threads = configs.get("max_threads")
    runtime_multiplier = configs.get("runtime_multiplier")
    molecule_cache = configs.get("molecule_cache")

    if path.isfile(molecule_cache):
        molecules = load_molecule(molecule_cache)

        for key, value in molecules.fingerprint_dict.values():

            print(key)
            print(value)


        for key, value in molecules.descriptor_2D_dict.values():

            print(key)
            print(value)


        for key, value in molecules.descriptor_3D_dict.values():

            print(key)
            print(value)


