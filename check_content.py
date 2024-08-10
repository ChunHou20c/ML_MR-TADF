# use to inspect the pickled data

from descriptor_generator.descriptor_generator import Molecule_Aggregate
import json
import pickle
import pandas as pd
from functools import reduce

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

with open("config.json", "r") as f:
    configs = json.load(f)

    molecule_cache = configs.get("molecule_cache")
    print(molecule_cache)

    if path.isfile(molecule_cache):
        molecules = load_molecule(molecule_cache)

        dataframes_list = [ i.drop_duplicates(subset="Name", keep="last") for i in molecules.fingerprint_dict.values()]
        # print(dataframes_list)
        merged_df = dataframes_list[0]

        for df in dataframes_list[1:]:
            merged_df = merged_df.merge(df, on='Name', how='outer')

        dataframes_list = [ i.drop_duplicates(subset="Name", keep="last") for i in molecules.descriptor_2D_dict.values()]
        
        for df in dataframes_list:
            merged_df = merged_df.merge(df, on='Name', how='outer')

        dataframes_list = [ i.drop_duplicates(subset="Name", keep="last") for i in molecules.descriptor_3D_dict.values()]
        
        for df in dataframes_list:
            merged_df = merged_df.merge(df, on='Name', how='outer')

        merged_df = merged_df.merge(molecules.rdkit_descriptors, on='Name', how='outer')

        merged_df = merged_df.dropna(axis=1)
        merged_df.to_csv("inspect.csv", index=False)



