# This class is mean to be u into a single object to ease processing
import os
from pathlib import Path

import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import Descriptors

class Molecule_Aggregate:
    """
    example
    ---------------
    path = "/home/chunhou/Dev/python/padelpy_descriptor_generation/data/molecules/optimized/all_molecules/"
    molecules = Molecule_Aggregate.from_path(path)
    molecules.check_partial_charge()
    
    savepath = "/home/chunhou/Dev/python/padelpy_descriptor_generation/data/molecules/optimized/addHs/"
    molecules.optimize()
    #molecules.to_mol_files(savepath)
    molecules.to_single_file("all_molecules.sdf")
    descriptors = molecules.generate_rdkit_descriptor()
    descriptors.to_csv("all_descriptor.csv", index=False)

    print(molecules)
    """

    def __init__(self, molecule_dictionary: dict[str, Chem.rdchem.Mol]):
        self.molecules = molecule_dictionary

    @classmethod
    def from_path(cls, path):
        """
        path: the path all the molecules located
        return: Molecule_Aggregator

        """
        
        molecule_dictionary = {}
        def extract_integer(filename):
            return int(filename.split('-')[0])

        for file in sorted(os.listdir(path), key = extract_integer):

            #get the full filename with path
            filename = os.path.join(path, file)
            molecule = Chem.MolFromMolFile(filename, removeHs=False)

            # get the basename of molecule without file extension
            key = Path(filename).stem
            molecule_dictionary.update({key: molecule})

        return cls(molecule_dictionary)  

    def __str__(self):
        return f"Molecule_Aggregator: length = {len(self.molecules)} \nkeys = {self.molecules.keys()}"

    def optimize(self):
        """Calling this function will optimize the geometry of the molecules using rdkit optimization"""


        # I didn't use iterator because I can't mutate the value in the dictionary using iterator
        for key in self.molecules.keys():
            print(key)

            try:
                AllChem.EmbedMolecule(self.molecules[key])
                AllChem.MMFFOptimizeMolecule(self.molecules[key])

            except Exception as E:
                print(E)
                print(key)

    def to_single_file(self, filename:str):
        """This method is useful for generating descriptor using padelpy"""

        with Chem.SDWriter(filename) as w:

            for key, molecule in self.molecules.items():
                molecule.SetProp("_Name",key)
                w.write(molecule)

    def generate_rdkit_descriptor(self)->pd.DataFrame:
        """
        generate rdkit descriptor and return a pandas.DataFrame
        """

        descriptors=[x[0] for x in Descriptors._descList]
        Descriptor_calculator = MoleculeDescriptors.MolecularDescriptorCalculator(descriptors)

        def calc_descriptor(key, molecule):
            results = Descriptor_calculator.CalcDescriptors(molecule)
            dictionary = dict(zip(descriptors, results))
            dictionary["Name"] = key
            return dictionary
        
        result_array = map(calc_descriptor, self.molecules.keys(), self.molecules.values())
        df = pd.DataFrame(result_array)

        # switch name to first column
        first_column = df.pop('Name')
        df.insert(0, 'Name', first_column)
        return df
    
    def check_partial_charge(self):
        
        for molecule in self.molecules.values():

            print(Chem.rdPartialCharges.ComputeGasteigerCharges(molecule, throwOnParamFailure=True))

    def to_mol_files(self, path:str):
        """
        Save the molecules in individual mol file
        dir: directory to save the files
        """

        for name, molecule in self.molecules.items():

            filename = os.path.join(path, f'{name}.mol')
            Chem.MolToMolFile(molecule, filename)


if __name__ == "__main__":
    pass


