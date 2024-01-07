# This class is mean to be u into a single object to ease processing
import os
from pathlib import Path
from typing import Iterable

import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import Descriptors

import glob
from os import path, remove
from padelpy import padeldescriptor
import pandas as pd

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

    def __init__(self, molecule_dictionary: dict[str, Chem.rdchem.Mol], padelpy_threads:int = 5):
        self.molecules = molecule_dictionary


        padelpy_metadata_path = path.join(path.dirname(__file__), 'padelpy_metadata')
        fingerprint_xml_files = glob.glob(f"{padelpy_metadata_path}/fingerprint_descriptors/*.xml")
        fingerprint_xml_files.sort()

        fingerprint_namelist = [x.replace(f"{padelpy_metadata_path}/fingerprint_descriptors/", "").replace(".xml","") for x in fingerprint_xml_files]
        self.savepath = "./cache.csv"
        self.fingerprint_dict = dict(zip(fingerprint_namelist, fingerprint_xml_files))
        self.descriptor_dict:dict[str, pd.DataFrame] = {}
        self.padelpy_threads = padelpy_threads

    @classmethod
    def from_path(cls, path, padelpy_threads:int = 5):
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

        return cls(molecule_dictionary, padelpy_threads)  

    def __str__(self):

        descriptor_summary = [(key, str(len(value))) for (key, value) in self.descriptor_dict.items()]

        summary_string = ""
        for key, length in descriptor_summary:
            summary_string += f"{key}, length {length}\n"


        return f"Molecule_Aggregator: length = {len(self.molecules)} \nkeys = {list(self.molecules.keys())} \n\nPadelpy Descriptors:\n{summary_string}"

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

    def to_single_file(self, filename:str, key_list:Iterable[str]):
        """This method is useful for generating descriptor using padelpy, generates a sdf file that contain all the molecules"""

        with Chem.SDWriter(filename) as w:

            for key, molecule in [(key, value) for (key, value) in self.molecules.items() if key in key_list]:
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

    def generate_padelpy_fingerprint(self, key_list:Iterable[str], max_run_time:int = 100, regenerate:bool = False):

        molecule_filename = "temp.sdf"
        self.to_single_file("temp.sdf", key_list)


        for key, filename in self.fingerprint_dict.items():
            try:

                padeldescriptor(
                    mol_dir=molecule_filename, 
                    maxruntime = max_run_time,
                    descriptortypes=filename,
                    d_file=self.savepath,
                    detectaromaticity=True,
                    standardizenitro=True,
                    standardizetautomers=True,
                    threads=self.padelpy_threads,
                    removesalt=True,
                    fingerprints=True)

            except Exception as E:

                print(E)
                continue

            if path.isfile(self.savepath):
                result = pd.read_csv(self.savepath)

                if regenerate:

                    current_dataframe = self.descriptor_dict[key]
                    self.descriptor_dict[key] = pd.concat([current_dataframe[current_dataframe.isnull().any(axis=1)==False], result], ignore_index=True)
                    print(result)
                else:

                    self.descriptor_dict[key] = result
                remove(self.savepath)

        if path.isfile(molecule_filename):
            remove(molecule_filename)

    def check_padelpy_descriptor_empty_list(self)->dict[str, list[str]]:
            
        status_dict:dict[str,list[str]] = {}
        for key, value in self.descriptor_dict.items():
            
            molecule_list_contain_null = value[value.isnull().any(axis=1)==True]['Name'].to_list()
            molecule_list_contain_null.sort()
            status_dict[key] = molecule_list_contain_null

        return status_dict



if __name__ == "__main__":
    pass


