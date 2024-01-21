# This class is mean to be u into a single object to ease processing
import os
from pathlib import Path
from typing import Iterable

import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
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
        
        # fingerprint
        fingerprint_xml_files = glob.glob(f"{padelpy_metadata_path}/fingerprint_descriptors/*.xml")
        fingerprint_xml_files.sort()
        fingerprint_namelist = [x.replace(f"{padelpy_metadata_path}/fingerprint_descriptors/", "").replace(".xml","") for x in fingerprint_xml_files]
        self.fingerprint_xml_dict = dict(zip(fingerprint_namelist, fingerprint_xml_files))

        # 2D descriptros
        descriptor_2D_xml_files = glob.glob(f"{padelpy_metadata_path}/bidimensional_descriptors/*.xml")
        descriptor_2D_xml_files.sort()
        descriptor_2D_namelist = [x.replace(f"{padelpy_metadata_path}/bidimensional_descriptors/", "").replace(".xml","") for x in descriptor_2D_xml_files]
        self.descriptor_2D_xml_dict = dict(zip(descriptor_2D_namelist, descriptor_2D_xml_files))

        # 2D descriptros
        descriptor_3D_xml_files = glob.glob(f"{padelpy_metadata_path}/tridimensional_descriptors/*.xml")
        descriptor_3D_xml_files.sort()
        descriptor_3D_namelist = [x.replace(f"{padelpy_metadata_path}/tridimensional_descriptors/", "").replace(".xml","") for x in descriptor_3D_xml_files]
        self.descriptor_3D_xml_dict = dict(zip(descriptor_3D_namelist, descriptor_3D_xml_files))

        self.savepath = "./cache.csv"

        self.fingerprint_dict:dict[str, pd.DataFrame] = {}
        self.descriptor_2D_dict:dict[str, pd.DataFrame] = {}
        self.descriptor_3D_dict:dict[str, pd.DataFrame] = {}
        self.rdkit_descriptors:pd.DataFrame = pd.DataFrame()
        
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

        fingerprint_summary = [(key, str(len(value))) for (key, value) in self.fingerprint_dict.items()]
        descriptor_2D_summary = [(key, str(len(value))) for (key, value) in self.descriptor_2D_dict.items()]
        descriptor_3D_summary = [(key, str(len(value))) for (key, value) in self.descriptor_3D_dict.items()]

        summary_string = ""
        for key, length in fingerprint_summary:
            summary_string += f"{key}, length {length}\n"

        for key, length in descriptor_2D_summary:
            summary_string += f"{key}, length {length}\n"

        for key, length in descriptor_3D_summary:
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

    def generate_rdkit_descriptor(self)->None:
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
        self.rdkit_descriptors = df

    def get_rdkit_descriptor(self)->pd.DataFrame:

        return self.rdkit_descriptors
    
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

    def to_image(self, path:str):
        """
        save the molecules to png file
        dir: directory to save the files
        """
        
        for name, molecule in self.molecules.items():

            filename = os.path.join(path, f'{name}.png')
            Draw.MolToFile(molecule, filename, size=(720,720), fitImage=False, imageType='png')

    def generate_padelpy_fingerprint(self, key_list:Iterable[str], max_run_time:int = 100, regenerate:bool = False):

        molecule_filename = "temp.sdf"
        self.to_single_file("temp.sdf", key_list)


        for key, filename in self.fingerprint_xml_dict.items():
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

                    if key in self.fingerprint_dict:
                        current_dataframe = self.fingerprint_dict[key]
                        self.fingerprint_dict[key] = pd.concat([current_dataframe[current_dataframe.isnull().any(axis=1)==False], result], ignore_index=True)
                    else:
                        continue
                else:

                    self.fingerprint_dict[key] = result
                remove(self.savepath)

        if path.isfile(molecule_filename):
            remove(molecule_filename)

    def check_padelpy_fingerprint_empty_list(self)->dict[str, list[str]]:
            
        status_dict:dict[str,list[str]] = {}
        for key, value in self.fingerprint_dict.items():
            
            molecule_list_contain_null = value[value.isnull().any(axis=1)==True]['Name'].to_list()
            molecule_list_contain_null.sort()
            status_dict[key] = molecule_list_contain_null

        return status_dict


    def generate_padelpy_2D_descriptor(self, key_list:Iterable[str], max_run_time:int = 100, regenerate:bool = False):

        molecule_filename = "temp.sdf"
        self.to_single_file("temp.sdf", key_list)


        for key, filename in self.descriptor_2D_xml_dict.items():

            try:
                padeldescriptor(
                    mol_dir=molecule_filename, 
                    maxruntime = max_run_time,
                    descriptortypes= filename,
                    d_file=self.savepath,
                    detectaromaticity=True,
                    standardizenitro=True,
                    standardizetautomers=True,
                    threads=self.padelpy_threads,
                    removesalt=True,
                    fingerprints=True,
                    d_2d=True)

            except Exception as E:

                print(E)
                continue

            if path.isfile(self.savepath):
                result = pd.read_csv(self.savepath)

                if regenerate:

                    if key in self.fingerprint_dict:
                        current_dataframe = self.fingerprint_dict[key]
                        self.fingerprint_dict[key] = pd.concat([current_dataframe[current_dataframe.isnull().any(axis=1)==False], result], ignore_index=True)
                    else:
                        continue
                else:

                    self.fingerprint_dict[key] = result
                remove(self.savepath)

        if path.isfile(molecule_filename):
            remove(molecule_filename)

    def check_padelpy_descriptor_2D_empty_list(self)->dict[str, list[str]]:
            
        status_dict:dict[str,list[str]] = {}
        for key, value in self.descriptor_2D_dict.items():
            
            molecule_list_contain_null = value[value.isnull().any(axis=1)==True]['Name'].to_list()
            molecule_list_contain_null.sort()
            status_dict[key] = molecule_list_contain_null

        return status_dict

    def generate_padelpy_3D_descriptor(self, key_list:Iterable[str], max_run_time:int = 100, regenerate:bool = False):

        molecule_filename = "temp.sdf"
        self.to_single_file("temp.sdf", key_list)


        for key, filename in self.descriptor_3D_xml_dict.items():

            try:

                padeldescriptor(
                    mol_dir=molecule_filename,
                    maxruntime = max_run_time,
                    descriptortypes= filename,
                    d_file=self.savepath,
                    detectaromaticity=True,
                    standardizenitro=True,
                    standardizetautomers=True,
                    threads=self.padelpy_threads,
                    removesalt=True,
                    fingerprints=True,
                    convert3d=True,
                    d_3d=True) 

            except Exception as E:

                print(E)
                continue

            if path.isfile(self.savepath):
                result = pd.read_csv(self.savepath)

                if regenerate:

                    if key in self.fingerprint_dict:
                        current_dataframe = self.fingerprint_dict[key]
                        self.fingerprint_dict[key] = pd.concat([current_dataframe[current_dataframe.isnull().any(axis=1)==False], result], ignore_index=True)
                    else:
                        continue
                else:

                    self.fingerprint_dict[key] = result
                remove(self.savepath)

        if path.isfile(molecule_filename):
            remove(molecule_filename)

    def check_padelpy_descriptor_3D_empty_list(self)->dict[str, list[str]]:
            
        status_dict:dict[str,list[str]] = {}
        for key, value in self.descriptor_3D_dict.items():
            
            molecule_list_contain_null = value[value.isnull().any(axis=1)==True]['Name'].to_list()
            molecule_list_contain_null.sort()
            status_dict[key] = molecule_list_contain_null

        return status_dict
if __name__ == "__main__":
    pass


