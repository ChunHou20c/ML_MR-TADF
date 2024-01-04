# this module is mean to reduce the work to calculate padelpy descriptors
import glob
from os import path, remove
from padelpy import padeldescriptor
import pandas as pd

class Molecule_File_Aggregator:

    def __init__(self, filepath:str, savepath:str) -> None:
        self.molecules = filepath
        self.savepath = savepath
        self.descriptor_dict = {}
        

    def __str__(self):
        return f"Molecule_File_Aggregator: length = {len(self.descriptor_dict)} \nkeys = {self.descriptor_dict.keys()}"

    def gen_fingerPrint(self):

        padelpy_metadata_path = path.join(path.dirname(__file__), 'padelpy_metadata')
        fingerprint_xml_files = glob.glob(f"{padelpy_metadata_path}/fingerprint_descriptors/*.xml")
        fingerprint_xml_files.sort()

        fingerprint_namelist = [x.replace(f"{padelpy_metadata_path}/fingerprint_descriptors/", "").replace(".xml","") for x in fingerprint_xml_files]


        for file, key in zip(fingerprint_xml_files, fingerprint_namelist):
            try:

                padeldescriptor(
                    mol_dir=self.molecules, 
                    maxruntime = 100,
                    descriptortypes= file,
                    d_file=self.savepath,
                    detectaromaticity=True,
                    standardizenitro=True,
                    standardizetautomers=True,
                    threads=1,
                    removesalt=True,
                    fingerprints=True)

            except Exception as E:

                print(E)
                continue

            if path.isfile(self.savepath):
                result = pd.read_csv(self.savepath)
                self.descriptor_dict[key] = result
                remove(self.savepath)




if __name__ == "__main__":
    
    molecule_file = Molecule_File_Aggregator("test.sdf", "result.csv");
    molecule_file.gen_fingerPrint()
    print(molecule_file)
