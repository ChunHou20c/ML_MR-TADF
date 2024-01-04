# this module is mean to reduce the work to calculate padelpy descriptors
import glob
from os import path, PathLike
from padelpy import padeldescriptor

class Molecule_File_Aggregator:

    def __init__(self, filepath:str) -> None:
        self.molecules = filepath
        

    def gen_fingerPrint(self):

        padelpy_metadata_path = path.join(path.dirname(__file__), 'padelpy_metadata')
        fingerprint_xml_files = glob.glob(f"{padelpy_metadata_path}/fingerprint_descriptors/*.xml")
        fingerprint_xml_files.sort()

        fingerprint_namelist = [x.replace("fingerprint_descriptors/", "").replace(".xml","") for x in fingerprint_xml_files]

        try:

            padeldescriptor(
                mol_dir=self.molecules, 
                maxruntime = 100,
                descriptortypes= fingerprint_xml_files[0],
                d_file="result.csv",
                detectaromaticity=True,
                standardizenitro=True,
                standardizetautomers=True,
                threads=1,
                removesalt=True,
                fingerprints=True)

        except Exception as E:

            print(E)


if __name__ == "__main__":
    
    molecule_file = Molecule_File_Aggregator("molecule.mol");
    molecule_file.gen_fingerPrint()
