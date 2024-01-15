import os
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

def from_xyz_to_mol(from_path, to_path):
    """
    path: the path all the molecules located

    """
    
    def extract_integer(filename):
        return int(filename.split('-')[0])

    for file in sorted(os.listdir(from_path), key = extract_integer):

        #get the full filename with path
        filename = os.path.join(from_path, file)
        molecule = Chem.MolFromXYZFile(filename)
        mol = Chem.Mol(molecule)
        rdDetermineBonds.DetermineBonds(mol,charge=-1)

        # get the basename of molecule without file extension
        key = Path(filename).stem

        save_filename = os.path.join(to_path, f'{key}.mol')
        Chem.MolToMolFile(molecule, save_filename)

from_xyz_to_mol("/path1", "/path2")
