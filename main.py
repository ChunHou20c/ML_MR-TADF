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

    with open(filename, 'wb') as f:

        molecule = pickle.load(f)

        return molecule

with open("config.json") as f:
    configs = json.load(f)


    threads = configs.get("max_threads")
    runtime_multiplier = configs.get("runtime_multiplier")
    molecule_cache = configs.get("molecule_cache")

    if path.isfile(molecule_cache):
        molecules = load_molecule(molecule_cache)

    else:
        #on first calculation do whatever that need to be done
        molecule_dir = configs.get("molecule_dir")
        molecules = Molecule_Aggregate.from_path(molecule_dir, threads)
        molecules.optimize()
        molecules.generate_rdkit_descriptor()

        molecules.generate_padelpy_fingerprint(molecules.molecules.keys(), max_run_time = runtime_multiplier, regenerate = False)
        molecules.generate_padelpy_2D_descriptor(molecules.molecules.keys(), max_run_time = runtime_multiplier, regenerate = False)
        molecules.generate_padelpy_3D_descriptor(molecules.molecules.keys(), max_run_time = runtime_multiplier, regenerate = False)
        
        cache_molecule(molecules, molecule_cache)

    for i in range(1,6):
    # this part might be repeated to make padelpy descriptor complete, thought it might never able to be completed
    # everytime a caculation is finished it should cache the molecule again to prevent loss

        fingerprint_status = molecules.check_padelpy_fingerprint_empty_list()
        descriptor_2D_status = molecules.check_padelpy_descriptor_2D_empty_list()
        descriptor_3D_status = molecules.check_padelpy_descriptor_3D_empty_list()

        for (key, incomplete_list) in fingerprint_status.items():
            molecules.generate_padelpy_fingerprint(incomplete_list, max_run_time = runtime_multiplier * (10 ** i) , regenerate=True)
            cache_molecule(molecules, molecule_cache)

        for (key, incomplete_list) in descriptor_2D_status.items():
            molecules.generate_padelpy_2D_descriptor(incomplete_list, max_run_time = runtime_multiplier * (10 ** i), regenerate=True)
            cache_molecule(molecules, molecule_cache)

        for (key, incomplete_list) in descriptor_3D_status.items():
            molecules.generate_padelpy_3D_descriptor(incomplete_list, max_run_time = runtime_multiplier * (10 ** i), regenerate=True)
            cache_molecule(molecules, molecule_cache)

print(molecules)
#analyze the completeness of the descriptors 
print(f'fingerprint status: {molecules.check_padelpy_fingerprint_empty_list()} \n')
print(f'2D descriptor status: {molecules.check_padelpy_descriptor_2D_empty_list()} \n')
print(f'3D descriptor status: {molecules.check_padelpy_descriptor_3D_empty_list()} \n')
#descriptors = molecules.generate_rdkit_descriptor()
#print(descriptors)
