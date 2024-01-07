from descriptor_generator.descriptor_generator import Molecule_Aggregate

molecules = Molecule_Aggregate.from_path("/home/chunhou/Dev/python/padelpy_descriptor_generation/data/molecules/optimized/fixed_first_batch/", 5)
molecules.generate_padelpy_fingerprint(molecules.molecules.keys(), max_run_time = 100, regenerate = False)
print(molecules)
status = molecules.check_padelpy_descriptor_empty_list()
print(status)

for (key, incomplete_list) in status.items():
    molecules.generate_padelpy_fingerprint(incomplete_list, max_run_time = 1000, regenerate=True)

status = molecules.check_padelpy_descriptor_empty_list()
print(status)
#descriptors = molecules.generate_rdkit_descriptor()
#print(descriptors)
