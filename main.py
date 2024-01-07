from descriptor_generator.descriptor_generator import Molecule_Aggregate

molecules = Molecule_Aggregate.from_path("/home/chunhou/Dev/python/padelpy_descriptor_generation/data/molecules/optimized/fixed_first_batch/", 5)
molecules.generate_padelpy_fingerprint(molecules.molecules.keys())
print(molecules)
status = molecules.check_padelpy_descriptor_status()
print(status)
#descriptors = molecules.generate_rdkit_descriptor()
#print(descriptors)
