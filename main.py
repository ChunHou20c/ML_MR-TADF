from descriptor_generator.descriptor_generator import Molecule_Aggregate

molecules = Molecule_Aggregate.from_path("/home/chunhou/Dev/python/padelpy_descriptor_generation/data/molecules/optimized/fixed_first_batch/")
descriptors = molecules.generate_rdkit_descriptor()
print(descriptors)
