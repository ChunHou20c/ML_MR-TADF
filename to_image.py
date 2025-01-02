from descriptor_generator.descriptor_generator import Molecule_Aggregate

molecule_dir = "/home/chunhou/Dev/python/ML_MR-TADF/data/processed/"
molecules = Molecule_Aggregate.from_path(molecule_dir)
molecules.to_image("/home/chunhou/Dev/python/ML_MR-TADF/data/image/")
