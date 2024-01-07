from descriptor_generator.descriptor_generator import Molecule_Aggregate
import json

with open("config.json") as f:
    configs = json.load(f)

    molecule_dir = configs.get("molecule_dir")
    threads = configs.get("max_threads")
    runtime_multiplier = configs.get("runtime_multiplier")

    molecules = Molecule_Aggregate.from_path(molecule_dir, threads)

    molecules.generate_padelpy_fingerprint(molecules.molecules.keys(), max_run_time = runtime_multiplier, regenerate = False)
    molecules.generate_padelpy_2D_descriptor(molecules.molecules.keys(), max_run_time = runtime_multiplier, regenerate = False)
    molecules.generate_padelpy_3D_descriptor(molecules.molecules.keys(), max_run_time = runtime_multiplier, regenerate = False)
    print(molecules)
    fingerprint_status = molecules.check_padelpy_fingerprint_empty_list()
    print(fingerprint_status)

    descriptor_2D_status = molecules.check_padelpy_descriptor_2D_empty_list()
    print(descriptor_2D_status)

    descriptor_3D_status = molecules.check_padelpy_descriptor_3D_empty_list()
    print(descriptor_3D_status)
    for (key, incomplete_list) in fingerprint_status.items():
        molecules.generate_padelpy_fingerprint(incomplete_list, max_run_time = runtime_multiplier * 10, regenerate=True)

    status = molecules.check_padelpy_fingerprint_empty_list()
print(status)
#descriptors = molecules.generate_rdkit_descriptor()
#print(descriptors)
