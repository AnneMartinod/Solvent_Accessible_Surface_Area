import test_v2 as fonctions

rayons = {'N' : 1.5, 'O': 1.4, 'SG' : 1.85, 'CA' : 2,
          'C' : 1.5, 'CB' : 1.5, 'CD1' : 1.85, 'CD2' : 1.85, 'CE1' : 1.85, 'CE2' : 1.85, 
          'CZ' : 1.85, 'CG' : 2, 'CG1': 2, 'CG2' : 2, 'CD' : 2, 'CD1' : 2, 'CD2' : 2,
          'CE' : 2, 'CG' : 2}

rayon_H2O = 1.4

if __name__ == "__main__" : 
    atoms = fonctions.parsepdb("3I40.pdb")
    atoms = fonctions.generate_protein_spheres(atoms, n_points=92, probe_radius=rayon_H2O)

    print(f"Nombre total d'atomes traités : {len(atoms)}")
    print(f"Exemple pour l'atome 1 ({atoms[0]['name']}) :")
    print(f"  - Centre : {atoms[0]['coords']}")
    print(f"  - Rayon étendu : {atoms[0]['radius_extended']:.2f} Å")
    print(f"  - Nombre de points générés : {len(atoms[0]['sphere_points'])}")
    print(f"  - Coordonnées du 1er point : {atoms[0]['sphere_points'][0]}")

    #print(atoms[0])
    print(len(atoms))