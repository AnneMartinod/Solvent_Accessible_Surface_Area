import src.sphere as fonctions
import src.calcul as formule
from pathlib import Path
from tqdm import tqdm


if __name__ == "__main__" : 
    rayons = {'N' : 1.5, 'O': 1.4, 'S' : 1.85, 'CA' : 2,
              'C' : 1.5, 'CB' : 1.5, 'CD1' : 1.85, 'CD2' : 1.85, 'CE1' : 1.85, 'CE2' : 1.85, 
              'CZ' : 1.85, 'CG' : 2, 'CG1': 2, 'CG2' : 2, 'CD' : 2, 'CD1' : 2, 'CD2' : 2,
              'CE' : 2}
    
    rayon_H2O = 1.4
    nb_points = 92

    liste_molecule = list(Path("data").glob('*.pdb'))             

    with open("results/resultats.txt", "w") as f:
        f.write("")

    for mol in tqdm(liste_molecule):
        atoms = fonctions.parsepdb(mol)
        atoms = fonctions.generate_protein_spheres(atoms, n_points=nb_points, probe_radius=rayon_H2O, rayons_dict=rayons)
    
        formule.accessibilite_point(atoms, nb_points, 8)
        surface = formule.surface_accessible(atoms, nb_points)

        with open("results/resultats.txt", "a") as f:
            f.write(f"Pour {mol.stem} : \n")
            f.write(f"Pourcentage de points accessibles : {surface[0]:.2f} %\n")
            f.write(f"Surface totale accessible de la protéine : {surface[1]:.2f} Å²\n\n")
