"""Calcul de la surface accessible au solvant (SASA) pour des structures PDB."""

import sphere as fonctions
import calcul as formule
from pathlib import Path
import time


if __name__ == "__main__" : 
    debut = time.time()

    rayons = {'N' : 1.5, 'O': 1.4, 'S' : 1.85, 'CA' : 2,
              'C' : 1.5, 'CB' : 1.5, 'CD1' : 1.85, 'CD2' : 1.85, 'CE1' : 1.85, 'CE2' : 1.85, 
              'CZ' : 1.85, 'CG' : 2, 'CG1': 2, 'CG2' : 2, 'CD' : 2,
              'CE' : 2}
    
    rayon_probe = 1.4
    #nb_points = 100

    liste_molecule = list(Path("data").glob('*.pdb'))     

    if 'nb_points' not in globals():
        nb_points = 92
        dossier_resultats = Path("results/default_results")
    else:
        dossier_resultats = Path(f"results/results_{nb_points}")

    dossier_resultats.mkdir(exist_ok=True, parents=True)
    

    with open(dossier_resultats / "results.txt", "w") as f:
        f.write(f"Avec {nb_points} points par sphère : \n\n")

        for mol in liste_molecule:
            debut_iter = time.time()
            print(f"Molécule {mol.stem}")
            
            atoms = fonctions.parsepdb(mol)
            atoms = fonctions.generate_protein_spheres(atoms, n_points=nb_points, probe_radius=rayon_probe, rayons_dict=rayons)
        
            formule.accessibilite_point(atoms, nb_points)

            surface = formule.surface_accessible(atoms, nb_points)
   
            f.write(f"Pour {mol.stem} : \n")
            f.write(f"Nombre total d'atomes : {len(atoms)} \n")
            f.write(f"Nombre total de points accessibles : {surface[2]} \n")
            f.write(f"Pourcentage de surface accessible : {surface[0]:.2f} %\n")
            f.write(f"Surface totale accessible de la protéine : {surface[1]:.2f} Å²\n")

            fin_iter = time.time()
            duree_iter = fin_iter - debut_iter
            print(f"Temps d'exécution : {duree_iter:.2f} secondes")
            f.write(f"Temps d'exécution : {duree_iter:.2f} secondes\n\n")

        fin = time.time()
        duree = fin - debut
        print(f"\nTemps d'exécution total : {duree:.2f} secondes")

        f.write(f"Temps d'exécution total : {duree:.2f} secondes")

    
