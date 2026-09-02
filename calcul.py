import math
import test_v2 as fonctions

rayons = {'N' : 1.5, 'O': 1.4, 'SG' : 1.85, 'CA' : 2,
          'C' : 1.5, 'CB' : 1.5, 'CD1' : 1.85, 'CD2' : 1.85, 'CE1' : 1.85, 'CE2' : 1.85, 
          'CZ' : 1.85, 'CG' : 2, 'CG1': 2, 'CG2' : 2, 'CD' : 2, 'CD1' : 2, 'CD2' : 2,
          'CE' : 2, 'CG' : 2}

rayon_H2O = 1.4
nb_points = 92
SEUIL = 8

def distance(atom1, atom2): 
    return math.sqrt((atom1[0]-atom2[0])**2 + (atom1[1]-atom2[1])**2 + (atom1[2] - atom2[2])**2)

def near_neighboor(atom1, atom2):
    if distance(atom1, atom2) < SEUIL : 
        return True
    return False

if __name__ == "__main__" : 
    atoms = fonctions.parsepdb("3I40.pdb")
    atoms = fonctions.generate_protein_spheres(atoms, n_points=nb_points, probe_radius=rayon_H2O)

    for i in range(len(atoms)): 
        atoms[i]['nb_points_accessibles'] = 0
        for j in range(nb_points):
            point_obstrus = False
            for k in range(len(atoms)):
                if i == k : continue
                else : 
                    d = distance(atoms[i]['sphere_points'][j], atoms[k]['coords'])
                    if d < atoms[k]['radius_extended'] :
                        point_obstrus = True
                        break
            if point_obstrus == False : 
                atoms[i]['nb_points_accessibles'] += 1

    for atom in atoms : 
        print(f"{atom['res_name']} {atom['name']} : {atom['nb_points_accessibles']}")

                


# pour i de 0 à nb_atoms :
#     nb_points_accessibles[i] = 0
#     pour k de 0 à 92 :                      # <-- k à l'extérieur maintenant
#         point_cache = faux
#         pour j de 0 à nb_atoms :            # <-- j à l'intérieur
#             si i == j : passe
#             d = distance_euclidienne(point[k] de sphère i, centre[j])
#             si d < rayon_etendu[j] :
#                 point_cache = vrai
#                 arrêter la boucle j          
#         si point_cache == faux :
#             nb_points_accessibles[i] += 1    # incrémenté DANS la boucle k, un point à la fois