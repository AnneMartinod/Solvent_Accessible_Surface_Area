import math
import matplotlib.pyplot as plt

# rayons = {'N' : 1.5, 'NE2' : 1.5, 'ND1' : 1.5, 'ND2' : 1.5, 'NE' : 1.5, 'NH1':1.5, 
#           'NH2' : 1.5, 'NZ' : 1.5, 'O': 1.4, 'OXT' : 1.4, 'OE1' : 1.4, 'OE2' : 1.4,
#           'OD1' : 1.4, 'OG' : 1.4, 'OG1' : 1.4, 'OH' : 1.4, 'SG' : 1.85, 'CA' : 2,
#           'C' : 1.5, 'CB' : 1.5, 'CD1' : 1.85, 'CD2' : 1.85, 'CE1' : 1.85, 'CE2' : 1.85, 
#           'CZ' : 1.85, 'CG' : 2, 'CG1': 2, 'CG2' : 2, 'CD' : 2, 'CD1' : 2, 'CD2' : 2,
#           'CE' : 2, 'CG' : 2}

rayons = {'N' : 1.5, 'O': 1.4, 'S' : 1.85, 'CA' : 2,
          'C' : 1.5, 'CB' : 1.5, 'CD1' : 1.85, 'CD2' : 1.85, 'CE1' : 1.85, 'CE2' : 1.85, 
          'CZ' : 1.85, 'CG' : 2, 'CG1': 2, 'CG2' : 2, 'CD' : 2, 'CD1' : 2, 'CD2' : 2,
          'CE' : 2, 'CG' : 2}

rayon_H2O = 1.4

def parsepdb(file_name) : 
    atoms = []
    with open(file_name, "r") as f : 
        for line in f :
            if line.startswith("ATOM"):
                atom_id = int(line[6:11].strip())
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain_id = line[21].strip()
                res_seq = int(line[22:26].strip())

                # Coordonnées (x, y, z) en Angströms
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())

                if atom_name[0] != 'H':
                    atoms.append({
                        "id": atom_id,
                        "name": atom_name,
                        "res_name": res_name,
                        "chain": chain_id,
                        "res_id": res_seq,
                        "coords": (x, y, z),
                    })
    return atoms

def saff_kuijlaars_sphere(n_points, radius, center=(0.0, 0.0, 0.0)):
    """
    Génère n_points répartis uniformément sur une sphère de rayon 'radius'
    centrée sur 'center' (cx, cy, cz) selon l'algorithme de Saff & Kuijlaars (1997).
    """
    points = []
    cx, cy, cz = center
    phi = 0.0

    for k in range(1, n_points + 1):
        # h va de -1 à +1
        h = -1.0 + 2.0 * (k - 1) / (n_points - 1)
        theta = math.acos(h)

        if k == 1 or k == n_points:
            phi = 0.0
        else:
            phi = (phi + 3.6 / math.sqrt(n_points * (1.0 - h * h))) % (2.0 * math.pi)

        # Coordonnées cartésiennes décalées au centre de l'atome
        x = cx + radius * math.sin(theta) * math.cos(phi)
        y = cy + radius * math.sin(theta) * math.sin(phi)
        z = cz + radius * math.cos(theta)

        points.append((x, y, z))

    return points


def get_vdw_radius(atom_name, rayons_dict):
    """Récupère le rayon VdW selon le nom de l'atome."""
    if atom_name in rayons_dict:
        return rayons_dict[atom_name]
    else:
        return rayons_dict[atom_name[0]]
    # if atom_name in rayons_dict:
    #     return rayons_dict[atom_name]
    # # Fallback sur l'élément (ex: 'C', 'N', 'O', 'S')
    # elem = atom_name[0]
    # return rayons_dict.get(elem, 1.70) # a verifier si mon dictionnaire est complet au début


def generate_protein_spheres(atoms, n_points, probe_radius):
    """
    Associe à chaque atome son nuage de points sur sa sphère étendue.
    """
    for atom in atoms:
        r_vdw = get_vdw_radius(atom["name"], rayons)
        r_etendu = r_vdw + probe_radius
        
        # Stockage du rayon étendu et de la liste de coordonnées des points
        atom["radius_extended"] = r_etendu
        atom["sphere_points"] = saff_kuijlaars_sphere(
            n_points=n_points,
            radius=r_etendu,
            center=atom["coords"]
        )
    return atoms

if __name__ == "__main__":
    atoms = parsepdb("3I40.pdb")
    atoms = generate_protein_spheres(atoms, n_points=92, probe_radius=rayon_H2O)

    print(f"Nombre total d'atomes traités : {len(atoms)}")
    print(f"Exemple pour l'atome 1 ({atoms[0]['name']}) :")
    print(f"  - Centre : {atoms[0]['coords']}")
    print(f"  - Rayon étendu : {atoms[0]['radius_extended']:.2f} Å")
    print(f"  - Nombre de points générés : {len(atoms[0]['sphere_points'])}")
    print(f"  - Coordonnées du 1er point : {atoms[0]['sphere_points'][0]}")

    #print(atoms[0])
    print(len(atoms))



    # def afficher_sphere(atom):
    #     """
    #     Affiche en 3D le centre d'un atome et son nuage de points
    #     sur la sphère étendue.
    #     """
    #     fig = plt.figure(figsize=(7, 7))
    #     ax = fig.add_subplot(111, projection='3d')

    #     # Coordonnées du centre de l'atome
    #     cx, cy, cz = atom["coords"]
    #     ax.scatter(cx, cy, cz, color="red", s=80, label="Centre atome")

    #     # Coordonnées de tous les points de la sphère
    #     xs = [p[0] for p in atom["sphere_points"]]
    #     ys = [p[1] for p in atom["sphere_points"]]
    #     zs = [p[2] for p in atom["sphere_points"]]
    #     ax.scatter(xs, ys, zs, color="blue", s=15, alpha=0.6, label="Points sphère")

    #     ax.set_xlabel("X (Å)")
    #     ax.set_ylabel("Y (Å)")
    #     ax.set_zlabel("Z (Å)")
    #     ax.set_title(f"Sphère de {atom['name']} (rayon étendu = {atom['radius_extended']:.2f} Å)")
    #     ax.legend()

    #     # Pour que la sphère ne soit pas déformée visuellement
    #     ax.set_box_aspect([1, 1, 1])

    #     plt.show()


    # # Utilisation
    # afficher_sphere(atoms[2])
# atoms = parsepdb("3I40.pdb")
# print(atoms[0])

# for atom in atoms : 
#     print(atom)

# print(len(atoms))

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