import math
import matplotlib.pyplot as plt

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


def generate_protein_spheres(atoms, n_points, probe_radius, rayons_dict):
    """
    Associe à chaque atome son nuage de points sur sa sphère étendue.
    """
    for atom in atoms:
        r_vdw = get_vdw_radius(atom["name"], rayons_dict)
        r_etendu = r_vdw + probe_radius
        
        # Stockage du rayon étendu et de la liste de coordonnées des points
        atom["radius_extended"] = r_etendu
        atom["sphere_points"] = saff_kuijlaars_sphere(
            n_points=n_points,
            radius=r_etendu,
            center=atom["coords"]
        )
    return atoms
