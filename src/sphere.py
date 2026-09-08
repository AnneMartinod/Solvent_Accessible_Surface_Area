"""Extraction des données d'un fichier PDB et génération des points sur les sphères
atomiques."""

import math

def parsepdb(file_name) : 
    """ Extraction des informations des atomes d'un fichier PDB en ignorant atome d'hydrogene.

    Arguments : 
        file_name (str): nom du fichier PDB.
    
    Returns : 
    liste de dictionnaires représentant chaque atome.
        Clé : 
            "id" (int) : numéro de l'atome.
            "name" (str): nom de l'atome.
            "res_name" (str): nom du résidu.
            "chain" (str): chaine protéique.
            "res_id" (int): numéro du résidu.
            "coords" (tuple[float, float, float]): Coordonnées atomiques (x, y, z).
    """

    atoms = []
    with open(file_name, "r") as f : 
        for line in f :
            if line.startswith("ATOM"):
                atom_id = int(line[6:11].strip())
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain_id = line[21].strip()
                res_seq = int(line[22:26].strip())

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
                        "coords": (x, y, z)
                    })
    return atoms

def saff_kuijlaars_sphere(radius, n_points, center):
    """
    Génère n_points répartis uniformément sur une sphère selon l'algorithme 
    de Saff & Kuijlaars (1997).

    Arguments : 
        n_points (int) : nombre de points à placer sur la sphère.
        radius (float) : rayon de la sphère.
        center (tuple) : centre de la sphère.

    Returns : 
        list(tuple): chaque tuple contient les coordonnées x, y, z d'un point de la sphère.
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
    """Récupère le rayon VdW associé à un atome.
    
    Arguments : 
        atom_name (str) : nom de l'atome.
        rayons_dict (dict) : dictionnaire avec pour clé les noms d'atome et valeur
            leurs rayons.
    
    Returns : 
        float : Rayon de van der Waals de l'atome recherché.
    """
    if atom_name in rayons_dict:
        return rayons_dict[atom_name]
    else:
        return rayons_dict[atom_name[0]]
    


def generate_protein_spheres(atoms, probe_radius, rayons_dict, n_points):
    """
    Associe à chaque atome son nuage de points sur sa sphère étendue.

    Arguments : 
        atoms (list[dict]) : liste de dictionnaires représentant chaque atome.
        n_points (int) : nombre de points à placer sur la sphère.
        probe_radius (float) : rayon de la sonde.
        rayons_dict (dict) : dictionnaire avec pour clé les noms d'atomes et pour valeurs
            leurs rayons.

    Returns : 
        liste[dict] : liste modifiées en place avec clés ajoutées :  
            "radius_extended" (float) : rayon étendu de l'atome.
            "sphere_points" (list[tuple]): Coordonnées des points de la sphère. 
    """
    for atom in atoms:
        r_vdw = get_vdw_radius(atom["name"], rayons_dict)
        r_etendu = r_vdw + probe_radius
        
        # Stockage du rayon étendu et de la liste de coordonnées des points
        atom["radius_extended"] = r_etendu
        atom["sphere_points"] = saff_kuijlaars_sphere(
            radius=r_etendu,
            n_points=n_points,
            center=atom["coords"]
        )
    return atoms
