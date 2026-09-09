"""Calcul de l'accessibilité au solvant d'une protéine."""

import math

def distance(atom1, atom2): 
    """Calcul de la distance euclidienne entre deux atomes.

    Arguments : 
        atom1 (tuple[float, float, float]): Coordonnées (x, y, z) du premier atome.
        atom2 (tuple[float, float, float]): Coordonnées (x, y, z) du deuxième atome.

    Returns : 
        float : distance entre les deux points.
    """
    return math.sqrt((atom1[0]-atom2[0])**2 + 
                    (atom1[1]-atom2[1])**2 + 
                    (atom1[2] - atom2[2])**2)

def near_neighbor(atoms):
    """Identifie les atomes voisins dont les sphères étendues se chevauchent.

    Arguments : 
        atoms (list[dict]) : liste des atomes contenant au moins les clés "coords"
            et "radius_extended".

    Returns : 
        dict[int, list[int]] : dictionnaire associant l'indice d'un atome aux indices 
            de ses voisins. 
    """
    dico = {}
    for i in range(len(atoms)):
        dico[i] = []
        for j in range(len(atoms)):
            if i == j :
                continue
            else : 
                if distance(atoms[i]['coords'], atoms[j]['coords']) < atoms[i]['radius_extended']+atoms[j]['radius_extended']:
                    dico[i].append(j)
    return dico

def accessibilite_point(atoms, nbpoints):
    """Calcule l'accessibilité au solvant de chaque atome en testant ses points de surface.

    Arguments : 
        atoms (list[dict]) : liste des atomes contenant "coords", "radius_extended"
            et "sphere_points".
        nbpoints (int) : nombre de points de la sphère. 

    Returns : 
        list[dict] : liste modifiée en place avec les clés ajoutées : 
            'nb_points_accessibles' (int) : nombre de points exposés au solvant.
            '% points accessibles' (float) : pourcentage de points accessibles.
    """
    voisins = near_neighbor(atoms)
    for i in range(len(atoms)): 
        atoms[i]['nb_points_accessibles'] = 0
        for j in range(nbpoints):
            point_obstrus = False
            for k in voisins[i]:
                d = distance(atoms[i]['sphere_points'][j], atoms[k]['coords'])
                if d < atoms[k]['radius_extended'] :
                    point_obstrus = True
                    break
            if point_obstrus == False : 
                atoms[i]['nb_points_accessibles'] += 1
        atoms[i]['% points accessibles'] = (atoms[i]['nb_points_accessibles']/nbpoints)*100

    return atoms

def surface_accessible(atoms, nbpoints):
    """Calcule la surface accessible au solvant totale (SASA), le pourcentage de points
    accessibles et le nombre total de points accessibles.

    Arguments : 
        atoms (list[dict]) : Liste des atomes contenant au moins les clés "radius_extended"
            et "nb_points_accessibles".
        nbpoints (int) : nombre de points de la sphère.

    Returns : 
        tuple[float, float, int] : tuple contenant : 
            pourcentage de surface accessibles (float).
            surface totale accessible (float).
            nombre total de points accessibles (int). 
    """
    points_accessibles = 0
    surface_totale_accessible = 0
    surface_maximum = 0
    for atom in atoms : 
        surface = 4 * math.pi * atom['radius_extended']**2
        surface_maximum += surface
        surface_1pt = surface / nbpoints
        points_accessibles += atom['nb_points_accessibles']
        surface_accessible = atom['nb_points_accessibles'] * surface_1pt
        surface_totale_accessible += surface_accessible

    perc_surface_access = (surface_totale_accessible/surface_maximum)*100

    return perc_surface_access, surface_totale_accessible, points_accessibles


