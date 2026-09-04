import math
import src.sphere as fonctions


def distance(atom1, atom2): 
    return math.sqrt((atom1[0]-atom2[0])**2 + (atom1[1]-atom2[1])**2 + (atom1[2] - atom2[2])**2)

def near_neighboor(atoms, seuil):
    dico = {}
    for i in range(len(atoms)):
        dico[i] = []
        for j in range(len(atoms)):
            if i == j :
                continue
            else : 
                if distance(atoms[i]['coords'], atoms[j]['coords']) < seuil:
                    dico[i].append(j)
    return dico

def accessibilite_point(atoms, nbpoints, seuil):
    voisins = near_neighboor(atoms, seuil)
    for i in range(len(atoms)): 
            atoms[i]['nb_points_accessibles'] = 0
            for j in range(nbpoints):
                point_obstrus = False
                for k in voisins[i]:
                    #if i == k : continue
                    #else : 
                    d = distance(atoms[i]['sphere_points'][j], atoms[k]['coords'])
                    if d < atoms[k]['radius_extended'] :
                        point_obstrus = True
                        break
                if point_obstrus == False : 
                    atoms[i]['nb_points_accessibles'] += 1
            atoms[i]['% points accessibles'] = (atoms[i]['nb_points_accessibles']/nbpoints)*100

def surface_accessible(atoms, nbpoints):
    points_accessibles = 0
    surface_totale_accessible = 0
    for atom in atoms : 
        surface = 4 * math.pi * atom['radius_extended']**2
        surface_1pt = surface / nbpoints
        points_accessibles += atom['nb_points_accessibles']
        surface_accessible = atom['nb_points_accessibles'] * surface_1pt
        surface_totale_accessible += surface_accessible

    perc_pt_access = points_accessibles/(nbpoints*len(atoms)) * 100
    return perc_pt_access, surface_totale_accessible

