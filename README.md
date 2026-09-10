# Solvent_Accessible_Surface_Area

L'objectif de ce projet est de calculer la surface accessible au solvant d'une protéine. Il repose sur l'article : Shrake, A; Rupley, JA. (1973). "Environment and exposure to solvent of protein atoms. Lysozyme and insulin". J Mol Biol 79 (2): 351–71. doi: 10.1016/0022-2836(73)90011-9

## Principe

Pour chaque atome de la protéine :
1. Une sphère "étendue" est définie autour de l'atome, dont le rayon correspond au rayon de van der Waals de l'atome **+** le rayon de la molécule de solvant (sonde, typiquement une molécule d'eau, rayon = 1.4 Å).
2. Un ensemble de points est réparti uniformément sur cette sphère grâce à l'algorithme de **Saff & Kuijlaars (1997)**.
3. Pour chaque point, on vérifie s'il est masqué par la sphère étendue d'un atome voisin. Si aucun voisin ne le recouvre, le point est considéré comme **accessible au solvant**.
4. Le pourcentage de points accessibles par atome, combiné à la surface de la sphère, permet d'estimer la surface accessible totale de la protéine.

## Structure du projet

```
.
├── data/                   # Fichiers pdb (possibilité d'en ajouter)
│   ├── file1.pdb
│   ├── file2.pdb
│   └── ...
├── results/                # Résultats générés (créé automatiquement)
│   ├── default_results     # Répertoire créé avec nombre de points par défaut
│   │   └── results.txt
│   └── results_[N]         # Répertoire créé en utilisant N points
│       └── results.txt
└── src/
    ├── sphere.py           # Parsing PDB + génération des points de sphère
    ├── calcul.py           # Calculs de distance, voisinage et accessibilité
    └── prog_princ.py       # Script principal

```

# Description des modules 

### `src/sphere.py`
- **`parsepdb(file_name)`** : lit un fichier PDB et extrait les atomes (identifiant, nom, résidu, chaîne, coordonnées), en ignorant les atomes d'hydrogène.
- **`saff_kuijlaars_sphere(radius, n_points, center)`** : génère `n_points` répartis uniformément sur une sphère de rayon et centre donnés.
- **`get_vdw_radius(atom_name, rayons_dict)`** : renvoie le rayon de van der Waals d'un atome à partir d'un dictionnaire de rayons.
- **`generate_protein_spheres(atoms, probe_radius, rayons_dict, n_points)`** : calcule le rayon étendu (van der Waals + sonde) de chaque atome et génère les points de sa sphère.

### `src/calcul.py`
- **`distance(atom1, atom2)`** : distance euclidienne entre deux points 3D.
- **`near_neighbor(atoms)`** : identifie, pour chaque atome, les atomes voisins dont les sphères étendues se chevauchent (optimisation du calcul).
- **`accessibilite_point(atoms, nbpoints)`** : détermine, pour chaque atome, le nombre et le pourcentage de points de sa sphère qui sont accessibles au solvant.
- **`surface_accessible(atoms, nbpoints)`** : calcule la SASA totale de la protéine, ainsi que le pourcentage de surface et le nombre total de points accessibles.

### `prog_princ.py`
Script principal qui :
1. Définit les rayons de van der Waals par type d'atome.
2. Parcourt tous les fichiers `.pdb` du dossier `data/`.
3. Calcule la SASA de chaque structure.
4. Écrit les résultats dans `results/[chemin_vers_résultats]/results.txt`.

## Prérequis

- Python v.3.12.2

## Exécution

1. Placer vous à la racine du projet.
2. Lancer le script principal :
```bash
python src/prog_princ.py
```
3. Les résultats sont enregistrés dans `results/[chemin_vers_résultats]/results.txt`, avec :
   - une ligne d'en-tête indiquant le nombre de points utilisés par sphère ;
   - pour chaque structure : le nombre total de points accessibles, le pourcentage de points accessibles, la surface totale accessible (en Å²) et son temps d'exécution ;
   - le temps d'exécution total nécessaire pour traiter toutes les molécules.


### Exemple de fichier de sortie

```
Avec 92 points par sphère : 

Pour 168L : 
Nombre total d'atomes : 6445 
Nombre total de points accessibles : 29702 
Pourcentage de surface accessible : 5.09 %
Surface totale accessible de la protéine : 38751.07 Å²
Temps d'exécution : 13.33 secondes

Pour 1CRN : 
Nombre total d'atomes : 327 
Nombre total de points accessibles : 2266 
Pourcentage de surface accessible : 7.84 %
Surface totale accessible de la protéine : 3014.73 Å²
Temps d'exécution : 0.16 secondes

[...]

Temps d'exécution total : 17.86 secondes
```

## Paramètres modifiables

Dans `prog_princ.py` :
- **`rayons`** : dictionnaire des rayons de van der Waals par type d'atome (Å).
- **`rayon_probe`** : rayon de la sonde (molécule de solvant), 1.4 Å par défaut.
- **`nb_points`** : nombre de points générés par sphère atomique (92 par défaut). Plus ce nombre est élevé, plus la précision augmente, au prix d'un temps de calcul plus long. Vous pouvez décommentez la ligne 18 pour modifier le nombre de points.

## Limites 

- Les atomes dont le nom n'est ni présent dans `rayons_dict` ni couvert par son premier caractère lèveront une `KeyError` dans `get_vdw_radius`.
- Le calcul de voisinage (`near_neighbor`) est de complexité O(n²) en nombre d'atomes, ce qui peut être lent sur de grosses structures.
