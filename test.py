rayons = {'N' : 1.5, 'O': 1.4, 'S' : 1.85, 'CG' : 1.85, 'CD1' : 1.85, 'CD2' : 1.85,
           'CE1' : 1.85, 'CE2' : 1.85, 'CE3' : 1.85, 'CZ' : 1.85, 'CZ2' : 1.85,
            'CZ3' : 1.85, 'CH2' : 1.85, 'CA' : 1.5 , "CB" : 1.5, 'C' : 2.0}

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

                # # Symbole de l'élément (ex: C, N, O, S)
                # element = line[76:78].strip()
                # if not element:  # Si la colonne est vide, déduire du nom de l'atome
                #     element = atom_name[0]

                atoms.append({
                    "id": atom_id,
                    "name": atom_name,
                    "res_name": res_name,
                    "chain": chain_id,
                    "res_id": res_seq,
                    "coords": (x, y, z),
                    # "element": element
                })
    return atoms

atoms = parsepdb("3I40.pdb")
print(atoms[0])

for atom in atoms : 
    print(atom)
   #print(atom["name"], ":", atom["coords"])

print(len(atoms))
# print(rayons[atoms[2]['name']])

# for atom in atoms:
#     if atom['name'].startswith("S") : 
#         print(rayons["S"])