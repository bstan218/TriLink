from Parser import ParsedFile
#create a bond class that has 
class Bond:
    def __init__(self, elements) -> None:
        self.elements = [eval(i) for i in elements]
        self.connections = self.elements[0:2]
        self.multiplicity = self.elements[2]

def IntBonds(molecule):
    bonds = []
    for bond in molecule.bonds:
        b = Bond(bond)
        bonds.append(b)
    return bonds

class Atom:
    def __init__(self, index) -> None:
        self.index = index
        self.bonds = {}

def Int(molecule):
    atoms = []
    for i in range(1,1+len(molecule.atoms)):
        atoms.append(Atom(i))

    bondlist = IntBonds(molecule)
    for bond in bondlist:
        # add bond to atom.bond_list, where atom is determined by the bond.atom1 and bond.atom2
        atom1 = atoms[int(bond.connections[0])-1]
        atom2 = atoms[int(bond.connections[1])-1]
        atom1.bonds[atom2.index] = bond.multiplicity
        atom2.bonds[atom1.index] = bond.multiplicity

    return atoms, bondlist

#Find all useful rotatable bonds
#useful rotatable bonds defined as a bond with:
# 1) multiplicity of 1
# 2) each atom in bond has at least one additonal connection
# 3) rotating bond doesn't unlink any other bonds in the structure
    #How do we determine this?
    #need to follow all bonds along "bondpath" from one atom in tested bond
    #if any bonds along bondpath reach back to  the other atom in tested bond, then the tested bond is invalid.

def findbond(bondlist, i1, i2):
    found = None
    for bond in bondlist:
        if i1 in bond.connections and i2 in bond.connections:
            found = bond

    return found

def provebond(bond, atomlist) -> bool:
    if bond.multiplicity > 1:
        return False

    atom1 = atomlist[int(bond.connections[0])-1]
    atom2 = atomlist[int(bond.connections[1])-1]
    if not len(atom1.bonds) > 1 or not len(atom2.bonds) > 1: #checks if both atoms have at least one additional connection
        return False
    
    #determines if given bond is able to rotate without unlinking any other bonds in structure
    def rotateable(atom, testatom):
        
        def helper_rotateable(atom, testatom, visited, validity):
            visited.append(atom.index)
            #print(visited)
            #print(validity)
            bonds = atom.bonds
            for a in bonds.keys():
                if a == testatom.index:
                    validity = False
                    break
                if not a in visited:
                    #print(a)
                    validity = helper_rotateable(atomlist[a-1], testatom, visited, validity)
                    
            return validity

        for i in atom.bonds.keys():
            if i != testatom.index:
                if not helper_rotateable(atomlist[i-1], testatom, [atom.index, testatom.index], True):
                    return False
        
        return True
  
    if not rotateable(atom1, atom2):
        return False
    
    return True


def main():
    inf = 'Tridentate\SNN\ABOBIE_Pt_1_3_NNS.mol'
    mol = ParsedFile(inf)
    atoms, bondlist = Int(mol)

    for bond in bondlist:
        if provebond(bond,atoms):
            print(bond.connections)


    #follow bondpath






######################################################
if __name__ == '__main__':
    main()
