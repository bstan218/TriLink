import statistics as s
import math as m

class coordsobj:
    def __init__(self, coords) -> None:
        self.coords = coords

    def assign_identities(self, identities):
        self.identities = identities
    
    def assign_ignorelist(self, ignorelist):
        self.ignorelist = ignorelist
    
    def __repr__(self) -> str:
        pass

class Bond:
    def __init__(self, elements) -> None:
        self.elements = [eval(i) for i in elements]
        self.connections = self.elements[0:2]
        self.multiplicity = self.elements[2]

class Atom:
    def __init__(self, index) -> None:
        self.index = index
        self.bonds = {}

class ParsedFile:
    def __init__(self, filepath, template=False, coreidentity=None, rotamer_coords=None):
        self.filepath = filepath
        self.filetype = filepath.split('.')[-1]
        self.filename = filepath.split('/')[-1].strip(f'.{self.filetype}')       
        self.rotamer_coords=rotamer_coords
        self.missing_data = False

        self.isrotamer=False
        if rotamer_coords:
            self.isrotamer=True

        if self.filetype == 'mol':
            self.parse_mol()
        elif self.filetype == 'com':
            self.parse_com()
        elif self.filetype == 'xyz':
            self.parse_xyz()
        else:
            raise TypeError(f'file type .{self.filetype} not parseable')
        
        if self.filetype != 'xyz':
            if len(self.connection_atoms) < 3:
                self.missing_data = True
            else:
                self.find_connection_coords()
                self.find_core_coordinate()
                self.center_molecule()
                self.istemplate = template
                if template:
                    self.coreidentity = coreidentity
                    self.find_core()
    
    def parse_xyz(self):
        with open(self.filepath, 'r') as inf:
            fi = []
            for line in inf:
                fi.append(line.strip('\n').strip('\t').strip(' '))
        line_lst=[]
        for i in fi:
            line_lst.append(i.split(None))
        
        self.atoms = []
        for line in line_lst[2:]:
            self.atoms.append(line)
               
    
    def parse_mol(self):
        with open(self.filepath, 'r') as inf:
            fi = []
            for line in inf:
                fi.append(line.strip('\n').strip('\t').strip(' '))
        line_lst=[]
        for i in fi:
            line_lst.append(i.split(None))

        self.connection_atoms = list(map(int, line_lst[0][0].split(',')))
        self.dentisity = len(self.connection_atoms)

        self.charge_lst = line_lst[1]
        self.charge = 0
        if len(self.charge_lst) >= 3:
            self.charge = round(float(self.charge_lst[2].strip(',')))
        
        self.atoms = []
        self.bonds = []

        for line in line_lst:
            if len(line) == 16:
                self.atoms.append(line)
            if len(line) == 7:
                self.bonds.append(line[0:3])
            
        self.coordinates = []
        self.atomtypes = []
        for atom in self.atoms:
            self.coordinates.append([float(atom[0]),float(atom[1]),float(atom[2])])
            self.atomtypes.append(atom[3])

    def parse_com(self):
        with open(self.filepath, 'r') as inf:
            fi = []
            for line in inf:
                fi.append(line.strip('\n').strip('\t').strip(' '))
        line_lst=[]
        for i in fi:
            line_lst.append(i)
        #print(line_lst)
        #parse line list using empty lists as separators between header, connection list, coordinates, and optional footer.
        self.header = []
        self.title = []
        atomsincludingcharge = []
        self.sections = [self.header, self.title, atomsincludingcharge]
        i = 0
        for line in line_lst:
            #if line is  not '', append it to current list.
            #if line is '', increment i IF current list has something
            if line != '':
                self.sections[i].append(line)
            elif self.sections[i] != []:
                i += 1
                if i >= len(self.sections):
                    self.sections.append([])

        self.atoms = []
        self.coordinates = []
        self.atomtypes = []
        for atom in atomsincludingcharge:
            lst = (atom.split(None))
            if len(lst) < 4:
                self.charge, self.spinmultiplicity = lst[0:2]
                self.charge = int(self.charge)
                self.spinmultiplicity = int(self.spinmultiplicity)
            else:
                self.atoms.append(lst)
                self.coordinates.append([float(lst[1]), float(lst[2]), float(lst[3])])
                self.atomtypes.append(lst[0])
        self.connection_atoms = [len(self.atoms), len(self.atoms)-1, len(self.atoms)-2]
        #reformat coordinates

    def find_connection_coords(self):
        c1, c2, c3 = self.connection_atoms[0:3]
        l1 = list(map(float,self.coordinates[c1-1][0:3]))
        l2 = list(map(float, self.coordinates[c2-1][0:3]))
        l3 = list(map(float, self.coordinates[c3-1][0:3]))
        l1.append(c1)
        l2.append(c2)
        l3.append(c3)
        self.connection_coordinates = [l1,l2,l3]
    
    def find_core_coordinate(self):
        #find which coords are furthest apart
        c1, c2, c3 = self.connection_coordinates[0:3]
        coordinates_indexes = dict(zip([tuple(c) for c in self.connection_coordinates], self.connection_atoms))

        def find_dif(c1,c2):
            xdif = abs(c1[0]-c2[0])
            ydif = abs(c1[1]-c2[1])
            zdif = abs(c1[2]-c2[2])
            dif = m.sqrt(xdif**2+ydif**2+zdif**2)
            return dif
        
        d1 = find_dif(c1,c2)
        d2 = find_dif(c1,c3)
        d3 = find_dif(c2,c3)
        if d1 > d2 and d1 > d3:
            f = [c1,c2,c3]
        elif d2 > d3:
            f = [c1,c3,c2]
        else:
            f = [c2,c3,c1]

        self.sideconnections = f[0:2]
        self.centerconnection = f[2]
        self.sideconnections_indexes = [coordinates_indexes[tuple(c)] for c in self.sideconnections]
        self.centerconnection_index = coordinates_indexes[tuple(self.centerconnection)]

        #find midpoint between those two points
        self.midpoint =  [s.mean([f[0][0],f[1][0]]), s.mean([f[0][1],f[1][1]]), s.mean([f[0][2],f[1][2]])]

    def get_proper_connection_type(self):
        self.find_core_coordinate()
        s1, s2 = [self.atomtypes[i-1] for i in self.sideconnections_indexes]
        center_atomtype = self.atomtypes[self.centerconnection_index-1]

        if s1 < s2:
            return s1 + center_atomtype + s2
        else:
            return s2 + center_atomtype + s1
        
        

    def center_molecule(self):
        self.centered_atoms = []
        for point in self.coordinates:
            newX = round((float(point[0])-self.midpoint[0]),12)
            newY = round((float(point[1])-self.midpoint[1]),12)
            newZ = round((float(point[2])-self.midpoint[2]),12)
            self.centered_atoms.append([newX,newY,newZ])
    
    def find_core(self):
        for i in range(0,len(self.atoms)):
            if self.coreidentity in self.atoms[i]:
                self.coreindex = i
                return
            
    def get_atom_types_set(self):
        retset = set(self.atomtypes)
        return retset
        
    def make_coords(self, coords):
        return_coords = coordsobj(coords)
        return_coords.assign_identities(self.atomtypes)
        return_coords.assign_ignorelist(self.connection_atoms)
        return return_coords
    
    def set_rotated_coordinates(self, coords):
        self.rotated_coordinates = self.make_coords(coords)
    
    def set_alt_rotated_coordinates(self, coords):
        self.alt_rotated_coordinates = self.make_coords(coords)

    def set_optimized_coordinates(self, coords):
        self.optimized_coordinates = self.make_coords(coords)

    def set_rocked_coordinates(self,coords):
        self.rocked_coordinates = self.make_coords(coords)
    
    def init_bonds(self):
        self.bond_atoms = []
        for i in range(1,1+len(self.atoms)):
            self.bond_atoms.append(Atom(i))

        self.bondlist = []
        for bond in self.bonds:
            b = Bond(bond)
            self.bondlist.append(b)

        for bond in self.bondlist:
            # add bond to atom.bonds, where atom is determined by the bond.atom1 and bond.atom2
            atom1 = self.bond_atoms[int(bond.connections[0])-1]
            atom2 = self.bond_atoms[int(bond.connections[1])-1]
            atom1.bonds[atom2.index] = bond.multiplicity
            atom2.bonds[atom1.index] = bond.multiplicity

    def provebond(self, bond) -> bool:
        atomlist = self.bond_atoms
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
    
    def find_useful_bonds(self):
        self.init_bonds()
        self.useful_bonds = []
        for bond in self.bondlist:
            if self.provebond(bond):
                self.useful_bonds.append(bond.connections)

class TLCode:
    def __init__(self) -> None:
        self.codes = ["corelesslist", "dummy code"]

    def replace_code(self, line, ligand=None, template=None):
        code = self.check_for_code(line)
        if code == "corelesslist":
            atom_list = ligand.get_atom_types_set() | template.get_atom_types_set()
            atom_list.remove(template.coreidentity)
            line = line.replace("t[corelesslist]", ' '.join(atom_list))
        
        return line

    def check_for_code(self, line):
        for code in self.codes:
            if f"t[{code}]" in line:
                return code
        return None
    

        #takes in line of code, determines if the code needs 



if __name__ == '__main__':
    xyzfi = ParsedFile("xtbopt.xyz")

    