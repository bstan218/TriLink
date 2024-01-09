import re
import RotateMolecule
import statistics as s
import math as m

class MolFile:
    def __init__(self, filepath, template=False, coreidentity=None):
        self.filepath = filepath
        self.filename = filepath.split('/')[-1].strip('.mol')
        #if re.match(r'/', self.filename):
            #self.filename = re.match                
        self.parse()
        self.find_connection_coords()
        self.find_core_coordinate()
        self.center_molecule()
        self.istemplate = template
        self.rotated_coordinates= None
        self.alt_rotated_coordinates = None
        if template:
            self.coreidentity = coreidentity
            self.find_core()
    
    def parse(self):
        with open(self.filepath, 'r') as inf:
            fi = []
            for line in inf:
                fi.append(line.strip('\n').strip('\t').strip(' '))
        line_lst=[]
        for i in fi:
            line_lst.append(i.split(None))
        self.connection_atoms = None
        try: 
            line_lst[0][0]
        except:
            pass
        else:
            if re.match(r'\d+,\d+,\d+',line_lst[0][0]):
                self.connection_atoms = list(map(int, line_lst[0][0].split(',')))
        self.charge_lst = line_lst[1]
        
        self.atoms = []
        for line in line_lst:
            if len(line) == 16:
                self.atoms.append(line)
        
        self.coordinates = []
        for atom in self.atoms:
            self.coordinates.append([float(atom[0]),float(atom[1]),float(atom[2])])

    def find_connection_coords(self):
        c1, c2, c3 = self.connection_atoms[0:3]
        l1 = list(map(float,self.atoms[c1-1][0:3]))
        l2 = list(map(float, self.atoms[c2-1][0:3]))
        l3 = list(map(float, self.atoms[c3-1][0:3]))
        l1.append(c1)
        l2.append(c2)
        l3.append(c3)
        self.connection_coordinates = [l1,l2,l3]
    
    def find_core_coordinate(self):
        #find which coords are furthest apart
        c1, c2, c3 = self.connection_coordinates[0:3]
        
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

        #find midpoint between those two points
        self.midpoint =  [s.mean([f[0][0],f[1][0]]), s.mean([f[0][1],f[1][1]]), s.mean([f[0][2],f[1][2]])]

    def center_molecule(self):
        self.centered_atoms = []
        for point in self.atoms:
            newX = round((float(point[0])-self.midpoint[0]),12)
            newY = round((float(point[1])-self.midpoint[1]),12)
            newZ = round((float(point[2])-self.midpoint[2]),12)
            self.centered_atoms.append([newX,newY,newZ])
    
    def find_core(self):
        for i in range(0,len(self.atoms)):
            if self.coreidentity in self.atoms[i]:
                self.coreindex = i+1
                return


    



if __name__ == '__main__':
    fi = "ligands/CCC-I.mol"
    lig = MolFile(fi)
    #print(lig.filename)