import numpy as np
import RotateMolecule

class BondRadii:
    def __init__(self) -> None:
        self.radii = {'H':0.37, 'C':0.77, 'N':0.7, 'O':0.66, 'F':0.64, 'Al':1.43, 'Si':1.77,
                      'P':1.1,'S':1.04,'Cl':0.99}
    
    def find_bond_length(self, a1, a2):
        a1d = 1.5
        if a1 in self.radii.keys():
            a1d = self.radii[a1]
        a2d = 1.5
        if a2 in self.radii.keys():
            a2d = self.radii[a2]
        
        return (a1d + a2d)

def compare_coords(coords1, coords2) -> list:
    exceptions_lst = []
    br = BondRadii()
    
    for i in range(0,len(coords1.coords)):
        for j in range(0,len(coords2.coords)):
            if i==33 and j==10:
                pass
            #check length between ligand[i] and template[j]
            xdif = abs(coords1.coords[i][0]- coords2.coords[j][0])
            ydif = abs(coords1.coords[i][1]- coords2.coords[j][1])
            zdif = abs(coords1.coords[i][2]- coords2.coords[j][2])
            dif = np.sqrt(xdif**2+ydif**2+zdif**2)
            
            threshold = -0.3 + br.find_bond_length(coords1.identities[i], coords2.identities[j])
            
            if dif < threshold:
                exceptions_lst.append([i+1,j+1,dif,threshold-dif])

    return exceptions_lst

def validate(coords1, coords2) -> bool:
    exceptions_lst = compare_coords(coords1, coords2)

    #remove any exceptions that include connection atoms
    to_remove = []
    for e in exceptions_lst:
        if e[0] in coords1.ignorelist or e[1] in coords2.ignorelist:
            to_remove.append(e)
    for e in to_remove:
        exceptions_lst.remove(e)
    
    validity = False
    #valid case
    if exceptions_lst == []:
        validity = True
    
    interference_val = 0
    for e in exceptions_lst:
        interference_val += e[3]

    #invalid case
    return validity, exceptions_lst, interference_val

def altornot(template_lst, ligand):
    norm_validity = True
    alt_validity = True
    norm_interference = 0
    alt_interference = 0
    for template in template_lst:
       
        validity1, invalidations1, interference1 = validate(ligand.rotated_coordinates, template.rotated_coordinates)
        validity2, invalidations2, interference2 = validate(ligand.rotated_coordinates, template.alt_rotated_coordinates)
        if interference1:                
            norm_interference += interference1
        if interference2:
            alt_interference += interference2
        if not validity1:
            norm_validity = False
        if not validity2:
            alt_validity = False
    
    #determine which initial rotation is more optimal
    #base case: both rotations valid
    validity = True #reflects if there are no exceptions in the set of coords given by alt
    alt = False
    if not norm_validity: #first rotation invalid
        alt = True
        
        if not alt_validity: #second rotation invalid
            validity = False
            #if second rotation is invalid, need to compare which coords have higher interference
            if norm_interference < alt_interference:
                alt = False
    
    return validity, alt

def rock_optimize_molecule(ligand, template) -> bool:
    #perform rock rotate on
    lowest_interference = validate(ligand.rotated_coordinates, template.optimized_coordinates)[2]
    placeholder_coords=None

    for t in range(5, 60, 1):
        for theta in range(-t, t+1, 2*t):
            if theta == 24:
                pass
            rocked_coords = RotateMolecule.rock_rotate(template, theta*np.pi/180)
            template.set_rocked_coordinates(rocked_coords)
            interference = validate(ligand.rotated_coordinates, template.rocked_coordinates)[2]
            if not interference:
                template.optimized_coordinates.coords = rocked_coords
                return True
            elif interference < lowest_interference:
                lowest_interference = interference
                placeholder_coords = rocked_coords
    if placeholder_coords:
        template.optimized_coordinates.coords = placeholder_coords
    return False

def twist_optimize_molecule(ligand, template):
    pass

def optimize_molecule(ligand, template):
    optimized = rock_optimize_molecule(ligand, template)
    if not optimized:
        optimized = twist_optimize_molecule(ligand, template)
    if not optimized:
        return False
    return True

def main_optimization(ligand, template_lst):
    ligand.set_rotated_coordinates(RotateMolecule.init_rotate_molecule(ligand)[0])
    
    either_alt_valid, alt = altornot(template_lst, ligand)

    #determine which templates need optimization
    fully_optimized = True
    invalidations = []
    for template in template_lst:

        template.set_optimized_coordinates(template.rotated_coordinates.coords)
        if alt:
            template.set_optimized_coordinates(template.alt_rotated_coordinates.coords)

        if not either_alt_valid: #continues if neither alternatives are already valid
            individual_validity = validate(ligand.rotated_coordinates, template.optimized_coordinates)[0]
            if not individual_validity:
                optimized = optimize_molecule(ligand, template)
                if not optimized:
                    fully_optimized = False
                    #create post-optimization invalidations list
                    for i in validate(ligand.rotated_coordinates, template.optimized_coordinates)[1]:
                        invalidations.append([template.filename,i])
    
    return fully_optimized, invalidations

def rotamer_optimization():
    pass