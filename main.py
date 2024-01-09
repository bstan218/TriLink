from Parser import MolFile
import RotateMolecule as rm
import os
import configparser
import numpy as np
import time

def rotate_molecule(ligand):
    
    def helper(sideindex):
        #define inputs for first rotation
        p11 = ligand.centered_atoms[sideindex]
        p1mag = np.linalg.norm(p11)
        p12 = np.array([p1mag,0,0])
        v1 = rm.FindRotationVector(p11,p12)
        t1 = rm.FindTheta(p11,p12)
        #first rotation
        first_rotated_atoms = []
        for atom in ligand.centered_atoms:
            first_rotated_atoms.append(rm.RotatePoint(atom,v1,t1))
        
        #define inputs for second rotation
        centerindex = ligand.centerconnection[3]-1
        p21= first_rotated_atoms[centerindex]
        p2mag = np.linalg.norm(p21)
        p22 = [p21[0], 0.0, np.sqrt(p2mag**2-p21[0]**2)]
        v2 = rm.FindRotationVector(p21,p22)
        t2 = rm.FindTheta(p21,p22)
        #second rotation
        final_rotated_atoms = []
        for atom in first_rotated_atoms:
            final_rotated_atoms.append(rm.RotatePoint(atom, v2,t2))
        
        return final_rotated_atoms

    sideindex = ligand.sideconnections[0][3]-1
    altsideindex = ligand.sideconnections[1][3]-1
    rotated_atoms = helper(sideindex)
    alt_rotated_atoms = helper(altsideindex)

    return rotated_atoms, alt_rotated_atoms

def trifilter(ligand, template, threshold, alt=False) -> bool:
    ligand_coords = ligand.rotated_coordinates
    if alt:
        ligand_coords = ligand.alt_rotated_coordinates
    template_coords = template.coordinates

    below_threshold = []
    min_dif = 10
    for i in range(0,len(ligand_coords)):
        for j in range(0,len(template_coords)):
            #check length between ligand[i] and template[j]
            if i+1 not in ligand.connection_atoms and j+1 not in template.connection_atoms:
                xdif = abs(ligand_coords[i][0]- template_coords[j][0])
                ydif = abs(ligand_coords[i][1]- template_coords[j][1])
                zdif = abs(ligand_coords[i][2]- template_coords[j][2])
                dif = np.sqrt(xdif**2+ydif**2+zdif**2)
                if dif < threshold:
                    below_threshold.append(f'Distance between Ligand Atom {i+1} and Template Atom {j+1}: {dif}')
                if dif < min_dif:
                    min_dif = dif


    if below_threshold == []:
        return True, None, min_dif
    return False, below_threshold, min_dif

def create_com_file(template_file, ligand_file, outdir, alt, headers, TS):
    tempname = template_file.filename
    outfname = f'{tempname}.com'
    non_ts_head, ts_head, ts_tail = headers[0:3]
    with open(f'{outdir}/{outfname}', 'w') as outf:
        if TS:
            with open(ts_head,'r') as header:
                for line in header:
                    outf.write(line)
        else:
            with open(non_ts_head,'r') as header:
                for line in header:
                    outf.write(f'{line}')
        outf.write(f' \n')
        outf.write(f'0 1\n')
        
        for i in range(0,len(template_file.rotated_coordinates)):
            if i+1 not in template_file.connection_atoms: 
                atomtype = template_file.atoms[i][3]
                x, y, z = template_file.rotated_coordinates[i][0:3]
                outf.write(f'{atomtype}              {x}   {y}   {z}\n')
        
        lig_coords = ligand_file.rotated_coordinates
        if alt:
            lig_coords = ligand_file.alt_rotated_coordinates
        for i in range(0,len(lig_coords)):
            atomtype = ligand_file.atoms[i][3]
            x, y, z = lig_coords[i][0:3]
            outf.write(f'{atomtype}              {x}   {y}   {z}\n')   

        outf.write('\n')


        if TS:
            #append bonds
            for connection in ligand_file.connection_atoms:
                newconnectionindex = int(connection) + len(template_file.rotated_coordinates) - 3
                adj = 0
                for i in template_file.connection_atoms:
                    if i < template_file.coreindex:
                        adj += 1
                newcoreindex = template_file.coreindex - adj
                outf.write(f'B {newcoreindex} {newconnectionindex} F\n')

            outf.write('\n')
            
            with open(ts_tail,'r') as footer:
                for line in footer:
                    outf.write(line)     
        
        
        outf.write('\n')

def Init(configFile):
    config = configparser.ConfigParser()
    config.read(configFile)
    LigandLib = config['IODir']['LigandLibDir']
    TemplateDir = config['IODir']['TemplateDir']
    HeadingsDir = config['ComFile']['HeadingsDir']
    non_ts_head = config['ComFile']['non_ts_heading']
    ts_head = config['ComFile']['ts_heading']
    ts_tail = config['ComFile']['ts_tail']
    threshold = float(config['Other']['Threshold'])
    headings = [f'{HeadingsDir}{non_ts_head}', f'{HeadingsDir}{ts_head}', f'{HeadingsDir}{ts_tail}']
    filtertemplatefile = config['Other']['filtertemplate']
    coreidentity = config['Other']['CoreIdentity']
    return LigandLib, TemplateDir, headings, threshold, filtertemplatefile, coreidentity

def main():
    #assign directories and files, make directories
    configFile = "config.ini"
    LigandLib, TemplateDir, headings, threshold, filtertemplatefile, coreidentity = Init(configFile)
    d = time.localtime()
    outdir = f'Output_{d[1]}{d[2]}{d[0]}{d[3]}{d[4]}{d[5]})'
    validligandsdir = f'{outdir}/valid_ligands/'
    invalidligandsdir = f'{outdir}/invalid_ligands/'
    validcomdir = f'{outdir}/com_files/'
    os.mkdir(outdir)
    os.mkdir(validligandsdir)
    os.mkdir(invalidligandsdir)
    os.mkdir(validcomdir)
    filtertemplate = MolFile(f'{TemplateDir}{filtertemplatefile}')

    for lig in os.listdir(LigandLib):
        if lig.endswith('.mol'):
            lig_filepath = f'{LigandLib}{lig}'
            ligand = MolFile(lig_filepath)
            ligand.rotated_coordinates, ligand.alt_rotated_coordinates = rotate_molecule(ligand)
            
            validity1, invalidations1, min_proximity1 = trifilter(ligand, filtertemplate, threshold)
            validity2, invalidations2, min_proximity2 = trifilter(ligand, filtertemplate, threshold, alt=True)
            #base case: both rotations valid
            validity = True
            alt = False
            invalidations = invalidations1
            if not validity1: #first rotation invalid
                alt = True
                invalidations = invalidations2
                if not validity2: #second rotation invalid
                    validity = False
            
            liganddir = invalidligandsdir
            if validity:
                liganddir = validligandsdir

                #create .com directory for ligand
                outdir = (f'{validcomdir}/{ligand.filename}')
                os.mkdir(outdir)
                #create .com files
                for template in os.listdir(TemplateDir):
                    if template != filtertemplatefile:
                        temp_filepath = f'{TemplateDir}{template}'
                        template = MolFile(temp_filepath, template=True, coreidentity=coreidentity)
                        TS = False
                        if template.filename.startswith('TS'):
                            TS = True
                        template.rotated_coordinates = rotate_molecule(template)[0]

                        create_com_file(template,ligand,outdir,alt,headings,TS)
            

            with open(f'{liganddir}{ligand.filename}.mol', 'w') as outf:
                with open(ligand.filepath, 'r') as inf:
                    for line in inf:
                        outf.write(line)
            if not validity:
                with open(f'{liganddir}{ligand.filename} Exceptions.txt', 'w') as outf:
                    for exception in invalidations:
                        outf.write(f'{exception}\n')

if __name__ == "__main__":
    main()