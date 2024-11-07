from Parser import ParsedFile, TLCode
import os
import configparser
import time
import OptimizeValidate
import RotateMolecule
import runxtb
import subprocess

def create_xyz_file(template_file, ligand_file, outdir, temporaryfile=False):
    templatename = template_file.filename
    outfname = f'{templatename}.xyz'
    outfpath = f'{outdir}/{outfname}'

    if temporaryfile:
        outfname = 'temp.xyz'
        outfpath = outfname
    
    
    with open(outfpath, 'w') as outf:

        temp_coords = template_file.optimized_coordinates.coords  
        lig_coords = ligand_file.rotated_coordinates.coords

        outf.write(f'{len(temp_coords)+len(lig_coords)-3}\n')
        outf.write('\n')
     
        for i in range(0,len(temp_coords)):
            if i+1 not in template_file.connection_atoms: 
                atomtype = template_file.atomtypes[i]
                x, y, z = temp_coords[i][0:3]
                outf.write(f'{atomtype}      {x}   {y}   {z}\n')
        
        for i in range(0,len(lig_coords)):
            atomtype = ligand_file.atomtypes[i]
            x, y, z = lig_coords[i][0:3]
            outf.write(f'{atomtype}      {x}   {y}   {z}\n')   

        outf.write('\n')

def create_com_file(template_file, ligand_file, outdir, xtbopt = None):
    tempname = template_file.filename
    outfname = f'{tempname}.gjf'


    with open(f'{outdir}/{outfname}', 'w') as outf:
        for line in template_file.header:
            outf.write(f'{line}\n')
        
        outf.write(f'\n')
        outf.write(f'{template_file.title[0]}\n')
        outf.write(f'\n')


        outf.write(f'{ligand_file.charge + template_file.charge} {template_file.spinmultiplicity}\n')
        
        if not xtbopt:
            temp_coords = template_file.optimized_coordinates.coords       
            for i in range(0,len(temp_coords)):
                if i+1 not in template_file.connection_atoms: 
                    atomtype = template_file.atomtypes[i]
                    x, y, z = temp_coords[i][0:3]
                    outf.write(f'{atomtype}              {x}   {y}   {z}\n')
            
            lig_coords = ligand_file.rotated_coordinates.coords
            for i in range(0,len(lig_coords)):
                atomtype = ligand_file.atomtypes[i]
                x, y, z = lig_coords[i][0:3]
                outf.write(f'{atomtype}              {x}   {y}   {z}\n')

        else:
            xtb_coords = xtbopt.atoms
            for i in range(0, len(xtb_coords)):
                atomtype = xtb_coords[i][0]
                x, y, z = xtb_coords[i][1:4]
                outf.write(f'{atomtype}            {x}     {y}     {z}\n')

        outf.write('\n')

        TLCodefinder = TLCode()
        for section in template_file.sections[3:]:
            for line in section:
                if TLCodefinder.check_for_code(line) != None:
                    line = TLCodefinder.replace_code(line, ligand_file, template_file)
                outf.write(f'{line}\n')
            
            outf.write('\n')

def Init(configFile):
    config = configparser.ConfigParser()
    config.read(configFile)
    LigandLib = config['IODir']['LigandLibDir']
    TemplateDir = config['IODir']['TemplateDir']
    coreidentity = config['Other']['CoreIdentity']
    xtbopt = config['Other']['xtbopt']
    runscript = config['Other']['runscript']
    script = config['Other']['script']
    if xtbopt == "True":
        xtbopt = True
    else:
        xtbopt = False
    if runscript == "True":
        runscript = True
    else:
        runscript = False    
    return LigandLib, TemplateDir, coreidentity, xtbopt, runscript, script



def main():
    #assign directories and files, make directories
    configFile = "config.ini"
    LigandLib, TemplateDir, coreidentity, xtbopt, runscript, script = Init(configFile)
    d = time.localtime()
    OutDirLabel = LigandLib.replace('/','_')
    xtblabel=""
    if xtbopt:
        xtblabel="xtb_"
    outdir = f'Output_{OutDirLabel}{xtblabel}{d[1]}{d[2]}{d[0]}.{d[3]}{d[4]}{d[5]}'
    validligdir = f'{outdir}/valid_ligands/'
    invalidligdir = f'{outdir}/invalid_ligands/'
    os.mkdir(outdir)
    os.mkdir(validligdir)
    os.mkdir(invalidligdir)

    #initialize template list
    template_lst = []
    for template in os.listdir(TemplateDir):
        temp_filepath = f'{TemplateDir}{template}'
        template = ParsedFile(temp_filepath, template=True, coreidentity=coreidentity)
        nonaltcoords, altcoords = RotateMolecule.init_rotate_molecule(template)
        template.set_rotated_coordinates(nonaltcoords)
        template.set_alt_rotated_coordinates(altcoords)
        template_lst.append(template)
    
    #main program loop
    for lig in os.listdir(LigandLib):
        print(lig)
        lig_filepath = f'{LigandLib}{lig}'
        ligand = ParsedFile(lig_filepath)
        if ligand.missing_data:
            continue
        
        valid, invalidations = OptimizeValidate.main_optimization(ligand, template_lst)

        #rotamer optimization goes here
        if not valid:
            pass
            #OptimizeValidate.rotamer_optimization(ligand, template_lst)

        #create .com directory for ligand
        outdir = (f'{invalidligdir}{ligand.filename}/')
        if valid:
            outdir = (f'{validligdir}{ligand.filename}/')
        
        os.mkdir(outdir)
        #create .com files
        for template in template_lst:
            if valid and xtbopt:
                create_xyz_file(template,ligand,outdir,temporaryfile=True)
                #xtb optimization
                #define connection atoms
                freezeatoms = [f"1-{len(template.atoms)-3}"]
                for atom in ligand.connection_atoms:
                    freezeatoms.append(str(len(template.atoms)+int(atom)-3))

                runxtb.xtb_optimize('temp.xyz', freezeatoms)
                xtbopt = ParsedFile('xtbopt.xyz')
                create_com_file(template,ligand,outdir,xtbopt)
            else:
                create_com_file(template,ligand,outdir)
                        
        
        if invalidations:
            with open(f'{outdir}invalidations.txt', 'w') as outf:
                for i in invalidations:
                    outf.write(f'{i}\n')
    
    if runscript:
        #for folder in 
        for dir in os.listdir(validligdir):
            os.chdir(f"./{validligdir}{dir}")
            subprocess.run(script.split(), shell=True)
    

if __name__ == "__main__":
    main()
    print('operation completed')