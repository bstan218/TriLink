from Parser import ParsedFile, TLCode
import os
import sys
import argparse
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

def initial_config():
    CONFIG_FILE = "config.ini"
    parser = configparser.ConfigParser()
    parser.read(CONFIG_FILE)

    configurations = {}
    
    configurations['ligands'] = parser['IODir']['LigandLibDir']
    configurations['templates'] = parser['IODir']['TemplateDir']
    configurations['core'] = parser['Other']['CoreIdentity']
    configurations['xtb'] = parser['Other']['xtbopt'] == 'True'
    configurations['runscript'] = parser['Other']['runscript'] == 'True'
    configurations['script'] = parser['Other']['script']

    return configurations

def trilink(**kwargs):
    #initialize template list
    template_dir = kwargs['templates']
    ligand_dir = kwargs['ligands']
    core_identity = kwargs['core']
    invalid_dir = kwargs['invalid']
    valid_dir = kwargs['valid']
    runscript = kwargs['runscript']
    script = kwargs['script']
    

    template_lst = []
    for template in os.listdir(template_dir):
        temp_filepath = f'{template_dir}{template}'
        template = ParsedFile(temp_filepath, template=True, coreidentity=core_identity)
        nonaltcoords, altcoords = RotateMolecule.init_rotate_molecule(template)
        template.set_rotated_coordinates(nonaltcoords)
        template.set_alt_rotated_coordinates(altcoords)
        template_lst.append(template)
    
    #main program loop
    for lig in os.listdir(ligand_dir):
        print(lig)
        lig_filepath = f'{ligand_dir}{lig}'
        ligand = ParsedFile(lig_filepath)
        if ligand.missing_data:
            continue
        
        valid, invalidations = OptimizeValidate.main_optimization(ligand, template_lst)

        #rotamer optimization goes here
        if not valid:
            pass
            #OptimizeValidate.rotamer_optimization(ligand, template_lst)

        #create .com directory for ligand
        outdir = (f'{invalid_dir}{ligand.filename}/')
        if valid:
            outdir = (f'{valid_dir}{ligand.filename}/')
        
        os.mkdir(outdir)
        #create .com files
        for template in template_lst:
            if valid and xtbopt:
                create_xyz_file(template,ligand,outdir,temporaryfile=True)
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
        for dir in os.listdir(valid_dir):
            os.chdir(f"./{valid_dir}{dir}")
            subprocess.run(script.split(), shell=True)


def main():
    params = initial_config()
    d = time.localtime()
    out_dir_label = params['ligands'].replace('/','_')
    xtblabel=""
    if params['xtb']:
        xtblabel="xtb_"

    parser = argparse.ArgumentParser()
    parser.add_argument('-output', '-o', help='Directory where the valid and invalid directories are created', default=None)
    parser.add_argument('-valid', '-v', help='Directory where valid structures are saved', default=None)
    parser.add_argument('-invalid', '-i', help='Directory where invalid structures are saved', default=None)
    args = parser.parse_args()

    if args.output:
        outdir = args.output
    else:
        outdir = f'Output_{out_dir_label}{xtblabel}{d[1]}{d[2]}{d[0]}.{d[3]}{d[4]}{d[5]}'

    if args.valid:
        validligdir = args.valid
    else:
        validligdir = f'{outdir}/valid_ligands/'
    params['valid'] = validligdir

    if args.invalid:
        invalidligdir = args.invalid
    else:
        invalidligdir = f'{outdir}/invalid_ligands/'
    params['invalid'] = invalidligdir

    os.makedirs(validligdir, exist_ok=True)
    os.makedirs(invalidligdir, exist_ok=True)

    trilink(params)
    

if __name__ == "__main__":
    main()
    print('operation completed')