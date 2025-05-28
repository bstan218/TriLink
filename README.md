## Overview
TriLink is a tool that assembles transition metal pincer complexes using ligands from the ReaLigands database. It is designed to automate the creation of large numbers of gaussian input files.

TriLink attempts to assemble a pincer complex between each provided ligand and each provided template. Before assembly, TriLink preforms a rudimentary filter algorithm that marks whether each ligand can be successfully mounted on every provided template. Unsuccessful ligands are considered invalid. For each valid structure, a gaussian input file is created. The contents of the gaussian input are determined by template files. 

## Dependencies
- numpy
- xtb (linux only)

## Input
- Tridentate Ligands (directory)
- Templates (directory)
	- contains transition and ground state templates (.com file format)
- config.ini - contains file locations and optimization settings

### Tridentate Ligand Directory  
Every ligand needs to be from the ReaLigands database OR be in the exact same format as these ligands, with charge and connection atoms documented in the same way. Don't have multiple levels of directory nesting. 

### Creating Template Files  
Each template has two purposes:  
1. Provides instructions for gaussian input files output by TriLink.  
2. Gives an initial positions for ligand connection points around the metal catalyst.  

To create template files, follow these steps:
1. For each step of the catalytic cycle, hand mount a ligand and run a gaussian optimization.
2. Replace each gaussian input file's geometry (xyz coordinates) with the optimized geometry.
3. In a molecule viewer like Gaussview or Chemcraft, open the edited gaussian input file. Remove all atoms belonging to the ligand except for those in direct contact with the metal (there should be 3 for tridentate ligands). For each of these three, convert its atom type to hydrogen. These three hydrogens are dummy atoms telling TriLink where to mount ligands.
4. Update the gaussian input file geometry again to match this ligand-stripped geometry. 
5. In the updated geometry, identify which 3 lines correspond to the dummy hydrogens. Move these three lines to the end of the geometry. For example, if the geometry has 30 lines and the dummy hydrogens are on lines 5, 8, and 12, cut these three lines and paste them so they are on lines 28, 29, and 30.
6. Make sure each file has .com file extension. Copy all template files to the template directory.  
- IMPORTANT: - In order to ensure consistent assembly across all template files, the order in which the 3 dummy atoms are listed is important. The first listed dummy atom should be on the same side of each template. 
- Ensure there are 3 dummy atoms representing the ligand connection atoms.
- Dummy atoms must be the last three atoms in the file.

### Config.ini File  
the config.ini file provides input parameters for trilink. Some parameters are work in progress and can be ignored. You will need to edit the following parameters:  
- LigandLibDir - path to tridentate ligand directory
- TemplateDir - path to template directory
- CoreIdentity - atomic identity of metal center, abbreviated to match its abbreviation in the template geometry.
- xtbopt (True or False) - True if wanting to run xtb optimization on valid systems, false if not. HIGHLY RECOMMENDED to run this optimization. **Xtb is only available on linux systems.**

## Output
- Valid ligands
- Invalid ligands with explanation for invalidity
- Pincer complexes assembled from valid ligands (.com file format)


## 
