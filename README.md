## Overview
TriLink is a tool that assembles transition metal pincer complexes using ligands from the ReaLigands database. It is designed to automate the creation of large numbers of gaussian input files.

TriLink attempts to assemble a pincer complex between each provided ligand and each provided template. Before assembly, TriLink preforms a rudimentary filter algorithm that marks whether each ligand can be successfully mounted on every provided template. Unsuccessful ligands are considered invalid. For each valid structure, a gaussian input file is created. The contents of the gaussian input are determined by template files. 

## Dependencies
- numpy
- xtb (linux only)

## Input
- Tridentate Ligands (folder)
- Templates (folder)
	- contains transition and ground states (.com file format)
	- also contains a filter template (.mol file format)
- config.ini - contains file locations and optimization settings

## Output
- Valid ligands
- Invalid ligands with explanation for invalidity
- Pincer complexes assembled from valid ligands (.com file format)

## Creating valid template files
- Ensure there are 3 dummy atoms representing the ligand connection atoms.
- Dummy atoms must be the last three atoms in the file.
- In order to ensure consistent assembly across all template files, the order in which the 3 dummy atoms are listed is important. The first listed dummy atom should be on the same side of each template.

## Ligand Directory  
Every ligand needs to be from the ReaLigands database OR be in the exact same format as these ligands, with charge and connection atoms documented in the same way. Don't have multiple levels of directory nesting. 
