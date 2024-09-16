import subprocess
import os



#takes in a .xyz file, returns an optimized .mol
def xtb_optimize(xyzfile, freezeatoms=None):
    if freezeatoms:
        with open("freeze.inp", 'w') as outf:
            outf.write(f"$fix\n  atoms: {','.join(freezeatoms)}\n$end")
    
    os.environ['OMP_STACKSIZE'] = '1G'
    runscript = ["xtb", f'{xyzfile}', "--opt", "[crude]"]
    if freezeatoms:
        runscript.append("--input")
        runscript.append("freeze.inp")
    subprocess.run(runscript)










if __name__ == "__main__":
    outfdir = ""