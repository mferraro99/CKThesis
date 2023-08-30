# This is a script for the conversion of the ΔH0f(0 K) to ΔH0f(298.15 K)
# The correction is made as suggested by "Thermochemistry in Gaussian" - Joseph W. Ochterski, Ph.D.
# The estimated ΔH0f(0 K) is obtained using the CBH approach at the highest level available

# The correction scheme is defined as:

# ΔH0f(298.15 K) = ΔH0f(0 K) + Hcorr_molecule - sum(xHcorr_atom_i)
# Where the Hcorr_molecule is determine from the log file of level1: Hcorr_molecule = Hcorr_log - ZPE
# x is the number of atoms (i) in the molecule and Hcorr_atom_i is the correction associated with atom_i

import os
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import PeriodicTable

# Define the directory in which the correction data are stored ΔH0f(0 K)
directory_corrections = "/g100/home/userexternal/mferrar2/python_codes/cbh/H_correction_referencespecies.dat"

# Open and save all the corrections information (element-numerical correction)
# 'correction' is a dictionary, using the atomic symbols as keys; it contains the the number of ith atom in the molecule and their numerical correction to ΔH0f(0 K)
correction = {}
with open(directory_corrections) as file:
    lines = file.readlines()
    for line in lines:
        try:
            splitted_line = line.replace("\n","").split()
            correction[splitted_line[0]]=[0,float(splitted_line[1])]
        except:
            pass

# Define the correction associated with the molecule, obtaining the informations from log file and zpe stored in me_files

# Grep the value of Hcorr+zpe from log file
try:
    Hcorr_err = False
    zpe_err = False
    os.system("grep 'Thermal correction to Enthalpy' 'geoms/reac1_l1.log' > Hcorrection.tmp")

    with open("Hcorrection.tmp") as hcorrection:
        line = hcorrection.readlines()
        line = line[0].split()
        Hcorr_plus_zpe = float(line[-1])
        Hcorr_err = True

    # Extract the value of zero point energy
    with open("me_files/reac1_zpe.me") as zpe_file:
        line = zpe_file.readlines()
        line = line[0].strip("\n").replace(" ","")
        zpe = float(line)
        zpe_err = True
except: # If the data are not recovered succesfully, an error is thrown out
    os.chdir("thermo")
    if not Hcorr_err:
        text_error = "Error: check 'reac1_l1.log' for thermal correction to enthalpy; can't compute ΔH0f(298.15 K)"
    elif not zpe_err:
        text_error = "Error: check 'reac1_zpe.me' for Zero Point Energy; can't compute ΔH0f(298.15 K)"
    error = open("DH298K_failed.out","w")
    error.write(text_error)
    error.close()
    sys.exit() # Stop the ΔH0f(298 K) estimation process if no data are found

# Calculate the correction associated with the molecule
Hcorr_molecule = (Hcorr_plus_zpe - zpe) * 627.509 # [kcal/mol]

# Calculate the correction associated with the atoms contained in the molecule

# Extract the InChI contained in 'data/name.dat'
try:
    with open("data/name.dat") as file:
        line = file.readlines()
        inchi = line[0]
except:
    os.chdir("thermo")
    text_error = "Error: check 'name.dat' for InChI identifier; can't compute ΔH0f(298.15 K)"
    error = open("DH298K_failed.out","w")
    error.write(text_error)
    error.close()
    sys.exit() # Stop the ΔH0f(298 K) estimation process if no InChI is found

mol = Chem.MolFromInchi(inchi)
mol = Chem.AddHs(mol)

# Evaluate how many atoms of every type there are in the molecule
for atom_i in range(mol.GetNumAtoms()): correction[mol.GetAtomWithIdx(atom_i).GetSymbol()][0] += 1

# Calculate the estimated correction associated with the atoms
Hcorr_atoms = 0
for key in correction: Hcorr_atoms = Hcorr_atoms + correction[key][0]*correction[key][1] # [kcal/mol]

check_level = False
# Extract the ΔH0f(0 K)
try: # Check if HL estimation is present
    with open("./me_files/reac1_DH0K.me") as hldh0:
        line = hldh0.readlines()
        dh0k =  float(line[0].strip("\n").replace(" ","")) #[kcal/mol]
        check_level = "HL"
except:
    try: # Check if level1 estimation is present
        with open("./output/reac1_DH0K.out") as dh0:
            line = dh0.readlines()
            dh0k=  float(line[0].strip("\n").replace(" ","")) #[kcal/mol]
            check_level = "l1"
    except: # No ΔH0f(0 K) is found, an error is thrown out
        os.chdir("thermo")
        text_error = "Error: no estimated ΔH0f(0 K) present in 'me_files' or 'output'; can't compute ΔH0f(298.15 K)"
        error = open("DH298K_failed.out","w")
        error.write(text_error)
        error.close()
        sys.exit() # Stop the ΔH0f(298 K) estimation process if no ΔH0f(0 K) is found

dh298 = dh0k + Hcorr_molecule - Hcorr_atoms

if check_level == "HL":
    text_ref = "HL estimated ΔH0f(298.15 K): "+str(dh298)+" [kcal/mol]\n\nReference species:\nElement | H0(298K)-H0(0K) [kcal/mol] | Number of atoms\n"
else:
    text_ref = "Level1 estimated ΔH0f(298.15 K): "+str(dh298)+" [kcal/mol]\n\nReference species:\nElement | H0(298K)-H0(0K) [kcal/mol] | Number of atoms\n"

for key in correction:
    if correction[key][0] != 0:
        if len(key) == 1:
            text_ref = text_ref + "   "+key+"                "+str(correction[key][1])+"                   "+str(correction[key][0])+"\n"

        else:
            text_ref = text_ref + "   "+key+"               "+str(correction[key][1])+"                   "+str(correction[key][0])+"\n"

if check_level == "HL":
    os.chdir("me_files")
    file_thermo = open("reac1_DH298K.me","w")
    file_thermo.write("  "+str(dh298))
    file_thermo.close()
    os.chdir("..")
else:
    os.chdir("output")
    file_thermo = open("reac1_DH298K.out","w")
    file_thermo.write("  "+str(dh298))
    file_thermo.close()
    os.chdir("..")

os.chdir("thermo")
file_thermo = open("DH298K.out","w")
file_thermo.write(text_ref)
file_thermo.close()
os.chdir("..")

os.system("rm -rf Hcorrection.tmp")
