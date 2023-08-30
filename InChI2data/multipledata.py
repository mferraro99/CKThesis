# This is a script for the automatic generation of a list of data directory

import os
from rdkit import Chem
from rdkit.Chem import AllChem

list_inchi,final_list = [],[]
with open("inchi_list.dat") as file:
    lines = file.readlines()
    for inchi in lines: list_inchi.append(inchi.strip())
    for inchi in range(len(list_inchi)):
        list_inchi[inchi] = list_inchi[inchi].split(" ")
        final_list.append(list_inchi[inchi][0])

list_smiles,final_smiles = [],[]
with open("smiles_list.dat") as file:
    lines = file.readlines()
    for inchi in lines: list_smiles.append(inchi.strip())
    for inchi in range(len(list_smiles)):
        list_smiles[inchi] = list_smiles[inchi].split(" ")
        final_smiles.append(list_smiles[inchi][0])

for i in range(len(final_list)):
    try:
        num_str = str(i+1)
        if i<=8: num_with_zero = num_str.zfill(3 + len(num_str))
        elif i<=98: num_with_zero = num_str.zfill(2 + len(num_str))
        elif i<=1000: num_with_zero = num_str.zfill(1 + len(num_str))
        print(num_with_zero)
        os.mkdir(str(num_with_zero))
        os.chdir(str(num_with_zero))

        text = open("inchi_file.dat","w")
        text.write(final_list[i])
        text.close()

        os.system("makedata")
        os.remove("inchi_file.dat")

        os.chdir("data")
        os.system("rm -rf name.dat")

        smiles = final_smiles[i]
        inchi = final_list[i]
        string = inchi+"\n"+smiles
        text = open("name.dat","w")
        text.write(string)
        text.close()
        os.chdir("../..")
    except:
        os.chdir("..")
        num_str = str(i+1)
        if i<=8: num_with_zero = num_str.zfill(3 + len(num_str))
        elif i<=98: num_with_zero = num_str.zfill(2 + len(num_str))
        elif i<=1000: num_with_zero = num_str.zfill(1 + len(num_str))
        print("ERROR: can't create "+num_with_zero+" data directory")
