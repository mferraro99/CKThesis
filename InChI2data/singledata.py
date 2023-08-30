# This is a script for the automatic generation of a directory "data" suitable for es2k, starting from a text file containing the InChI string of
# the desired molecule
# The es2k simulation can be runned directly
# The order of definition of the atoms follows the rules of x2z, written by Yuri Georgievskii, rearranged in order to prioritize the definition of
# heavy atoms (supported C,O,N,Cl,F,S -> Not working so well for silicated compounds)

# !!! RDKIT MUST BE ACTIVATED, EITHER USING AUTOMECH OR SEPARATED ENVIRONMENT !!!
# The environment in Cineca is rdkit-env

# This directory is the source of theory.dat, estoktp.dat, me_head.dat and hl_molpro.dat -> Must be specified
directory_getfile = "/g100/home/userexternal/mferrar2/python_codes/es2k_datafile" # NO / AT THE END!

import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

clean_data = []
with open('inchi_file.dat') as file:
    stringWithMolData = file.readlines() #Translate TextIOWrapper in a string
    for inchi in stringWithMolData: clean_data.append(inchi.strip())# Cycle to clean the input data and work better with the matrix
    inchi = list(filter("".__ne__,clean_data))
    inchi = inchi[0]

m = Chem.inchi.MolFromInchi(inchi)
smiles = Chem.MolToSmiles(m)
m_with_h = Chem.AddHs(m)
Chem.AllChem.EmbedMolecule(m_with_h)

xyz_separated = []

xyz_rdkit = Chem.MolToMolBlock(m_with_h).split("\n")

for i_file in xyz_rdkit: xyz_separated.append(list(filter("".__ne__,i_file.strip().split(" "))))

xyz_separated = list(filter(None,xyz_separated))
num_atoms = int(xyz_separated[1][0])

clear_xyz = ["    "+str(num_atoms)+"\n"]

for i in range(2,2+num_atoms): clear_xyz.append(xyz_separated[i][3]+"  "+xyz_separated[i][0]+"  "+xyz_separated[i][1]+"  "+xyz_separated[i][2])

merged = "\n".join(clear_xyz)
text = open(r"matrix.xyz","w")
text.write(merged)
text.close()

os.system("x2dat matrix.xyz gaussian_matrix")
clean_data = []
with open('gaussian_matrix') as gmat:
    file = gmat.readlines() #Translate TextIOWrapper in a string
    for i_file in file: # Cycle to clean the input data and work better with the matrix
        clean_data.append(i_file.strip())
    clean_data_start = list(filter("".__ne__,clean_data))

counter_first_atom = 0 # We find the position of the first atom
for i in range(len(clean_data_start)):
    if clean_data_start[i].find("atomlabel") != -1:
        pos_first_atom = counter_first_atom + 2
        break
    counter_first_atom = counter_first_atom + 1

for i in range(len(clean_data_start)):
    if clean_data_start[i].find("intcoor") != -1:
        pos_last_atom = i-1
        break

counter_heavy_atoms = 0
for j in range(pos_first_atom,pos_last_atom+1):
    if clean_data_start[j][0] == "C" or clean_data_start[j][0] == "O" or clean_data_start[j][0] == "N" or clean_data_start[j][0] == "S" or clean_data_start[j][0] == "F": counter_heavy_atoms = counter_heavy_atoms + 1

switch_pos = 3
if counter_heavy_atoms > 2 and clean_data_start[pos_first_atom+2][0] == "H":
    wrong_xyz = open(r"matrix.xyz","r+")
    file1 = wrong_xyz.readlines()
    while clean_data_start[pos_first_atom+2][0] == "H": #This means that we need to change our xyz
        clean_xyz = []

        for pos_first_h in range(len(file1)):
            if file1[pos_first_h].find("H") != -1: break
        heavy_atoms = pos_first_h-2

        file1[2],file1[switch_pos] = file1[switch_pos],file1[2]
        switch_pos = switch_pos + 1

        if switch_pos == heavy_atoms + 2:
            switch_pos = 3

        merged = "".join(file1)
        wrong_xyz.seek(0)
        wrong_xyz.write(merged)
        wrong_xyz.truncate()

        os.system("x2dat matrix.xyz gaussian_matrix")
        clean_data = []
        with open('gaussian_matrix') as gmat:
            file = gmat.readlines() #Translate TextIOWrapper in a string
            for i_file in file: clean_data.append(i_file.strip())# Cycle to clean the input data and work better with the matrix
            clean_data_start = list(filter("".__ne__,clean_data))

os.system("rm -rf x2z_out.dat")

clean_data = []
with open('gaussian_matrix') as gmat:
    file = gmat.readlines() #Translate TextIOWrapper in a string
    for i_file in file: clean_data.append(i_file.strip())# Cycle to clean the input data and work better with the matrix
    clean_data_start = list(filter("".__ne__,clean_data))

counter_first_atom = 0 # We find the position of the first atom
for i in range(len(clean_data_start)):
    if clean_data_start[i].find("atomlabel") != -1:
        pos_first_atom = counter_first_atom + 2
        break
    counter_first_atom = counter_first_atom + 1

#Find the atom number
clean_atom_num = []
clean_atom_num = clean_data_start[pos_first_atom-3].split(" ")

atom_num = int(clean_atom_num[1])
atom_num_nodummy = int(clean_atom_num[0])
atom_dummy = int(clean_atom_num[1]) - int(clean_atom_num[0])

#Create the starting atom matrix and first guess matrix
start_matrix = []
start_first_guess = []
for i in range(pos_first_atom,pos_first_atom+atom_num): start_matrix.append(list(filter("".__ne__,clean_data_start[i].split(" "))))

counter_first_dist = 0 # We find the position of the first guessed distance
for i in range(len(clean_data_start)):
    if clean_data_start[i].find("intcoor") != -1:
        pos_first_distance = counter_first_dist + 1
        break
    counter_first_dist = counter_first_dist + 1

#Create a list of heavy atoms and hydrogens, ordered
heavy_atom_list = []
hydrogen_list = []

for i in range(len(start_matrix)):
    if start_matrix[i][0].find("C") != -1 or start_matrix[i][0].find("O") != -1 or start_matrix[i][0].find("N") != -1 or start_matrix[i][0].find("X") != -1 or start_matrix[i][0].find("S") != -1 or start_matrix[i][0].find("F") != -1: heavy_atom_list.append(start_matrix[i])
    else: hydrogen_list.append(start_matrix[i])
ordered_matrix = heavy_atom_list + hydrogen_list

#Create a void twin matrix to change the info contained in the ordered matrix
if len(ordered_matrix) == 2: twin_matrix = [[''],['','','']]
else:
    twin_matrix = [[''],['','',''],['','','','','']]
    if len(ordered_matrix) > 3:
        for i in range(len(ordered_matrix)-3): twin_matrix.append(['','','','','','',''])

#Create a list of first defined atoms
change_num_due_to_dummy = 0
for i in range(len(twin_matrix)):
    if ordered_matrix[i][0][0] == "X": change_num_due_to_dummy = change_num_due_to_dummy + 1
    if ordered_matrix[i][0][0] != "X" and  ordered_matrix[i][0][1].isdigit(): twin_matrix[i][0] = ordered_matrix[i][0][0]+str(i+1-change_num_due_to_dummy)
    elif ordered_matrix[i][0][0] != "X" and not ordered_matrix[i][0][1].isdigit(): twin_matrix[i][0] = ordered_matrix[i][0][0] + ordered_matrix[i][0][1] +str(i+1-change_num_due_to_dummy)
    elif ordered_matrix[i][0][0] == "X" and atom_dummy < 2: twin_matrix[i][0] = ordered_matrix[i][0][0]
    else: twin_matrix[i][0] = ordered_matrix[i][0]

for i in range(len(ordered_matrix)):
    if ordered_matrix[i][0] != twin_matrix[i][0] and ordered_matrix[i][0][0] != "X": #Only if we had a change in nomenclature
        for j in range(i,len(ordered_matrix)):
            for k in range(1,len(ordered_matrix[j])-1,2):
                if twin_matrix[j][k] == '' and ordered_matrix[j][k] == ordered_matrix[i][0]:
                    twin_matrix[j][k] = twin_matrix[i][0]

if atom_dummy < 2:
    for i in range(1,len(ordered_matrix)):
        for j in range(1,len(ordered_matrix[i])-1,2):
            if twin_matrix[i][j] == '' and ordered_matrix[i][j][0] != "X": twin_matrix[i][j] = ordered_matrix[i][j]
            elif twin_matrix[i][j] == '' and ordered_matrix[i][j][0] == "X": twin_matrix[i][j] = "X"

    #Copy the unchanged nomenclature in the twin matrix and assign the name of distances, angles and dh -> Is digit takes into account atoms with more than one letter as symbol (suck as Cl, Si)
    change_num_due_to_dummy = 0
    for i in range(1,len(ordered_matrix)):
        if twin_matrix[i][0] == "X": change_num_due_to_dummy = change_num_due_to_dummy + 1
        if len(twin_matrix[i]) == 3:
            if twin_matrix[i][0] != "X":
                if twin_matrix[i][0][1].isdigit() and twin_matrix[i][1][1].isdigit(): twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][1][0] + str(i+1-change_num_due_to_dummy)
                elif not twin_matrix[i][0][1].isdigit() and twin_matrix[i][1][1].isdigit(): twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][0][1] + twin_matrix[i][1][0] + str(i+1-change_num_due_to_dummy)
                elif twin_matrix[i][0][1].isdigit() and not twin_matrix[i][1][1].isdigit(): twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][1][0] + twin_matrix[i][1][1] + str(i+1-change_num_due_to_dummy)
                elif not twin_matrix[i][0][1].isdigit() and not twin_matrix[i][1][1].isdigit(): twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][0][1] + twin_matrix[i][1][0] + twin_matrix[i][1][1] + str(i+1-change_num_due_to_dummy)
            else: twin_matrix[i][2] = ordered_matrix[i][2]
        elif len(twin_matrix[i]) == 5 or len(twin_matrix[i]) == 7:
            if twin_matrix[i][0] != "X" and len(twin_matrix[i][1]) > 1 and len(twin_matrix[i][3]) > 1:
                if twin_matrix[i][0][1].isdigit() and twin_matrix[i][1][1].isdigit() and twin_matrix[i][3][1].isdigit():
                    twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][1][0] + str(i+1-change_num_due_to_dummy)
                    twin_matrix[i][4] = twin_matrix[i][0][0] + twin_matrix[i][1][0] + twin_matrix[i][3][0] + str(i+1-change_num_due_to_dummy)
                elif not twin_matrix[i][0][1].isdigit() and twin_matrix[i][1][1].isdigit() and twin_matrix[i][3][1].isdigit():
                    twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][0][1] + twin_matrix[i][1][0] + str(i+1-change_num_due_to_dummy)
                    twin_matrix[i][4] = twin_matrix[i][0][0] + twin_matrix[i][0][1] + twin_matrix[i][1][0] + twin_matrix[i][3][0] + str(i+1-change_num_due_to_dummy)
                elif twin_matrix[i][0][1].isdigit() and not twin_matrix[i][1][1].isdigit() and twin_matrix[i][3][1].isdigit():
                    twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][1][0] + twin_matrix[i][1][1] + str(i+1-change_num_due_to_dummy)
                    twin_matrix[i][4] = twin_matrix[i][0][0] + twin_matrix[i][1][0] + twin_matrix[i][1][1] + twin_matrix[i][3][0] + str(i+1-change_num_due_to_dummy)
                elif twin_matrix[i][0][1].isdigit() and twin_matrix[i][1][1].isdigit() and not twin_matrix[i][3][1].isdigit():
                    twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][1][0] + str(i+1-change_num_due_to_dummy)
                    twin_matrix[i][4] = twin_matrix[i][0][0] + twin_matrix[i][1][0] + twin_matrix[i][3][0] + twin_matrix[i][3][1] + str(i+1-change_num_due_to_dummy)
                elif not twin_matrix[i][0][1].isdigit() and not twin_matrix[i][1][1].isdigit() and twin_matrix[i][3][1].isdigit():
                    twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][0][1] + twin_matrix[i][1][0] + twin_matrix[i][1][1] + str(i+1-change_num_due_to_dummy)
                    twin_matrix[i][4] = twin_matrix[i][0][0] + twin_matrix[i][0][1] + twin_matrix[i][1][0] + twin_matrix[i][1][1] + twin_matrix[i][3][0] + str(i+1-change_num_due_to_dummy)
                elif not twin_matrix[i][0][1].isdigit() and  twin_matrix[i][1][1].isdigit() and not twin_matrix[i][3][1].isdigit():
                    twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][0][1] + twin_matrix[i][1][0] + str(i+1-change_num_due_to_dummy)
                    twin_matrix[i][4] = twin_matrix[i][0][0] + twin_matrix[i][0][1] + twin_matrix[i][1][0] + twin_matrix[i][3][0] + twin_matrix[i][3][1] + str(i+1-change_num_due_to_dummy)
                elif twin_matrix[i][0][1].isdigit() and not twin_matrix[i][1][1].isdigit() and not twin_matrix[i][3][1].isdigit():
                    twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][1][0] + twin_matrix[i][1][1] + str(i+1-change_num_due_to_dummy)
                    twin_matrix[i][4] = twin_matrix[i][0][0] + twin_matrix[i][1][0] + twin_matrix[i][1][1] + twin_matrix[i][3][0] + twin_matrix[i][3][1] + str(i+1-change_num_due_to_dummy)
                elif not twin_matrix[i][0][1].isdigit() and not twin_matrix[i][1][1].isdigit() and not twin_matrix[i][3][1].isdigit():
                    twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][0][1] + twin_matrix[i][1][0] + twin_matrix[i][1][1] + str(i+1-change_num_due_to_dummy)
                    twin_matrix[i][4] = twin_matrix[i][0][0] + twin_matrix[i][0][1] + twin_matrix[i][1][0] + twin_matrix[i][1][1] + twin_matrix[i][3][0] + twin_matrix[i][3][1] + str(i+1-change_num_due_to_dummy)
            elif twin_matrix[i][0] != "X" and len(twin_matrix[i][1]) > 1 and len(twin_matrix[i][3]) == 1:
                if twin_matrix[i][0][1].isdigit() and twin_matrix[i][1][1].isdigit():
                    twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][1][0] + str(i+1-change_num_due_to_dummy)
                    twin_matrix[i][4] = twin_matrix[i][0][0] + twin_matrix[i][1][0] + twin_matrix[i][3][0] + str(i+1-change_num_due_to_dummy)
                elif not twin_matrix[i][0][1].isdigit() and twin_matrix[i][1][1].isdigit():
                    twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][0][1] + twin_matrix[i][1][0] + str(i+1-change_num_due_to_dummy)
                    twin_matrix[i][4] = twin_matrix[i][0][0] + twin_matrix[i][0][1] + twin_matrix[i][1][0] + twin_matrix[i][3][0] + str(i+1-change_num_due_to_dummy)
                elif twin_matrix[i][0][1].isdigit() and not twin_matrix[i][1][1].isdigit():
                    twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][1][0] + twin_matrix[i][1][1] + str(i+1-change_num_due_to_dummy)
                    twin_matrix[i][4] = twin_matrix[i][0][0] + twin_matrix[i][1][0] + twin_matrix[i][1][1] + twin_matrix[i][3][0] + str(i+1-change_num_due_to_dummy)
                elif not twin_matrix[i][0][1].isdigit() and not twin_matrix[i][1][1].isdigit():
                    twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][0][1] + twin_matrix[i][1][0] + twin_matrix[i][1][1] + str(i+1-change_num_due_to_dummy)
                    twin_matrix[i][4] = twin_matrix[i][0][0] + twin_matrix[i][0][1] + twin_matrix[i][1][0] + twin_matrix[i][1][1] + twin_matrix[i][3][0] + str(i+1-change_num_due_to_dummy)
            elif twin_matrix[i][0] != "X" and len(twin_matrix[i][1]) == 1 and len(twin_matrix[i][3]) > 1:
                if twin_matrix[i][0][1].isdigit() and twin_matrix[i][1][3].isdigit():
                    twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][1][0] + str(i+1-change_num_due_to_dummy)
                    twin_matrix[i][4] = twin_matrix[i][0][0] + twin_matrix[i][1][0] + twin_matrix[i][3][0] + str(i+1-change_num_due_to_dummy)
                elif not twin_matrix[i][0][1].isdigit() and twin_matrix[i][1][3].isdigit():
                    twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][0][1] + twin_matrix[i][1][0] + str(i+1-change_num_due_to_dummy)
                    twin_matrix[i][4] = twin_matrix[i][0][0] + twin_matrix[i][0][1] + twin_matrix[i][1][0] + twin_matrix[i][3][0] + str(i+1-change_num_due_to_dummy)
                elif twin_matrix[i][0][1].isdigit() and not twin_matrix[i][1][3].isdigit():
                    twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][1][0] + twin_matrix[i][1][1] + str(i+1-change_num_due_to_dummy)
                    twin_matrix[i][4] = twin_matrix[i][0][0] + twin_matrix[i][1][0] + twin_matrix[i][1][1] + twin_matrix[i][3][0] + twin_matrix[i][3][1] + str(i+1-change_num_due_to_dummy)
                elif not twin_matrix[i][0][1].isdigit() and not twin_matrix[i][1][3].isdigit():
                    twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][0][1] + twin_matrix[i][1][0] + twin_matrix[i][1][1] + str(i+1-change_num_due_to_dummy)
                    twin_matrix[i][4] = twin_matrix[i][0][0] + twin_matrix[i][0][1] + twin_matrix[i][1][0] + twin_matrix[i][1][1] + twin_matrix[i][3][0] + twin_matrix[i][3][1] + str(i+1-change_num_due_to_dummy)
            else:
                twin_matrix[i][2] = ordered_matrix[i][2]
                twin_matrix[i][4] = ordered_matrix[i][4]
        if len(twin_matrix[i]) == 7:
            if twin_matrix[i][0] != "X":
                twin_matrix[i][6] = "dih"+str(i+1-change_num_due_to_dummy)
            else:
                twin_matrix[i][6] = ordered_matrix[i][6]

    minus_dummy = 0
    #Create a list of first guess distances
    for i in range(1,len(twin_matrix)):
        if twin_matrix[i][0][0] == "X":
            if i == 1: minus_dummy = minus_dummy + 1
            elif i == 2: minus_dummy = minus_dummy + 2
            else: minus_dummy = minus_dummy + 3

    if atom_num>2:
        for i in range(pos_first_distance,pos_first_distance+3*atom_num-6-minus_dummy): start_first_guess.append(list(filter("".__ne__,clean_data_start[i].split(" "))))
    else:
        for i in range(pos_first_distance,pos_first_distance+3*atom_num-5-minus_dummy): start_first_guess.append(list(filter("".__ne__,clean_data_start[i].split(" "))))
    #Create a void twin list of first guesses to change the nomenclature

    twin_list = []
    for i in range(3*atom_num-6): twin_list.append(['',''])
    if atom_num == 2: twin_list = ['','']
    ### Renaiming of the first guess list

    if atom_num == 2: twin_list[0] = twin_matrix[1][2]
    else:
        for i in range(1,atom_num): twin_list[i-1][0] = twin_matrix[i][2]
        for i in range(atom_num-1,2*atom_num-3): twin_list[i][0] = twin_matrix[i-atom_num+3][4]
        for i in range(2*atom_num-3,len(twin_list)): twin_list[i][0] = "dih"+str(i-2*atom_num+7-minus_dummy)
        twin_list_final = []
        for j in range(len(twin_list)-1,-1,-1):
            if twin_list[j][0].find("C") == -1 and twin_list[j][0].find("S") == -1 and twin_list[j][0].find("F") == -1 and twin_list[j][0].find("O") == -1 and twin_list[j][0].find("N") == -1 and twin_list[j][0].find("H") == -1 and twin_list[j][0].find("dih") == -1: twin_list[j].pop()
        twin_list_final = []
        for k in range(len(twin_list)):
            if len(twin_list[k]) == 2: twin_list_final.append(twin_list[k])
        twin_list = twin_list_final
    if atom_num == 2: twin_list[1] = start_first_guess[0][1]
    else:
        for i in range(1,len(twin_matrix)):
            for j in range(len(twin_list)):
                if i >= 1 and twin_matrix[i][2] == twin_list[j][0]:
                    for k in range(len(start_first_guess)):
                        if start_first_guess[k][0] == ordered_matrix[i][2]:
                            twin_list[j][1] = start_first_guess[k][1]
                            break
                if i >= 2 and twin_matrix[i][4] == twin_list[j][0]:
                    for k in range(len(start_first_guess)):
                        if start_first_guess[k][0] == ordered_matrix[i][4]:
                            twin_list[j][1] = start_first_guess[k][1]
                            break
                if i >= 3 and twin_matrix[i][6] == twin_list[j][0]:
                    for k in range(len(start_first_guess)):
                        if start_first_guess[k][0] == ordered_matrix[i][6]:
                            twin_list[j][1] = start_first_guess[k][1]
                            break
    ###

    # Change matrix and list based on the position of X (if present)
    if twin_matrix[1][0] == "X":
        list_without_dummy = []
        for i in range(2,len(twin_matrix)):
            if twin_matrix[i][3] == "X":
                for j in range(len(twin_list)):
                    if twin_list[j][0] == twin_matrix[i][4]:
                        twin_matrix[i][4] = "90."
                        twin_list[j].pop()
                break
        for i in range(3,len(twin_matrix)):
            if twin_matrix[i][5] == "X":
                for j in range(len(twin_list)):
                    if twin_list[j][0] == twin_matrix[i][6]:
                        twin_matrix[i][6] = "180."
                        twin_list[j].pop()
                break
        for i in range(len(twin_list)):
            if len(twin_list[i]) != 1: list_without_dummy.append(twin_list[i])
        twin_list = list_without_dummy
    if atom_num > 2 and twin_matrix[2][0] == "X":
        list_without_dummy = []
        for i in range(3,len(twin_matrix)):
            if twin_matrix[i][5] == "X":
                for j in range(len(twin_list)):
                    if twin_list[j][0] == twin_matrix[i][6]:
                        twin_matrix[i][6] = "180."
                        twin_list[j].pop()
                break
        for i in range(len(twin_list)):
            if len(twin_list[i]) != 1: list_without_dummy.append(twin_list[i])
        twin_list = list_without_dummy
else:
    for i in range(1,len(ordered_matrix)):
        for j in range(1,len(ordered_matrix[i])-1,2):
            if twin_matrix[i][j] == '': twin_matrix[i][j] = ordered_matrix[i][j]

    change_num_due_to_dummy = 0
    for i in range(1,len(ordered_matrix)):
        if twin_matrix[i][0][0] == "X": change_num_due_to_dummy = change_num_due_to_dummy + 1
        if len(twin_matrix[i]) == 3:
            if twin_matrix[i][0][0] != "X":
                if twin_matrix[i][0][1].isdigit() and twin_matrix[i][1][1].isdigit(): twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][1][0] + str(i+1-change_num_due_to_dummy)
                elif not twin_matrix[i][0][1].isdigit() and twin_matrix[i][1][1].isdigit(): twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][0][1] + twin_matrix[i][1][0] + str(i+1-change_num_due_to_dummy)
                elif twin_matrix[i][0][1].isdigit() and not twin_matrix[i][1][1].isdigit(): twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][1][0] + twin_matrix[i][1][1] + str(i+1-change_num_due_to_dummy)
                elif not twin_matrix[i][0][1].isdigit() and not twin_matrix[i][1][1].isdigit(): twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][0][1] + twin_matrix[i][1][0] + twin_matrix[i][1][1] + str(i+1-change_num_due_to_dummy)
            else: twin_matrix[i][2] = ordered_matrix[i][2]
        elif len(twin_matrix[i]) == 5 or len(twin_matrix[i]) == 7:
            if twin_matrix[i][0][0] != "X" and len(twin_matrix[i][1]) > 1 and len(twin_matrix[i][3]) > 1:
                if twin_matrix[i][0][1].isdigit() and twin_matrix[i][1][1].isdigit() and twin_matrix[i][3][1].isdigit():
                    twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][1][0] + str(i+1-change_num_due_to_dummy)
                    twin_matrix[i][4] = twin_matrix[i][0][0] + twin_matrix[i][1][0] + twin_matrix[i][3][0] + str(i+1-change_num_due_to_dummy)
                elif not twin_matrix[i][0][1].isdigit() and twin_matrix[i][1][1].isdigit() and twin_matrix[i][3][1].isdigit():
                    twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][0][1] + twin_matrix[i][1][0] + str(i+1-change_num_due_to_dummy)
                    twin_matrix[i][4] = twin_matrix[i][0][0] + twin_matrix[i][0][1] + twin_matrix[i][1][0] + twin_matrix[i][3][0] + str(i+1-change_num_due_to_dummy)
                elif twin_matrix[i][0][1].isdigit() and not twin_matrix[i][1][1].isdigit() and twin_matrix[i][3][1].isdigit():
                    twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][1][0] + twin_matrix[i][1][1] + str(i+1-change_num_due_to_dummy)
                    twin_matrix[i][4] = twin_matrix[i][0][0] + twin_matrix[i][1][0] + twin_matrix[i][1][1] + twin_matrix[i][3][0] + str(i+1-change_num_due_to_dummy)
                elif twin_matrix[i][0][1].isdigit() and twin_matrix[i][1][1].isdigit() and not twin_matrix[i][3][1].isdigit():
                    twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][1][0] + str(i+1-change_num_due_to_dummy)
                    twin_matrix[i][4] = twin_matrix[i][0][0] + twin_matrix[i][1][0] + twin_matrix[i][3][0] + twin_matrix[i][3][1] + str(i+1-change_num_due_to_dummy)
                elif not twin_matrix[i][0][1].isdigit() and not twin_matrix[i][1][1].isdigit() and twin_matrix[i][3][1].isdigit():
                    twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][0][1] + twin_matrix[i][1][0] + twin_matrix[i][1][1] + str(i+1-change_num_due_to_dummy)
                    twin_matrix[i][4] = twin_matrix[i][0][0] + twin_matrix[i][0][1] + twin_matrix[i][1][0] + twin_matrix[i][1][1] + twin_matrix[i][3][0] + str(i+1-change_num_due_to_dummy)
                elif not twin_matrix[i][0][1].isdigit() and  twin_matrix[i][1][1].isdigit() and not twin_matrix[i][3][1].isdigit():
                    twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][0][1] + twin_matrix[i][1][0] + str(i+1-change_num_due_to_dummy)
                    twin_matrix[i][4] = twin_matrix[i][0][0] + twin_matrix[i][0][1] + twin_matrix[i][1][0] + twin_matrix[i][3][0] + twin_matrix[i][3][1] + str(i+1-change_num_due_to_dummy)
                elif twin_matrix[i][0][1].isdigit() and not twin_matrix[i][1][1].isdigit() and not twin_matrix[i][3][1].isdigit():
                    twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][1][0] + twin_matrix[i][1][1] + str(i+1-change_num_due_to_dummy)
                    twin_matrix[i][4] = twin_matrix[i][0][0] + twin_matrix[i][1][0] + twin_matrix[i][1][1] + twin_matrix[i][3][0] + twin_matrix[i][3][1] + str(i+1-change_num_due_to_dummy)
                elif not twin_matrix[i][0][1].isdigit() and not twin_matrix[i][1][1].isdigit() and not twin_matrix[i][3][1].isdigit():
                    twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][0][1] + twin_matrix[i][1][0] + twin_matrix[i][1][1] + str(i+1-change_num_due_to_dummy)
                    twin_matrix[i][4] = twin_matrix[i][0][0] + twin_matrix[i][0][1] + twin_matrix[i][1][0] + twin_matrix[i][1][1] + twin_matrix[i][3][0] + twin_matrix[i][3][1] + str(i+1-change_num_due_to_dummy)
            elif twin_matrix[i][0] != "X" and len(twin_matrix[i][1]) > 1 and len(twin_matrix[i][3]) == 1:
                if twin_matrix[i][0][1].isdigit() and twin_matrix[i][1][1].isdigit():
                    twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][1][0] + str(i+1-change_num_due_to_dummy)
                    twin_matrix[i][4] = twin_matrix[i][0][0] + twin_matrix[i][1][0] + twin_matrix[i][3][0] + str(i+1-change_num_due_to_dummy)
                elif not twin_matrix[i][0][1].isdigit() and twin_matrix[i][1][1].isdigit():
                    twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][0][1] + twin_matrix[i][1][0] + str(i+1-change_num_due_to_dummy)
                    twin_matrix[i][4] = twin_matrix[i][0][0] + twin_matrix[i][0][1] + twin_matrix[i][1][0] + twin_matrix[i][3][0] + str(i+1-change_num_due_to_dummy)
                elif twin_matrix[i][0][1].isdigit() and not twin_matrix[i][1][1].isdigit():
                    twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][1][0] + twin_matrix[i][1][1] + str(i+1-change_num_due_to_dummy)
                    twin_matrix[i][4] = twin_matrix[i][0][0] + twin_matrix[i][1][0] + twin_matrix[i][1][1] + twin_matrix[i][3][0] + str(i+1-change_num_due_to_dummy)
                elif not twin_matrix[i][0][1].isdigit() and not twin_matrix[i][1][1].isdigit():
                    twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][0][1] + twin_matrix[i][1][0] + twin_matrix[i][1][1] + str(i+1-change_num_due_to_dummy)
                    twin_matrix[i][4] = twin_matrix[i][0][0] + twin_matrix[i][0][1] + twin_matrix[i][1][0] + twin_matrix[i][1][1] + twin_matrix[i][3][0] + str(i+1-change_num_due_to_dummy)
            elif twin_matrix[i][0] != "X" and len(twin_matrix[i][1]) == 1 and len(twin_matrix[i][3]) > 1:
                if twin_matrix[i][0][1].isdigit() and twin_matrix[i][1][3].isdigit():
                    twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][1][0] + str(i+1-change_num_due_to_dummy)
                    twin_matrix[i][4] = twin_matrix[i][0][0] + twin_matrix[i][1][0] + twin_matrix[i][3][0] + str(i+1-change_num_due_to_dummy)
                elif not twin_matrix[i][0][1].isdigit() and twin_matrix[i][1][3].isdigit():
                    twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][0][1] + twin_matrix[i][1][0] + str(i+1-change_num_due_to_dummy)
                    twin_matrix[i][4] = twin_matrix[i][0][0] + twin_matrix[i][0][1] + twin_matrix[i][1][0] + twin_matrix[i][3][0] + str(i+1-change_num_due_to_dummy)
                elif twin_matrix[i][0][1].isdigit() and not twin_matrix[i][1][3].isdigit():
                    twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][1][0] + twin_matrix[i][1][1] + str(i+1-change_num_due_to_dummy)
                    twin_matrix[i][4] = twin_matrix[i][0][0] + twin_matrix[i][1][0] + twin_matrix[i][1][1] + twin_matrix[i][3][0] + twin_matrix[i][3][1] + str(i+1-change_num_due_to_dummy)
                elif not twin_matrix[i][0][1].isdigit() and not twin_matrix[i][1][3].isdigit():
                    twin_matrix[i][2] = twin_matrix[i][0][0] + twin_matrix[i][0][1] + twin_matrix[i][1][0] + twin_matrix[i][1][1] + str(i+1-change_num_due_to_dummy)
                    twin_matrix[i][4] = twin_matrix[i][0][0] + twin_matrix[i][0][1] + twin_matrix[i][1][0] + twin_matrix[i][1][1] + twin_matrix[i][3][0] + twin_matrix[i][3][1] + str(i+1-change_num_due_to_dummy)
            else:
                twin_matrix[i][2] = ordered_matrix[i][2]
                twin_matrix[i][4] = ordered_matrix[i][4]
        if len(twin_matrix[i]) == 7:
            if twin_matrix[i][0][0] != "X": twin_matrix[i][6] = "dih"+str(i+1-change_num_due_to_dummy)
            else: twin_matrix[i][6] = ordered_matrix[i][6]

    minus_dummy = 0
    #Create a list of first guess distances
    for i in range(1,len(twin_matrix)):
        if twin_matrix[i][0][0] == "X":
            if i == 1: minus_dummy = minus_dummy + 1
            elif i == 2: minus_dummy = minus_dummy + 2
            else: minus_dummy = minus_dummy + 3

    if atom_num>2:
        for i in range(pos_first_distance,pos_first_distance+3*atom_num-6-minus_dummy): start_first_guess.append(list(filter("".__ne__,clean_data_start[i].split(" "))))
    else:
        for i in range(pos_first_distance,pos_first_distance+3*atom_num-5-minus_dummy): start_first_guess.append(list(filter("".__ne__,clean_data_start[i].split(" "))))

    twin_list = []
    twin_list_R = []
    twin_list_A = []
    twin_list_D = []
    dh_change_dummy = 0
    for i in range(atom_num-1): twin_list_R.append(['',''])
    for i in range(atom_num-2): twin_list_A.append(['',''])
    for i in range(atom_num-3): twin_list_D.append(['',''])

    if atom_num == 2: twin_list[0] = twin_matrix[1][2]
    else:
        for i in range(1,len(twin_matrix)):
            if twin_matrix[i][0][0] == "X": dh_change_dummy = dh_change_dummy + 1
            twin_list_R[i-1][0] = twin_matrix[i][2]
            if i > 1: twin_list_A[i-2][0] = twin_matrix[i][4]
            if atom_num > 3 and i > 2:
                if not ordered_matrix[i][6][0].isdigit(): twin_list_D[i-3][0] = "dih"+str(i+1-dh_change_dummy)

        for j in range(len(twin_list_R)-1,-1,-1):
            if twin_list_R[j][0].find("C") == -1 and twin_list_R[j][0].find("S") == -1 and twin_list_R[j][0].find("F") == -1 and twin_list_R[j][0].find("O") == -1 and twin_list_R[j][0].find("N") == -1 and twin_list_R[j][0].find("H") == -1: twin_list_R[j].pop()
        for j in range(len(twin_list_A)-1,-1,-1):
            if twin_list_A[j][0].find("C") == -1 and twin_list_A[j][0].find("S") == -1 and twin_list_A[j][0].find("F") == -1 and twin_list_A[j][0].find("O") == -1 and twin_list_A[j][0].find("N") == -1 and twin_list_A[j][0].find("H") == -1: twin_list_A[j].pop()
        for j in range(len(twin_list_D)-1,-1,-1):
             if twin_list_D[j][0].find("dih") == -1: twin_list_D[j].pop()

        twin_list= twin_list_R + twin_list_A + twin_list_D
        twin_list_final = []
        for k in range(len(twin_list)):
            if len(twin_list[k]) == 2: twin_list_final.append(twin_list[k])
        twin_list = twin_list_final

    if atom_num == 2: twin_list[1] = start_first_guess[0][1]
    else:
        for i in range(1,len(twin_matrix)):
            for j in range(len(twin_list)):
                if i >= 1 and twin_matrix[i][2] == twin_list[j][0]:
                    for k in range(len(start_first_guess)):
                        if start_first_guess[k][0] == ordered_matrix[i][2]:
                            twin_list[j][1] = start_first_guess[k][1]
                            break
                if i >= 2 and twin_matrix[i][4] == twin_list[j][0]:
                    for k in range(len(start_first_guess)):
                        if start_first_guess[k][0] == ordered_matrix[i][4]:
                            twin_list[j][1] = start_first_guess[k][1]
                            break
                if i >= 3 and twin_matrix[i][6] == twin_list[j][0]:
                    for k in range(len(start_first_guess)):
                        if start_first_guess[k][0] == ordered_matrix[i][6]:
                            twin_list[j][1] = start_first_guess[k][1]
                            break

    # Change matrix and list based on the position of X (if present)
    if twin_matrix[1][0][0] == "X":
        list_without_dummy = []
        change_dh = True
        for i in range(2,len(twin_matrix)):
            if twin_matrix[i][3][0] == "X":
                for j in range(len(twin_list)):
                    if twin_list[j][0] == twin_matrix[i][4]:
                        twin_matrix[i][4] = "90."
                        twin_list[j].pop()
                break
        for i in range(3,len(twin_matrix)):
            if twin_matrix[i][5][0] == "X":
                for j in range(len(twin_list)):
                    if twin_list[j][0] == twin_matrix[i][6]:
                        twin_matrix[i][6] = "180."
                        twin_list[j].pop()
                        change_dh = False
                break

        if change_dh:
            for i in range(3,len(twin_matrix)):
                if twin_matrix[i][3][0] == "X" and not twin_matrix[i][4][0].isdigit():
                    for j in range(len(twin_list)):
                        if twin_list[j][0] == twin_matrix[i][6]:
                            twin_matrix[i][6] = "180."
                            twin_list[j].pop()
                    break
        for i in range(len(twin_list)):
            if len(twin_list[i]) != 1: list_without_dummy.append(twin_list[i])
        twin_list = list_without_dummy
    if atom_num > 2 and twin_matrix[2][0][0] == "X":
        list_without_dummy = []
        for i in range(3,len(twin_matrix)):
            if twin_matrix[i][5][0] == "X":
                for j in range(len(twin_list)):
                    if twin_list[j][0] == twin_matrix[i][6]:
                        twin_matrix[i][6] = "180."
                        twin_list[j].pop()
                break
        for i in range(len(twin_list)):
            if len(twin_list[i]) != 1: list_without_dummy.append(twin_list[i])
        twin_list = list_without_dummy
    for i in range(len(twin_matrix)):
        if twin_matrix[i][0][0] == "X":
            for j in range(i,len(twin_matrix)):
                for k in range(1,len(twin_matrix[j])-1,2):
                    if twin_matrix[j][k] == twin_matrix[i][0] and i == 1: twin_matrix[j][k] = "X"+str(i)
                    elif twin_matrix[j][k] == twin_matrix[i][0] and i > 1: twin_matrix[j][k] = "X"+str(i-1)
            if i == 1:  twin_matrix[i][0] = "X"+str(i)
            else: twin_matrix[i][0] = "X"+str(i-1)

for i in range(len(clean_data_start)):
    if clean_data_start[i].find("ntau") != -1: num_tau = int(clean_data_start[i+1])
    if clean_data_start[i].find("nhind") != -1 and clean_data_start[i].find("namehind,hindmn,hindmx,nhindsteps") == -1: num_hr = int(clean_data_start[i+1])
    if clean_data_start[i].find("taumn,taumx sampling interval") != -1: start_tau = int(i+1)
    if clean_data_start[i].find("namehind,hindmn,hindmx,nhindsteps") != -1:
        start_hr = int(i+1)
        break

tau=[]
for i in range(start_tau,start_tau+num_tau): tau.append(list(filter("".__ne__,clean_data_start[i].split(" "))))
hr = []
for i in range(start_hr,start_hr+num_hr): hr.append(list(filter("".__ne__,clean_data_start[i].split(" "))))

for i in range(3,len(ordered_matrix)):
    for j in range(len(tau)):
        if tau[j][0] == ordered_matrix[i][6]: tau[j][0] = twin_matrix[i][6]
    for k in range(len(hr)):
        if hr[k][0] == ordered_matrix[i][6]: hr[k][0] = twin_matrix[i][6]

for i in range(len(hr)):
    for j in range(2*atom_num_nodummy-3,len(twin_list)):
        if hr[i][0] == twin_list[j][0]: twin_list[j].pop()

dh_list = twin_list[2*atom_num-3:]
twin_list = twin_list[:2*atom_num-3]

for i in range(len(dh_list)):
    if len(dh_list[i]) != 1: twin_list.append(dh_list[i])

#From upper to lower case
if atom_num == 2:
    twin_matrix[0][0] = twin_matrix[0][0].lower()
    twin_matrix[1][0] = twin_matrix[1][0].lower()
    twin_matrix[1][1] = twin_matrix[1][1].lower()
    twin_matrix[1][2] = twin_matrix[1][2].lower()
    twin_list[0] = twin_list[0].lower()
else:
    for i in range(len(twin_matrix)):
        for j in range(len(twin_matrix[i])): twin_matrix[i][j] = twin_matrix[i][j].lower()
    for i in range(len(twin_list)): twin_list[i][0] = twin_list[i][0].lower()

#Bubble-sorting the list of tau and hr
for i in range(len(tau)-2):
    for j in range(i+1,len(tau)-1):
        if int(tau[j][0][3:]) < int(tau[i][0][3:]): tau[j],tau[i]= tau[i],tau[j]
for i in range(len(hr)-2):
    for j in range(i+1,len(hr)-1):
        if int(hr[j][0][3:]) < int(hr[i][0][3:]): hr[j],hr[i]= hr[i],hr[j]

tau_merged = ""
for i in range(len(tau)): tau_merged = tau_merged+"  ".join(tau[i])+"\n"
hr_merged = ""

for i in range(len(hr)): hr_merged = hr_merged+"  ".join(hr[i])+"\n"

matrix_merged = ""
for i in range(len(twin_matrix)):
    merged_atom = ""
    for j in range(len(twin_matrix[i])):
        if len(twin_matrix[i][j]) == 1: merged_atom = merged_atom + twin_matrix[i][j] + "      "
        elif len(twin_matrix[i][j]) == 2: merged_atom = merged_atom + twin_matrix[i][j] + "     "
        elif len(twin_matrix[i][j]) == 3: merged_atom = merged_atom + twin_matrix[i][j] + "    "
        elif len(twin_matrix[i][j]) == 4: merged_atom = merged_atom + twin_matrix[i][j] + "   "
        elif len(twin_matrix[i][j]) == 5: merged_atom = merged_atom + twin_matrix[i][j] + "  "
        else: merged_atom = merged_atom + twin_matrix[i][j] + " "
    matrix_merged = matrix_merged+merged_atom+"\n"

list_merged = ""
if atom_num == 2: list_merged = twin_list[0] + "   " + twin_list[1]+"\n"
else:
    for i in range(len(twin_list)):
        merged_data = ""
        for j in range(len(twin_list[i])):
            if len(twin_list[i][j]) == 3: merged_data = merged_data + twin_list[i][j] + "    "
            elif len(twin_list[i][j]) == 4: merged_data = merged_data + twin_list[i][j] + "   "
            elif len(twin_list[i][j]) == 5: merged_data = merged_data + twin_list[i][j] + "  "
            else: merged_data = merged_data + twin_list[i][j] + " "
        list_merged = list_merged+merged_data+"\n"

sample = clean_data_start[1].split(" ")
clean_data_start[1] = " ".join(sample)

merged ="\n".join(clean_data_start[:2])+"\n\nntau\n"+str(len(tau))+"\n--> name and sampling interval\n" + tau_merged + "\nnhind\n"+str(len(hr))+"\n-->namehind,hindmn,hindmx,nhindsteps\n"+ hr_merged + "\nnatom natomt ilin\n" + \
        clean_data_start[pos_first_atom-3] + "\n" + "\ncharge  spin  atomlabel\n" + clean_data_start[len(tau)+len(hr)+11]+ "\n"+matrix_merged + \
            "\nintcoor\n"+ list_merged + "\n" + "\n".join(clean_data_start[len(clean_data_start)-6:len(clean_data_start)-4]) + "\n\n" + \
                "\n".join(clean_data_start[len(clean_data_start)-4:len(clean_data_start)-1]) + "\n\nend"
text = open(r"final_gaussian.dat","w")
text.write(merged)
text.close()

os.mkdir("data")
os.system("cp final_gaussian.dat ./data/reac1.dat")
os.system("rm -rf final_gaussian.dat")
os.system("rm -rf matrix.xyz")
os.system("rm -rf gaussian_matrix")

os.chdir("data")
os.system("cp -r "+directory_getfile+"/* .") #Directory con moduli da aggiornare se cambia fonte
names = inchi + "\n" + smiles
text = open(r"name.dat","w")
text.write(names)
text.close()

exact_mol_wt = Descriptors.ExactMolWt(m)

clean = []
with open('me_head.dat') as file:
    lines = file.readlines() #Translate TextIOWrapper in a string
    for line in lines: clean.append(line.strip())# Cycle to clean the input data and work better with the matrix

for i in range(pos_first_atom,pos_first_atom+atom_num): start_matrix.append(list(filter("".__ne__,clean_data_start[i].split(" "))))
splitted_lines = []
for line in clean: splitted_lines.append(line.split(" "))

final_file = ""
for line in splitted_lines:
    if line[0] == "Masses[amu]": line[23] = str(exact_mol_wt)
    string = ""
    for element in line:
        if element == "": string = string + " "
        else: string = string +" "+ element
    final_file = final_file + string + "\n"

os.system("rm -f me_head.dat")
text = open("me_head.dat","w")
text.write(final_file)
text.close()

os.chdir("..")
