os.system("rm -rf bondenergies.out")
os.system("rm -rf bondenergies_hl.out")
os.system("rm -rf failed_original")
# Set how many fragments have been computed
with open(r"nfrags.out") as file:
    line = file.readlines()
    line = line[0].strip("\n").replace(" ","")
    ncomputedfrags = float(line)

# Check if the original molecule has been succesfully computed
os.chdir("originalmol")
check_l1 = False
try:
    checkHL,en,hlen,zpe = collectdata()
    check_l1 = True
    os.chdir("..")
except:
    os.chdir("..")
    text_fail = "Error: original molecule not computed correctly, can't compute bond energies"
    failed_originalmol = open("failed_original","w")
    failed_originalmol.write(text_fail)
    failed_originalmol.close()
    sys.exit()

# If the try-except is passed, at least level1 has been computed

bond_l1,bond_hl = [],[]
for i in range(int(ncomputedfrags)):

    checkHL1,checkHL2 = False,False
    # Collect and check which levels have been computed for frags 1 and 2
    os.chdir("frags_"+str(i+1)+"/fragment1")
    check_l1_frag1 = False
    try:
        checkHL1,en1,hlen1,zpe1 = collectdata()
        check_l1_frag1 = True
        os.chdir("..")
    except:
        os.chdir("..")

    os.chdir("fragment2")
    check_l1_frag2 = False
    try:
        checkHL2,en2,hlen2,zpe2 = collectdata()
        check_l1_frag2 = True
        os.chdir("../..")
    except:
        os.chdir("../..")

    if checkHL: # If the HL of the original molecule has been computed, the HL estimations can be made
        if checkHL1 and checkHL2:
            bond_hl.append(str((hlen1+zpe1+hlen2+zpe2-hlen-zpe)*627.509))
        elif not checkHL1 and checkHL2:
            bond_hl.append("failed 1")
        elif checkHL1 and not checkHL2:
            bond_hl.append("failed 2")
        else:
            bond_hl.append("failed 1 2")
    else:
        bond_hl = "Error: original molecule not computed with HL protocol"

    if check_l1: # If the HL of the original molecule has been computed, the HL estimations can be made

        if check_l1_frag1 and check_l1_frag2:
            bond_l1.append(str((en1+zpe1+en2+zpe2-en-zpe)*627.509))
        elif not check_l1_frag1 and check_l1_frag2:
            bond_l1.append("failed 1")
        elif check_l1_frag1 and not check_l1_frag2:
            bond_l1.append("failed 2")
        else:
            bond_l1.append("failed 1 2")
    else:
        bond_l1 = "Error: original molecule not computed at level1"

if not "Error" in bond_l1:
    text_l1 = ""
    for i in range(len(bond_l1)):
        text_l1 = text_l1 + str(i+1)+ " " + bond_l1[i] + "\n"
else:
    text_l1 = "Error: original molecule not computed at level1"

if not "Error" in bond_hl:
    text_hl = ""
    for i in range(len(bond_hl)):
        text_hl = text_hl + str(i+1)+ " " + bond_hl[i] + "\n"
else:
    text_hl = "Error: original molecule not computed with HL protocol"

level1 = open("bondenergies.out","w")
level1.write(text_l1)
level1.close()

hl = open("bondenergies_hl.out","w")
hl.write(text_hl)
hl.close()
[mferrar2@login03 frag_generation_splitted]$ clear
[mferrar2@login03 frag_generation_splitted]$ cat collect_frags.py
# Code to collect and calculate all the bond energy computed by molecule fragmentation

import os
import sys

####### Collect data (level1 or level1 and HL if computer) #######
def collectdata():
    with open("me_files/reac1_zpe.me") as zpe_file:
        line = zpe_file.readlines()
        line = line[0].strip("\n").replace(" ","")
        zpe = float(line)

    HL_computed = False
    try:
        with open("me_files/reac1_en.me") as en:
            line = en.readlines()
            line = line[0].strip("\n").replace(" ","")
            energyHL = float(line)
        with open("geoms/reac1_l1.xyz") as en:
            line = en.readlines()
            line = line[1].strip("\n").replace(" ","")
            energy = float(line)
        HL_computed = True
        return [True,energy,energyHL,zpe]
    except:
        with open("geoms/reac1_l1.xyz") as en:
            line = en.readlines()
            line = line[1].strip("\n").replace(" ","")
            energy = float(line)
        return [False,energy,0,zpe]
###################

os.system("rm -rf bondenergies.out")
os.system("rm -rf bondenergies_hl.out")
os.system("rm -rf failed_original")
# Set how many fragments have been computed
with open(r"nfrags.out") as file:
    line = file.readlines()
    line = line[0].strip("\n").replace(" ","")
    ncomputedfrags = float(line)

# Check if the original molecule has been succesfully computed
os.chdir("originalmol")
check_l1 = False
try:
    checkHL,en,hlen,zpe = collectdata()
    check_l1 = True
    os.chdir("..")
except:
    os.chdir("..")
    text_fail = "Error: original molecule not computed correctly, can't compute bond energies"
    failed_originalmol = open("failed_original","w")
    failed_originalmol.write(text_fail)
    failed_originalmol.close()
    sys.exit()

# If the try-except is passed, at least level1 has been computed

bond_l1,bond_hl = [],[]
for i in range(int(ncomputedfrags)):

    checkHL1,checkHL2 = False,False
    # Collect and check which levels have been computed for frags 1 and 2
    os.chdir("frags_"+str(i+1)+"/fragment1")
    check_l1_frag1 = False
    try:
        checkHL1,en1,hlen1,zpe1 = collectdata()
        check_l1_frag1 = True
        os.chdir("..")
    except:
        os.chdir("..")

    os.chdir("fragment2")
    check_l1_frag2 = False
    try:
        checkHL2,en2,hlen2,zpe2 = collectdata()
        check_l1_frag2 = True
        os.chdir("../..")
    except:
        os.chdir("../..")

    if checkHL: # If the HL of the original molecule has been computed, the HL estimations can be made
        if checkHL1 and checkHL2:
            bond_hl.append(str((hlen1+zpe1+hlen2+zpe2-hlen-zpe)*627.509))
        elif not checkHL1 and checkHL2:
            bond_hl.append("failed 1")
        elif checkHL1 and not checkHL2:
            bond_hl.append("failed 2")
        else:
            bond_hl.append("failed 1 2")
    else:
        bond_hl = "Error: original molecule not computed with HL protocol"

    if check_l1: # If the HL of the original molecule has been computed, the HL estimations can be made

        if check_l1_frag1 and check_l1_frag2:
            bond_l1.append(str((en1+zpe1+en2+zpe2-en-zpe)*627.509))
        elif not check_l1_frag1 and check_l1_frag2:
            bond_l1.append("failed 1")
        elif check_l1_frag1 and not check_l1_frag2:
            bond_l1.append("failed 2")
        else:
            bond_l1.append("failed 1 2")
    else:
        bond_l1 = "Error: original molecule not computed at level1"

if not "Error" in bond_l1:
    text_l1 = ""
    for i in range(len(bond_l1)):
        text_l1 = text_l1 + str(i+1)+ " " + bond_l1[i] + "\n"
else:
    text_l1 = "Error: original molecule not computed at level1"

if not "Error" in bond_hl:
    text_hl = ""
    for i in range(len(bond_hl)):
        text_hl = text_hl + str(i+1)+ " " + bond_hl[i] + "\n"
else:
    text_hl = "Error: original molecule not computed with HL protocol"

level1 = open("bondenergies.out","w")
level1.write(text_l1)
level1.close()

hl = open("bondenergies_hl.out","w")
hl.write(text_hl)
hl.close()
