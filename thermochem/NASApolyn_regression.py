# This is a code for the automatic estimation of the NASA coefficient in the CHEM-KIN format
# They are composed of two sets of seven coefficients each; the first set describes the thermodynamic properties
# at high temperature (usually higher then 1000 K), while the second set describes the thermodynamic properties
# at low temperature.
# The input of this code consists of:
# - A series of cp estimated with the thermo subroutine, spaced with a relatively small dT (10 K usually)
# - Two values (one at high T and one at low T) of entropy (to obtain the low and high entropy coefficient)
# - The value of DH0f(0 K) estimated with the CBH code
# *Two regressions in series are performed with each set of coefficients with "numpy.polyfit" function, in order
#  to have an estimation of the first five coefficients of each set.
#  Then a second series of regression is performed using "optimize.curve_fit" function, with Trust Region Reflective method,
#  obtaining the definitive set of the first five coefficients, for each temperature interval.
# *The sixth coefficient of each interval is obtained using the CBH estimated DH0f:
# -For the low T interval, the DH0f(0 K) is directly corrected to give DH0f(298.15 K), and the sixth coefficient is obtained
#  using the explicit formula
# -For the high T interval, the estimated DH0f(@HighT) is obatined exploiting the relation between the enthalpy and the specific heat:
#  DH0f(@HighT) = DH0f(298.15 K) + integral(298.15 K -> Tsplit) [ cp_lowT ] + integral(Tsplit -> HighT) [ cp_highT ]
#  since the coefficients of cp_lowT and cp_highT have been estimated before, using the double regression
# *The seventh coefficient for each interval is obtained using the explicit expression, having two S data (low and high T)
#  obtained using the thermo subroutine.
# Once all the coefficients have been estimated they are saved in a CHEM-KIN format, in the thermo directory.
# The 15th number appearing in the CHEM-KIN format is the DH0f(298.15 K), which is added as bonus information in the most modern format of CHEM-KIN
# The CHEM-KIN format of the NASA polynomials provides also the molecular formula, the phase of the molecule (G, gas, standard for this code),
# the low T boundary, the high T boundary and the split temperature, all in this order in the first row of the CHEM-KIN.
# The second row is composed by the high T coefficients from 1 to 5; the third row is composed by the sixth and seventh high T coefficients and
# by the low T coefficients from 1 to 3; the fourth row is composed by the low T coefficients from 4 to 7 and from the DH0f(298.15)
# The DH0f(298.15) is expressed in [cal/mol]

# cp data provided in [cal/mol/K] and S data provided in [cal/mol]
# CHEM-KIN NASA format will give parameters for the estimation of cp, S and H @T in consistent units [cal], [mol], [K]

# The standard form of the NASA polynomials is as follow (ex for ethane):
'''
C2H6                                        G   200.000  3500.000  1860.000    1
 1.07141839e+00 2.16759611e-02-1.00212824e-05 2.21316648e-09-1.89921174e-13    2
-9.26431727e+03 1.51090904e+01 1.22401699e+00 1.86466552e-02-3.89501801e-06    3
-1.99855390e-09 7.61794965e-13-9.08368307e+03 1.52122522e+01-1.57528560e+04    4
'''
# The functional dependance of cp, H and S with respect to T is as follow:

# cp/R = a0 + a1*T + a2*T^2 + a3*T^3 + a4*T^4

# H/(RT) = a0 + 1/2*a1*T + 1/3*a2*T^2 + 1/4*a3*T^3 + 1/5*a4*T^4 + a5/T

# S/R = a0*ln(T) + a1*T + 1/2*a2*T^2 + 1/3*a3*T^3 + 1/4*a4*T^4 + a6

'''
The code needs a working rdkit environment and the InChI identificator in the directory ./data/name.dat in order to compute the molecular formula,
which is mandatory to have a working canonical CHEM-KIN NASA format
'''

import os
import sys
import numpy
import scipy
from scipy import optimize
from rdkit import Chem
from rdkit.Chem import AllChem

# Constants definition
R = 1.987 # Universal gas constant [cal/mol/K]
kb = 1.380649 * 10**(-23) # Boltzmann constant [J/K]
Na = 6.022 * 10**23 # Avogadro number [1/mol]

# Define the phase of parameters estimation -> By default is considered gas phase
phase_molecule = "G"

########################## FUNCTIONS DEFINITION ##########################
# Functional dependance of cp
def calculate_cp(T,a1,a2,a3,a4,a5):
    cp = (a1 + a2*T + a3*T**2 + a4*T**3 + a5*T**4)
    return cp

# The most suitable split temperature is automatically established by the code, using as first guess the middle data temperature provided by the cp data file
# Check if the join of low and high temperature interval is consistent (no change in derivative)
def check_joinTinterval(cp_splitT,cp_splitT_minus1,cp_splitT_minus2,cp_spliT_plus1,cp_spliT_plus2,):
    if cp_splitT_minus2<=cp_splitT_minus1 and cp_splitT_minus1<=cp_splitT and cp_splitT<=cp_spliT_plus1 and cp_spliT_plus1<=cp_spliT_plus2: return True
    elif cp_splitT_minus2>=cp_splitT_minus1 and cp_splitT_minus1>=cp_splitT and cp_splitT>=cp_spliT_plus1 and cp_spliT_plus1>=cp_spliT_plus2: return True
    else: return False
##########################################################################

# Define the molecular formula exploiting rdkit function
try: # Define the InChI from data/name.dat
    with open("./data/name.dat") as name:
        line = name.readlines()
        inchi = line[0].strip("\n").replace(" ","")
except: # If no InChI is found, the CHEM-KIN format is not computed because it's not valid
    os.chdir("thermo")
    text_error = "Error: no InChI found in 'data/name.dat'; can't compute molecular formula."
    error = open("nasa_polyn.out","w")
    error.write(text_error)
    error.close()
    sys.exit() # Stop the regression process if no InChI is found in data/name.dat

# Define the mol file and convert it into the molecular formula
mol = Chem.MolFromInchi(inchi)
molecular_formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
'''
# Take the value of DH0f(0 K) either from me_files/reac1_dh0f_0K.me (for HL) or output/reac1_dh0f_0K.out (for level1)
# The most precise value is taken
try: # Check if HL estimation is present
    with open("./me_files/reac1_dh0f_0K.me") as hldh0:
        line = hldh0.readlines()
        dh0k =  float(line[0].strip("\n").replace(" ","")) #[kcal/mol]
except:
    try: # Check if level1 estimation is present
        with open("./output/reac1_dh0f_0K.out") as dh0:
            line = dh0.readlines()
            dh0k =  float(line[0].strip("\n").replace(" ","")) #[kcal/mol]
    except: # No ΔH0f(0 K) is found, an error is thrown out
        os.chdir("thermo")
        text_error = "Error: no estimated ΔH0f(0 K) present in 'me_files' or 'output', can't compute NASA polynomials"
        error = open("nasa_polyn.out","w")
        error.write(text_error)
        error.close()
        sys.exit() # Stop the regression process if no ΔH0f(0 K) is found

# Extract the values of Eel and ZPE to correct ΔH0f(0 K) -> ΔH0f(298.15 K)
try: # The high level has been computed -> We can estimate DH0 at level1 and also at HL
    with open("./me_files/reac1_en.me") as en:
        line = en.readlines()
        line = line[0].strip("\n").replace(" ","")
        energy = float(line)
except:
    try:
        with open("geoms/reac1_l1.xyz") as en:
            line = en.readlines()
            line = line[1].strip("\n").replace(" ","")
            energy = float(line)
    except:
        os.chdir("thermo")
        text_error = "Error: no electronic energy present in 'me_files' or 'output', can't compute NASA polynomials"
        error = open("nasa_polyn.out","w")
        error.write(text_error)
        error.close()
        sys.exit() # Stop the regression process if no Eel is found

try:
    with open("me_files/reac1_zpe.me") as zpe_file:
        line = zpe_file.readlines()
        line = line[0].strip("\n").replace(" ","")
        zpe = float(line)
except:
    os.chdir("thermo")
    text_error = "Error: no zpe present in 'me_files', can't compute NASA polynomials"
    error = open("nasa_polyn.out","w")
    error.write(text_error)
    error.close()
    sys.close() # Stop the regression process of no zpe is found
'''
# Extract the value of ΔH0f(298.15 K) after correction of ΔH0f(0 K)

try:
    # Try to extract HL estimation
    with open("./me_files/reac1_DH298K.me") as file:
        line = file.readlines()
        dh298 =  float(line[0].strip("\n").replace(" ","")) #[kcal/mol]

except:
    try: # Try to extract level1 estimation
        with open("./output/reac1_DH298K.out") as file:
            line = file.readlines()
            dh298 =  float(line[0].strip("\n").replace(" ","")) #[kcal/mol]
    except:
        os.chdir("thermo")
        text_error = "Error: no ΔH0f(298.15 K) found. Can't estimate NASA polynomials"
        error = open("nasa_polyn.out","w")
        error.write(text_error)
        error.close()
        sys.exit() # Stop the regression process if no cp data are found

# Open the files containing cp and S data
raw_data,cp_data = [],[]

try: # Extract the cp data
    with open("./thermo/cp.out") as data:
        file = data.readlines()
    for i in file: raw_data.append(i.strip())
    intermediate_data = list(filter("".__ne__,raw_data))
    for j in range(len(intermediate_data)):
        cp_data.append(list(filter("".__ne__,intermediate_data[j].split(" "))))
except: # If no cp file is found, throw an error
    os.chdir("thermo")
    text_error = "Error: no cp data found."
    error = open("nasa_polyn.out","w")
    error.write(text_error)
    error.close()
    sys.exit() # Stop the regression process if no cp data are found

raw_data,S_data = [],[]
try: # Extract the S data
    with open("./thermo/S.out") as data:
        file = data.readlines()
    for i in file: raw_data.append(i.strip())
    intermediate_data = list(filter("".__ne__,raw_data))
    for j in range(len(intermediate_data)):
        S_data.append(list(filter("".__ne__,intermediate_data[j].split(" "))))
except: # If no S file is found, throw an error
    os.chdir("thermo")
    text_error = "Error: no S data found."
    error = open("nasa_polyn.out","w")
    error.write(text_error)
    error.close()
    sys.exit() # Stop the regression process if no S data are found

# Check if only two S data are present
if len(S_data) != 2:
    os.chdir("thermo")
    text_error = "Error: more then two S data found. Need exactly two."
    error = open("nasa_polyn.out","w")
    error.write(text_error)
    error.close()
    sys.exit() # Stop the regression process if more then 2 S data are found

# A preliminary regression is made, in order to check in which direction the automatic regression process should move (higher or lower than the first guess temperature)
low_T = []
high_T = []

low_cp = []
high_cp = []

# Find the position of the middle temperature ad choose two temperatures (one higher and one lower) to evaluate in which direction we need to move
if len(cp_data) % 2 == 0: # Even number of cp data
    T_low = cp_data[int(len(cp_data)/2-5)][0]
    T_med = cp_data[int(len(cp_data)/2)][0]
    T_high = cp_data[int(len(cp_data)/2+5)][0]
else: # Odd number of cp data
    T_low = cp_data[int((len(cp_data)+1)/2-5)][0]
    T_med = cp_data[int((len(cp_data)+1)/2)][0]
    T_high = cp_data[int((len(cp_data)+1)/2+5)][0]

# Save the values of R2 at the three different temperatures to evaluate in which direction (lower or higher than middle T) the final regression should be made
R2_cumulative_check = []
T_splitrange = [T_low,T_med,T_high]

for T in T_splitrange:
    #Create a set of low T and high T data, accoarding to T split range
    for j in range(len(cp_data)):

        # Save low and high temperature cp data
        if cp_data[j][1]<=T:
            low_T.append(float(cp_data[j][0]))
            low_cp.append(float(cp_data[j][1]))
        else:
            high_T.append(float(cp_data[j][0]))
            high_cp.append(float(cp_data[j][1]))

    #Fit a1->a5 with cp data on cp at low T
    coeff_lowT = numpy.polyfit(low_T,low_cp,4)

    #Fit a1->a5 with cp data on cp at high T
    coeff_highT = numpy.polyfit(high_T,high_cp,4)

    # Low T R2 evaluation
    cpR_mean_lowT = numpy.mean(low_cp)/len(low_cp)

    SS_res_cpR_lowT = 0
    SS_tot_cpR_lowT = 0

    # R2 and sigma2 of cp regression calculation
    for i in range(len(low_cp)):
        cpR_lowT_model = calculate_cp(low_T[i],coeff_lowT[4],coeff_lowT[3],coeff_lowT[2],coeff_lowT[1],coeff_lowT[0]) # [cal/mol/K]

        SS_res_cpR_lowT = SS_res_cpR_lowT + (low_cp[i]-cpR_lowT_model)*(low_cp[i]-cpR_lowT_model)
        SS_tot_cpR_lowT = SS_tot_cpR_lowT + (low_cp[i]-cpR_mean_lowT)*(low_cp[i]-cpR_mean_lowT)

    R2_cpR_lowT = 1-SS_res_cpR_lowT/SS_tot_cpR_lowT

    # High T R2 evaluation
    cpR_mean_highT = numpy.mean(high_cp)/len(high_cp)

    SS_res_cpR_highT = 0
    SS_tot_cpR_highT = 0

    # R2 and sigma2 of cp regression calculation
    for i in range(len(high_cp)):
        cpR_highT_model = calculate_cp(high_T[i],coeff_highT[4],coeff_highT[3],coeff_highT[2],coeff_highT[1],coeff_highT[0]) # [cal/mol/K]

        SS_res_cpR_highT = SS_res_cpR_highT + (high_cp[i]-cpR_highT_model)*(high_cp[i]-cpR_highT_model)
        SS_tot_cpR_highT = SS_tot_cpR_highT + (high_cp[i]-cpR_mean_highT)*(high_cp[i]-cpR_mean_highT)

    R2_cpR_highT = 1-SS_res_cpR_highT/SS_tot_cpR_highT

    R2_cumulative_check.append(R2_cpR_highT+R2_cpR_lowT)

# The direction of the recursive regression process is established on the basis of the R2 results
if R2_cumulative_check[2] > R2_cumulative_check[1]: where_to_go = "higher"
if R2_cumulative_check[0] > R2_cumulative_check[1]: where_to_go = "lower"

check_join = True

# The middle position is defined, while the direction has been already established by the previous preliminary regression
if len(cp_data) % 2 == 0: position = int(len(cp_data)/2)
else: position = int((len(cp_data)+1)/2)

while check_join:
    low_T = []
    high_T = []

    low_cp = []
    high_cp = []

    # Create the set of low and high T data
    for j in range(position+10):
        low_T.append(float(cp_data[j][0]))
        low_cp.append(float(cp_data[j][1]))
    for k in range(position-10,len(cp_data)):
        high_T.append(float(cp_data[k][0]))
        high_cp.append(float(cp_data[k][1]))

    #Fit a1->a5 with data on cp at low T
    coeff_lowT = numpy.polyfit(low_T,low_cp,4)
    fg_lowT = [coeff_lowT[4], coeff_lowT[3],coeff_lowT[2],coeff_lowT[1],coeff_lowT[0]]
    coeff_lowT , cov1_lowT = optimize.curve_fit(calculate_cp, xdata = low_T, ydata = low_cp, p0 = fg_lowT, method = "trf")

    #Fit a1->a5 with data on cp at high T
    coeff_highT = numpy.polyfit(high_T,high_cp,4)
    fg_highT = [coeff_highT[4], coeff_highT[3],coeff_highT[2],coeff_highT[1],coeff_highT[0]]
    coeff_highT , cov1_lowT = optimize.curve_fit(calculate_cp, xdata = high_T, ydata = high_cp, p0 = fg_highT, method = "trf")

    # Five values of cp are calculated (around the actual T considered) in order to establish if the temperature is correct for joining low nad high T intervals
    cp_splitT_minus2 = calculate_cp(float(cp_data[position-2][0]),coeff_lowT[4],coeff_lowT[3],coeff_lowT[2],coeff_lowT[1],coeff_lowT[0])
    cp_splitT_minus1 = calculate_cp(float(cp_data[position-1][0]),coeff_lowT[4],coeff_lowT[3],coeff_lowT[2],coeff_lowT[1],coeff_lowT[0])
    cp_splitT = calculate_cp(float(cp_data[position][0]),coeff_lowT[4],coeff_lowT[3],coeff_lowT[2],coeff_lowT[1],coeff_lowT[0])
    cp_spliT_plus1 = calculate_cp(float(cp_data[position+1][0]),coeff_highT[4],coeff_highT[3],coeff_highT[2],coeff_highT[1],coeff_highT[0])
    cp_spliT_plus2 = calculate_cp(float(cp_data[position+2][0]),coeff_highT[4],coeff_highT[3],coeff_highT[2],coeff_highT[1],coeff_highT[0])

    check_join = check_joinTinterval(cp_splitT,cp_splitT_minus1,cp_splitT_minus2,cp_spliT_plus1,cp_spliT_plus2)

    if check_join == False: # If the join is satisfactory, the regression process is stopped
        where_to_go = "Done"
        Tsplit = float(cp_data[position][0])

    # If the regression process should continue, the position of the split temperature is shifted according to the preliminary R2 evaluation
    if where_to_go == "higher": position = position + 1
    elif where_to_go == "lower": position = position - 1

# Save the two sets of coefficients calculated
a_lowT = [coeff_lowT[0]/R,coeff_lowT[1]/R,coeff_lowT[2]/R,coeff_lowT[3]/R,coeff_lowT[4]/R,"",""]
a_highT = [coeff_highT[0]/R,coeff_highT[1]/R,coeff_highT[2]/R,coeff_highT[3]/R,coeff_highT[4]/R,"",""]

# Calculate the seventh coefficient with explicit formula for S, both low and high temperature
TS_low,Slow = float(S_data[0][0]),float(S_data[0][1])
TS_high,Shigh =float(S_data[1][0]),float(S_data[1][1])

a_lowT[6] = str(Slow/R - a_lowT[0]*numpy.log(TS_low) - a_lowT[1]*TS_low - 1/2 * a_lowT[2]*TS_low*TS_low - 1/3 * a_lowT[3]*TS_low*TS_low*TS_low - 1/4 * a_lowT[4]*TS_low*TS_low*TS_low*TS_low)
a_highT[6]= str(Shigh/R - a_highT[0]*numpy.log(TS_high) - a_highT[1]*TS_high - 1/2 * a_highT[2]*TS_high*TS_high - 1/3 * a_highT[3]*TS_high*TS_high*TS_high - 1/4 * a_highT[4]*TS_high*TS_high*TS_high*TS_high)

# Calculate the sixth coefficient with explicit formula for H, both low and high temperature
TH_low,Hlow = 298.15,dh298*1000

# To obtain TH_high, the relation between cp and H is exploited
TH_high = float(S_data[1][0])
Hhigh = Hlow + R*(a_lowT[0]*(Tsplit-298.15) + 1/2*a_lowT[1]*(Tsplit**2-298.15**2) + 1/3*a_lowT[2]*(Tsplit**3-298.15**3) + 1/4*a_lowT[3]*(Tsplit**4-298.15**4) + 1/5*a_lowT[4]*(Tsplit**5-298.15**5))+\
        R*(a_highT[0]*(TH_high-Tsplit) + 1/2*a_highT[1]*(TH_high**2-Tsplit**2) + 1/3*a_highT[2]*(TH_high**3-Tsplit**3) + 1/4*a_highT[3]*(TH_high**4-Tsplit**4) + 1/5*a_highT[4]*(TH_high**5-Tsplit**5))

a_lowT[5] = str(Hlow/R - a_lowT[0]*TH_low - 1/2*a_lowT[1]*TH_low**2 - 1/3 * a_lowT[2]*TH_low**3 - 1/4 * a_lowT[3]*TH_low**4 - 1/5 * a_lowT[4]*TH_low**5)
a_highT[5] = str(Hhigh/R - a_highT[0]*TH_high - 1/2*a_highT[1]*TH_high**2 - 1/3 * a_highT[2]*TH_high**3 - 1/4 * a_highT[3]*TH_high**4 - 1/5 * a_highT[4]*TH_high**5)

coeff_lowT_final,coeff_highT_final = a_lowT,a_highT

# Formatting the coefficients in CHEMKIN format
highT_coeff_arranged = []
for i in range(len(coeff_highT_final)):
    formatted_number = "{:.8e}".format(float(coeff_highT_final[i]))
    if float(coeff_highT_final[i]) >=0.:
        highT_coeff_arranged.append(" "+str(formatted_number))
    else:
        highT_coeff_arranged.append(str(formatted_number))

lowT_coeff_arranged = []
for i in range(len(coeff_lowT_final)):
    formatted_number = "{:.8e}".format(float(coeff_lowT_final[i]))
    if float(coeff_lowT_final[i]) >=0.:
        lowT_coeff_arranged.append(" "+str(formatted_number))
    else:
        lowT_coeff_arranged.append(str(formatted_number))

# Formatting the first line
first_line = [" "]*80
first_line[len(first_line)-1] = "1"
first_line[44] = phase_molecule
for j in range(len(molecular_formula)): first_line[j] = molecular_formula[j]

formatted_lowT = str("{:.3f}".format(low_T[0]))
formatted_highT = str("{:.3f}".format(high_T[len(high_T)-1]))
formatted_splitT = str("{:.3f}".format(float(Tsplit)))

pos = 0
for i in range(48,48+len(formatted_lowT)):
    first_line[i] = formatted_lowT[pos]
    pos = pos + 1

pos = 0
for i in range(57,57+len(formatted_highT)):
    first_line[i] = formatted_highT[pos]
    pos = pos + 1

pos = 0
for i in range(67,67+len(formatted_splitT)):
    first_line[i] = formatted_splitT[pos]
    pos = pos + 1

formatted_dh0 = "{:.8e}".format(float(dh298))
if dh298 >=0.:
    dh0_chemkin = " "+formatted_dh0
else:
    dh0_chemkin = formatted_dh0

first_line_merged = "".join(map(str,first_line))
merged = (first_line_merged+"\n"+ \
          highT_coeff_arranged[0]+highT_coeff_arranged[1]+highT_coeff_arranged[2]+highT_coeff_arranged[3]+highT_coeff_arranged[4]+"    2"+"\n"+ \
          highT_coeff_arranged[5]+highT_coeff_arranged[6]+lowT_coeff_arranged[0]+lowT_coeff_arranged[1]+lowT_coeff_arranged[2]+"    3"+"\n"+ \
          lowT_coeff_arranged[3]+lowT_coeff_arranged[4]+lowT_coeff_arranged[5]+lowT_coeff_arranged[6]+dh0_chemkin+"    4")+"\n"

os.chdir("thermo")
text = open("nasa_polyn.out","w")
text.write(merged)
text.close()
