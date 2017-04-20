import matplotlib.pyplot as plt
import numpy as np
import os

#########################################################################################################################
#  2016 AeroDyn standalone performance code interface, Robynne Murray, NREL
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
#  limitations under the License.
#########################################################################################################################

NumCases=14      #Specify the number of wind speed cases being modelled

TSR_curve=[0]*NumCases      #Allocate arrays
Cp_curve=[0]*NumCases
Ct_curve=[0]*NumCases

#########################################################################################################################
# Run the AeroDyn model #################################################################################################

print "  "
print "Running AeroDyn model....."
print "  "

os.system('AeroDyn_Driver_Win32.exe       NREL_MHK_driver.dvr')      #Execute the AeroDyn driver

## Obtain performance results from output files #########################################################################
#these depend on what you want to output from the code, this is an example of plotting the power and thrust coefficients

for i in range (1,NumCases):        #Open the output files in a loop depending on number of cases

    outputfile= np.loadtxt('NREL_MHK_driver.'+ str(i) + '.out',  skiprows=9, unpack=True, dtype={'names': ('time', 'speed', 'TSR','power', 'thrust', 'Cp', 'Ct'),'formats': ('f8', 'f8', 'f8', 'f8', 'f8', 'f8' ,'f8')})
    time= np.array(outputfile[0], copy=True)
    speed= np.array(outputfile[1], copy=True)
    TSR= np.array(outputfile[2], copy=True)
    power= np.array(outputfile[3], copy=True)
    thrust= np.array(outputfile[4], copy=True)
    Cp= np.array(outputfile[5], copy=True)
    Ct= np.array(outputfile[6], copy=True)

    TSR=np.mean(TSR)        #Take the mean at each case number to compare rotor power and thrust for each TSR
    Cp=np.mean(Cp)
    Ct=np.mean(Ct)

    TSR_curve[i] = TSR      #Save the mean values at each case number in an array
    Cp_curve[i] = Cp
    Ct_curve[i] = Ct

TSR_curve.remove(0)     #Remove the zeros in the first position on each array
Cp_curve.remove(0)
Ct_curve.remove(0)


# Experimental data to used as a comparison #############################################################################
# NRELS814 blade towing tank tests 1 m/s experimental data (Doman et al. 2015)

TSR_1ms = [5.154  ,  5.2481  ,  4.1788  ,  5.3368  ,  5.6008  ,  4.4512 ,   5.8118  ,  6.3936 ,   6.4383    ,6.47  ,  2.6807,
           3.8156  ,  2.9098  ,  4.0923   , 3.2547   , 3.395  ,  3.5304  ,  4.9203  ,  5.1858  ,  5.0863   , 5.1953  ,  5.461,
           4.2857  ,  5.69   , 5.9083  ,  4.6959   , 5.1415  ,  5.4787  ,  6.3389   , 6.5125   , 6.0076]

Cp_1ms = [0.197098098  ,  0.183734715  ,  0.276288124   , 0.172545535   , 0.139383428   , 0.263988854   , 0.103488738,
          0.012983471 ,   0.005088612 ,- 0.004697768  ,  0.083281111   , 0.285510447   , 0.115685477   , 0.277107375,
          0.235418338  ,  0.271346802  ,  0.285640201 ,   0.220800135  ,  0.190994755  ,  0.207992389 ,   0.194035496,
          0.156783173  ,  0.273998112  ,  0.12827253   , 0.092727219   , 0.24109398   , 0.194976516   , 0.159702106,
          0.028713836 , 0.004507393   , 0.081687159]

CT_1ms = [0.442395982  ,  0.439659786 ,   0.447107987  ,  0.437968211 ,   0.433157916  ,  0.452024199 ,   0.422243535,
          0.402088862 ,   0.402276456 ,   0.397961199 ,   0.297915076 ,   0.441426551  ,  0.319362235,    0.446535708,
          0.39207852  ,  0.417446943  ,  0.428406264  ,  0.444680865 ,   0.439730264   , 0.444171494  ,  0.444610215,
          0.434856675  ,  0.450929759  ,  0.42892612  ,  0.419905371  ,  0.445418213  ,  0.442165114 ,   0.437500535,
          0.405959554  ,  0.395807402  ,  0.418948347]

## Plots %###############################################################################################################
fig = plt.figure(1)
ax = fig.add_subplot(1, 1, 1)
ax.plot(TSR_1ms, Cp_1ms, 'bs', TSR_curve, Cp_curve, '-g',TSR_1ms, CT_1ms, 'gs', TSR_curve, Ct_curve, '-r')
ax.set_ylabel('Cp, Ct')
ax.set_xlabel('TSR')
plt.show()

# Save the figure:
fig.savefig('Cp and Ct vs TSR, Strathclyde turbine.pdf')




