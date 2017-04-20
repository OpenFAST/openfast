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

NumCases=16    #Specify the number of wind speed cases being modelled

TSR_curve=[0]*NumCases      #Allocate arrays
Cp_curve=[0]*NumCases
Ct_curve=[0]*NumCases

#########################################################################################################################
# Run the AeroDyn model #################################################################################################

print "  "
print "Running AeroDyn model....."
print "  "

os.system('AeroDyn_Driver_Win32.exe       Bahaj_MHK_driver.dvr')      #Execute the AeroDyn driver

## Obtain performance results from output files ##########################################################################
#these depend on what you want to output from the code, this is an example of plotting the power and thrust coefficients

for i in range (1,NumCases):        #Open the output files in a loop depending on number of cases

    outputfile= np.loadtxt('Bahaj_MHK_driver.'+ str(i) + '.out',  skiprows=9, unpack=True, dtype={'names': ('time', 'speed', 'TSR','power', 'thrust', 'Cp', 'Ct'),'formats': ('f8', 'f8', 'f8', 'f8', 'f8', 'f8' ,'f8')})
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

    TSR_curve[i] = TSR       #Save the mean values at each case number in an array
    Cp_curve[i] = Cp
    Ct_curve[i] = Ct

TSR_curve.remove(0)      #Remove the zeros in the first position on each array
Cp_curve.remove(0)
Ct_curve.remove(0)

# Experimental data to used as a comparison #############################################################################
# NACA63815 blade cavitation tunnel tests 1.73 m/s experimental data (Bahaj et al. 2007)
Cp_dataTSR173 = [4.17061600000000, 4.42338100000000 ,4.66034800000000, 4.89731400000000 ,5.13428100000000 ,5.37124800000000, 5.37124800000000 ,5.59241700000000, 5.84518200000000, 6.08214800000000 ,6.30331800000000 ,6.54028400000000, 6.77725100000000 ,7.01421800000000 ,7.21958900000000, 7.44075800000000, 7.69352300000000 ,8.15165900000000 ,8.60979500000000 ,9.06793000000000 ,9.54186400000000 ,9.87361800000000]
Cp_data173 = [0.413793000000000, 0.430885000000000, 0.437181000000000, 0.446177000000000 ,0.445277000000000, 0.454273000000000 ,0.457871000000000 ,0.452474000000000 ,0.454273000000000 ,0.452474000000000 ,0.449775000000000 ,0.447976000000000 ,0.441679000000000 ,0.435382000000000, 0.429085000000000, 0.412894000000000, 0.409295000000000, 0.389505000000000 ,0.361619000000000, 0.331034000000000, 0.287856000000000 ,0.258171000000000]
Ct_dataTSR173 = [4.18495300000000, 4.43573700000000 ,4.67084600000000, 4.90595600000000, 5.14106600000000 ,5.39185000000000, 5.37617600000000 ,5.61128500000000, 5.86206900000000 ,5.87774300000000 ,6.09717900000000 ,6.09717900000000 ,6.31661400000000 ,6.55172400000000 ,6.78683400000000 ,7.02194400000000 ,7.25705300000000, 7.46081500000000 ,7.71159900000000 ,8.16614400000000 ,8.63636400000000,9.09090900000000, 9.54545500000000 ,9.89028200000000]
Ct_data173 = [0.644510000000000, 0.672997000000000, 0.701484000000000 ,0.726409000000000, 0.740653000000000, 0.753116000000000 ,0.772700000000000, 0.778042000000000 ,0.795846000000000 ,0.802967000000000, 0.811869000000000 ,0.817211000000000 ,0.820772000000000 ,0.836795000000000, 0.851039000000000 ,0.859941000000000, 0.874184000000000, 0.884866000000000 ,0.890208000000000, 0.902671000000000, 0.915134000000000 ,0.925816000000000, 0.925816000000000 ,0.938279000000000]

## Plots ###############################################################################################################
fig = plt.figure(1)
ax = fig.add_subplot(1, 1, 1)
ax.plot(Cp_dataTSR173, Cp_data173, 'bs', TSR_curve, Cp_curve, '-g',Ct_dataTSR173, Ct_data173, 'gs', TSR_curve, Ct_curve, '-r')
ax.set_ylabel('Cp, Ct')
ax.set_xlabel('TSR')
plt.show()

# Save the figure:
fig.savefig('Cp and Ct vs TSR, Bahaj at al. turbine.pdf')




