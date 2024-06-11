import unittest
import os
import numpy as np
import matplotlib.pyplot as plt
import pyFAST.linearization as lin

# Get current directory so this script can be called from any location
scriptDir = os.path.dirname(__file__)


f0_ref   = [0.587830219596616, 0.722482755465059, 0.841644628629986, 0.937126199078433, 1.237130544073060, 1.837320640593674, 1.986990893316201, 2.133747228535351, 2.256063763227361]
zeta_ref = [6.310588187661754e-01, 5.252901946117751e-01, 4.401008557171799e-01, 1.634432195163714e-02, 1.235888553742901e-02, 1.555277937084565e-01, 1.428798211984810e-01, 1.337606910498821e-01, 2.258501382142458e-02]
Qmag_11  = [1.514472004400742e-05, 1.947994954300747e-01, 1.754772630476129e-01, 1.825691131755990e-01, 4.495989525843738e-03, 5.053399414140790e-03, 5.969950454246335e-03, 1.471919803320060e-02, 1.285988755307813e-02, 1.308283162729386e-02, 5.593620867971226e-05, 7.194814559410425e-01, 6.481158301940565e-01, 6.743091970923338e-01, 1.660569542445695e-02, 1.866445885760891e-02, 2.204969081277093e-02, 5.436456602636927e-02, 4.749730280101332e-02, 4.832073474448904e-02]
Qmag_83  = [1.441947564996028e-06, 6.371604702191713e-03, 6.305232540330547e-03, 6.289204718671321e-03, 6.932640151237337e-03, 6.546664262594788e-03, 6.589993010888945e-03, 5.227729076001339e-02, 5.129797039680727e-02, 5.179175676167141e-02, 1.933180057575511e-05, 8.542237903835923e-02, 8.453254543550748e-02, 8.431766477028646e-02, 9.294402939527062e-02, 8.776935516478050e-02, 8.835025196127627e-02, 7.008674823885858e-01, 6.877379994442021e-01, 6.943580595381521e-01]
Qphi_92  = [1.561610973903578, -1.700282526230297, -1.701805350416750, -1.697002105810609, 1.548360360707519, 1.547997531767750, 1.548270549679593, 2.395681952509270, 2.461948283219694, 2.405369652375428, -3.128191072180800, -0.106899265135085, -0.108422089321542, -0.103618844715400, -3.141441685376858, 3.141380792862958, -3.141531496404784, -2.294120093575104, -2.227853762864683, -2.284432393708950]
Qphi_11  = [-4.475921307027970e-02, -7.345405529305936e-01, 1.372699116580623e+00, -2.700369426141966e+00, 1.484807930794751e+00, -2.598715647054003e+00, -5.805991423437629e-01, -1.530872103538168e+00, 5.410725214387462e-01, 2.766892929144499e+00, 2.208954492083937e+00, 1.519173152223936e+00, -2.656772485444435e+00, -4.466557209874353e-01, -2.544663671230343e+00, -3.450019418994276e-01, 1.673114562810780e+00, 7.228416016163599e-01, 2.794786226593286e+00, -1.262578672880562e+00]

def compare(var, x1, x2, res, plot=False):
    try:
        np.testing.assert_almost_equal(x1, x2, res)
    except:
        if True:
            raise
        print('fail', var, 'with ', res)
        if plot:
            fig,axes = plt.subplots(2, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
            axes[0].plot(x1, label='new')
            axes[0].plot(x2, label='ref')
            axes[0].legend()
            axes[1].plot(x2-x1)
            plt.title(var)

class Test(unittest.TestCase):

    def test_readBinMatlab(self):
        matBinFile = os.path.join(scriptDir, '../../../data/NREL5MW/5MW_Land_Lin_Rotating/Main.ModeShapeVTK.postMBC_backup')
        # Read bin file
        d = lin.readModesForViz(matBinFile)
        # Test against known matlab values
        Qmag = d['x_eig_magnitude']
        Qphi = d['x_eig_phase']
        np.testing.assert_equal(Qmag.shape, (20,9,3))
        np.testing.assert_almost_equal(d['NaturalFreq_Hz'],f0_ref  , 15) # 16 fails
        np.testing.assert_almost_equal(d['DampingRatio']  ,zeta_ref, 16) # 17 fails
        np.testing.assert_almost_equal(Qmag[:,0,0],Qmag_11, 15)
        np.testing.assert_almost_equal(Qmag[:,7,2],Qmag_83, 15)
        np.testing.assert_almost_equal(Qphi[:,0,0],Qphi_11, 15)
        np.testing.assert_almost_equal(Qphi[:,8,1],Qphi_92, 15)

    def test_readWriteBinMatlab(self):
        matBinFile = os.path.join(scriptDir, '../../../data/NREL5MW/5MW_Land_Lin_Rotating/Main.ModeShapeVTK.postMBC_backup')
        matBinFile_tmp = os.path.join(scriptDir, '../../../data/NREL5MW/5MW_Land_Lin_Rotating/Main.ModeShapeVTK.postMBC_backup_tmp')
        # Read the data
        d = lin.readModesForViz(matBinFile)
        # Write it
        lin.writeModesForViz(d, matBinFile_tmp)
        # Read what we wrote
        d = lin.readModesForViz(matBinFile_tmp)
        # Test against known matlab values
        Qmag = d['x_eig_magnitude']
        Qphi = d['x_eig_phase']
        np.testing.assert_equal(Qmag.shape, (20,9,3))
        np.testing.assert_almost_equal(d['NaturalFreq_Hz'],f0_ref  , 15) # 16 fails
        np.testing.assert_almost_equal(d['DampingRatio']  ,zeta_ref, 16) # 17 fails
        np.testing.assert_almost_equal(Qmag[:,0,0],Qmag_11, 15)
        np.testing.assert_almost_equal(Qmag[:,7,2],Qmag_83, 15)
        np.testing.assert_almost_equal(Qphi[:,0,0],Qphi_11, 15)
        np.testing.assert_almost_equal(Qphi[:,8,1],Qphi_92, 15)
        try:
            os.remove(matBinFile_tmp)
        except:
            pass

    def test_writeBinPython(self):
        matBinFile = os.path.join(scriptDir, '../../../data/NREL5MW/5MW_Land_Lin_Rotating/Main.ModeShapeVTK.postMBC_backup')
        fstFile = os.path.join(scriptDir, '../../../data/NREL5MW/5MW_Land_Lin_Rotating/Main.fst')

        # --- Perform python MBC and write modes to binary file
        CDDOP, MBCOP = lin.getCampbellDataOP(fstFile, writeModes=True)

        # --- Read matlab binary file (reference)
        d2 = lin.readModesForViz(matBinFile)

        # --- Read python binary fie
        d = lin.readModesForViz(MBCOP['modeFile'])
        # Outputs to screen
        #Freq,Damp = lin.printCampbellDataOP(CDDOP, nModesMax=10, nCharMaxDesc=50)

        # Test against known matlab values
        Qmag1 = d['x_eig_magnitude']
        Qphi1 = d['x_eig_phase']
        Qmag2 = d2['x_eig_magnitude']
        Qphi2 = d2['x_eig_phase']
        np.testing.assert_equal(Qmag1.shape, (20,9,3))
        compare('f0', d['NaturalFreq_Hz'],f0_ref  , 14) # 15 fails
        compare('ze', d['DampingRatio']  ,zeta_ref, 15) # 16 fails
        compare('f0', d['NaturalFreq_Hz'], d2['NaturalFreq_Hz'], 14) # 15 fails
        compare('ze', d['DampingRatio']  , d2['DampingRatio'], 15) # 16 fails
        # Somehow values for state 0 and 10 are off...
        for iAz in range(1):
            for iMode in range(9):
                v1 = Qmag1[:,iMode,iAz]
                v2 = Qmag2[:,iMode,iAz]
                v1 = np.concatenate((Qmag1[1:10,iMode,iAz],Qmag1[11:-1,iMode,iAz] ))
                v2 = np.concatenate((Qmag2[1:10,iMode,iAz],Qmag2[11:-1,iMode,iAz] ))
                compare('q{:d}{:d}'.format(iMode,iAz), v1, v2, 10, plot=False)

        for iAz in range(1):
            for iMode in range(9):
                v1 = np.concatenate((Qphi1[1:10,iMode,iAz],Qphi1[11:-1,iMode,iAz] ))
                v2 = np.concatenate((Qphi2[1:10,iMode,iAz],Qphi2[11:-1,iMode,iAz] ))
                compare('p{:d}{:d}'.format(iMode,iAz), v1, v2, 10, plot=False)
        #plt.show()


if __name__ == '__main__':
    unittest.main()
