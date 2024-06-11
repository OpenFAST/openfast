"""

"""
import numpy as np
import scipy.linalg as scp
import os
from pyFAST.linearization.linfile import get_Mats


# --------------------------------------------------------------------------------}
# --- Utils 
# --------------------------------------------------------------------------------{
def get_tt_inverse(sin_col, cos_col):

    c1 = cos_col[0];
    c2 = cos_col[1];
    c3 = cos_col[2];
    
    s1 = sin_col[0];
    s2 = sin_col[1];
    s3 = sin_col[2];
    
    ttv = [ [c2*s3 - s2*c3,  c3*s1 - s3*c1, c1*s2 - s1*c2],
            [  s2 - s3 ,       s3 - s1,       s1 - s2 ],
            [  c3 - c2 ,       c1 - c3,       c2 - c1 ] ]
    ttv = ttv/(1.5*np.sqrt(3));

    return ttv

def get_new_seq(rot_triplet,ntot):
#  rot_triplet is size n x 3
    #print('rot_triplet ', len(rot_triplet))
    rot_triplet=np.array(rot_triplet,dtype=int)
    if (rot_triplet.size==0):
        nRotTriplets=0
        nb=0
        #return np.array([]),0,0;
    else:
        nRotTriplets,nb = rot_triplet.shape;

    #print(nRotTriplets,nb,rot_triplet.flatten())
    
    if (nb != 3 and nRotTriplets != 0 and ntot!= 0):
        print('**ERROR: the number of column vectors in the rotating triplet must equal 3, the num of blades');
        new_seq = np.range(1,ntot)
    else:
        non_rotating = np.ones(ntot,dtype=int);
        #print(non_rotating)
        non_rotating[rot_triplet.flatten()] = 0; # if they are rotating, set them false;
        a=np.array(np.nonzero(non_rotating)).flatten()
        b=(rot_triplet.reshape(nRotTriplets*nb, 1)).flatten()
        new_seq = np.concatenate((a,b));

        #print(new_seq)
    return new_seq,nRotTriplets,nb

def eiganalysis(A, ndof2=None, ndof1=None):
    """ """
    mbc={}
    m, ns = A.shape;
    if(m!=ns):
        raise Exception('**ERROR: the state-space matrix is not a square matrix.');

    if ndof2 is None:
        # Assume that all DOF are second order
        ndof1 = 0;
        ndof2 = int(ns/2)

        if np.mod(ns,2) != 0:
            raise Exception('**ERROR: the input matrix is not of even order.');
    elif ndof1 is None:
        # Assume that first order states are "reminder" 
        ndof1 = ns - 2*ndof2;
        if ndof1 < 0:
            raise Exception('**ERROR: ndof2 must be no larger than half the dimension of the state-space matrix.');      
    else:
        if ns != int(2*ndof2 + ndof1):
            raise Exception('**ERROR: the dimension of the state-space matrix must equal 2*ndof2 + ndof1.');

    ndof = ndof2 + ndof1;

    # Matlab code - KEEP ME
    #[origEigenVects, origEvals] = eig(A,'vector'); %,'nobalance'
    #positiveImagEvals = find( imag(origEvals) > 0);
    #mbc.Evals             = origEvals(positiveImagEvals);
    #mbc.EigenVects        = origEigenVects([1:ndof2  (ndof2*2+1):ns],positiveImagEvals); % save q2 and q1, throw away q2_dot
    #    EigenVects_save   = origEigenVects(:,positiveImagEvals); % save these for VTK visualization;
    #real_Evals = real(mbc.Evals);
    #imag_Evals = imag(mbc.Evals);
    #mbc.NaturalFrequencies = sqrt( real_Evals.^2 + imag_Evals.^2 );

    origEvals, origEigenVects = np.linalg.eig(A); #,'nobalance'
    # errorInSolution = norm(A * mbc.EigenVects - mbc.EigenVects* diag(mbc.EigenVals) )
    # these eigenvalues aren't sorted, so we just take the ones with
    # positive imaginary parts to get the pairs for modes with damping < 1:
    positiveImagEvals = np.argwhere( origEvals.imag > 0.0);
    
    mbc['Evals']             = origEvals[positiveImagEvals];
    row=np.array(list(range(0,ndof2))+list(range(ndof2*2+1,ns)))
    col=positiveImagEvals

    mbc['EigenVects']=origEigenVects[row,col].transpose(); # save q2 and q1, throw away q2_dot
    EigenVects_save=origEigenVects[:,positiveImagEvals]; # save these for VTK visualization;

    real_Evals = mbc['Evals'].real;
    imag_Evals = mbc['Evals'].imag;

    mbc['NaturalFrequencies'] = np.sqrt( real_Evals**2 + imag_Evals**2 );
    mbc['DampRatios'] = -real_Evals/mbc['NaturalFrequencies']
    mbc['DampedFrequencies']  = imag_Evals;

    mbc['NumRigidBodyModes'] = ndof - len(positiveImagEvals);

    mbc['NaturalFreqs_Hz'] = mbc['NaturalFrequencies']/(2.0*np.pi)
    mbc['DampedFreqs_Hz']  = mbc['DampedFrequencies']/(2.0*np.pi);
    mbc['MagnitudeModes']  = np.abs(mbc['EigenVects']);
    mbc['PhaseModes_deg']  = np.angle(mbc['EigenVects'])*180.0/np.pi;
    return mbc, EigenVects_save[:,:,0]

# --------------------------------------------------------------------------------}
# --- Main function 
# --------------------------------------------------------------------------------{
def fx_mbc3(FileNames, verbose=True, starSub=None, removeStatesPattern=None, removeTwrAzimuth=False):
    """ 
    Perform MBC2 funciton based on a list of lin files.
    NOTE: variable names and data structure match MATLAB implementation.

    INPUTS:
     - FileNames: list of lin files for a given operating point

     - starSub: if None, raise an error if `****` are present
                otherwise replace *** with `starSub` (e.g. 0)
                see FASTLinearizationFile. 
     - removeStatesPattern: remove states matching a giving description pattern.
               e.g:  'tower|Drivetrain'  or '^AD'
               see FASTLinearizationFile. 
     - removeTwrAzimuth: if False do nothing
                otherwise discard lin files where azimuth in [60, 180, 300]+/-4deg (close to tower). 

    NOTE: unlike the matlab function, fx_mbc3 does not write the modes for VTK visualization
          Instead use the wrapper function def getCDDOP from pyFAST.linearization.tools

    Original contribution by: Srinivasa B. Ramisett, ramisettisrinivas@yahoo.com, http://ramisetti.github.io
    """

    MBC={}
    matData, _ = get_Mats(FileNames, verbose=verbose, starSub=starSub, removeStatesPattern=removeStatesPattern, removeTwrAzimuth=removeTwrAzimuth)

    # print('matData[Omega] ', matData['Omega'])
    # print('matData[OmegaDot] ', matData['OmegaDot'])
    
    MBC['DescStates'] = matData['DescStates'] # save this in the MBC type for possible campbell_diagram processing later 
    MBC['ndof2'] = matData['ndof2']
    MBC['ndof1'] = matData['ndof1']
    MBC['RotSpeed_rpm'] = np.mean(matData['Omega'])*(30/np.pi); #rad/s to rpm
    MBC['WindSpeed'] = np.mean(matData['WindSpeed']) # NOTE: might be NaN for old files
        
    # print('RotSpeed_rpm ',MBC['RotSpeed_rpm'])
    # print('ndof1 ', MBC['ndof1'])
    # print('ndof2 ', MBC['ndof2'])
    # print(matData['RotTripletIndicesStates2'])
    
    # TODO infer number of blades from ElastoDyn
    # TODO differentiate between number of blades, and the "nb" for MBC
    #  nb = 3; % number of blades required for MBC3 
    # ---------- Multi-Blade-Coordinate transformation -------------------------------------------
    new_seq_dof2, dummy, nb  = get_new_seq(matData['RotTripletIndicesStates2'],matData['ndof2']); # these are the first ndof2 states (not "first time derivative" states); these values are used to calculate matrix transformations
    new_seq_dof1, dummy, nb2 = get_new_seq(matData['RotTripletIndicesStates1'],matData['ndof1']); # these are the first-order ndof1 states; these values are used to calculate matrix transformations

    nb = max(nb,nb2);
    if (nb==0):
        #print('*** fx_mbc3: no states were found, so assuming turbine has 3 blades. ***')
        #nb = 3 # TODO, somehow in the past, we assumed 3 blades.
        print('*** fx_mbc3: no states were found. Setting number of blades to 0. Skipping MBC3 ***')


    new_seq_states=np.concatenate((new_seq_dof2, new_seq_dof2+matData['ndof2']))
    if new_seq_dof1.size!=0:
        new_seq_states=np.concatenate((new_seq_states,new_seq_dof1+matData['NumStates2']))
    
    #new_seq_states = [new_seq_dof2;  new_seq_dof2+matData['ndof2'];  new_seq_dof1+matData['NumStates2']]; # combine the second-order states, including "first time derivatives", with first-order states (assumes ordering of displacements and velocities in state matrices); these values are used to calculate matrix transformations 
    # second-order NonRotating q2, second-order Rotating q2, 
    # second-order NonRotating q2_dot, second-order Rotating q2_dot, 
    # first-order NonRotating q1, first-order Rotating q1


    if nb == 3:
        MBC['performedTransformation'] = True;

        if matData['n_RotTripletStates2'] + matData['n_RotTripletStates1'] < 1:
            print('*** There are no rotating states. MBC transformation, therefore, cannot be performed.');

        # perhaps just warn and perform eigenanalysis anyway?
        if (matData['n_RotTripletStates2']*nb > matData['ndof2']):
            print('**ERROR: the rotating second-order dof exceeds the total num of second-order dof');
        elif (matData['n_RotTripletStates1']*nb > matData['ndof1']):
            print('**ERROR: the rotating first-order dof exceeds the total num of first-order dof');

        new_seq_inp,dummy,dummy = get_new_seq(matData['RotTripletIndicesCntrlInpt'],matData['NumInputs']);
        new_seq_out,dummy,dummy = get_new_seq(matData['RotTripletIndicesOutput'],matData['NumOutputs']);

        n_FixFrameStates2 = matData['ndof2']      - matData['n_RotTripletStates2']*nb;  # fixed-frame second-order dof
        n_FixFrameStates1 = matData['ndof1']      - matData['n_RotTripletStates1']*nb;  # fixed-frame first-order dof
        n_FixFrameInputs  = matData['NumInputs']  - matData['n_RotTripletInputs']*nb;   # fixed-frame control inputs
        n_FixFrameOutputs = matData['NumOutputs'] - matData['n_RotTripletOutputs']*nb;  # fixed-frame outputs

        #print(n_FixFrameOutputs,n_FixFrameInputs, n_FixFrameStates1, n_FixFrameStates2)

        if ( len(matData['Omega']) != matData['NAzimStep']):
            print('**ERROR: the size of Omega vector must equal matData.NAzimStep, the num of azimuth steps')
        if ( len(matData['OmegaDot']) != matData['NAzimStep']):
            print('**ERROR: the size of OmegaDot vector must equal matData.NAzimStep, the num of azimuth steps');


        nLin = matData['A'].shape[-1]
        MBC['A'] = np.zeros(matData['A'].shape)
        MBC['B'] = np.zeros((len(new_seq_states),len(new_seq_inp),matData['NAzimStep']))
        if 'C' in matData.keys():
            MBC['C']=np.zeros(matData['C'].shape)
        else:
            MBC['C']=np.zeros((0,0,nLin))
        if 'D' in matData.keys():
            MBC['D']=np.zeros(matData['D'].shape)
        else:
            MBC['D']=np.zeros((0,0,nLin))
        
        # print('new_seq_inp ',new_seq_inp)
        # print('new_seq_out ',new_seq_out)
        # print('new_seq_states ', new_seq_states)
            
        # begin azimuth loop 
        for iaz in reversed(range(matData['NAzimStep'])):
            #(loop backwards so we don't reallocate memory each time [i.e. variables with iaz index aren't getting larger each time])

            temp=np.arange(nb)
            # compute azimuth positions of blades:
            az = matData['Azimuth'][iaz]*np.pi/180.0 + 2*np.pi/nb* temp ; # Eq. 1, azimuth in radians

            # get rotor speed squared
            OmegaSquared = matData['Omega'][iaz]**2;

            #print(OmegaSquared)

            # compute transformation matrices
            cos_col = np.cos(az);
            sin_col = np.sin(az);

            tt=np.column_stack((np.ones(3),cos_col,sin_col))  # Eq. 9, t_tilde
            ttv = get_tt_inverse(sin_col, cos_col);     # inverse of tt (computed analytically in function below)
            tt2 = np.column_stack((np.zeros(3), -sin_col,  cos_col))     # Eq. 16 a, t_tilde_2
            tt3 = np.column_stack((np.zeros(3), -cos_col, -sin_col))     # Eq. 16 b, t_tilde_3
            
            #---
            T1 = np.eye(n_FixFrameStates2);                # Eq. 11 for second-order states only
            #print('B ',T1, n_FixFrameStates2, matData['n_RotTripletStates2'])
            for ii in range(matData['n_RotTripletStates2']):
                T1 = scp.block_diag(T1,tt)
            
            T1v = np.eye(n_FixFrameStates2);               # inverse of T1
            for ii in  range(matData['n_RotTripletStates2']):
                T1v = scp.block_diag(T1v, ttv);

            T2 = np.zeros([n_FixFrameStates2,n_FixFrameStates2]);              # Eq. 14  for second-order states only
            for ii in range(matData['n_RotTripletStates2']):
                T2 = scp.block_diag(T2, tt2);

            #print('T1, T1v, T2 ',T1.shape, T1v.shape, T2.shape)
            #---    
            T1q = np.eye(n_FixFrameStates1);               # Eq. 11 for first-order states (eq. 8 in MBC3 Update document)
            for ii in range(matData['n_RotTripletStates1']):
                T1q = scp.block_diag(T1q, tt);

            T1qv = np.eye(n_FixFrameStates1);              # inverse of T1q
            for ii in range(matData['n_RotTripletStates1']):
                T1qv = scp.block_diag(T1qv, ttv);

            T2q = np.zeros([n_FixFrameStates1,n_FixFrameStates1]);             # Eq. 14 for first-order states (eq.  9 in MBC3 Update document)
            for ii in range(matData['n_RotTripletStates1']):
                T2q = scp.block_diag(T2q, tt2);

            #print('T1q, T1qv, T2q ',T1q.shape, T1qv.shape, T2q.shape)
            #     T1qc = np.eye(matData.NumHDInputs);            # inverse of T1q

            #---
            T3 = np.zeros([n_FixFrameStates2,n_FixFrameStates2]);              # Eq. 15
            for ii in range(matData['n_RotTripletStates2']):
                T3 = scp.block_diag(T3, tt3);

            #---
            T1c = np.eye(n_FixFrameInputs);                # Eq. 21
            for ii in range(matData['n_RotTripletInputs']):
                T1c = scp.block_diag(T1c, tt)

            T1ov = np.eye(n_FixFrameOutputs);              # inverse of Tlo (Eq. 23)
            for ii in range(matData['n_RotTripletOutputs']):
                T1ov = scp.block_diag(T1ov, ttv);

            #print('T3, T1c, T1ov ',T3.shape, T1c.shape, T1ov.shape, matData['A'].shape)
            # mbc transformation of first-order matrices
            #  if ( MBC.EqnsOrder == 1 ) # activate later

            #print('Before ',T1c)
    
            if 'A' in matData:
                #A =  matData['A'][:,:,iaz]
                #A_reordered =  A[np.ix_(new_seq_states, new_seq_states)]
                # Eq. 29
                L1=np.concatenate((T1, np.zeros([matData['ndof2'],matData['ndof2']]), np.zeros([matData['ndof2'], matData['ndof1']])), axis=1)
                L2=np.concatenate((matData['Omega'][iaz]*T2,T1,np.zeros([matData['ndof2'], matData['ndof1']])), axis=1)
                L3=np.concatenate((np.zeros([matData['ndof1'], matData['ndof2']]),np.zeros([matData['ndof1'], matData['ndof2']]), T1q), axis=1)            
                L=np.matmul(matData['A'][new_seq_states[:,None],new_seq_states,iaz], np.concatenate((L1,L2,L3),axis=0))

                R1=np.concatenate((matData['Omega'][iaz]*T2, np.zeros([matData['ndof2'],matData['ndof2']]), np.zeros([matData['ndof2'], matData['ndof1']])), axis=1)
                R2=np.concatenate((OmegaSquared*T3 + matData['OmegaDot'][iaz]*T2,  2*matData['Omega'][iaz]*T2, np.zeros([matData['ndof2'], matData['ndof1']])),axis=1)
                R3=np.concatenate((np.zeros([matData['ndof1'], matData['ndof2']]), np.zeros([matData['ndof1'], matData['ndof2']]),  matData['Omega'][iaz]*T2q), axis=1)
        
                R=np.concatenate((R1,R2,R3),axis=0)

                MBC['A'][new_seq_states[:,None],new_seq_states,iaz]=np.matmul(scp.block_diag(T1v, T1v, T1qv),(L-R))

                # ffname='AAA'+str(iaz)+'.txt'
                # with open(ffname, "a") as f:
                #     np.savetxt(f,MBC['A'][:,:,iaz],fmt='%5.4f')
                #     f.write('\n')

            if 'B' in matData:
                # Eq. 30
                MBC['B'][new_seq_states[:,None],new_seq_inp,iaz]=np.matmul(np.matmul(scp.block_diag(T1v, T1v, T1qv), matData['B'][new_seq_states[:,None],new_seq_inp,iaz]),T1c)
            
                # ffname='BBB'+str(iaz)+'.txt'
                # with open(ffname, "a") as f:
                #     np.savetxt(f,MBC['B'][:,:,iaz],fmt='%5.4f')
                #     f.write('\n')

            if 'C' in matData:
                # Eq. 31

                L1=np.concatenate((T1, np.zeros([matData['ndof2'],matData['ndof2']]), np.zeros([matData['ndof2'], matData['ndof1']])),axis=1)
                L2=np.concatenate((matData['Omega'][iaz]*T2, T1, np.zeros([matData['ndof2'], matData['ndof1']])), axis=1)
                L3=np.concatenate((np.zeros([matData['ndof1'], matData['ndof2']]), np.zeros([matData['ndof1'], matData['ndof2']]), T1q), axis=1)
            
                MBC['C'][new_seq_out[:,None], new_seq_states,iaz]=np.matmul(np.matmul(T1ov,matData['C'][new_seq_out[:,None],new_seq_states,iaz]),np.concatenate((L1,L2,L3),axis=0))

                # ffname='CCC'+str(iaz)+'.txt'
                # with open(ffname, "a") as f:
                #     np.savetxt(f,MBC['C'][:,:,iaz],fmt='%5.4f')
                #     f.write('\n')

            if 'D' in matData:
               # Eq. 32
                MBC['D'][new_seq_out[:,None],new_seq_inp,iaz] = np.matmul(np.matmul(T1ov,matData['D'][new_seq_out[:,None],new_seq_inp,iaz]), T1c)

                # ffname='DDD'+str(iaz)+'.txt'
                # with open(ffname, "a") as f:
                #     np.savetxt(f,MBC['D'][:,:,iaz],fmt='%5.4f')
                #     f.write('\n')

        # end   # end of azimuth loop
    else:
        print(' fx_mbc3 WARNING: Number of blades is ', str(nb), ' not 3. MBC transformation was not performed.')
        MBC['performedTransformation'] = False;
    
        # initialize matrices
        if 'A' in matData:
            MBC['A'] = matData['A'] # initalize matrix
        if 'B' in matData:
            MBC['B'] = matData['B'] # initalize matrix
        if 'C' in matData:
            MBC['C'] = matData['C'] # initalize matrix
        if 'D' in matData:
            MBC['D'] = matData['D'] # initalize matrix

    # ------------- Eigensolution and Azimuth Averages -------------------------
    if 'A' in MBC:
        MBC['AvgA'] = np.mean(MBC['A'],axis=2); # azimuth-average of azimuth-dependent MBC.A matrices
        MBC['eigSol'], EigenVects_save = eiganalysis(MBC['AvgA'], matData['ndof2'], matData['ndof1']);
        MBC['EigenVects_save'] = EigenVects_save
        MBC['nb'] = nb

    if 'B' in MBC:
        MBC['AvgB'] = np.mean(MBC['B'],axis=2); # azimuth-average of azimuth-dependent MBC.B matrices
        # ffname='BBB_avg'+'.txt'
        # with open(ffname, "a") as f:
        #     np.savetxt(f,MBC['AvgB'],fmt='%5.4f')
        #     f.write('\n')

    if 'C' in MBC:
        MBC['AvgC'] = np.mean(MBC['C'],axis=2); # azimuth-average of azimuth-dependent MBC.C matrices
        # ffname='CCC_avg'+'.txt'
        # with open(ffname, "a") as f:
        #     np.savetxt(f,MBC['AvgC'],fmt='%5.4f')
        #     f.write('\n')

    if 'D' in MBC:
        MBC['AvgD'] = np.mean(MBC['D'],axis=2); # azimuth-average of azimuth-dependent MBC.D matrices
        # ffname='DDD_avg'+'.txt'
        # with open(ffname, "a") as f:
        #     np.savetxt(f,MBC['AvgD'],fmt='%5.4f')
        #     f.write('\n')

    return MBC, matData


#%% ------------------------------------------------------------------------
def formatModesForViz(MBC, matData, nb, EigenVects_save):
    """ 
    get data required for VTK visualization:
        MBC.eigSol.EigenVects_save(:,SortedFreqIndx)       
    """
    nAzimuth = len(matData['Azimuth'])
    nStates, nModes = EigenVects_save.shape

    SortedFreqIndx = np.argsort((MBC['eigSol']['NaturalFreqs_Hz']).flatten(),kind="heapsort")

    #put these in order of natural frequency:
    VTK = dict()
    VTK['NaturalFreq_Hz'] = MBC['eigSol']['NaturalFreqs_Hz'][  SortedFreqIndx] # nModes
    VTK['DampedFreq_Hz']  = MBC['eigSol']['DampedFreqs_Hz'][   SortedFreqIndx] # nModes
    VTK['DampingRatio']   = MBC['eigSol']['DampRatios'][       SortedFreqIndx] # nModes
    x_eig                 =             EigenVects_save[:,     SortedFreqIndx] # nStates x nModes
    # Adopt a convention such that the real part of the first state is positive (arbitrary)
    S=np.sign(np.real(x_eig[0,:]))
    x_eig = S * x_eig
    VTK['x_eig'] = np.tile(x_eig, nAzimuth).reshape((x_eig.shape[0],x_eig.shape[1],nAzimuth))

    if MBC['performedTransformation']:
        # inverse MBC3 (Eq. 4, to move from collective, sine, cosine back to blade 1, blade 2, blade 3):
        dof1_offset = MBC['ndof2']*2
        for iaz, azimuth in enumerate(matData['Azimuth']):
            # MBC3 transformation matrices
            az = azimuth*np.pi/180.0 + 2*np.pi/nb* np.arange(nb)# % Eq. 1, azimuth in radians
            az = az.reshape((-1,1)) # column vector
            tt = np.column_stack( (np.ones((nb,1)), np.cos(az), np.sin(az))) #% Eq. 9, t_tilde
            # MBC on second order states
            I3_2nd = np.array(matData['RotTripletIndicesStates2']).astype(int)
            for i2 in range(I3_2nd.shape[0]):
                i3x   = I3_2nd[i2,:]
                i3xdot= I3_2nd[i2,:]+MBC['ndof2']
                VTK['x_eig'][i3x   , :, iaz] = tt.dot(x_eig[i3x   ,:])
                VTK['x_eig'][i3xdot, :, iaz] = tt.dot(x_eig[i3xdot,:])
            # MBC on first order states
            I3_1st = np.array(matData['RotTripletIndicesStates1']).astype(int)
            for i1 in range(I3_1st.shape[0]):
                  i3x = I3_1st[i1] + dof1_offset
                  VTK['x_eig'][i3x, :, iaz] = tt.dot(x_eig[i3x,:])
    # put this in order states are stored in FAST
    I = matData['StateOrderingIndx'] 
    VTK['x_desc'] = np.array(MBC['DescStates'])[I]
    VTK['x_eig']  = VTK['x_eig'][I,:,:]             # nStates x nModes x nAzimuth
    VTK['x_eig_magnitude'] = np.abs(  VTK['x_eig']) # nStates x nModes x nAzimuth
    VTK['x_eig_phase']     = np.angle(VTK['x_eig']) # nStates x nModes x nAzimuth
    return VTK




if __name__=='__main__':
    pass

    # FileNames=['5MW_Land_ModeShapes-1.fst', '5MW_Land_ModeShapes-2.fst', '5MW_Land_ModeShapes-3.fst', '5MW_Land_ModeShapes-6.fst', '5MW_Land_ModeShapes-7.fst'];
    #FileNames=['5MW_Land_BD_Linear-1.fst', '5MW_Land_BD_Linear-2.fst', '5MW_Land_BD_Linear-3.fst', '5MW_Land_BD_Linear-6.fst', '5MW_Land_BD_Linear-7.fst'];

    #FileNames=['5MW_Land_BD_Linear-1.fst'];

    #FileNames=['DLC-1.1/5MW_Land_BD_Linear-7.1.lin', 'DLC-1.1/5MW_Land_BD_Linear-7.2.lin']
    #FileNames=['/Users/sramiset/Desktop/OpenFAST/5MW_Land_BD_Linear/5MW_Land_BD_Linear-1.1.lin','/Users/sramiset/Desktop/OpenFAST/5MW_Land_BD_Linear/5MW_Land_BD_Linear-1.2.lin']
    # CampbellData=runMBC(FileNames)
    # print('Preparing campbell diagram data!');
    # # TO DO read x-axis for wind speed or rotor speed from csv file
    # #op_csv=pd.read_csv('input.csv', sep=',')
    # OP=[2,4,6,8,10]

    # modeID_table,modesDesc=IdentifyModes(CampbellData)

    # #print(modesDesc)

    # nModes=modeID_table.shape[0]
    # nRuns=modeID_table.shape[1]
    # cols=[item[0] for item in list(modesDesc.values())]
    # #cols.append('1P');cols.append('3P');cols.append('6P')
    # #cols.append('9P');cols.append('12P')
    # frequency=pd.DataFrame(np.nan, index=np.arange(nRuns), columns=cols)
    # dampratio=pd.DataFrame(np.nan, index=np.arange(nRuns), columns=cols)
    # FreqPlotData=np.zeros((nRuns,nModes))
    # DampPlotData=np.zeros((nRuns,nModes))
    # for i in range(nRuns):
    #     for modeID in range(len(modesDesc)): # list of modes we want to identify
    #         idx=int(modeID_table[modeID,i])
    #         FreqPlotData[i,modeID]=CampbellData[i]['Modes'][idx]['NaturalFreq_Hz']
    #         DampPlotData[i,modeID]=CampbellData[i]['Modes'][idx]['DampingRatio']
    #         #print(i,modeID,modesDesc[modeID][0],FreqPlotData[i,modeID])
    #     frequency.iloc[i,:]=FreqPlotData[i,:]
    #     dampratio.iloc[i,:]=DampPlotData[i,:]
        
    # for i in range(len(OP)):
    #     # for 15 DOF
    #     frequency.index.values[i]=OP[i]
    #     dampratio.index.values[i]=OP[i]

    # # import openpyxl
    # # xfile = openpyxl.load_workbook('/Users/sramiset/Desktop/OpenFAST/mbc3_py/CampbellDiagram_Template.xlsx')
        
    # pCD.plotCampbellData(OP,frequency,dampratio)

    # frequency['1P']=np.nan
    # frequency['3P']=np.nan
    # frequency['6P']=np.nan
    # frequency['9P']=np.nan
    # frequency['12P']=np.nan

    # print(nRuns)
    # for i in range(nRuns):
    #     # for 1P,3P,6P,9P,and 12P harmonics
    #     tmp=OP[i]/60.0
    #     print(i,tmp)
    #     LZ=15
    #     frequency.iloc[i,LZ]=tmp
    #     frequency.iloc[i,LZ+1]=3*tmp
    #     frequency.iloc[i,LZ+2]=6*tmp
    #     frequency.iloc[i,LZ+3]=9*tmp
    #     frequency.iloc[i,LZ+4]=12*tmp
    # print(frequency)
    # frequency.transpose().to_excel(r'CampbellData.xlsx')
