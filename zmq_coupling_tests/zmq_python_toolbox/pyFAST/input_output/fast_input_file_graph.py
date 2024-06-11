import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Local 
# try:
from pyFAST.input_output.tools.graph import *
# except ImportError:
#from welib.FEM.graph import *


# --------------------------------------------------------------------------------}
# --- Wrapper to convert a "fast" input file dictionary into a graph
# --------------------------------------------------------------------------------{
def fastToGraph(data, **kwargs):
    if 'BeamProp' in data.keys():
        return subdynToGraph(data, **kwargs)
    
    if 'SmplProp' in data.keys():
        return hydrodynToGraph(data, **kwargs)

    if 'DOF2Nodes' in data.keys():
        return subdynSumToGraph(data, **kwargs)

    raise NotImplementedError('Graph for object with keys: {}'.format(data.keys()))

# --------------------------------------------------------------------------------}
# --- SubDyn
# --------------------------------------------------------------------------------{
def subdynToGraph(sd, propToNodes=False, propToElem=False):
    """
    sd: dict-like object as returned by weio

    -propToNodes: if True, the element properties are also transferred to the nodes for convenience.
                 NOTE: this is not the default because a same node can have two different diameters in SubDyn (it's by element)
    """
    type2Color=[
            (0.1,0.1,0.1), # Watchout based on background
            (0.753,0.561,0.05),  # 1 Beam
            (0.541,0.753,0.05),  # 2 Cable
            (0.753,0.05,0.204),  # 3 Rigid
            (0.918,0.702,0.125), # 3 Rigid
        ]

    Graph = GraphModel()
    # --- Properties
    if 'BeamProp' in sd.keys():
        BProps     = sd['BeamProp']
        Graph.addNodePropertySet('Beam')
        for ip,P in enumerate(BProps):
            prop= NodeProperty(ID=P[0], E=P[1], G=P[2], rho=P[3], D=P[4], t=P[5] )
            Graph.addNodeProperty('Beam',prop)

    if 'CableProp' in sd.keys():
        CProps     = sd['CableProp']
        Graph.addNodePropertySet('Cable')
        for ip,P in enumerate(CProps):
            Chan = -1 if len(P)<5 else P[4]
            prop= NodeProperty(ID=P[0], EA=P[1], rho=P[2], T0=P[3], Chan=Chan)
            Graph.addNodeProperty('Cable',prop)

    if 'RigidProp' in sd.keys():
        RProps     = sd['RigidProp']
        Graph.addNodePropertySet('Rigid')
        for ip,P in enumerate(RProps):
            prop= NodeProperty(ID=P[0], rho=P[1])
            Graph.addNodeProperty('Rigid',prop)

    # --- Nodes and DOFs
    Nodes     = sd['Joints']
    for iNode,N in enumerate(Nodes):
        Type= 1 if len(N)<=4 else N[4]
        node = Node(ID=N[0], x=N[1], y=N[2], z=N[3], Type=Type)
        Graph.addNode(node)

    # --- Elements
    Members  = sd['Members'].astype(int)
    PropSets = ['Beam','Cable','Rigid']
    for ie,E in enumerate(Members):
        Type=1 if len(E)==5 else E[5]
        #elem= Element(E[0], E[1:3], propset=PropSets[Type-1], propIDs=E[3:5])
        elem= Element(E[0], E[1:3], Type=PropSets[Type-1], propIDs=E[3:5], propset=PropSets[Type-1])
        elem.data['object']='cylinder'
        elem.data['color'] = type2Color[Type]
        Graph.addElement(elem)
        # Nodal prop data
        if propToNodes:
            # NOTE: this is disallowed by default because a same node can have two different diameters in SubDyn (it's by element)
            Graph.setElementNodalProp(elem, propset=PropSets[Type-1], propIDs=E[3:5])
        if propToElem:
            Graph.setElementNodalPropToElem(elem) # TODO, this shouldn't be needed

    # --- Concentrated Masses (in global coordinates), node data
    for iC, CM in enumerate(sd['ConcentratedMasses']):
        #CMJointID, JMass, JMXX, JMYY, JMZZ, JMXY, JMXZ, JMYZ, MCGX, MCGY, MCGZ   
        nodeID = CM[0]
        n = Graph.getNode(nodeID)
        M66 = np.zeros((6,6))
        if len(CM)==11:
            m = CM[1]
            x, y ,z = (CM[8], CM[9], CM[10])
            Jxx = CM[2]; Jyy = CM[3]; Jzz = CM[4]
            Jxy = CM[5]; Jxz = CM[6]; Jyz = CM[7];
        else:
            raise NotImplementedError('TODO legacy')
            m = CM[1]
            Jxx = CM[2]; Jyy = CM[3]; Jzz = CM[4]
            Jxy = 0; Jyz =0; Jzz = 0; x,y,z=0,0,0
        M66[0, :] =[   m     ,   0     ,   0     ,   0                 ,  z*m                , -y*m                 ]
        M66[1, :] =[   0     ,   m     ,   0     , -z*m                ,   0                 ,  x*m                 ]
        M66[2, :] =[   0     ,   0     ,   m     ,  y*m                , -x*m                ,   0                  ]
        M66[3, :] =[   0     , -z*m    ,  y*m    , Jxx + m*(y**2+z**2) , Jxy - m*x*y         , Jxz  - m*x*z         ]
        M66[4, :] =[  z*m    ,   0     , -x*m    , Jxy - m*x*y         , Jyy + m*(x**2+z**2) , Jyz  - m*y*z         ]
        M66[5, :] =[ -y*m    , x*m     ,   0     , Jxz - m*x*z         , Jyz - m*y*z         , Jzz  + m*(x**2+y**2) ]
        n.setData({'addedMassMatrix':M66})

    # Nodal data
    for iN,N in enumerate(sd['InterfaceJoints']):
        nodeID   = int(N[0])
        Graph.setNodalData(nodeID,IBC=N[1:])
    for iN,N in enumerate(sd['BaseJoints']):
        NN=[int(n) if i<7 else n for i,n in enumerate(N)]
        nodeID   = NN[0]
        Graph.setNodalData(nodeID,RBC=NN[1:])
    #     print('CMass')
    #     print(sd['ConcentratedMasses'])

    return Graph





# --------------------------------------------------------------------------------}
# --- HydroDyn
# --------------------------------------------------------------------------------{
def hydrodynToGraph(hd, propToNodes=False, propToElem=False):
    """
     hd: dict-like object as returned by weio

    -propToNodes: if True, the element properties are also transferred to the nodes for convenience.
                 NOTE: this is not the default because a same node can have two different diameters in SubDyn (it's by element)

    - propToElem: This might be due to misunderstanding of graph..
    """
    def type2Color(Pot):
        if Pot:
            return (0.753,0.05,0.204),  # Pot flow
        else:
            return (0.753,0.561,0.05),  # Morison


    Graph = GraphModel()

    # --- Properties
    if 'SectionProp' in hd.keys():
        # NOTE: setting it as element property since two memebrs may connect on the same node with different diameters/thicknesses
        Graph.addNodePropertySet('Section')
        for ip,P in enumerate(hd['SectionProp']):
            # PropSetID    PropD         PropThck
            prop= NodeProperty(ID=P[0], D=P[1], t=P[2])
            Graph.addNodeProperty('Section',prop)

    # --- Hydro Coefs - will be stored in AxCoefs, SimpleCoefs, DepthCoefs, MemberCoefs
    if 'AxCoefs' in hd.keys():
        Graph.addNodePropertySet('AxCoefs')
        for ip,P in enumerate(hd['AxCoefs']):
            prop= NodeProperty(ID=P[0], JAxCd=P[1], JAxCa=P[2], JAxCp=P[3])
            Graph.addNodeProperty('AxCoefs',prop)
    if 'SmplProp' in hd.keys():
        Graph.addNodePropertySet('SimpleCoefs')
        for ip,P in enumerate(hd['SmplProp']):
            #      SimplCd    SimplCdMG    SimplCa    SimplCaMG    SimplCp    SimplCpMG   SimplAxCd  SimplAxCdMG   SimplAxCa  SimplAxCaMG  SimplAxCp   SimplAxCpMG
            if len(P)==12:
                prop= NodeProperty(ID=ip+1, Cd=P[0], CdMG=P[1], Ca=P[2], CaMG=P[3], Cp=P[4], CpMG=P[5], AxCd=P[6], AxCdMG=P[7], AxCa=P[8], AxCaMG=P[9], AxCp=P[10], AxCpMG=P[11])
            elif len(P)==10:
                prop= NodeProperty(ID=ip+1, Cd=P[0], CdMG=P[1], Ca=P[2], CaMG=P[3], Cp=P[4], CpMG=P[5], AxCa=P[6], AxCaMG=P[7], AxCp=P[8], AxCpMG=P[9])
            else:
                raise NotImplementedError()
            Graph.addNodeProperty('SimpleCoefs',prop)
    if 'DpthProp' in hd.keys():
        Graph.addMiscPropertySet('DepthCoefs')
        for ip,P in enumerate(hd['DpthProp']):
            # Dpth      DpthCd   DpthCdMG   DpthCa   DpthCaMG       DpthCp   DpthCpMG   DpthAxCd   DpthAxCdMG   DpthAxCa   DpthAxCaMG   DpthAxCp   DpthAxCpMG
            prop= Property(ID=ip+1, Dpth=P[0], Cd=P[1], CdMG=P[2], Ca=P[3], CaMG=P[4], Cp=P[5], CpMG=P[6], AxCd=P[7], AxCdMG=P[8], AxCa=P[9], AxCaMG=P[10], AxCp=P[11], AxCpMG=P[12])
            Graph.addMiscProperty('DepthCoefs',prop)

    if 'MemberProp' in hd.keys():
        # Member-based hydro coefficinet
        Graph.addMiscPropertySet('MemberCoefs')
        for ip,P in enumerate(hd['MemberProp']):
            # MemberID    MemberCd1     MemberCd2    MemberCdMG1   MemberCdMG2    MemberCa1     MemberCa2    MemberCaMG1   MemberCaMG2    MemberCp1     MemberCp2    MemberCpMG1   MemberCpMG2   MemberAxCd1   MemberAxCd2  MemberAxCdMG1 MemberAxCdMG2  MemberAxCa1   MemberAxCa2  MemberAxCaMG1 MemberAxCaMG2  MemberAxCp1  MemberAxCp2   MemberAxCpMG1   MemberAxCpMG2
            prop = ElemProperty(ID=ip+1, MemberID=P[0], Cd1=P[1], Cd2=P[2], CdMG1=P[3], CdMG2=P[4], Ca1=P[5], Ca2=P[6], CaMG1=P[7], CaMG2=P[8], Cp1=P[9], Cp2=P[10], CpMG1=P[11], CpMG2=P[12], AxCd1=P[14], AxCd2=P[15], axCdMG1=P[16], axCdMG2=P[17], AxCa1=P[18], AxCa2=P[19], AxCaMG1=P[20], AxCaMG2=P[21], AxCp1=P[22], AxCp2=P[23])
            Graph.addMiscProperty('MemberCoefs',prop)
    # ---
    if 'FillGroups' in hd.keys():
        # Filled members
        Graph.addMiscPropertySet('FillGroups')
        print('>>> TODO Filled Groups')
        #for ip,P in enumerate(hd['FillGroups']):
        #    #                       FillNumM FillMList             FillFSLoc     FillDens
        #    raise NotImplementedError('hydroDynToGraph, Fill List might not be properly set, verify below')
        #    prop = MiscProperty(ID=ip+1, FillNumM=P[0], FillMList=P[1],  FillFSLoc=P[2], FillDens=P[3])
        #    Graph.addMiscProperty('FillGroups',prop)

    if 'MGProp' in hd.keys():
        # Marine Growth
        Graph.addMiscPropertySet('MG')
        for ip,P in enumerate(hd['MGProp']):
            # MGDpth     MGThck       MGDens
            # (m)        (m)         (kg/m^3)
            prop = Property(ID=ip+1, MGDpth=P[0], MGThck=P[1],  MGDens=P[2])
            Graph.addMiscProperty('FillGroups',prop)

    # --- Nodes
    Nodes     = hd['Joints']
    for iNode,N in enumerate(Nodes):
        node = Node(ID=N[0], x=N[1], y=N[2], z=N[3])
        Graph.addNode(node)
        Graph.setNodeNodalProp(node, 'AxCoefs', N[4])
   
    # --- Elements
    PropSets=['SimpleCoefs','DepthCoefs','MemberCoefs']
    Members   = hd['Members']
    for ie,E in enumerate(Members):
        # MemberID  MJointID1  MJointID2  MPropSetID1  MPropSetID2  MDivSize   MCoefMod  PropPot 
        EE   = E[:5].astype(int)
        Type = int(E[6]) # MCoefMod
        Pot  = E[7].lower()[0]=='t'
        elem= Element(ID=EE[0], nodeIDs=EE[1:3], propIDs=EE[3:5], propset='Section', CoefMod=PropSets[Type-1], DivSize=float(E[5]), Pot=Pot)
        elem.data['object']='cylinder'
        elem.data['color'] = type2Color(Pot)
        Graph.addElement(elem)
        # Nodal prop data NOTE: can't do that anymore for memebrs with different diameters at the same node
        if propToNodes:
            # NOTE: not by default because of feature with members with different diameters at the same node
            Graph.setElementNodalProp(elem, propset='Section', propIDs=EE[3:5])
        if propToElem:
            Graph.setElementNodalPropToElem(elem) # TODO, this shouldn't be needed

        if Type==1:
            # Simple
            Graph.setElementNodalProp(elem, propset='SimpleCoefs', propIDs=[1,1])
        else:
            print('>>> TODO type DepthCoefs and MemberCoefs')
            # NOTE: this is disallowed by default because a same node can have two different diameters in SubDyn (it's by element)
            #Graph.setElementNodalProp(elem, propset=PropSets[Type-1], propIDs=E[3:5])

    return Graph


# --------------------------------------------------------------------------------}
# --- SubDyn Summary file 
# --------------------------------------------------------------------------------{
def subdynSumToGraph(data, Graph=None):
    """ 
     data: dict-like object as returned by weio
    """
    type2Color=[
            (0.1,0.1,0.1), # Watchout based on background
            (0.753,0.561,0.05),  # 1 Beam
            (0.541,0.753,0.05),  # 2 Cable
            (0.753,0.05,0.204),  # 3 Rigid
            (0.918,0.702,0.125), # 3 Rigid
        ]

    #print(data.keys())
    DOF2Nodes = data['DOF2Nodes']
    nDOF      = data['nDOF_red']

    if Graph is None:
        Graph = GraphModel()

    # --- Nodes and DOFs
    Nodes = data['Nodes']
    for iNode,N in enumerate(Nodes):
        if len(N)==9: # Temporary fix
            #N[4]=np.float(N[4].split()[0])
            N=N.astype(np.float32)
        ID = int(N[0])
        nodeDOFs=DOF2Nodes[(DOF2Nodes[:,1]==ID),0] # NOTE: these were reindex to start at 0
        node = Node(ID=ID, x=N[1], y=N[2], z=N[3], Type=int(N[4]), DOFs=nodeDOFs)
        Graph.addNode(node)

    # --- Elements
    Elements = data['Elements']
    for ie,E in enumerate(Elements):
        nodeIDs=[int(E[1]),int(E[2])]
        #  shear_[-]       Ixx_[m^4]       Iyy_[m^4]       Jzz_[m^4]          T0_[N]
        D = np.sqrt(E[7]/np.pi)*4 # <<< Approximation basedon area TODO use I as well
        elem= Element(int(E[0]), nodeIDs, Type=int(E[5]), Area=E[7], rho=E[8], E=E[7], G=E[8], D=D)
        elem.data['object']='cylinder'
        elem.data['color'] = type2Color[int(E[5])]
        Graph.addElement(elem)

    #print(self.extent)
    #print(self.maxDimension)

    # --- Graph Modes
    # Very important sortDims should be None to respect order of nodes
    dispGy, posGy, InodesGy, dispCB, posCB, InodesCB = data.getModes(sortDim=None) 
    for iMode in range(dispGy.shape[2]):
        Graph.addMode(displ=dispGy[:,:,iMode],name='GY{:d}'.format(iMode+1), freq=1/(2*np.pi))

    for iMode in range(dispCB.shape[2]):
        Graph.addMode(displ=dispCB[:,:,iMode],name='CB{:d}'.format(iMode+1), freq=data['CB_frequencies'][iMode]) 

    #print(Graph.toJSON())


    return Graph



if __name__ == '__main__':
    from .fast_input_file import FASTInputFile

    filename='../../_data/Monopile/MT100_SD.dat'
    # filename='../../_data/Monopile/TetraSpar_SubDyn_v3.dat'

    sd = FASTInputFile(filename)
#     sd.write('OutMT.dat')
    Graph = sd.toGraph()
    Graph.divideElements(2)
    print(Graph)
    print(Graph.sortNodesBy('z'))
    # print(Graph.nodalDataFrame(sortBy='z'))
    print(Graph.points)
    print(Graph.connectivity)
    print(Graph)

# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib import collections  as mc
# from mpl_toolkits.mplot3d import Axes3D
# fig = plt.figure()
# ax = fig.add_subplot(1,2,1,projection='3d')
# 
# lines=Graph.toLines(output='coord')
# for l in lines:
# #     ax.add_line(l)
#     ax.plot(l[:,0],l[:,1],l[:,2])
# 
# ax.autoscale()
# ax.set_xlim([-40,40])
# ax.set_ylim([-40,40])
# ax.set_zlim([-40,40])
# # ax.margins(0.1)
# 
# plt.show()


