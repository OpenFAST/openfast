""" 
Basics Classes for a "geometrical" graph model: 
  - nodes have a position (x,y,z), and some data (taken from a list of properties)
  - elements (links) connect nodes, they contain some data (taken from a list of properties)

An ordering of Elements, Nodes, and Properties is present, but whenever possible,
the "ID" is used to identify them, instead of their index.


Nodes: 
   Node.ID:    unique ID (int) of the node. IDs never change.
   Node.x,y,z: coordinate of the nodes
   Node.data : dictionary of data stored at the node

Elements: 
   Elem.ID:      unique ID (int) of the element. IDs never change.
   Elem.nodeIDs: list of node IDs making up the element
   Elem.nodes  : list of nodes making the element (by reference) # NOTE: this has cross reference!
   Elem.nodeProps   : properties # Nodal properties. NOTE: cannot be transfered to node because of how SubDyn handles it..
   Elem.data   : dictionary of data stored at the element
   # Optional
   Elem.propset: string referring to the property set in the dictionary of properties
   Elem.propIDs: IDs used for the properties of this element at each node

NodePropertySets: dictionary of NodeProperties
   Node Property: 
      NProp.ID:  unique ID of the node proprety
      NProp.data: dictionary of data
          

ElemPropertySets: dictionary of ElemProperties

"""

import numpy as np
import pandas as pd


# --------------------------------------------------------------------------------}
# --- Node
# --------------------------------------------------------------------------------{
class Node(object):
    def __init__(self, ID, x, y, z=0, **kwargs):
        self.ID = int(ID)
        self.x  = x
        self.y  = y
        self.z  = z
        self.data  = kwargs

    def setData(self, data_dict):
        """ set or add data"""
        for k,v in data_dict.items():
            #if k in self.data.keys():
            #    print('Warning overriding key {} for node {}'.format(k,self.ID))
            self.data[k]=v

    def __repr__(self):
        s='<Node{:4d}> x:{:7.2f} y:{:7.2f} z:{:7.2f} {:}'.format(self.ID, self.x, self.y, self.z, self.data)
        return s

# --------------------------------------------------------------------------------}
# --- Properties  
# --------------------------------------------------------------------------------{
class Property(dict):
    def __init__(self, ID, data=None, **kwargs):
        """ 
        data is a dictionary
        """
        dict.__init__(self)
        self.ID= int(ID)
        self.update(kwargs)
        if data is not None:
            self.update(data)

    @property
    def data(self):
        return {k:v for k,v in self.items() if k!='ID'}

    def __repr__(self):
        s='<Prop{:4d}> {:}'.format(self.ID, self.data)
        return s

class NodeProperty(Property):
    def __init__(self, ID, data=None, **kwargs):
        Property.__init__(self, ID, data, **kwargs)
    def __repr__(self):
        s='<NPrp{:4d}> {:}'.format(self.ID, self.data)
        return s
    
class ElemProperty(Property):
    def __init__(self, ID, data=None, **kwargs):
        Property.__init__(self, ID, data, **kwargs)
    def __repr__(self):
        s='<EPrp{:4d}> {:}'.format(self.ID, self.data)
        return s


# --------------------------------------------------------------------------------}
# --- Elements 
# --------------------------------------------------------------------------------{
class Element(dict):
    def __init__(self, ID, nodeIDs, nodes=None, propset=None, propIDs=None, properties=None, **kwargs):
        """ 

        """
        self.ID      = int(ID)
        self.nodeIDs = nodeIDs
        self.propset = propset
        self.propIDs = propIDs
        self.data    = kwargs     # Nodal data
        self.nodes   = nodes      # Typically a trigger based on nodeIDs
        self.nodeProps= properties # Typically a trigger based on propIDs. Otherwise list of dictionaries
        if (self.propIDs is not None) and (self.propset is None):
            raise Exception('`propset` should be provided if `propIDs` are provided')
        if (self.propIDs is not None) and (self.propset is not None) and properties is not None:
            raise Exception('When providing `propset` & `propIDs`, properties should not be provided')
        if nodes is not None:
            if len(nodes)!=len(nodeIDs):
                raise Exception('List of nodes has different length than list of nodeIDs')
            for i, (ID,n) in enumerate(zip(nodeIDs,nodes)):
                if n.ID!=ID:
                    raise Exception('Node ID do not match {}/={} for node index {}'.format(n.ID,ID,i))
        
    @property
    def length(self):
        n1=self.nodes[0]
        n2=self.nodes[1]
        return np.sqrt((n1.x-n2.x)**2+(n1.y-n2.y)**2+(n1.z-n2.z)**2)

    def setData(self, data_dict):
        """ set or add data"""
        for k,v in data_dict.items():
            #if k in self.data.keys():
            #    print('Warning overriding key {} for node {}'.format(k,self.ID))
            self.data[k]=v

    def __repr__(self):
        s='<Elem{:4d}> NodeIDs: {} {}'.format(self.ID, self.nodeIDs, self.data)
        if self.propIDs is not None:
            s+=' {'+'propIDs:{} propset:{}'.format(self.propIDs, self.propset)+'}'
        if self.nodes is not None:
            s+=' l={:.2f}'.format(self.length)
        return s


# --------------------------------------------------------------------------------}
# --- Mode 
# --------------------------------------------------------------------------------{
class Mode(dict):
    def __init__(self, data, name, freq=1, **kwargs):
        dict.__init__(self)

        self['name']=name
        self['freq']=freq
        self['data']=data # displacements nNodes x 3 assuming a given sorting of nodes

    def __repr__(self):
        s='<Mode> name:{:4s} freq:{:} '.format(self['name'], self['freq'])
        return s

    def reSort(self,I):
        self['data']=self['data'][I,:]

# --------------------------------------------------------------------------------}
# --- Graph
# --------------------------------------------------------------------------------{
class GraphModel(object):
    def __init__(self, Elements=None, Nodes=None, NodePropertySets=None, ElemPropertySets=None, MiscPropertySets=None ):
        self.Elements         = Elements if Elements is not None else []
        self.Nodes            = Nodes    if Nodes    is not None else []
        self.NodePropertySets = NodePropertySets if NodePropertySets is not None else {}
        self.ElemPropertySets = ElemPropertySets if ElemPropertySets is not None else {}
        self.MiscPropertySets = MiscPropertySets if MiscPropertySets is not None else {}
        # Dynamics
        self.Modes   = []
        self.Motions = []
        # Optimization variables
        self._nodeIDs2Elements   = {} # dictionary with key NodeID and value list of ElementID
        self._nodeIDs2Elements   = {} # dictionary with key NodeID and value list of elements
        self._elementIDs2NodeIDs = {} # dictionary with key ElemID and value list of nodes IDs
        self._connectivity =[]# 

    def addNode(self,node):
        self.Nodes.append(node)

    def addElement(self,elem):
        # Giving nodes to element if these were not provided
        elem.nodes=[self.getNode(i) for i in elem.nodeIDs]
        # Giving props to element if these were not provided
        if elem.propIDs is not None:
            elem.nodeProps=[self.getNodeProperty(elem.propset, i) for i in elem.propIDs]
        self.Elements.append(elem)

    # --- Getters
    def getNode(self, nodeID):
        for n in self.Nodes:
            if n.ID==nodeID:
                return n
        raise KeyError('NodeID {} not found in Nodes'.format(nodeID))

    def getElement(self, elemID):
        for e in self.Elements:
            if e.ID==elemID:
                return e
        raise KeyError('ElemID {} not found in Elements'.format(elemID))

    def getNodeProperty(self, setname, propID):
        for p in self.NodePropertySets[setname]:
            if p.ID==propID:
                return p
        raise KeyError('PropID {} not found for Node propset {}'.format(propID,setname))

    def getElementProperty(self, setname, propID):
        for p in self.ElemPropertySets[setname]:
            if p.ID==propID:
                return p
        raise KeyError('PropID {} not found for Element propset {}'.format(propID,setname))

    def getMiscProperty(self, setname, propID):
        for p in self.MiscPropertySets[setname]:
            if p.ID==propID:
                return p
        raise KeyError('PropID {} not found for Misc propset {}'.format(propID,setname))

    # ---
    @property
    def nodeIDs2ElementIDs(self):
        """ Return list of elements IDs connected to each node"""
        if len(self._nodeIDs2ElementIDs) == 0:
            # Compute list of connected elements for each node
            self._nodeIDs2ElementIDs=dict()
            for i,n in enumerate(self.Nodes):
                self._nodeIDs2ElementIDs[n.ID] = [e.ID for e in self.Elements if n.ID in e.nodeIDs]
        return self._nodeIDs2ElementIDs

    @property
    def nodeIDs2Elements(self):
        """ Return list of elements connected to each node"""
        if len(self._nodeIDs2Elements) == 0:
            # Compute list of connected elements for each node
            self._nodeIDs2Elements
            for i,n in enumerate(self.Nodes):
                self._nodeIDs2Elements[n.ID] = [e for e in self.Elements if n.ID in e.nodeIDs]
        return self._nodeIDs2Elements


    @property
    def elementIDs2NodeIDs(self):
        """ returns """
        if len(self._elementIDs2NodeIDs) ==0:
            self._elementIDs2NodeIDs =dict()
            for e in self.Elements:
                self._elementIDs2NodeIDs[e.ID] = [n.ID for n in e.nodes] 
        return self._elementIDs2NodeIDs


    @property
    def connectivity(self):
        """ returns connectivity, assuming points are indexed starting at 0 
        NOTE: this is basically element2Nodes but reindexed
        """
        if len(self._connectivity) ==0:
            self._connectivity = [[self.Nodes.index(n)  for n in e.nodes] for e in self.Elements]
        return self._connectivity


    # --- Handling of (element/material) Properties
    def addElementPropertySet(self, setname):
        self.ElemPropertySets[setname]= []

    def addNodePropertySet(self, setname):
        self.NodePropertySets[setname]= []

    def addMiscPropertySet(self, setname):
        self.MiscPropertySets[setname]= []

    def addNodeProperty(self, setname, prop):
        if not isinstance(prop, NodeProperty):
            print(type(prop))
            raise Exception('Property needs to inherit from NodeProperty')
        self.PropertySets[setname].append(prop)

    def addNodeProperty(self, setname, prop):
        if not isinstance(prop, NodeProperty):
            print(type(prop))
            raise Exception('Property needs to inherit from NodeProperty')
        self.NodePropertySets[setname].append(prop)

    def addElementProperty(self, setname, prop):
        if not isinstance(prop, ElemProperty):
            print(type(prop))
            raise Exception('Property needs to inherit from ElementProperty')
        self.ElemPropertySets[setname].append(prop)

    def addMiscProperty(self, setname, prop):
        if not isinstance(prop, ElemProperty):
            print(type(prop))
            raise Exception('Property needs to inherit from Property')
        self.MiscPropertySets[setname].append(prop)

    # --- Data and node and element prop setters
    def setElementNodalProp(self, elem, propset, propIDs):
        """ 
        Set Nodal Properties to each node of an element
        """
        for node, pID in zip(elem.nodes, propIDs):
            node.setData(self.getNodeProperty(propset, pID).data)

    def setElementNodalPropToElem(self, elem):
        """ 
        Set Element Properties to an element
        TODO: this seems to be a hack. It should be automatic I think...
        """
        propset=elem.propset
        propIDs=elem.propIDs
        # USING PROPID 0!!!
        elem.setData(self.getNodeProperty(propset, propIDs[0]).data)
        # TODO average the two maybe..

    def setNodeNodalProp(self, node, propset, propID):
        """ 
        Set Nodal Properties to a node
        """
        node.setData(self.getNodeProperty(propset, propID).data)

    def setNodalData(self, nodeID, **data_dict):
        self.getNode(nodeID).setData(data_dict)

    def __repr__(self):
        s='<{} object> with keys:\n'.format(type(self).__name__)
        s+='- Nodes ({}):\n'.format(len(self.Nodes))
        s+='\n'.join(str(n) for n in self.Nodes)
        s+='\n- Elements ({}):\n'.format(len(self.Elements))
        s+='\n'.join(str(n) for n in self.Elements)
        s+='\n- NodePropertySets ({}):'.format(len(self.NodePropertySets))
        for k,v in self.NodePropertySets.items():
            s+='\n> {} ({}):\n'.format(k, len(v))
            s+='\n'.join(str(p) for p in v)
        s+='\n- ElementPropertySets ({}):'.format(len(self.ElemPropertySets))
        for k,v in self.ElemPropertySets.items():
            s+='\n> {} ({}):\n'.format(k, len(v))
            s+='\n'.join(str(p) for p in v)
        s+='\n- MiscPropertySets ({}):'.format(len(self.MiscPropertySets))
        for k,v in self.MiscPropertySets.items():
            s+='\n> {} ({}):\n'.format(k, len(v))
            s+='\n'.join(str(p) for p in v)
        s+='\n- Modes ({}):\n'.format(len(self.Modes))
        s+='\n'.join(str(m) for m in self.Modes)
        s+='\n- Motions ({}):'.format(len(self.Motions))
        for m in self.Motions:
            s+='\n> {}\n'.format({k:v for k,v in m.items() if not isintance(v,np.ndarray)})
        return s

    # --------------------------------------------------------------------------------}
    # --- Geometrical properties 
    # --------------------------------------------------------------------------------{
    @property
    def extent(self):
        xmax=np.max([node.x for node in self.Nodes])
        ymax=np.max([node.y for node in self.Nodes])
        zmax=np.max([node.z for node in self.Nodes])
        xmin=np.min([node.x for node in self.Nodes])
        ymin=np.min([node.y for node in self.Nodes])
        zmin=np.min([node.z for node in self.Nodes])
        return [xmin,ymin,zmin],[xmax,ymax,zmax],[xmax-xmin,ymax-ymin,zmax-zmin]

    @property
    def maxDimension(self):
        _,_,D=self.extent
        return np.max(D)

    @property
    def points(self):
        nNodes = len(self.Nodes)
        Points = np.zeros((nNodes,3))
        for i,n in enumerate(self.Nodes):
            Points[i,:]=(n.x, n.y, n.z)
        return Points

    def toLines(self, output='coord'):
        if output=='coord':
            lines = np.zeros((len(self.Elements), 2, 3)) # 
            for ie, e in enumerate(self.Elements):
                n1=e.nodes[0]
                n2=e.nodes[-1]
                lines[ie, 0, : ] = (n1.x, n1.y, n1.z)
                lines[ie, 1, : ] = (n2.x, n2.y, n2.z)
        elif output=='lines3d':
            import mpl_toolkits.mplot3d as plt3d
            lines=[]
            for ie, e in enumerate(self.Elements):
                n1=e.nodes[0]
                n2=e.nodes[-1]
                line = plt3d.art3d.Line3D((n1.x,n2.x), (n1.y,n2.y), (n1.z,n2.z))
                lines.append(line)
        else:
            raise NotImplementedError()

        return lines

    # --------------------------------------------------------------------------------}
    # --- Change of connectivity
    # --------------------------------------------------------------------------------{
    def connecticityHasChanged(self):
        self._nodeIDs2ElementIDs = dict()
        self._nodeIDs2Elements   = dict()
        self._elementIDs2NodeIDs = dict()
        self._connectivity=[]

    def updateConnectivity(self):
        for e in self.Elements:
            e.nodes=[self.getNode(i) for i in e.nodeIDs]

        for e in self.Elements:
            e.nodeProps = [self.getNodeProperty(e.propset, ID) for ID in e.propIDs]

        # Potentially call nodeIDs2ElementIDs etc


    def _divideElement(self, elemID, nPerElement, maxElemId, keysNotToCopy=[]):
        """ divide a given element by nPerElement (add nodes and elements to graph) """ 
        if len(self.Modes)>0:
            raise Exception('Cannot divide graph when mode data is present')
        if len(self.Motions)>0:
            raise Exception('Cannot divide graph when motion data is present')

        maxNodeId=np.max([n.ID for n in self.Nodes])
        e = self.getElement(elemID)
        newElems = []
        if len(e.nodes)==2:
            n1=e.nodes[0]
            n2=e.nodes[1]
            subNodes=[n1]
            for iSub in range(1,nPerElement):
                maxNodeId += 1
                #data_dict  = n1.data.copy()
                data_dict  = dict()
                fact       = float(iSub)/nPerElement
                # Interpolating position
                x          = n1.x*(1-fact)+n2.x*fact
                y          = n1.y*(1-fact)+n2.y*fact
                z          = n1.z*(1-fact)+n2.z*fact
                # Interpolating data (only if floats)
                for k,v in n1.data.items():
                    if k not in keysNotToCopy:
                        try:
                            data_dict[k] = n1.data[k]*(1-fact) + n2.data[k]*fact
                        except:
                            data_dict[k] = n1.data[k]
                ni = Node(maxNodeId, x, y, z, **data_dict)
                subNodes.append(ni)
                self.addNode(ni)
            subNodes+=[n2]
            e.nodes  =subNodes[0:2]
            e.nodeIDs=[e.ID for e in e.nodes]
            for i in range(1,nPerElement):
                maxElemId+=1
                elem_dict = e.data.copy()
                # Creating extra properties if necessary
                if e.propIDs is not None:
                    if all(e.propIDs==e.propIDs[0]):
                        # No need to create a new property
                        propIDs=e.propIDs
                        propset=e.propset
                    else:
                        raise NotImplementedError('Division of element with different properties on both ends. TODO add new property.')
                elem= Element(maxElemId, [subNodes[i].ID, subNodes[i+1].ID], propset=propset, propIDs=propIDs, **elem_dict )
                newElems.append(elem)
        return newElems


    def sortNodesBy(self,key):
        """ Sort nodes, will affect the connectivity, but node IDs remain the same"""

        # TODO, that's quite doable
        if len(self.Modes)>0:
            raise Exception('Cannot sort nodes when mode data is present')
        if len(self.Motions)>0:
            raise Exception('Cannot sort nodes when motion data is present')

        nNodes = len(self.Nodes)
        if key=='x':
            values=[n.x for n in self.Nodes]
        elif key=='y':
            values=[n.y for n in self.Nodes]
        elif key=='z':
            values=[n.z for n in self.Nodes]
        elif key=='ID':
            values=[n.ID for n in self.Nodes]
        else:
            values=[n[key] for n in self.Nodes]
        I= np.argsort(values)
        self.Nodes=[self.Nodes[i] for i in I]

        # Trigger, remove precomputed values related to connectivity:
        self.connecticityHasChanged()

        return self

    def divideElements(self, nPerElement, excludeDataKey='', excludeDataList=[], method='append', keysNotToCopy=[]):
        """ divide all elements by nPerElement (add nodes and elements to graph)

        - excludeDataKey: is provided, will exclude elements such that e.data[key] in `excludeDataList`

        - method: append or insert

        - keysNotToCopy: when duplicating node and element data, make sure not to duplicate data with these keys
                         For instance if a node that has a boundary condition, it should not be passed to the 
                         node that is created when dividing an element.

        Example: 
           to avoid dividing elements of `Type` 'Cable' or `Rigid`, call as follows:
             self.divideElements(3, excludeDataKey='Type', excludeDataList=['Cable','Rigid'] )

        """ 
        maxNodeId=np.max([n.ID for n in self.Nodes])
        maxElemId=np.max([e.ID for e in self.Elements])

        if nPerElement<=0:
            raise Exception('nPerElement should be more than 0')

        newElements=[]
        for ie in np.arange(len(self.Elements)): # cannot enumerate since length increases
            elemID = self.Elements[ie].ID
            if method=='insert':
                newElements+=[self.getElement(elemID)] # newElements contains
            if (len(excludeDataKey)>0 and self.Elements[ie].data[excludeDataKey] not in excludeDataList) or len(excludeDataKey)==0:
                elems = self._divideElement(elemID, nPerElement, maxElemId, keysNotToCopy)
                maxElemId+=len(elems)
                newElements+=elems
            else:
                print('Not dividing element with ID {}, based on key `{}` with value `{}`'.format(elemID, excludeDataKey,self.Elements[ie].data[excludeDataKey]))
        # Adding elements at the end
        if method=='append':
            pass
        elif method=='insert':
            self.Elements=[] # We clear all elements
        else:
            raise NotImplementedError('Element Insertions')

        for e in newElements:
            self.addElement(e)

        # Trigger, remove precomputed values related to connectivity:
        self.connecticityHasChanged()

        return self
                    
    # --------------------------------------------------------------------------------}
    # --- Dynamics
    # --------------------------------------------------------------------------------{
    def addMode(self,displ,name=None,freq=1):
        if name is None:
            name='Mode '+str(len(self.Modes))
        mode = Mode(data=displ, name=name, freq=freq)
        self.Modes.append(mode)


    # --------------------------------------------------------------------------------}
    # --- Ouputs / converters
    # --------------------------------------------------------------------------------{
    def nodalDataFrame(self, sortBy=None):
        """ return a DataFrame of all the nodal data """
        data=dict()
        nNodes=len(self.Nodes)
        for i,n in enumerate(self.Nodes):
            if i==0:
                data['ID'] = np.zeros(nNodes).astype(int)
                data['x']  = np.zeros(nNodes)
                data['y']  = np.zeros(nNodes)
                data['z']  = np.zeros(nNodes)

            data['ID'][i] = n.ID
            data['x'][i]  = n.x
            data['y'][i]  = n.y
            data['z'][i]  = n.z
            for k,v in n.data.items():
                if k not in data:
                    data[k] = np.zeros(nNodes)
                try:
                    data[k][i]=v
                except:
                    pass
        df = pd.DataFrame(data)
        # Sorting 
        if sortBy is not None:
            df.sort_values([sortBy],inplace=True,ascending=True)
            df.reset_index(drop=True,inplace=True) 
        return df


    def toJSON(self,outfile=None):
        d=dict();
        Points=self.points
        d['Connectivity'] = self.connectivity
        d['Nodes']        = Points.tolist()
        
        d['ElemProps']=list()
        for iElem,elem in enumerate(self.Elements):
            Shape = elem.data['shape'] if 'shape' in elem.data.keys() else 'cylinder'
            Type  = elem.data['Type'] if 'Type' in elem.data.keys() else 1
            try:
                Diam  = elem.D
            except:
                Diam  = elem.data['D'] if 'D' in elem.data.keys() else 1
            if Shape=='cylinder':
                d['ElemProps'].append({'shape':'cylinder','type':Type, 'Diam':Diam})
            else:
                raise NotImplementedError()


        d['Modes']=[
                {
                    'name': self.Modes[iMode]['name'],
                    'omega':self.Modes[iMode]['freq']*2*np.pi, #in [rad/s]
                    'Displ':self.Modes[iMode]['data'].tolist()
                }  for iMode,mode in enumerate(self.Modes)]
        d['groundLevel']=np.min(Points[:,2]) # TODO

        if outfile is not None:
            import json
            jsonFile=outfile
            with open(jsonFile, 'w', encoding='utf-8') as f:
                #f.write(to_json(d))
                try:
                    #f.write(unicode(json.dumps(d, ensure_ascii=False))) #, indent=2)
                    #f.write(json.dumps(d, ensure_ascii=False)) #, indent=2)
                    f.write(json.dumps(d, ensure_ascii=False))
                except:
                    print('>>> FAILED')
                    json.dump(d, f, indent=0) 
        return d

# 



INDENT = 3
SPACE = " "
NEWLINE = "\n"
# Changed basestring to str, and dict uses items() instead of iteritems().

def to_json(o, level=0):
  ret = ""
  if isinstance(o, dict):
    if level==0:
        ret += "{" + NEWLINE
        comma = ""
        for k, v in o.items():
          ret += comma
          comma = ",\n"
          ret += SPACE * INDENT * (level + 1)
          ret += '"' + str(k) + '":' + SPACE
          ret += to_json(v, level + 1)
        ret += NEWLINE + SPACE * INDENT * level + "}"
    else:
        ret += "{" 
        comma = ""
        for k, v in o.items():
          ret += comma
          comma = ",\n"
          ret += SPACE
          ret += '"' + str(k) + '":' + SPACE
          ret += to_json(v, level + 1)
        ret += "}"

  elif isinstance(o, str):
    ret += '"' + o + '"'
  elif isinstance(o, list):
    ret += "[" + ",".join([to_json(e, level + 1) for e in o]) + "]"
  # Tuples are interpreted as lists
  elif isinstance(o, tuple):
    ret += "[" + ",".join(to_json(e, level + 1) for e in o) + "]"
  elif isinstance(o, bool):
    ret += "true" if o else "false"
  elif isinstance(o, int):
    ret += str(o)
  elif isinstance(o, float):
    ret += '%.7g' % o
  elif isinstance(o, numpy.ndarray) and numpy.issubdtype(o.dtype, numpy.integer):
    ret += "[" + ','.join(map(str, o.flatten().tolist())) + "]"
  elif isinstance(o, numpy.ndarray) and numpy.issubdtype(o.dtype, numpy.inexact):
    ret += "[" + ','.join(map(lambda x: '%.7g' % x, o.flatten().tolist())) + "]"
  elif o is None:
    ret += 'null'
  else:
    raise TypeError("Unknown type '%s' for json serialization" % str(type(o)))
  return ret
