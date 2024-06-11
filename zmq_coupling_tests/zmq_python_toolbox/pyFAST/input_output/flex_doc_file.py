import numpy as np
import pandas as pd
import os
import re
try:
    from .file import File, WrongFormatError, BrokenFormatError
except:
    File = dict
    class WrongFormatError(Exception): pass
    class BrokenFormatError(Exception): pass

class FLEXDocFile(File):

    @staticmethod
    def defaultExtensions():
        return ['.out','doc']

    @staticmethod
    def formatName():
        return 'FLEX WaveKin file'

    def _read(self):
        with open(self.filename, 'r', errors="surrogateescape") as f:
            line1=f.readline().strip()
            if line1.find('#Program')!=0:
                raise WrongFormatError()
            lines=f.read().splitlines()+[line1]

        ArraysHeader=[
        '#Blade_ShapeFunction_DOF1_Shape',
        '#Blade_ShapeFunction_DOF2_Shape',
        '#Blade_ShapeFunction_DOF3_Shape',
        '#Blade_ShapeFunction_DOF4_Shape',
        '#Blade_ShapeFunction_DOF5_Shape',
        '#Blade_ShapeFunction_DOF6_Shape',
        '#Blade_ShapeFunction_DOF7_Shape',
        '#Blade_ShapeFunction_DOF8_Shape',
        '#Blade_ShapeFunction_DOF9_Shape',
        '#Blade_ShapeFunction_DOF10_Shape',
        '#Tower_SectionData',
        '#Tower_ShapeFunction_DOF1_Shape',
        '#Tower_ShapeFunction_DOF2_Shape',
        '#Tower_ShapeFunction_DOF3_Shape',
        '#Tower_ShapeFunction_DOF4_Shape',
        '#Tower_ShapeFunction_DOF5_Shape',
        '#Tower_ShapeFunction_DOF6_Shape',
        '#Tower_ShapeFunction_DOF7_Shape',
        '#Tower_ShapeFunction_DOF8_Shape',
        '#Foundation_SectionData',
        '#Foundation_ShapeFunction_DOF1_Shape',
        '#Foundation_ShapeFunction_DOF2_Shape',
        '#Global_Mode1_Shape',
        ]

        ArraysNoHeader=[
        '#Blade_IsolatedMassMatrix'
        '#Blade_IsolatedStiffnessMatrix',
        '#Blade_IsolatedDampingMatrix',
        '#Foundation_IsolatedMassMatrix',
        '#Foundation_IsolatedStiffnessMatrix',
        '#Foundation_IsolatedStiffnessMatrixCorrection',
        '#Foundation_IsolatedDampingMatrix',
        '#Tower_MassMatrix',
        '#Tower_IsolatedMassMatrix',
        '#Tower_IsolatedStiffnessMatrix',
        '#Tower_IsolatedStiffnessMatrixCorrection',
        '#Tower_IsolatedDampingMatrix',
        '#EVA_MassMatrix',
        '#EVA_StiffnessMatrix',
        '#EVA_DampingMatrix',
        '#EVA_Eigenvectors',
        '#EVA_Eigenfrequencies',
        '#EVA_Eigenvalues',
        '#EVA_Damping',
        '#EVA_LogDec',
        ]

        numeric_const_pattern = r'[-+]? (?: (?: \d* \. \d+ ) | (?: \d+ \.? ) )(?: [Ee] [+-]? \d+ ) ?'
        rx = re.compile(numeric_const_pattern, re.VERBOSE)

        i=0
        while i<len(lines):
            l=lines[i]
            if len(l.strip())==0:
                i+=1; continue
            if l.find('###')==0:
                i+=1; continue
            if l[0]=='#':
                sp=l.split()
                # --- Array with header
                if sp[0] in ArraysHeader:
                    array_lines=[]
                    i=i+1
                    header=lines[i]
                    i=i+1
                    while i<len(lines) and len(lines[i])>0 and lines[i][0]!='#':
                        array_lines.append(np.array([float(v) if v not in ['Fnd','Twr','COG'] else ['Fnd','Twr','COG'].index(v) for v in lines[i].split()]))
                        i=i+1
                    # --- Process array
                    M = np.array(array_lines)
                    # --- Process header
                    cols=header.split()
                    try:
                        ii=int(cols[0])
                        header=' '.join(cols[1:])
                    except:
                        pass
                    if header.find('[')<=0:
                        cols=header.split()
                    else:
                        header=header.replace('rough','rough_[-]')
                        header=header.replace('n ','n_[-] ')
                        spcol=header.split(']')
                        cols= [v.strip().replace(' ','_').replace('[','_[').replace('__','_').replace('__','_')+']' for v in spcol[:-1]]

                    if len(cols)!=M.shape[1]:
                        cols=['C{}'.format(j) for j in range(M.shape[1])]
                    # --- Store
                    keys = sp[0].split('_')
                    keys[0]=keys[0][1:]
                    if keys[0] not in self.keys():
                        self[keys[0]] =  dict()
                    subkey = '_'.join(keys[1:])
                    df = pd.DataFrame(data = M, columns = cols)
                    self[keys[0]][subkey] = df
                    continue
                # --- Array with no header
                elif sp[0] in ArraysNoHeader:
                    array_lines=[]
                    i=i+1
                    header=lines[i]
                    i=i+1
                    while i<len(lines) and len(lines[i])>0 and lines[i][0]!='#':
                        array_lines.append(np.array([float(v) for v in lines[i].split()]))
                        i=i+1
                    # --- Process array
                    M = np.array(array_lines)
                    # --- Store
                    keys = sp[0].split('_')
                    keys[0]=keys[0][1:]
                    if keys[0] not in self.keys():
                        self[keys[0]] =  dict()
                    subkey = '_'.join(keys[1:])
                    self[keys[0]][subkey] = M
                    continue
                else:
                    # --- Regular
                    keys = sp[0].split('_')
                    key=keys[0][1:]
                    subkey = '_'.join(keys[1:])
                    values= ' '.join(sp[1:])
                    try:
                        dat= np.array(rx.findall(values)).astype(float)
                        if len(dat)==1:
                            dat=dat[0]
                    except:
                        dat = values

                    if len(key.strip())>0:
                        if len(subkey)==0:
                            if key not in self.keys():
                                self[key] = dat
                            else:
                                print('>>> line',i,l)
                                print(self.keys())
                                raise Exception('Duplicate singleton key:',key)
                        else:
                            if key not in self.keys():
                                self[key] = dict()
                            self[key][subkey]=dat
            i+=1

            # --- Some adjustements
            try:
                df=self['Global']['Mode1_Shape']
                df=df.drop('C0',axis=1)
                df.columns=['H_[m]','U_FA_[-]','U_SS_[]']
                self['Global']['Mode1_Shape']=df
            except:
                pass

#     def _write(self):
#         with open(self.filename,'w') as f:
#             f.write(self.toString)

    def __repr__(self):
        s='<{} object> with keys:\n'.format(type(self).__name__)
        for k in self.keys():
            if type(self[k]) is dict:
                s+='{:15s}: dict with keys {}\n'.format(k, list(self[k].keys()))
            else:
                s+='{:15s} : {}\n'.format(k,self[k])
        return s

    def _toDataFrame(self):
        dfs={}
        for k,v in self.items():
            if type(v) is pd.DataFrame:
                dfs[k]=v
            #if type(v) is np.ndarray:
            #    if len(v.shape)>1:
            #        dfs[k]=v
            if type(self[k]) is dict:
                for k2,v2 in self[k].items():
                    if type(v2) is pd.DataFrame:
                        dfs[k+'_'+k2]=v2
        return dfs

