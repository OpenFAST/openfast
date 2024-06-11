""" 
Wrapper around wetb to read/write htc files.
TODO: rewrite of c2_def might not be obvious
"""
from .file import File

import numpy as np
import pandas as pd
import os

from .wetb.hawc2.htc_file import HTCFile
from .hawc2_st_file import HAWC2StFile

class HAWC2HTCFile(File):

    @staticmethod
    def defaultExtensions():
        return ['.htc']

    @staticmethod
    def formatName():
        return 'HAWC2 htc file'

    def _read(self):
        self.data = HTCFile(self.filename)

    def _write(self):
        self.data.save(self.filename)

    def __repr__(self):
        s='<{} object>\n'.format(type(self).__name__)
        s+='| Attributes:\n'
        s+='| - data: HTCFile with keys: {}`\n'.format(list(self.data.keys()))
        s+='| Derived attributes:\n'
        s+='| * bodyNames: {}\n'.format(self.bodyNames)
        s+='| * bodyDict: dict with keys bodyNames\n'
        s+='| Methods: bodyByName, bodyC2, setBodyC2\n'
        return s

    @property
    def bodyNames(self):
        return list(self.bodyDict.keys())

    @property
    def bodyDict(self):
        struct = self.data.new_htc_structure
        bodyKeys = [k for k in struct.keys() if k.startswith('main_body')]
        #print('>>> keys', struct.keys())
        bdDict={}
        for k in bodyKeys:
            bodyName = struct[k].name[0]
            bdDict[bodyName] = struct[k]
        return bdDict

    def bodyByName(self, bodyname):
        """ return body inputs given a body name"""
        bodyDict= self.bodyDict
        if bodyname not in bodyDict.keys():
            raise Exception('No body found with name {} in file {}'.format(bodyname,self.filename))
        return bodyDict[bodyname]

    def bodyC2(self, bdy):
        """ return body C2_def given body inputs"""
        try:
            nsec = bdy.c2_def.nsec[0]
        except:
            raise Exception('body has no c2_def section')
        val = np.array([bdy.c2_def[k].values[0:] for k in bdy.c2_def.keys() if k.startswith('sec')])
        val = val.reshape((-1,5)).astype(float)
        val = val[np.argsort(val[:,0]),:]
        return val

    def setBodyC2(self, bdy, val):
        """ set body C2_def given body inputs and new c2def"""
        # TODO different number of section
        nsec     = bdy.c2_def.nsec[0]
        nsec_new = val.shape[0]
        sec_keys = [k for k in bdy.c2_def.keys() if k.startswith('sec')]
        bdy.c2_def.nsec = nsec_new
        if nsec != nsec_new:
            if nsec_new<nsec:
                for k in sec_keys[nsec_new:]:
                    del bdy.c2_def.k
                    del bdy.c2_def.contents[k]
                sec_keys = sec_keys[:nsec_new]
                # we delete
            else:
                raise NotImplementedError('Setting c2_def with different number of sections')
        for i, k in enumerate(sec_keys):
            bdy.c2_def[k].values[0] = int(val[i][0])
            bdy.c2_def[k].values[1:] = val[i][1:]
        pass


    def _toDataFrame(self):
        dfs ={}
        j=0

        simdir  = os.path.dirname(self.filename)

        # --- C2 def
        bodyKeys = [k for k in self.data.new_htc_structure.keys() if k.startswith('main_body')]
        for k in bodyKeys:
            bdy = self.data.new_htc_structure[k]
            try:
                val = self.bodyC2(bdy)
            except:
                continue
            name = bdy.name[0]
            val = val[:,[3,1,2,4]]
            dfs[name+'_c2'] = pd.DataFrame(data=val, columns=['z_[m]', 'x_[m]', 'y_[m]','twist_[deg]'])
            # potentially open st files..
            if "timoschenko_input" in bdy:
                tim = bdy.timoschenko_input
                H2_stfile = os.path.join(simdir, tim.filename[0])
                if not os.path.exists(H2_stfile):
                    # Try with a parent directory..
                    H2_stfile = os.path.join(simdir, '../',tim.filename[0])

                if not os.path.exists(H2_stfile):
                    print('[WARN] st file referenced in htc file was not found for body {}.\nSt file: {}\nhtc file {}'.format(name, H2_stfile, self.filename))
                else:
                    dfs_st = HAWC2StFile(H2_stfile).toDataFrame(extraCols=False)
                    if 'set' in tim.keys():
                        mset = tim.set[0]
                        iset = tim.set[1]
                        sSet = '{}_{}'.format(mset,iset)
                    else:
                        sSet = '1_1'
                    if sSet not in dfs_st:
                        raise Exception('Set {} not found in file {}'.format(sSet, H2_stfile))
                    else:
                        dfs[name+'_st'] = dfs_st[sSet]
        # --- Potentially ae and pc

        return dfs
