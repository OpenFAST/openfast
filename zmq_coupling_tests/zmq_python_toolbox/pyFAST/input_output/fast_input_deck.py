import os
import numpy as np
import re
import pandas as pd

from .fast_input_file import FASTInputFile

__all__  = ['FASTInputDeck']
# --------------------------------------------------------------------------------}
# --- Full FAST input deck
# --------------------------------------------------------------------------------{
class FASTInputDeck(dict):
    """Container for input files that make up a FAST input deck"""

    def __init__(self, fullFstPath='', readlist=['all'], verbose=False):
        """Read FAST master file and read inputs for FAST modules

        INPUTS:
          - fullFstPath: 
          - readlist: list of module files to be read, or ['all'], modules are identified as follows:
                ['Fst','ED','AD','BD','BDbld','EDtwr','EDbld','ADbld','AF','AC','IW','HD','SrvD','SD','MD']
                where: 
                 AF: airfoil polars
                 AC: airfoil coordinates (if present)

        """
        # Sanity
        if type(verbose) is not bool: 
            raise Exception('`verbose` arguments needs to be a boolean')

        # Main Data
        self.inputFilesRead = {}
        self.filename = fullFstPath
        self.verbose  = verbose
        self.readlist = readlist
        if not type(self.readlist) is list:
            self.readlist=[readlist]
        if 'all' in self.readlist:
            self.readlist = ['Fst','ED','AD','BD','BDbld','EDtwr','EDbld','ADbld','AF','AC','IW','HD','SrvD','SD','MD']
        else:
            self.readlist = ['Fst']+self.readlist


        # --- Harmonization with AeroElasticSE
        self.FAST_ver       = 'OPENFAST'
        self.path2dll       = None   # Path to dll file

        self.fst_vt={}
        self.fst_vt['description']       = ''            
        self.fst_vt['Fst']               = None
        self.fst_vt['ElastoDyn']         = None
        self.fst_vt['ElastoDynBlade']    = None
        self.fst_vt['ElastoDynTower']    = None
        self.fst_vt['InflowWind']        = None
        self.fst_vt['AeroDyn14']         = None
        self.fst_vt['AeroDyn15']         = None
        self.fst_vt['AeroDynBlade']      = None 
        self.fst_vt['AeroDynTower']      = None
        self.fst_vt['AeroDynPolar']      = None
        self.fst_vt['ServoDyn']          = None
        self.fst_vt['DISCON_in']         = None
        self.fst_vt['HydroDyn']          = None
        self.fst_vt['MoorDyn']           = None
        self.fst_vt['SubDyn']            = None
        self.fst_vt['MAP']               = None
        self.fst_vt['BeamDyn']           = None
        self.fst_vt['BeamDynBlade']      = None # Small change of interface
        self.fst_vt['af_data']           = [] # Small change of interface
        self.fst_vt['ac_data']           = [] # TODO, how is it stored in WEIS?


        self.ADversion=''

        # Read all inputs files
        if len(fullFstPath)>0:
            self.read()

    @property
    def ED(self):
        ED = self.fst_vt['ElastoDyn']
        if ED is None:
            if 'ED' not in self.readlist:
                self.readlist.append('ED')
            if self.verbose:
                print('>>> Reading ED', self.ED_path)
            self.fst_vt['ElastoDyn'] = self._read(self.fst_vt['Fst']['EDFile'],'ED')
            return self.fst_vt['ElastoDyn']
        else:
            return ED


    def readAD(self, filename=None, readlist=None, verbose=False, key='AeroDyn15'):
        """ 
        readlist: 'AD','AF','AC'
        """
        if readlist is not None:
            readlist_bkp = self.readlist
            self.readlist=readlist
            if not type(self.readlist) is list:
                self.readlist=[readlist]
            if 'all' in self.readlist:
                self.readlist = ['Fst','ED','AD','BD','BDbld','EDtwr','EDbld','ADbld','AF','AC','IW','HD','SrvD','SD','MD']

        if filename is None:
            filename = self.fst_vt['Fst']['AeroFile']
            baseDir  = os.path.dirname(self.fst_vt['Fst']['AeroFile'])
        else:
            baseDir  = os.path.dirname(filename)

        self.verbose  = verbose

        self.fst_vt[key] = self._read(filename,'AD')

        if self.fst_vt[key] is not None:
            # Blades
            bld_file = os.path.join(baseDir, self.fst_vt[key]['ADBlFile(1)'])
            self.fst_vt['AeroDynBlade'] = self._read(bld_file,'ADbld')
            #self.fst_vt['AeroDynBlade'] = []
            #for i in range(3):
            #    bld_file = os.path.join(os.path.dirname(self.fst_vt['Fst']['AeroFile']), self.fst_vt[key]['ADBlFile({})'.format(i+1)])
            #    self.fst_vt['AeroDynBlade'].append(self._read(bld_file,'ADbld'))
            # Polars
            self.fst_vt['af_data']=[] # TODO add to "AeroDyn"
            for afi, af_filename in enumerate(self.fst_vt['AeroDyn15']['AFNames']):
                af_filename = os.path.join(baseDir,af_filename).replace('"','')
                try: 
                    polar = self._read(af_filename, 'AF')
                except:
                    polar=None
                    print('[FAIL] reading polar {}'.format(af_filename))
                self.fst_vt['af_data'].append(polar)
                if polar is not None:
                    coordFile = polar['NumCoords']
                    if isinstance(coordFile,str):
                        coordFile = coordFile.replace('"','')
                        baseDirCoord=os.path.dirname(af_filename)
                        if coordFile[0]=='@':
                            ac_filename = os.path.join(baseDirCoord,coordFile[1:])
                            coords = self._read(ac_filename, 'AC')
                            self.fst_vt['ac_data'].append(coords)

        # --- Backward compatibility
        self.AD  = self.fst_vt[key]
        self.ADversion='AD15' if key=='AeroDyn15' else 'AD14'

        if readlist is not None:
            self.readlist=readlist_bkp

    @property
    def FAST_InputFile(self):
        return os.path.basename(self.filename)   # FAST input file (ext=.fst)
    @property
    def FAST_directory(self):
        return os.path.dirname(self.filename)    # Path to fst directory files


    @property
    def inputFiles(self):
        files=[]
        files+=[self.ED_path, self.ED_twr_path, self.ED_bld_path]
        files+=[self.BD_path, self.BD_bld_path]
        files+=[self.SD_path]
        return [f for f in files if f not in self.unusedNames]

    def _relpath(self, k1, k2=None, k3=None):
        try:
            if k2 is None:
                return self.fst_vt['Fst'][k1].replace('"','')
            else:
                parent = os.path.dirname(self.fst_vt['Fst'][k1]).replace('"','')
                if type(k3)==list:
                    for k in k3:
                        if k in self.fst_vt[k2].keys():
                            child =  self.fst_vt[k2][k].replace('"','')
                else:
                    child =  self.fst_vt[k2][k3].replace('"','')
                return os.path.join(parent, child)
        except:
            return 'none'

    @property
    def ED_path(self): return self._fullpath(self._relpath('EDFile'))
    @property
    def SD_path(self): return self._fullpath(self._relpath('SubFile'))
    @property
    def BD_path(self): return self._fullpath(self._relpath('BDBldFile(1)'))
    @property
    def BD_bld_path(self): return self._fullpath(self._relpath('BDBldFile(1)','BeamDyn','BldFile'))
    @property
    def ED_twr_path(self): return self._fullpath(self._relpath('EDFile','ElastoDyn','TwrFile'))
    @property
    def ED_bld_path(self): return self._fullpath(self._relpath('EDFile','ElastoDyn',['BldFile(1)','BldFile1']))



    def _fullpath(self, relfilepath):
        relfilepath = relfilepath.replace('"','')
        basename = os.path.basename(relfilepath)
        if basename.lower() in self.unusedNames:
            return 'none'
        else:
            return os.path.join(self.FAST_directory, relfilepath)


    def read(self, filename=None):
        """ 
        Read all OpenFAST inputs files, based on the requested list of modules `readlist`
        """
        if filename is not None:
            self.filename = filename

        # Read main file (.fst, or .drv) and store into key "Fst"
        if self.verbose:
            print('Reading:', self.FAST_InputFile)
        self.fst_vt['Fst'] = self._read(self.FAST_InputFile, 'Fst')
        if self.fst_vt['Fst'] is None:
            raise Exception('Error reading main file {}'.format(self.filename))
        keys = self.fst_vt['Fst'].keys()

        # Detect driver or OpenFAST version
        if 'NumTurbines' in keys:
            self.version='AD_driver'
        elif 'DynamicSolve' in keys:
            self.version='BD_driver'
        elif 'InterpOrder' in self.fst_vt['Fst'].keys():
            self.version='OF2'
        else:
            self.version='F7'


        if self.version=='AD_driver':
            # ---- AD Driver
            # InflowWind
            if self.fst_vt['Fst']['CompInflow']>0:
                self.fst_vt['InflowWind'] = self._read(self.fst_vt['Fst']['InflowFile'],'IW')

            self.readAD(key='AeroDyn15')

        elif self.version=='BD_driver':
            # --- BD driver
            self.fst_vt['BeamDyn'] = self._read(self.fst_vt['Fst']['InputFile'],'BD')
            if self.fst_vt['BeamDyn'] is not None:
                # Blades
                bld_file = os.path.join(os.path.dirname(self.fst_vt['Fst']['InputFile']), self.fst_vt['BeamDyn']['BldFile'])
                print('bld_file', bld_file)
                self.fst_vt['BeamDynBlade']= self._read(bld_file,'BDbld')

            del self.fst_vt['af_data']
            del self.fst_vt['ac_data']

        elif self.version=='OF2':
            # ---- Regular OpenFAST file
            # ElastoDyn
            if 'EDFile' in self.fst_vt['Fst'].keys():
                self.fst_vt['ElastoDyn'] = self._read(self.fst_vt['Fst']['EDFile'],'ED')
                if self.fst_vt['ElastoDyn'] is not None:
                    twr_file = self.ED_twr_path
                    bld_file = self.ED_bld_path
                    self.fst_vt['ElastoDynTower'] = self._read(twr_file,'EDtwr')
                    self.fst_vt['ElastoDynBlade'] = self._read(bld_file,'EDbld')

            # InflowWind
            if self.fst_vt['Fst']['CompInflow']>0:
                self.fst_vt['InflowWind'] = self._read(self.fst_vt['Fst']['InflowFile'],'IW')

            # AeroDyn
            if self.fst_vt['Fst']['CompAero']>0:
                key = 'AeroDyn14' if self.fst_vt['Fst']['CompAero']==1 else 'AeroDyn15'
                self.readAD(key=key, readlist=self.readlist)

            # ServoDyn
            if self.fst_vt['Fst']['CompServo']>0:
                self.fst_vt['ServoDyn'] = self._read(self.fst_vt['Fst']['ServoFile'],'SrvD')
                # TODO Discon

            # HydroDyn
            if self.fst_vt['Fst']['CompHydro']== 1:
                self.fst_vt['HydroDyn'] = self._read(self.fst_vt['Fst']['HydroFile'],'HD')

            # SubDyn
            if self.fst_vt['Fst']['CompSub'] == 1:
                self.fst_vt['SubDyn'] = self._read(self.fst_vt['Fst']['SubFile'], 'SD')

            # Mooring
            if self.fst_vt['Fst']['CompMooring']==1:
                self.fst_vt['MAP'] = self._read(self.fst_vt['Fst']['MooringFile'],'MD')
            if self.fst_vt['Fst']['CompMooring']==2:
                self.fst_vt['MoorDyn'] = self._read(self.fst_vt['Fst']['MooringFile'],'MD')

            # BeamDyn
            if self.fst_vt['Fst']['CompElast'] == 2:
                self.fst_vt['BeamDyn'] = self._read(self.fst_vt['Fst']['BDBldFile(1)'],'BD')
                if self.fst_vt['BeamDyn'] is not None:
                    # Blades
                    bld_file = os.path.join(os.path.dirname(self.fst_vt['Fst']['BDBldFile(1)']), self.fst_vt['BeamDyn']['BldFile'])
                    self.fst_vt['BeamDynBlade']= self._read(bld_file,'BDbld')

        # --- Backward compatibility
        self.fst = self.fst_vt['Fst']
        self._ED  = self.fst_vt['ElastoDyn']
        if not hasattr(self,'AD'):
            self.AD = None
        if self.AD is not None:
            self.AD.Bld1 = self.fst_vt['AeroDynBlade']
            self.AD.AF  = self.fst_vt['af_data']
        self.IW  = self.fst_vt['InflowWind']
        self.BD  = self.fst_vt['BeamDyn']
        self.BDbld  = self.fst_vt['BeamDynBlade']
        self.SD  = self.fst_vt['SubDyn']

    @ property
    def unusedNames(self):
        return ['unused','nan','na','none']

    def _read(self, relfilepath, shortkey):
        """ read any openfast input """
        relfilepath =relfilepath.replace('"','')
        basename = os.path.basename(relfilepath)

        # Only read what the user requested to be read
        if shortkey not in self.readlist:
            if self.verbose:
                print('>>> Skipping ',shortkey)
            return None

        # Skip "unused" and "NA"
        if basename.lower() in self.unusedNames:
            if self.verbose:
                print('>>> Unused ',shortkey)
            return None

        # Attempt reading
        if relfilepath.startswith(self.FAST_directory):
            fullpath = relfilepath
        else:
            fullpath = os.path.join(self.FAST_directory, relfilepath)
        try:
            data = FASTInputFile(fullpath)
            if self.verbose:
                print('>>> Read: ',fullpath)
            self.inputFilesRead[shortkey] = fullpath
            return data
        except FileNotFoundError:
            print('[WARN] File not found '+fullpath)
            return None



    def write(self, filename=None, prefix='', suffix='', directory=None):
        """ Write a standardized input file deck"""
        if filename is None:
            filename=self.filename # Overwritting
        self.filename=filename
        if directory is None:
            directory = os.path.dirname(filename)
        else:
            # Making sure filename is within directory
            filename = os.path.join(directory, os.path.basename(filename))
        if not os.path.exists(directory):
            os.makedirs(directory)
  
        basename = os.path.splitext(os.path.basename(filename))[0]


        fst = self.fst_vt['Fst']


        if self.version=='AD_driver':
            raise NotImplementedError()

        elif self.version=='BD_driver':
            # --- BD driver
            filename_BD     = os.path.join(directory, prefix+'BD'+suffix+'.dat')
            filename_BD_bld = os.path.join(directory, prefix+'BD_bld'+suffix+'.dat')
            fst['InputFile'] = '"' + os.path.basename(filename_BD) + '"'
            fst.write(filename)
            BD = self.fst_vt['BeamDyn'] 
            BD['BldFile'] = '"'+os.path.basename(filename_BD_bld)+'"'
            self.fst_vt['BeamDynBlade'].write(filename_BD_bld)  # TODO TODO pick up the proper blade file!
            BD.write(filename_BD)

        elif self.version=='OF2':

            # Filenames
            filename_ED     = os.path.join(directory,prefix+'ED'+suffix+'.dat')      if fst['CompElast']>0   else 'none'
            filename_IW     = os.path.join(directory,prefix+'IW'+suffix+'.dat')      if fst['CompInflow']>0  else 'none'
            filename_BD     = os.path.join(directory,prefix+'BD'+suffix+'.dat')      if fst['CompElast']==2  else 'none'
            filename_AD     = os.path.join(directory,prefix+'AD'+suffix+'.dat')      if fst['CompAero']>0    else 'none'
            filename_HD     = os.path.join(directory,prefix+'HD'+suffix+'.dat')      if fst['CompHydro']>0   else 'none'
            filename_SD     = os.path.join(directory,prefix+'SD'+suffix+'.dat')      if fst['CompSub']>0     else 'none'
            filename_MD     = os.path.join(directory,prefix+'MD'+suffix+'.dat')      if fst['CompMooring']>0 else 'none'
            filename_SvD    = os.path.join(directory,prefix+'SvD'+suffix+'.dat')     if fst['CompServo']>0   else 'none'
            filename_Ice    = os.path.join(directory,prefix+'Ice'+suffix+'.dat')     if fst['CompIce']>0     else 'none'
            filename_ED_bld = os.path.join(directory,prefix+'ED_bld'+suffix+'.dat')  if fst['CompElast']>0   else 'none'
            filename_ED_twr = os.path.join(directory,prefix+'ED_twr'+suffix+'.dat')  if fst['CompElast']>0   else 'none'
            filename_BD_bld = os.path.join(directory,prefix+'BD_bld'+suffix+'.dat')  if fst['CompElast']>0   else 'none'
            # TODO AD Profiles and OLAF

            fst['EDFile']       = '"' + os.path.basename(filename_ED) + '"'
            fst['BDBldFile(1)'] = '"' + os.path.basename(filename_BD) + '"'
            fst['BDBldFile(2)'] = '"' + os.path.basename(filename_BD) + '"'
            fst['BDBldFile(3)'] = '"' + os.path.basename(filename_BD) + '"'
            fst['InflowFile']   = '"' + os.path.basename(filename_IW) + '"'
            fst['AeroFile']     = '"' + os.path.basename(filename_AD) + '"'
            fst['ServoFile']    = '"' + os.path.basename(filename_AD) + '"'
            fst['HydroFile']    = '"' + os.path.basename(filename_HD) + '"'
            fst['SubFile']      = '"' + os.path.basename(filename_SD) + '"'
            fst['MooringFile']  = '"' + os.path.basename(filename_MD) + '"'
            fst['IceFile']      = '"' + os.path.basename(filename_Ice)+ '"'
            fst.write(filename)


            ED =  self.fst_vt['ElastoDyn']
            if fst['CompElast']>0:
                ED['TwrFile'] = '"' + os.path.basename(filename_ED_twr)+ '"'
                self.fst_vt['ElastoDynTower'].write(filename_ED_twr)
            if fst['CompElast']==1:
                if 'BldFile1' in ED.keys():
                    ED['BldFile1'] = '"' + os.path.basename(filename_ED_bld)+ '"'
                    ED['BldFile2'] = '"' + os.path.basename(filename_ED_bld)+ '"'
                    ED['BldFile3'] = '"' + os.path.basename(filename_ED_bld)+ '"'
                else:
                    ED['BldFile(1)']   = '"' + os.path.basename(filename_ED_bld)+ '"'
                    ED['BldFile(2)']   = '"' + os.path.basename(filename_ED_bld)+ '"'
                    ED['BldFile(3)']   = '"' + os.path.basename(filename_ED_bld)+ '"'
                self.fst_vt['ElastoDynBlade'].write(filename_ED_bld)

            elif fst['CompElast']==2:
                BD = self.fst_vt['BeamDyn'] 
                BD['BldFile'] = '"'+os.path.basename(filename_BD_bld)+'"'
                self.fst_vt['BeamDynBlade'].write(filename_BD_bld)  # TODO TODO pick up the proper blade file!
                BD.write(filename_BD)
            ED.write(filename_ED)


            if fst['CompInflow']>0:
                self.fst_vt['InflowWind'].write(filename_IW)

            if fst['CompAero']>0:
                self.fst_vt['AeroDyn15'].write(filename_AD)
                # TODO other files

            if fst['CompServo']>0:
                self.fst_vt['ServoDyn'].write(filename_SvD)

            if fst['CompHydro']==1:
                self.fst_vt['HydroDyn'].write(filename_HD)

            if fst['CompSub']==1:
                self.fst_vt['SubDyn'].write(filename_SD)
            elif fst['CompSub']==2:
                raise NotImplementedError()

            if fst['CompMooring']==1:
                self.fst_vt['MAP'].write(filename_MD)
            if self.fst_vt['Fst']['CompMooring']==2:
                self.fst_vt['MoorDyn'].write(filename_MD)

        return filename



    def __repr__(self):
        s='<weio.FastInputDeck object>'+'\n'
        s+='filename   : '+self.filename+'\n'
        s+='readlist   : {}'.format(self.readlist)+'\n'
        s+='version    : '+self.version+'\n'
        s+='AD version : '+self.ADversion+'\n'
        s+='fst_vt     : dict{'+','.join([k for k,v in self.fst_vt.items() if v is not None])+'}\n'
        s+='inputFiles : {}\n'.format(self.inputFiles)
        s+='inputFilesRead : {}\n'.format(self.inputFilesRead)
        s+='\n'
        return s

if __name__ == "__main__":
    fst=FASTInputDeck('NREL5MW.fst')
    print(fst)
