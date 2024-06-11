# --- Making main readers available
from .csv_file import CSVFile
from .excel_file import ExcelFile
from .fast_input_deck import FASTInputDeck
from .fast_input_file import FASTInputFile
from .fast_linearization_file import FASTLinearizationFile
from .fast_output_file import FASTOutputFile
from .fast_summary_file import FASTSummaryFile
from .fast_wind_file import FASTWndFile
# from .bmodes_out_file         import BModesOutFile
from .hawc2_pc_file           import HAWC2PCFile
from .hawc2_ae_file           import HAWC2AEFile
from .hawc2_dat_file          import HAWC2DatFile
from .hawc2_htc_file          import HAWC2HTCFile
from .hawc2_st_file           import HAWC2StFile
# from .hawcstab2_pwr_file      import HAWCStab2PwrFile
# from .hawcstab2_ind_file      import HAWCStab2IndFile
# from .hawcstab2_cmb_file      import HAWCStab2CmbFile
from .mannbox_file            import MannBoxFile 
# from .flex_blade_file         import FLEXBladeFile
# from .flex_profile_file       import FLEXProfileFile
# from .flex_out_file           import FLEXOutFile
# from .flex_doc_file           import FLEXDocFile
# from .flex_wavekin_file       import FLEXWaveKinFile
# from .turbsim_ts_file         import TurbSimTSFile
from .turbsim_file            import TurbSimFile
# from .netcdf_file             import NetCDFFile
# from .tdms_file               import TDMSFile
# from .tecplot_file            import TecplotFile 
from .vtk_file                import VTKFile
# from .bladed_out_file         import BladedFile
# from .parquet_file            import ParquetFile
# from .cactus_file             import CactusFile
# from .rosco_performance_file  import ROSCOPerformanceFile
from .raawmat_file import RAAWMatFile


# --- Generic reader / fileformat detection
from .file  import File, WrongFormatError, BrokenFormatError, FileNotFoundError, EmptyFileError, OptionalImportError
from .file_formats  import FileFormat, isRightFormat
import sys
import os
import numpy as np

class FormatNotDetectedError(Exception):
    pass

class UserFormatImportError(Exception):
    pass


_FORMATS=None

def fileFormats(userpath=None, ignoreErrors=False, verbose=False):
    """ return list of fileformats supported by the library
    If userpath is provided, 

    OUTPUTS:
      if ignoreErrors is True:
          formats,  errors
      else:
          formats

    """
    global _FORMATS
    errors=[]
    if _FORMATS is not None:
        if ignoreErrors:
            return _FORMATS, errors
        else:
            return _FORMATS
    # --- Library formats
    from .fast_input_file         import FASTInputFile
    from .fast_output_file        import FASTOutputFile
    from .csv_file                import CSVFile
    from .fast_wind_file          import FASTWndFile
    from .fast_linearization_file import FASTLinearizationFile
    from .fast_summary_file       import FASTSummaryFile
    from .bmodes_out_file         import BModesOutFile
    from .hawc2_pc_file           import HAWC2PCFile
    from .hawc2_ae_file           import HAWC2AEFile
    from .hawc2_dat_file          import HAWC2DatFile
    from .hawc2_htc_file          import HAWC2HTCFile
    from .hawc2_st_file           import HAWC2StFile
    from .hawcstab2_pwr_file      import HAWCStab2PwrFile
    from .hawcstab2_ind_file      import HAWCStab2IndFile
    from .hawcstab2_cmb_file      import HAWCStab2CmbFile
    from .mannbox_file            import MannBoxFile 
    from .flex_blade_file         import FLEXBladeFile
    from .flex_profile_file       import FLEXProfileFile
    from .flex_out_file           import FLEXOutFile
    from .flex_doc_file           import FLEXDocFile
    from .flex_wavekin_file       import FLEXWaveKinFile
    from .excel_file              import ExcelFile
    from .turbsim_ts_file         import TurbSimTSFile
    from .turbsim_file            import TurbSimFile
    from .netcdf_file             import NetCDFFile
    from .tdms_file               import TDMSFile
    from .tecplot_file            import TecplotFile 
    from .vtk_file import VTKFile
    from .bladed_out_file         import BladedFile
    from .parquet_file            import ParquetFile
    from .pickle_file             import PickleFile        
    from .cactus_file             import CactusFile
    from .raawmat_file            import RAAWMatFile
    from .rosco_discon_file       import ROSCODISCONFile
    from .rosco_performance_file  import ROSCOPerformanceFile
    priorities = []
    formats = []
    def addFormat(priority, fmt):
        priorities.append(priority)
        formats.append(fmt)
    addFormat(0, FileFormat(CSVFile))
    addFormat(0, FileFormat(ExcelFile))
    addFormat(10, FileFormat(TecplotFile))
    addFormat(10, FileFormat(BladedFile))
    addFormat(20, FileFormat(FASTInputFile))
    addFormat(20, FileFormat(FASTOutputFile))
    addFormat(20, FileFormat(FASTWndFile))
    addFormat(20, FileFormat(FASTLinearizationFile))
    addFormat(20, FileFormat(FASTSummaryFile))
    addFormat(20, FileFormat(TurbSimTSFile))
    addFormat(20, FileFormat(TurbSimFile))
    addFormat(30, FileFormat(HAWC2DatFile))
    addFormat(30, FileFormat(HAWC2HTCFile))
    addFormat(30, FileFormat(HAWC2StFile))
    addFormat(30, FileFormat(HAWC2PCFile))
    addFormat(30, FileFormat(HAWC2AEFile))
    addFormat(30, FileFormat(HAWCStab2PwrFile))
    addFormat(30, FileFormat(HAWCStab2IndFile))
    addFormat(30, FileFormat(HAWCStab2CmbFile))
    addFormat(30, FileFormat(MannBoxFile))
    addFormat(40, FileFormat(FLEXBladeFile))
    addFormat(40, FileFormat(FLEXProfileFile))
    addFormat(40, FileFormat(FLEXOutFile))
    addFormat(40, FileFormat(FLEXWaveKinFile))
    addFormat(40, FileFormat(FLEXDocFile))
    addFormat(50, FileFormat(BModesOutFile))
    addFormat(50, FileFormat(ROSCODISCONFile))
    addFormat(50, FileFormat(ROSCOPerformanceFile))
    addFormat(60, FileFormat(NetCDFFile))
    addFormat(60, FileFormat(VTKFile))
    addFormat(60, FileFormat(TDMSFile))
    addFormat(60, FileFormat(ParquetFile))
    addFormat(60, FileFormat(PickleFile))
    addFormat(70, FileFormat(CactusFile))
    addFormat(70, FileFormat(RAAWMatFile))

    # --- User defined formats from user path
    UserClasses, UserPaths, UserModules, UserModuleNames, errors = userFileClasses(userpath, ignoreErrors, verbose=verbose)
    for cls, f in zip(UserClasses, UserPaths):
        try:
            ff = FileFormat(cls)
        except Exception as e:
            s='Error registering a user fileformat.\n\nThe module location was: {}\n\nThe class name was: {}\n\nMake sure the class has `defaultExtensions` and `formatName` as static methods.\n\nThe exception was:\n{}'.format(f, cls.__name__, e)
            if ignoreErrors:
                errors.append(s)
                continue
            else:
                raise UserFormatImportError(s)
        # Use class.priority 
        try:
            priority = cls.priority()
        except:
            priority=2
        addFormat(priority, ff)

    # --- Sort fileformats by priorities
    formats = np.asarray(formats)[np.argsort(priorities, kind='stable')]

    _FORMATS=formats
    if ignoreErrors:
        return formats, errors
    else:
        return formats



def userFileClasses(userpath=None, ignoreErrors=False, verbose=True):
    """ return list of user file class in UserData folder"""
    if userpath is None:
        dataDir = defaultUserDataDir()
        userpath = os.path.join(dataDir, 'weio')
    errors          = []
    UserClasses     = []
    UserPaths       = []
    UserModules     = []
    UserModuleNames = []
    if os.path.exists(userpath):
        if verbose:
            print('>>> Looking for user modules in folder:',userpath)
        import glob
        from importlib.machinery import SourceFileLoader
        import inspect
        pyfiles = glob.glob(os.path.join(userpath,'*.py'))
        # Loop through files, look for classes of the form ClassNameFile, 
        for f in pyfiles:
            if f in ['__init__.py']:
                continue
            mod_name = os.path.basename(os.path.splitext(f)[0])
            try:
                if verbose:
                    print('>>> Trying to load user module:',f)
                module = SourceFileLoader(mod_name,f).load_module()
            except Exception as e:
                s='Error importing a user module.\n\nThe module location was: {}\n\nTry importing this module to debug it.\n\nThe Exception was:\n{}'.format(f, e)
                if ignoreErrors:
                    errors.append(s)
                    continue
                else:
                    raise UserFormatImportError(s)
            found=False
            for name, obj in inspect.getmembers(module):
                if inspect.isclass(obj):
                    classname = obj.__name__.lower()
                    if classname!='file' and classname.find('file')>=0 and classname.find('error')<0:
                        if verbose:
                            print('    Found File class with name:',obj.__name__)
                        UserClasses.append(obj)
                        UserPaths.append(f)
                        UserModules.append(module)
                        UserModuleNames.append(mod_name)
                        found=True # allowing only one class per file for now..
                        break
            if not found:
                s='Error finding a class named "*File" in the user module.\n\nThe module location was: {}\n\nNo class containing the string "File" in its name was found.'.format(f)
                if ignoreErrors:
                    errors.append(s)
                else:
                    raise UserFormatImportError(s)
    return UserClasses, UserPaths, UserModules, UserModuleNames, errors


def defaultUserDataDir():
    """
    Returns a parent directory path
    where persistent application data can be stored.
    # linux: ~/.local/share
    # macOS: ~/Library/Application Support
    # windows: C:/Users/<USER>/AppData/Roaming
    """
    home = os.path.expanduser('~')
    ptfm = sys.platform
    if ptfm == "win32":
        return os.path.join(home , 'AppData','Roaming')
    elif ptfm.startswith("linux"):
        return os.path.join(home, '.local', 'share')
    elif ptfm == "darwin":
        return os.path.join(home, 'Library','Application Support')
    else:
        print('>>>>>>>>>>>>>>>>> Unknown Platform', sys.platform)
        return './UserData'



def detectFormat(filename, **kwargs):
    """ Detect the file formats by looping through the known list. 
        The method may simply try to open the file, if that's the case
        the read file is returned. """
    import os
    import re
    global _FORMATS
    if _FORMATS is None:
        formats=fileFormats()
    else:
        formats=_FORMATS
    ext = os.path.splitext(filename.lower())[1]
    detected = False
    i = 0 
    while not detected and i<len(formats):
        myformat = formats[i]
        if ext in myformat.extensions:
            extMatch = True
        else:
            # Try patterns if present
            extPatterns = [ef.replace('.',r'\.').replace('$',r'\$').replace('*','[.]*') for ef in myformat.extensions if '*' in ef]
            if len(extPatterns)>0:
                extPatMatch = [re.match(pat, ext) is not None for pat in extPatterns]
                extMatch = any(extPatMatch)
            else:
                extMatch = False
        if extMatch: # we have a match on the extension
            valid, F = isRightFormat(myformat, filename, **kwargs)
            if valid:
                #print('File detected as :',myformat)
                detected=True
                return myformat,F

        i += 1

    if not detected:
        raise FormatNotDetectedError('The file format could not be detected for the file: '+filename)

def read(filename, fileformat=None, **kwargs):
    F = None
    if not os.path.exists(filename):
        raise FileNotFoundError('weio cannot read the following file because it does not exist:\n   Inp. path: {}\n   Abs. path: {}'.format(filename, os.path.abspath(filename)))
    # Detecting format if necessary
    if fileformat is None:
        fileformat,F = detectFormat(filename, **kwargs)
    # Reading the file with the appropriate class if necessary
    if not isinstance(F, fileformat.constructor):
        F=fileformat.constructor(filename=filename)
    return F



