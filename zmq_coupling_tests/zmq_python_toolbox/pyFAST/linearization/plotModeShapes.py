from paraview.simple import * 
import os 
import sys
import glob
import numpy as np

if len(sys.argv) not in [3,4]:
    print('Error: Not enough argument provided.')
    print("""
usage:
    plotModeShapes ROOT_PATH STATEFILE [OUTPUTDIR]

where:
    ROOT_PATH: Path of the form VTKDIR\ROOTNAME
               where VTKDIR is the directory where vtp files are located 
               Modes and LinTimes will be selected using the pattern:
               ROOT_PATH.Mode*.LinTime*.vtp
    STATEFILE: State file used for visualization

""")
    sys.exit(-1)

# --- Input Parameters
RootPath    = os.path.normpath(sys.argv[1])
StateFile   = os.path.normpath(sys.argv[2])
if len(sys.argv)==4:
    OutputDir = os.path.normpath(sys.argv[3])
else:
    OutputDir='./'

# --- Script Parameters
Suffix = 'LinTime1.' # Depends on viz option VTKLinTim= 1 or 2
fps    =  3          # frames per second (rate to save in the .avi file)
nModes = 15           # number of modes to visualize

vFPS= np.linspace(5,1,nModes).astype(int)

# --- Derived params
parentDir   = os.path.dirname(RootPath)
rootSim     = os.path.basename(RootPath)
mainDirName = os.path.abspath(parentDir)

# --- Constants
StructureModule = 'ED'
BladeMesh       = "AD_Blade"

print('')
print('RootName   :',rootSim)
print('MainDirName:',mainDirName)
print('StateFile  :',StateFile)
print('')

for iMode in range(nModes):  # iMode starts at 0, so add 1
   rootMode    = rootSim+'.Mode{:d}.'.format(iMode+1)+Suffix
   absrootMode = os.path.join(mainDirName, rootMode)
   print('***' + absrootMode + '***')
   # nLinTimes = len(glob.glob(os.path.join(mainDirName,RootName)+'.Mode1.LinTime*.AD_Blade1..vtp'))
   
   # determine number of leading zeros in this mode shape
   nLeadingZeros = 0
   exists = False
   while (not exists) and nLeadingZeros < 6:
      nLeadingZeros = nLeadingZeros + 1
      txt = '{:0' + str(nLeadingZeros) + 'd}'
      fileLeadingZeros = txt.format(1)
      Blade1File = absrootMode + BladeMesh + '1.' + fileLeadingZeros + '.vtp'
      exists = os.path.isfile(Blade1File)

   #print(Blade1File)
   if not exists:
      print('  Could not find files to load.')
   else:
      LoadState(StateFile, LoadStateDataFileOptions='Choose File Names',
          DataDirectory=mainDirName,
          a5MW_Land_DLL_WTurbMode1LinTime1AD_Blade10FileName             =[absrootMode + BladeMesh + '1.' + fileLeadingZeros + '.vtp'],
          a5MW_Land_DLL_WTurbMode1LinTime1AD_Blade20FileName             =[absrootMode + BladeMesh + '2.' + fileLeadingZeros + '.vtp'],
          a5MW_Land_DLL_WTurbMode1LinTime1AD_Blade30FileName             =[absrootMode + BladeMesh + '3.' + fileLeadingZeros + '.vtp'],
          a5MW_Land_DLL_WTurbMode1LinTime1Blade1Surface0FileName         =[absrootMode + 'Blade1Surface.' + fileLeadingZeros + '.vtp'],
          a5MW_Land_DLL_WTurbMode1LinTime1Blade2Surface0FileName         =[absrootMode + 'Blade2Surface.' + fileLeadingZeros + '.vtp'],
          a5MW_Land_DLL_WTurbMode1LinTime1Blade3Surface0FileName         =[absrootMode + 'Blade3Surface.' + fileLeadingZeros + '.vtp'],
          a5MW_Land_DLL_WTurbMode1LinTime1ED_Hub0FileName                =[absrootMode + StructureModule + '_Hub.' + fileLeadingZeros + '.vtp'],
          a5MW_Land_DLL_WTurbMode1LinTime1ED_Nacelle0FileName            =[absrootMode + StructureModule + '_Nacelle.' + fileLeadingZeros + '.vtp'],
          a5MW_Land_DLL_WTurbMode1LinTime1ED_TowerLn2Mesh_motion0FileName=[absrootMode + StructureModule + '_TowerLn2Mesh_motion.' + fileLeadingZeros + '.vtp'],
          a5MW_Land_DLL_WTurbMode1LinTime1HubSurface0FileName            =[absrootMode + 'HubSurface.' + fileLeadingZeros + '.vtp'],
          a5MW_Land_DLL_WTurbMode1LinTime1NacelleSurface0FileName        =[absrootMode + 'NacelleSurface.' + fileLeadingZeros + '.vtp'],
          a5MW_Land_DLL_WTurbMode1LinTime1TowerSurface0FileName          =[absrootMode + 'TowerSurface.' + fileLeadingZeros + '.vtp']
      )
      ## find new sources
      # blade 1
      for iBlade in range(3):
         Blade = FindSource(rootMode + BladeMesh + str(iBlade+1) + '...vtp')
         SetActiveSource(Blade)
         ExtendFileSeries(Blade)
         Blade = FindSource(rootMode + 'Blade' + str(iBlade+1) + 'Surface...vtp')
         SetActiveSource(Blade)
         ExtendFileSeries(Blade)
      # Hub
      Hub = FindSource(rootMode + StructureModule + '_Hub...vtp')
      SetActiveSource(Hub)
      ExtendFileSeries(Hub)
      Hub = FindSource(rootMode + 'HubSurface...vtp')
      SetActiveSource(Hub)
      ExtendFileSeries(Hub)

      # nacelle
      Nacelle = FindSource(rootMode + StructureModule + '_Nacelle...vtp')
      SetActiveSource(Nacelle)
      ExtendFileSeries(Nacelle)
      Nacelle = FindSource(rootMode + 'NacelleSurface...vtp')
      SetActiveSource(Nacelle)
      ExtendFileSeries(Nacelle)
      
      # tower
      Tower = FindSource(rootMode + StructureModule + '_TowerLn2Mesh_motion...vtp')
      SetActiveSource(Tower)
      ExtendFileSeries(Tower)
      Tower = FindSource(rootMode + 'TowerSurface...vtp')
      SetActiveSource(Tower)
      ExtendFileSeries(Tower)

      #####
      SetActiveView(GetRenderView()) 
      #view = GetActiveView() 
      layout = GetLayout()
      
      animFile= os.path.join(OutputDir, rootSim+'.Mode{:d}.avi'.format(iMode+1))
      print('Saving animation... ',animFile, end='')
      WriteAnimation(animFile, viewOrLayout=layout, FrameRate=vFPS[iMode], ImageResolution=(1544,784), Compression=True)#  ImageResolution=(1544,784) 
#      SaveAnimation(rootMode + 'avi', viewOrLayout=layout, FrameRate=fps, ImageResolution=(1544,784) ) 
      # this .pvsm file defaults to (2734,1178) without ImageResolution arguments, resulting in a bunch of warnings
      # For some reason, ParaView is ignoring the FrameRate argument and always uses a value of 1.
      print(' Done.')


