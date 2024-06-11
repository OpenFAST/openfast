""" 
Read/Write Mann Box input file

"""
import os
import numpy as np
try:
    from .file import File, EmptyFileError, BrokenFormatError
except:
    EmptyFileError = type('EmptyFileError', (Exception,),{})
    BrokenFormatError = type('BrokenFormatError', (Exception,),{})
    File=dict

class MannBoxInputFile(File):
    def __init__(self, filename=None, **kwargs):
        self.filename = None
        # Default Init
        self['fieldDim']  = None
        self['nComp']     = None
        self['nx']        = None
        self['ny']        = None
        self['nz']        = None
        self['Lx']        = None
        self['Ly']        = None
        self['Lz']        = None
        self['type']      = None
        self['U']         = None
        self['z']         = None
        self['zNone']     = None
        self['spectrum']  = None
        # Basic
        self['alpha_eps'] = None
        self['L']         = None
        self['Gamma']     = None
        self['seed']      = -1  
        self['file_u']    = None
        self['file_v']    = None
        self['file_w']    = None
        for k,v in kwargs.items():
            print('>>> Setting',k,v)

        if filename:
            self.read(filename=filename)

    def read(self, filename):
        if filename:
            self.filename = filename
        if not self.filename:
            raise Exception('No filename provided')
        if not os.path.isfile(self.filename):
            raise OSError(2,'File not found:',self.filename)
        if os.stat(self.filename).st_size == 0:
            raise EmptyFileError('File is empty:',self.filename)

        with open(self.filename) as f:
            print('IO: Reading Mann input file: '+self.filename)
            self['fieldDim'] = int(f.readline())
            self['nComp']    = int(f.readline())
            if(self['nComp']==3):
                self['nx'] = int(f.readline())
                self['ny'] = int(f.readline())
                self['nz'] = int(f.readline())
                self['Lx'] = float(f.readline())
                self['Ly'] = float(f.readline())
                self['Lz'] = float(f.readline())
            else:
                raise Exception('nComp=2')
            self['type']=f.readline().strip()
            if(self['type']=='basic'):
                self['alpha_eps'] = float(f.readline())
                self['L']         = float(f.readline())
                self['Gamma']     = float(f.readline())
            else:
                raise Exception('not basic')
            self['seed']=int(f.readline())
            if(self['nComp']>=1):
                self['file_u']=f.readline().strip()
            if(self['nComp']>=2):
                self['file_v']=f.readline().strip()
            if(self['nComp']>=3):
                self['file_w']=f.readline().strip()
    
    def write(self, filename):
        """ Write mann box """
        if filename:
            self.filename = filename

        self.defaultFileNames(self.filename)

        with open(self.filename, 'w') as f:
            f.write(str(self['fieldDim'])+'\n')
            f.write(str(self['nComp'])+'\n')
            if(self['nComp']==3):
                f.write(str(self['nx'])+'\n')
                f.write(str(self['ny'])+'\n')
                f.write(str(self['nz'])+'\n')
                f.write(str(self['Lx'])+'\n')
                f.write(str(self['Ly'])+'\n')
                f.write(str(self['Lz'])+'\n')
            else:
                raise Exception('nComp=2')
            
            f.write(self['type']+'\n')
            if(self['type']=='basic'):
                f.write(str(self['alpha_eps'])+'\n');
                f.write(str(self['L'])+'\n');
                f.write(str(self['Gamma'])+'\n');
            else:
                raise Exception('not basic')
            
            f.write(str(self['seed'])+'\n');
            if(self['nComp']>=1):
                f.write(self['file_u']+'\n');
            if(self['nComp']>=2):
                f.write(self['file_v']+'\n');
            if(self['nComp']>=3):
                f.write(self['file_w']+'\n');

    def defaultFileNames(self, inpfile):
        base = os.path.splitext(os.path.basename(inpfile))[0]
        nx,ny,nz = self['nx'], self['ny'], self['nz']
        self['file_u'] = base+'_{}x{}x{}.u'.format(nx,ny,nz)
        self['file_v'] = base+'_{}x{}x{}.v'.format(nx,ny,nz)
        self['file_w'] = base+'_{}x{}x{}.w'.format(nx,ny,nz)

#     def set3Disotropic(self,nx,ny,nz,Lx,Ly,Lz,alpha_eps,L):
#         self.fieldDim=3;
#         self.NComp=3;
#         self.nx=nx;
#         self.ny=ny;
#         self.nz=nz;
#         self.Lx=Lx;
#         self.Ly=Ly;
#         self.Lz=Lz;
#         self.alpha_eps=alpha_eps;
#         self.L=L;
#         # Isotropic
#         self.Gamma=0.0; 
#         self.type='basic';
# 
#     def auto_file_names(self,folder=''):
#         if folder is None:
#             folder=''
#         import os
#         if(self.type=='basic' and self.Gamma==0):
#             if self.fieldDim==3:
#                 filenameBase='Isotropic3D_'+str(self.nx)+'_'+str(self.ny)+'_'+str(self.nz)+'_h'+str(int(100*float(self.Lx)/float(self.nx-1))/100.)
#                 self.filename = os.path.join(folder,filenameBase+'.maninp')
#                 self.file_u   = os.path.join(folder,filenameBase+'_u.dat')
#                 self.file_v   = os.path.join(folder,filenameBase+'_v.dat')
#                 self.file_w   = os.path.join(folder,filenameBase+'_w.dat')
# 

    def __repr__(self):
        s = '<{} object> with keys:\n'.format(type(self).__name__)
        for k,v in self.items():
            s += '{:15s}: {}\n'.format(k,v)
        return s



# 
#     def get_outfile_u(self):
#         return self.file_u
#     def get_outfile_v(self):
#         return self.file_v
#     def get_outfile_w(self):
#         return self.file_w
# 
# if __name__ == "__main__":
#     import sys
#     if len(sys.argv)>1:
#         print 'called with: ', int(sys.argv[1])
# 
#     nx=1024;
#     ny=256;
#     nz=256;
#     Lx=100;
#     Ly=10;
#     Lz=10;
#     alpha_eps=0.01
#     L=10;
# 
#     inp=MannInputFile();
#     inp.set3Disotropic(nx,ny,nz,Lx,Ly,Lz,alpha_eps,L)
#     inp.auto_file_names()
#     inp.write();
# #     inp2=MannInputFile(inp.filename);
# #     inp2.filename='out10.inp'
# #     inp2.write()
# 
