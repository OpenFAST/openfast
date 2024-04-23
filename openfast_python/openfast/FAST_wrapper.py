import os
import subprocess
import platform
import time

class FAST_wrapper(object):

    def __init__(self, **kwargs):

        self.FAST_exe = None   # Path to executable
        self.FAST_InputFile = None   # FAST input file (ext=.fst)
        self.FAST_directory = None   # Path to fst directory files

        # Optional population class attributes from key word arguments
        for k, w in kwargs.items():
            try:
                setattr(self, k, w)
            except Exception:
                pass

        super(FAST_wrapper, self).__init__()

    def execute(self):

        self.input_file = os.path.join(self.FAST_directory, self.FAST_InputFile)

        try:
            if platform.system()!='Windows' and self.FAST_exe[-4:]=='.exe':
                self.FAST_exe = self.FAST_exe[:-4]
        except Exception:
            pass
        exec_str = []
        exec_str.append(self.FAST_exe)
        exec_str.append(self.FAST_InputFile)

        olddir = os.getcwd()
        os.chdir(self.FAST_directory)

        run_idx = 0
        start = time.time()
        while run_idx < 2:
            try:
                subprocess.run(exec_str, check=True)
                failed = False
                run_idx = 2
            except subprocess.CalledProcessError as e:
                if e.returncode > 1 and run_idx < 1: # This probably failed because of a temporary library access issue, retry
                    print('Error loading OpenFAST libraries, retrying.')
                    failed = False
                    run_idx += 1
                else: # Bad OpenFAST inputs, or we've already retried
                    print('OpenFAST Failed: {}'.format(e))
                    failed = True
                    run_idx = 2
            except Exception as e:
                print('OpenFAST Failed: {}'.format(e))
                failed = True
                run_idx = 2
                
        runtime = time.time() - start
        print('Runtime: \t{} = {:<6.2f}s'.format(self.FAST_InputFile, runtime))

        os.chdir(olddir)

        return failed
