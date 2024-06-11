'''
Created on 20/01/2014

See documentation of HTCFile below

'''
# from wetb.utils.process_exec import pexec
# from wetb.hawc2.hawc2_pbs_file import HAWC2PBSFile
# import jinja2
# from wetb.utils.cluster_tools.os_path import fixcase, abspath, pjoin

from collections import OrderedDict
from .htc_contents import HTCContents, HTCSection, HTCLine
from .htc_extensions import HTCDefaults, HTCExtensions
import os

# --- cluster_tools/os_path
def fmt_path(path):
    return path.lower().replace("\\", "/")

def repl(path):
    return path.replace("\\", "/")

def abspath(path):
    return repl(os.path.abspath(path))

def relpath(path, start=None):
    return repl(os.path.relpath(path, start))

def realpath(path):
    return repl(os.path.realpath(path))

def pjoin(*path):
    return repl(os.path.join(*path))

def fixcase(path):
    path = realpath(str(path)).replace("\\", "/")
    p, rest = os.path.splitdrive(path)
    p += "/"
    for f in rest[1:].split("/"):
        f_lst = [f_ for f_ in os.listdir(p) if f_.lower() == f.lower()]
        if len(f_lst) > 1:
            # use the case sensitive match
            f_lst = [f_ for f_ in f_lst if f_ == f]
        if len(f_lst) == 0:
            raise IOError("'%s' not found in '%s'" % (f, p))
        # Use matched folder
        p = pjoin(p, f_lst[0])
    return p
# --- end os_path

class HTCFile(HTCContents, HTCDefaults, HTCExtensions):
    """Wrapper for HTC files

    Examples:
    ---------
    >>> htcfile = HTCFile('htc/test.htc')
    >>> htcfile.wind.wsp = 10
    >>> htcfile.save()

    #---------------------------------------------
    >>> htc = HTCFile(filename=None, modelpath=None) # create minimal htcfile

    #Add section
    >>> htc.add_section('hydro')

    #Add subsection
    >>> htc.hydro.add_section("hydro_element")

    #Set values
    >>> htc.hydro.hydro_element.wave_breaking = [2, 6.28, 1] # or
    >>> htc.hydro.hydro_element.wave_breaking = 2, 6.28, 1

    #Set comments
    >>> htc.hydro.hydro_element.wave_breaking.comments = "This is a comment"

    #Access section
    >>> hydro_element = htc.hydro.hydro_element #or
    >>> hydro_element = htc['hydro.hydro_element'] # or
    >>> hydro_element = htc['hydro/hydro_element'] # or
    >>> print (hydro_element.wave_breaking) #string represenation
    wave_breaking    2 6.28 1;    This is a comment
    >>> print (hydro_element.wave_breaking.name_) # command
    wave_breaking
    >>> print (hydro_element.wave_breaking.values) # values
    [2, 6.28, 1
    >>> print (hydro_element.wave_breaking.comments) # comments
    This is a comment
    >>> print (hydro_element.wave_breaking[0]) # first value
    2

    #Delete element
    htc.simulation.logfile.delete()
    #or
    del htc.simulation.logfile #Delete logfile line. Raise keyerror if not exists

    """

    filename = None
    jinja_tags = {}
    htc_inputfiles = []
    level = 0
    modelpath = "../"
    initial_comments = None

    def __init__(self, filename=None, modelpath=None, jinja_tags={}):
        """
        Parameters
        ---------
        filename : str
            Absolute filename of htc file
        modelpath : str
            Model path relative to htc file
        """
        if filename is not None:
            try:
                filename = fixcase(abspath(filename))
                with self.open(str(filename)):
                    pass
            except Exception:
                pass

            self.filename = filename

        self.jinja_tags = jinja_tags
        self.modelpath = modelpath or self.auto_detect_modelpath()

        if filename and self.modelpath != "unknown" and not os.path.isabs(self.modelpath):
            drive, p = os.path.splitdrive(os.path.join(os.path.dirname(str(self.filename)), self.modelpath))
            self.modelpath = os.path.join(drive, os.path.splitdrive(os.path.realpath(p))[1]).replace("\\", "/")
        if self.modelpath != 'unknown' and self.modelpath[-1] != '/':
            self.modelpath += "/"

        self.load()

    def auto_detect_modelpath(self):
        if self.filename is None:
            return "../"

        #print (["../"*i for i in range(3)])
        import numpy as np
        input_files = HTCFile(self.filename, 'unknown').input_files()
        if len(input_files) == 1:  # only input file is the htc file
            return "../"
        rel_input_files = [f for f in input_files if not os.path.isabs(f)]

        def isfile_case_insensitive(f):
            try:
                f = fixcase(f)  # raises exception if not existing
                return os.path.isfile(f)
            except IOError:
                return False
        found = ([np.sum([isfile_case_insensitive(os.path.join(os.path.dirname(self.filename), "../" * i, f))
                          for f in rel_input_files]) for i in range(4)])

        if max(found) > 0:
            relpath = "../" * np.argmax(found)
            return abspath(pjoin(os.path.dirname(self.filename), relpath))
        else:
            print("Modelpath cannot be autodetected for '%s'.\nInput files not found near htc file" % self.filename)
            return 'unknown'

    def load(self):
        self.contents = OrderedDict()
        self.initial_comments = []
        self.htc_inputfiles = []
        if self.filename is None:
            lines = self.empty_htc.split("\n")
        else:
            lines = self.readlines(self.filename)

        lines = [l.strip() for l in lines]

        #lines = copy(self.lines)
        while lines:
            if lines[0].startswith(";"):
                self.initial_comments.append(lines.pop(0).strip() + "\n")
            elif lines[0].lower().startswith("begin"):
                self._add_contents(HTCSection.from_lines(lines))
            else:
                line = HTCLine.from_lines(lines)
                if line.name_ == "exit":
                    break
                self._add_contents(line)

    def readfilelines(self, filename):
        with self.open(self.unix_path(os.path.abspath(filename.replace('\\', '/'))), encoding='cp1252') as fid:
            txt = fid.read()
        if txt[:10].encode().startswith(b'\xc3\xaf\xc2\xbb\xc2\xbf'):
            txt = txt[3:]
        if self.jinja_tags:
            template = jinja2.Template(txt)
            txt = template.render(**self.jinja_tags)
        return txt.replace("\r", "").split("\n")

    def readlines(self, filename):
        if filename != self.filename:  # self.filename may be changed by set_name/save. Added it when needed instead
            self.htc_inputfiles.append(filename)
        htc_lines = []
        lines = self.readfilelines(filename)
        for l in lines:
            if l.lower().lstrip().startswith('continue_in_file'):
                filename = l.lstrip().split(";")[0][len("continue_in_file"):].strip().lower()

                if self.modelpath == 'unknown':
                    p = os.path.dirname(self.filename)
                    try:
                        lu = [os.path.isfile(os.path.abspath(os.path.join(p, "../" * i, filename.replace("\\", "/"))))
                              for i in range(4)].index(True)
                        filename = os.path.join(p, "../" * lu, filename)
                    except ValueError:
                        print('[FAIL] Cannot continue in file: {}'.format(filename))
                        filename = None
                else:
                    filename = os.path.join(self.modelpath, filename)
                    if not os.path.isfile(filename):
                        print('[FAIL] Cannot continue in file: {}'.format(filename))
                        filename=None
                if filename is not None:
                    #print('[INFO] Continuing in file: {}'.format(filename))
                    for line in self.readlines(filename):
                        if line.lstrip().lower().startswith('exit'):
                            break
                        htc_lines.append(line)
            else:
                htc_lines.append(l)
        return htc_lines

    def __setitem__(self, key, value):
        self.contents[key] = value

    def __str__(self):
        self.contents  # load
        return "".join(self.initial_comments + [c.__str__(1) for c in self] + ["exit;"])

    def save(self, filename=None):
        """Saves the htc object to an htc file.

        Args:
            filename (str, optional): Specifies the filename of the htc file to be saved. 
            If the value is none, the filename attribute of the object will be used as the filename. 
            Defaults to None.
        """
        self.contents  # load if not loaded
        if filename is None:
            filename = self.filename
        else:
            self.filename = filename
        # exist_ok does not exist in Python27
        if not os.path.exists(os.path.dirname(filename)) and os.path.dirname(filename) != "":
            os.makedirs(os.path.dirname(filename))  # , exist_ok=True)
        with self.open(filename, 'w', encoding='cp1252') as fid:
            fid.write(str(self))

    def set_name(self, name, subfolder=''):
        """Sets the base filename of the simulation files. 

        Args:
            name (str): Specifies name of the log file, dat file (for animation), hdf5 file (for visualization) and htc file.
            subfolder (str, optional): Specifies the name of a subfolder to place the files in. 
                If the value is an empty string, no subfolders will be created. 
                Defaults to ''.

        Returns:
            None
        """
        # if os.path.isabs(folder) is False and os.path.relpath(folder).startswith("htc" + os.path.sep):
        self.contents  # load if not loaded

        def fmt_folder(folder, subfolder): return "./" + \
            os.path.relpath(os.path.join(folder, subfolder)).replace("\\", "/")

        self.filename = os.path.abspath(os.path.join(self.modelpath, fmt_folder(
            'htc', subfolder), "%s.htc" % name)).replace("\\", "/")
        if 'simulation' in self and 'logfile' in self.simulation:
            self.simulation.logfile = os.path.join(fmt_folder('log', subfolder), "%s.log" % name).replace("\\", "/")
            if 'animation' in self.simulation:
                self.simulation.animation = os.path.join(fmt_folder(
                    'animation', subfolder), "%s.dat" % name).replace("\\", "/")
            if 'visualization' in self.simulation:
                f = os.path.join(fmt_folder('visualization', subfolder), "%s.hdf5" % name).replace("\\", "/")
                self.simulation.visualization[0] = f
        elif 'test_structure' in self and 'logfile' in self.test_structure:  # hawc2aero
            self.test_structure.logfile = os.path.join(fmt_folder('log', subfolder), "%s.log" % name).replace("\\", "/")
        if 'output' in self:
            self.output.filename = os.path.join(fmt_folder('res', subfolder), "%s" % name).replace("\\", "/")

    def set_time(self, start=None, stop=None, step=None):
        self.contents  # load if not loaded
        if stop is not None:
            self.simulation.time_stop = stop
        else:
            stop = self.simulation.time_stop[0]
        if step is not None:
            self.simulation.newmark.deltat = step
        if start is not None:
            self.output.time = start, stop
            if "wind" in self:  # and self.wind.turb_format[0] > 0:
                self.wind.scale_time_start = start

    def expected_simulation_time(self):
        return 600

    def pbs_file(self, hawc2_path, hawc2_cmd, queue='workq', walltime=None,
                 input_files=None, output_files=None, copy_turb=(True, True)):
        walltime = walltime or self.expected_simulation_time() * 2
        if len(copy_turb) == 1:
            copy_turb_fwd, copy_turb_back = copy_turb, copy_turb
        else:
            copy_turb_fwd, copy_turb_back = copy_turb

        input_files = input_files or self.input_files()
        if copy_turb_fwd:
            input_files += [f for f in self.turbulence_files() if os.path.isfile(f)]

        output_files = output_files or self.output_files()
        if copy_turb_back:
            output_files += self.turbulence_files()

        return HAWC2PBSFile(hawc2_path, hawc2_cmd, self.filename, self.modelpath,
                            input_files, output_files,
                            queue, walltime)

    def input_files(self):
        self.contents  # load if not loaded
        if self.modelpath == "unknown":
            files = [str(f).replace("\\", "/") for f in [self.filename] + self.htc_inputfiles]
        else:
            files = [os.path.abspath(str(f)).replace("\\", "/") for f in [self.filename] + self.htc_inputfiles]
        if 'new_htc_structure' in self:
            for mb in [self.new_htc_structure[mb]
                       for mb in self.new_htc_structure.keys() if mb.startswith('main_body')]:
                if "timoschenko_input" in mb:
                    files.append(mb.timoschenko_input.filename[0])
                files.append(mb.get('external_bladedata_dll', [None, None, None])[2])
        if 'aero' in self:
            files.append(self.aero.ae_filename[0])
            files.append(self.aero.pc_filename[0])
            files.append(self.aero.get('external_bladedata_dll', [None, None, None])[2])
            files.append(self.aero.get('output_profile_coef_filename', [None])[0])
            if 'dynstall_ateflap' in self.aero:
                files.append(self.aero.dynstall_ateflap.get('flap', [None] * 3)[2])
            if 'bemwake_method' in self.aero:
                files.append(self.aero.bemwake_method.get('a-ct-filename', [None] * 3)[0])
        for dll in [self.dll[dll] for dll in self.get('dll', {}).keys() if 'filename' in self.dll[dll]]:
            files.append(dll.filename[0])
            f, ext = os.path.splitext(dll.filename[0])
            files.append(f + "_64" + ext)
        if 'wind' in self:
            files.append(self.wind.get('user_defined_shear', [None])[0])
            files.append(self.wind.get('user_defined_shear_turbulence', [None])[0])
            files.append(self.wind.get('met_mast_wind', [None])[0])
        if 'wakes' in self:
            files.append(self.wind.get('use_specific_deficit_file', [None])[0])
            files.append(self.wind.get('write_ct_cq_file', [None])[0])
            files.append(self.wind.get('write_final_deficits', [None])[0])
        if 'hydro' in self:
            if 'water_properties' in self.hydro:
                files.append(self.hydro.water_properties.get('water_kinematics_dll', [None])[0])
                files.append(self.hydro.water_properties.get('water_kinematics_dll', [None, None])[1])
        if 'soil' in self:
            if 'soil_element' in self.soil:
                files.append(self.soil.soil_element.get('datafile', [None])[0])
        try:
            dtu_we_controller = self.dll.get_subsection_by_name('dtu_we_controller')
            theta_min = dtu_we_controller.init.constant__5[1]
            if theta_min >= 90:
                files.append(os.path.join(os.path.dirname(
                    dtu_we_controller.filename[0]), "wpdata.%d" % theta_min).replace("\\", "/"))
        except Exception:
            pass

        try:
            files.append(self.force.dll.dll[0])
        except Exception:
            pass

        def fix_path_case(f):
            if os.path.isabs(f):
                return self.unix_path(f)
            elif self.modelpath != "unknown":
                try:
                    return "./" + os.path.relpath(self.unix_path(os.path.join(self.modelpath, f)),
                                                  self.modelpath).replace("\\", "/")
                except IOError:
                    return f
            else:
                return f
        return [fix_path_case(f) for f in set(files) if f]

    def output_files(self):
        self.contents  # load if not loaded
        files = []
        for k, index in [('simulation/logfile', 0),
                         ('simulation/animation', 0),
                         ('simulation/visualization', 0),
                         ('new_htc_structure/beam_output_file_name', 0),
                         ('new_htc_structure/body_output_file_name', 0),
                         ('new_htc_structure/struct_inertia_output_file_name', 0),
                         ('new_htc_structure/body_eigenanalysis_file_name', 0),
                         ('new_htc_structure/constraint_output_file_name', 0),
                         ('wind/turb_export/filename_u', 0),
                         ('wind/turb_export/filename_v', 0),
                         ('wind/turb_export/filename_w', 0)]:
            line = self.get(k)
            if line:
                files.append(line[index])
        if 'new_htc_structure' in self:
            if 'system_eigenanalysis' in self.new_htc_structure:
                f = self.new_htc_structure.system_eigenanalysis[0]
                files.append(f)
                files.append(os.path.join(os.path.dirname(f), 'mode*.dat').replace("\\", "/"))
            if 'structure_eigenanalysis_file_name' in self.new_htc_structure:
                f = self.new_htc_structure.structure_eigenanalysis_file_name[0]
                files.append(f)
                files.append(os.path.join(os.path.dirname(f), 'mode*.dat').replace("\\", "/"))
        files.extend(self.res_file_lst())

        for key in [k for k in self.contents.keys() if k.startswith("output_at_time")]:
            files.append(self[key]['filename'][0] + ".dat")
        return [f.lower() for f in files if f]

    def turbulence_files(self):
        self.contents  # load if not loaded
        if 'wind' not in self.contents.keys() or self.wind.turb_format[0] == 0:
            return []
        elif self.wind.turb_format[0] == 1:
            files = [self.get('wind.mann.filename_%s' % comp, [None])[0] for comp in ['u', 'v', 'w']]
        elif self.wind.turb_format[0] == 2:
            files = [self.get('wind.flex.filename_%s' % comp, [None])[0] for comp in ['u', 'v', 'w']]
        return [f for f in files if f]

    def res_file_lst(self):
        self.contents  # load if not loaded
        res = []
        for output in [self[k] for k in self.keys()
                       if self[k].name_.startswith("output") and not self[k].name_.startswith("output_at_time")]:
            dataformat = output.get('data_format', 'hawc_ascii')
            res_filename = output.filename[0]
            if dataformat[0] == "gtsdf" or dataformat[0] == "gtsdf64":
                res.append(res_filename + ".hdf5")
            elif dataformat[0] == "flex_int":
                res.extend([res_filename + ".int", os.path.join(os.path.dirname(res_filename), 'sensor')])
            else:
                res.extend([res_filename + ".sel", res_filename + ".dat"])
        return res

    def _simulate(self, exe, skip_if_up_to_date=False):
        self.contents  # load if not loaded
        if skip_if_up_to_date:
            from os.path import isfile, getmtime, isabs
            res_file = os.path.join(self.modelpath, self.res_file_lst()[0])
            htc_file = os.path.join(self.modelpath, self.filename)
            if isabs(exe):
                exe_file = exe
            else:
                exe_file = os.path.join(self.modelpath, exe)
            #print (from_unix(getmtime(res_file)), from_unix(getmtime(htc_file)))
            if (isfile(htc_file) and isfile(res_file) and isfile(exe_file) and
                str(HTCFile(htc_file)) == str(self) and
                    getmtime(res_file) > getmtime(htc_file) and getmtime(res_file) > getmtime(exe_file)):
                if "".join(self.readfilelines(htc_file)) == str(self):
                    return

        self.save()
        htcfile = os.path.relpath(self.filename, self.modelpath)
        assert any([os.path.isfile(os.path.join(f, exe)) for f in [''] + os.environ['PATH'].split(";")]), exe
        return pexec([exe, htcfile], self.modelpath)

    def simulate(self, exe, skip_if_up_to_date=False):
        errorcode, stdout, stderr, cmd = self._simulate(exe, skip_if_up_to_date)
        if ('simulation' in self.keys() and "logfile" in self.simulation and
                os.path.isfile(os.path.join(self.modelpath, self.simulation.logfile[0]))):
            with self.open(os.path.join(self.modelpath, self.simulation.logfile[0])) as fid:
                log = fid.read()
        else:
            log = "%s\n%s" % (str(stdout), str(stderr))

        if errorcode or 'Elapsed time' not in log:
            log_lines = log.split("\n")
            error_lines = [i for i, l in enumerate(log_lines) if 'error' in l.lower()]
            if error_lines:
                import numpy as np
                line_i = np.r_[np.array([error_lines + i for i in np.arange(-3, 4)]).flatten(),
                               np.arange(-5, 0) + len(log_lines)]
                line_i = sorted(np.unique(np.maximum(np.minimum(line_i, len(log_lines) - 1), 0)))

                lines = ["%04d %s" % (i, log_lines[i]) for i in line_i]
                for jump in np.where(np.diff(line_i) > 1)[0]:
                    lines.insert(jump, "...")

                error_log = "\n".join(lines)
            else:
                error_log = log
            raise Exception("\nError code: %s\nstdout:\n%s\n--------------\nstderr:\n%s\n--------------\nlog:\n%s\n--------------\ncmd:\n%s" %
                            (errorcode, str(stdout), str(stderr), error_log, cmd))
        return str(stdout) + str(stderr), log

    def simulate_hawc2stab2(self, exe):
        errorcode, stdout, stderr, cmd = self._simulate(exe, skip_if_up_to_date=False)

        if errorcode:
            raise Exception("\nstdout:\n%s\n--------------\nstderr:\n%s\n--------------\ncmd:\n%s" %
                            (str(stdout), str(stderr), cmd))
        return str(stdout) + str(stderr)

    def deltat(self):
        return self.simulation.newmark.deltat[0]

    def compare(self, other):
        if isinstance(other, str):
            other = HTCFile(other)
        return HTCContents.compare(self, other)

    @property
    def open(self):
        return open

    def unix_path(self, filename):
        filename = os.path.realpath(str(filename)).replace("\\", "/")
        ufn, rest = os.path.splitdrive(filename)
        ufn += "/"
        for f in rest[1:].split("/"):
            f_lst = [f_ for f_ in os.listdir(ufn) if f_.lower() == f.lower()]
            if len(f_lst) > 1:
                # use the case sensitive match
                f_lst = [f_ for f_ in f_lst if f_ == f]
            if len(f_lst) == 0:
                raise IOError("'%s' not found in '%s'" % (f, ufn))
            else:  # one match found
                ufn = os.path.join(ufn, f_lst[0])
        return ufn.replace("\\", "/")


#
#     def get_body(self, name):
#         lst = [b for b in self.new_htc_structure if b.name_=="main_body" and b.name[0]==name]
#         if len(lst)==1:
#             return lst[0]
#         else:
#             if len(lst)==0:
#                 raise ValueError("Body '%s' not found"%name)
#             else:
#                 raise NotImplementedError()
#

class H2aeroHTCFile(HTCFile):
    def __init__(self, filename=None, modelpath=None):
        HTCFile.__init__(self, filename=filename, modelpath=modelpath)

    @property
    def simulation(self):
        return self.test_structure

    def set_time(self, start=None, stop=None, step=None):
        if stop is not None:
            self.test_structure.time_stop = stop
        else:
            stop = self.simulation.time_stop[0]
        if step is not None:
            self.test_structure.deltat = step
        if start is not None:
            self.output.time = start, stop
            if "wind" in self and self.wind.turb_format[0] > 0:
                self.wind.scale_time_start = start


class SSH_HTCFile(HTCFile):
    def __init__(self, ssh, filename=None, modelpath=None):
        object.__setattr__(self, 'ssh', ssh)
        HTCFile.__init__(self, filename=filename, modelpath=modelpath)

    @property
    def open(self):
        return self.ssh.open

    def unix_path(self, filename):
        rel_filename = os.path.relpath(filename, self.modelpath).replace("\\", "/")
        _, out, _ = self.ssh.execute("find -ipath ./%s" % rel_filename, cwd=self.modelpath)
        out = out.strip()
        if out == "":
            raise IOError("'%s' not found in '%s'" % (rel_filename, self.modelpath))
        elif "\n" in out:
            raise IOError("Multiple '%s' found in '%s' (due to case senitivity)" % (rel_filename, self.modelpath))
        else:
            drive, path = os.path.splitdrive(os.path.join(self.modelpath, out))
            path = os.path.realpath(path).replace("\\", "/")
            return os.path.join(drive, os.path.splitdrive(path)[1])


if "__main__" == __name__:
    f = HTCFile(r"C:/Work/BAR-Local/Hawc2ToBeamDyn/sim.htc", ".")
    print(f.input_files())
#     f.save(r"C:\mmpe\HAWC2\models\DTU10MWRef6.0\htc\DTU_10MW_RWT_power_curve.htc")
#
#     f = HTCFile(r"C:\mmpe\HAWC2\models\DTU10MWRef6.0\htc\DTU_10MW_RWT.htc", "../")
#     f.set_time = 0, 1, .1
#     print(f.simulate(r"C:\mmpe\HAWC2\bin\HAWC2_12.8\hawc2mb.exe"))
#
#     f.save(r"C:\mmpe\HAWC2\models\DTU10MWRef6.0\htc\DTU_10MW_RWT.htc")
