import os
import numpy as np


class AEFile(object):

    """Read and write the HAWC2 AE (aerodynamic blade layout) file

    examples
    --------
    >>> aefile = AEFile(r"tests/test_files/NREL_5MW_ae.txt")
    >>> aefile.thickness(36) # Interpolated thickness at radius 36
    23.78048780487805
    >>> aefile.chord(36) # Interpolated chord at radius 36
    3.673
    >>> aefile.pc_set_nr(36) # pc set number at radius 36
    1
    >>> ae= AEFile()
        ae.add_set(radius=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
                   chord=[1.1, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1],
                   thickness=[100.0, 100.0, 90.0, 80.0, 70.0, 60.0, 50.0, 40.0, 30.0, 20.0, 10.0],
                   pc_set_id=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    >>> str(ae)
    1 r[m]           Chord[m]    T/C[%]  Set no.
    1 11
    0.00000000000000000e+00   1.10000000000000009e+00   1.00000000000000000e+02     1
    1.00000000000000006e-01   1.00000000000000000e+00   1.00000000000000000e+02     1
    2.00000000000000011e-01   9.00000000000000022e-01   9.00000000000000000e+01     1
    2.99999999999999989e-01   8.00000000000000044e-01   8.00000000000000000e+01     1
    4.00000000000000022e-01   6.99999999999999956e-01   7.00000000000000000e+01     1
    5.00000000000000000e-01   5.99999999999999978e-01   6.00000000000000000e+01     1
    5.99999999999999978e-01   5.00000000000000000e-01   5.00000000000000000e+01     1
    6.99999999999999956e-01   4.00000000000000022e-01   4.00000000000000000e+01     1
    8.00000000000000044e-01   2.99999999999999989e-01   3.00000000000000000e+01     1
    9.00000000000000022e-01   2.00000000000000011e-01   2.00000000000000000e+01     1
    1.00000000000000000e+00   1.00000000000000006e-01   1.00000000000000000e+01     1
    """

    cols = ['radius', 'chord', 'relative_thickness', 'setnr']

    def __init__(self, filename=None):
        self.ae_sets = {}
        if filename is not None:
            self._read_file(filename)

    def _value(self, radius, column, set_nr=1):
        ae_data = self.ae_sets[set_nr]
        if radius is None:
            return ae_data[:, column]
        else:
            return np.interp(radius, ae_data[:, 0], ae_data[:, column])

    def chord(self, radius=None, set_nr=1):
        return self._value(radius, 1, set_nr)

    def thickness(self, radius=None, set_nr=1):
        return self._value(radius, 2, set_nr)

    def radius_ae(self, radius=None, set_nr=1):
        radii = self.ae_sets[set_nr][:, 0]
        if radius:
            return radii[np.argmin(np.abs(radii - radius))]
        else:
            return radii

    def pc_set_nr(self, radius, set_nr=1):
        ae_data = self.ae_sets[set_nr]
        index = np.searchsorted(ae_data[:, 0], radius)
        index = max(1, index)
        # --- Emmanuel's addition
        maxRad = np.max(ae_data[:,0])
        index2 = np.argmin(np.abs(ae_data[:, 0]-radius))
        if abs(ae_data[index2,0]-radius)<1e-4*maxRad:
            # We are very close to an ae location, we use this set
            return ae_data[index2, 3]
        # Otherwise we look at index before or after
        setnrs = ae_data[index - 1:index + 1, 3]
        if setnrs[0] != setnrs[-1]:
            print('[WARN] AE file, at radius {}, should return a set between {}. Using first one.'.format(radius,setnrs))
        return setnrs[0]

    def add_set(self, radius, chord, thickness, pc_set_id, set_id=None):
        '''This method will add another set to the ae data'''
        if set_id is None:
            set_id = 1
            while set_id in self.ae_sets.keys():
                set_id += 1
        self.ae_sets[set_id] = np.array([radius, chord, thickness, pc_set_id]).T
        return set_id

    def __str__(self):
        '''This method will create a string that is formatted like an ae file with the data in this class'''
        n_sets = len(self.ae_sets)
        retval = str(n_sets) + ' r[m]           Chord[m]    T/C[%]  Set no.\n'
        for st_idx, st in self.ae_sets.items():
            retval += str(st_idx) + ' ' + str(len(st)) + '\n'
            for line in st:
                retval += '%25.17e %25.17e %25.17e %5d\n' % (line[0], line[1], line[2], line[3])
        return retval

    def save(self, filename):
        if not os.path.isdir(os.path.dirname(filename)):
            # fails if dirname is empty string
            if len(os.path.dirname(filename)) > 0:
                os.makedirs(os.path.dirname(filename))
        with open(filename, 'w') as fid:
            fid.write(str(self))

    def _read_file(self, filename):
        ''' This method will read in the ae data from a HAWC2 ae file'''
        with open(filename) as fid:
            lines = fid.readlines()
        nsets = int(lines[0].split()[0])
        lptr = 1
        self.ae_sets = {}
        for _ in range(1, nsets + 1):
            set_nr, n_rows = [int(v) for v in lines[lptr].split()[:2]]
            lptr += 1
            data = np.array([[float(v) for v in l.split()[:4]] for l in lines[lptr:lptr + n_rows]])
            self.ae_sets[set_nr] = data
            lptr += n_rows


def main():
    if __name__ == "__main__":
        ae = AEFile(os.path.dirname(__file__) + "/tests/test_files/NREL_5MW_ae.txt")
        print(ae.radius_ae(36))
        print(ae.thickness())
        print(ae.chord(36))
        print(ae.pc_set_nr(36))
        ae.add_set(radius=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
                   chord=[1.1, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1],
                   thickness=[100.0, 100.0, 90.0, 80.0, 70.0, 60.0, 50.0, 40.0, 30.0, 20.0, 10.0],
                   pc_set_id=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
        print(str(ae))


main()
