
# Input output file readers

This package contains readers and writers for typical files used in OpenFAST simulations:
- FAST input file  (`.fst, .dat, .txt`), `Class: FASTInputFile, file: fast_input_file.py`. 
- FAST output file (`.out, .outb, .elev`), `Class: FASTOutputFile, file: fast_output_file.py`. 
- FAST linearization file (`.lin`), `Class: FASTLinearizationFile, file: fast_linearization_file.py`. 
- FAST summary file (`.sum.yaml`), `Class: FASTSummaryFile, file: fast_summary_file.py`. 
- TurbSim binary file (`.bts`), `Class: TurbSimFile, file: turbsim_file.py`. 
- CSV file (`.csv, .dat, .txt`), `Class: CSVFile, file: csv_file.py`. 


## Main architecture and interface

A separate python file and class is used for each file format.
The classes inherit from the standard `File` class, present in the file `file.py`.

The object returned by each class is (or behaves as) a dictionary.
The main methods are:
- `object = class()` : create an instance of a file object
- `object = class(filename)`: create an instance of a file object, and read a given file
- `object.read(filename)`: read the given file
- `object.write(filename)`: write the object to a file (may overwrite)
- `object.toDataFrame()`: attempts to convert object to a pandas DataFrame

Additional methods may be present depending on the file format.


## Examples
Examples scripts are found in this [folder](examples). 
Below are simple examples to get started:


Read an AeroDyn file, modifies some values and write the modified file:
```python
from pyFAST.input_output import FASTInputFile
filename = 'AeroDyn.dat'
f = FASTInputFile(filename)
f['TwrAero'] = True
f['AirDens'] = 1.225
f.write('AeroDyn_Changed.dat')
```

Read an OpenFAST binary output file and convert it to a pandas DataFrame
```python
from pyFAST.input_output import FASTOutputFile
df = FASTOutputFile('5MW.outb').toDataFrame()
time  = df['Time_[s]']
Omega = df['RotSpeed_[rpm]']
```

Read a TurbSim binary file
```python 
from pyFAST.input_output import TurbSimFile
ts = TurbSimFile('Turb.bts')
print(ts.keys())
print(ts['u'].shape)  
```
For more examples on how to manipulate TurbSim files see
see [examples](examples/Example_TurbSimBox.py)

TurbSim input file processing (file modification, run and result export) 
see [examples](examples/Example_TurbSim_Processing.py).
