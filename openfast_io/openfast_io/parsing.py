"""Low-level parsing utilities for OpenFAST input files.

Extracted from FAST_reader.py to enable use by both the legacy monolithic
reader and the new per-module IO classes without circular imports.

All functions preserve their original behavior exactly.
"""
import os
import re


def readline_filterComments(f):
    """
    Filter out comments and empty lines from a file

    Args:
    f: file handle

    Returns:
    line: next line in the file that is not a comment or empty
    """
    read = True
    while read:
        line = f.readline().strip()
        if len(line)>0:
            if line[0] != '!':
                read = False
    return line

def readline_ignoreComments(f, char = '#'): # see line 64 in NWTC_IO.f90
    """
    returns line before comment character

    Args:
    f: file handle
    char: comment character

    Returns:
    line: content of next line in the file before comment character
    """

    line = f.readline().strip().split(char)

    return line[0]

def read_array(f,len,split_val=None,array_type=str):
    """
    Read an array of values from a line in a file

    Args:
    f: file handle
    len: number of values to read
    split_val: value to stop reading at
    array_type: type of values to return

    Returns:
    arr: list of values read from the file line with the specified type
    """


    strings = re.split(',| ',f.readline().strip())
    while '' in strings:    # remove empties
        strings.remove('')

    if len is None and split_val is None:
        raise Exception('Must have len or split_val to use read_array')
    
    if len is not None:
        arr = strings[:len]    # select len strings
    else:
        arr =  []
        for s in strings:
            if s != split_val:
                arr.append(s)
            else:
                break

    if array_type==str:
        arr = [ar.replace('"','') for ar in arr]  # remove quotes and commas
    elif array_type==float:
        arr = [float_read(ar) for ar in arr]
    elif array_type==int:
        arr = [int_read(ar) for ar in arr]
    elif array_type==bool:
        arr = [bool_read(ar) for ar in arr]
    else:
        raise Exception(f"read_array with type {str(array_type)} not currently supported")

    return arr

def fix_path(name):
    """ 
    split a path, then reconstruct it using os.path.join 
    
    Args:
    name: path to fix

    Returns:
    new: reconstructed path
    """
    name = re.split("\\|/", name)
    new = name[0]
    for i in range(1,len(name)):
        new = os.path.join(new, name[i])
    return new

def bool_read(text):
    """
    Read a boolean value from a string
    
    Args:
    text: string to read

    Returns:
    True if the string is 'true', False otherwise
    """
    if 'default' in text.lower():
        return str(text)
    else:
        text = text.lower()
        if text == 'true' or text == 't':
            return True
        else:
            return False

def float_read(text):
    """
    Read a float value from a string, with error handling for 'default' values

    Args:
    text: string to read

    Returns:
    float value if the string can be converted, string otherwise
    """
    if 'default' in text.lower():
        return str(text)
    else:
        try:
            return float(text)
        except:
            return str(text)

def int_read(text):
    """
    Read an integer value from a string, with error handling for 'default' values

    Args:
    text: string to read

    Returns:
    int value if the string can be converted, string otherwise
    """
    if 'default' in text.lower():
        return str(text)
    else:
        try:
            return int(text)
        except:
            return str(text)

def quoted_read(text):
    """
    Read a quoted value from a string (i.e. a value between quotes)

    Args:
    text: string to read

    Returns:
    quoted value if the string is quoted, unquoted value otherwise

    """
    if '"' in text:
        return text.split('"')[1]
    elif "'" in text:
        return text.split("'")[1]
    else:
        return text


def fmt_field(val, min_width: int = 28) -> str:
    """Format a single value field for OpenFAST input file lines.

    Returns a left-justified string that is guaranteed to have at least 2
    trailing spaces regardless of how long the string representation of *val*
    is.  This prevents the common bug where a value whose string exactly fills
    the format width concatenates with the parameter name on the right.

    Usage::

        f.write(fmt_field(dvr['SomeVec']) + 'ParamName  - description\\n')

    Args:
        val:       The value to format.  May be any type — it is coerced to
                   ``str``.
        min_width: Minimum field width (not including the mandatory 2-space
                   separator).  Defaults to 28.

    Returns:
        ``str(val).ljust(max(len(str(val)) + 2, min_width))``
    """
    s = str(val)
    width = max(len(s) + 2, min_width)
    return s.ljust(width)
