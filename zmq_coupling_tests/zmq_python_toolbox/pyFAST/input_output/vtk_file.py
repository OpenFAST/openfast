""" 
Read/Write VTK files

Part of weio library: https://github.com/ebranlard/weio

"""
import pandas as pd
import numpy as np
import numpy
import os
from functools import reduce
import collections

try:
    from .file import File, EmptyFileError, WrongFormatError, BrokenFormatError
except ImportError:
    EmptyFileError   = type('EmptyFileError', (Exception,),{})
    WrongFormatError =  type('WrongFormatError', (Exception,),{})
    BrokenFormatError =  type('BrokenFormatError', (Exception,),{})
    File=dict

class VTKFile(File):
    """ 
    Read/write a VTK file (.vtk). 

    Main attributes for grids:
    ---------
    - xp_grid, yp_grid, zp_grid: vectors of points locations
    - point_data_grid: dictionary containing data at the grid points

    Main attributes for mesh:
    ---------
    - points
    - point_data
    - cells
    - cell_data

    Main methods
    ------------
    - read, write

    Examples
    --------
        vtk = VTKFile('DisXZ1.vtk')
        x  = vtk.x_grid
        z  = vtk.z_grid
        Ux = vtk.point_data_grid['DisXZ'][:,0,:,0]
    
    """
    @staticmethod
    def defaultExtensions():
        return ['.vtk','.vtp']

    @staticmethod
    def formatName():
        return 'VTK file'

    def __init__(self,filename=None,**kwargs):
        self.filename = None
        # For regular grid
        self.xp_grid=None  # location of points
        self.yp_grid=None
        self.zp_grid=None
        self.point_data_grid = None

        # Main Data
        self.points = None  
        self.field_data = {}
        self.point_data = {} 
        self.dataset = {} 



        # Data for reading only
        self.cell_data_raw = {}
        self.c = None     # Cell
        self.ct = None    # CellTypes
        self.active = None
        self.is_ascii = False
        self.split = []
        self.num_items = 0
        self.section = None

        # Propagate read
        if filename:
            self.read(filename=filename,**kwargs)


    def read(self, filename=None):
        """ read a VTK file """
        if filename:
            self.filename = filename
        if not self.filename:
            raise Exception('No filename provided')
        if not os.path.isfile(self.filename):
            raise OSError(2,'File not found:',self.filename)
        if os.stat(self.filename).st_size == 0:
            raise EmptyFileError('File is empty:',self.filename)

        with open(filename, "rb") as f:
            # initialize output data
            # skip header and title
            f.readline()
            f.readline()

            data_type = f.readline().decode("utf-8").strip().upper()
            if data_type not in ["ASCII", "BINARY"]:
                raise ReadError('Unknown VTK data type ',data_type)
            self.is_ascii = data_type == "ASCII"

            while True:
                line = f.readline().decode("utf-8")
                if not line:
                    # EOF
                    break

                line = line.strip()
                if len(line) == 0:
                    continue

                self.split = line.split()
                self.section = self.split[0].upper()

                if self.section in vtk_sections:
                    _read_section(f, self)
                else:
                    _read_subsection(f, self)

        # --- Postpro
        _check_mesh(self) # generate points if needed
        cells, cell_data = translate_cells(self.c, self.ct, self.cell_data_raw)
        self.cells     = cells
        self.cell_data = cell_data

        if self.dataset['type']=='STRUCTURED_POINTS':
            self.point_data_grid = {}
            # We provide point_data_grid, corresponds to point_data but reshaped
            for k,PD in self.point_data.items():
                # NOTE: tested foe len(y)=1, len(z)=1
                self.point_data_grid[k]=PD.reshape(len(self.xp_grid), len(self.yp_grid), len(self.zp_grid),PD.shape[1], order='F')


    def write(self, filename=None, binary=True):
        """ 
        Write to unstructured grid
        TODO structured grid
        """

        if filename:
            self.filename = filename
        if not self.filename:
            raise Exception('No filename provided')


        def pad(array):
            return np.pad(array, ((0, 0), (0, 1)), "constant")

        if self.points.shape[1] == 2:
            points = pad(self.points)
        else:
            points = self.points

        if self.point_data:
            for name, values in self.point_data.items():
                if len(values.shape) == 2 and values.shape[1] == 2:
                    self.point_data[name] = pad(values)

        for name, data in self.cell_data.items():
            for k, values in enumerate(data):
                if len(values.shape) == 2 and values.shape[1] == 2:
                    data[k] = pad(data[k])

        with open(filename, "wb") as f:
            f.write(b"# vtk DataFile Version 4.2\n")
            f.write("written \n".encode("utf-8"))
            f.write(("BINARY\n" if binary else "ASCII\n").encode("utf-8"))
            f.write(b"DATASET UNSTRUCTURED_GRID\n")

            # write points and cells
            _write_points(f, points, binary)
            _write_cells(f, self.cells, binary)

            # write point data
            if self.point_data:
                num_points = self.points.shape[0]
                f.write("POINT_DATA {}\n".format(num_points).encode("utf-8"))
                _write_field_data(f, self.point_data, binary)

            # write cell data
            if self.cell_data:
                total_num_cells = sum(len(c.data) for c in self.cells)
                f.write("CELL_DATA {}\n".format(total_num_cells).encode("utf-8"))
                _write_field_data(f, self.cell_data, binary)


    def __repr__(self):
        s='<{} object> with keys:\n'.format(type(self).__name__)
        for k,v in self.items():
            s+=' - {}: {}\n'.format(k,v)
        return s


    def __repr__(self):
        """ print function """
        def show_grid(v,s):
            if v is None:
                return
            if len(v)==0:
                return
            if len(v)==1:
                lines.append('- {}: [{}], n: {}'.format(s,v[0],len(v)))
            else:
                lines.append('- {}: [{} ... {}],  dx: {}, n: {}'.format(s,v[0],v[-1],v[1]-v[0],len(v)))

        lines = ['<{} object> with attributes:'.format(type(self).__name__)]
        show_grid(self.xp_grid, 'xp_grid')
        show_grid(self.yp_grid, 'yp_grid')
        show_grid(self.zp_grid, 'zp_grid')

        if self.point_data_grid:
            lines.append('- point_data_grid:')
            for k,v in self.point_data_grid.items():
                lines.append('  "{}" : {}'.format(k,v.shape))

        lines.append('- points {}'.format(len(self.points)))
        if len(self.cells) > 0:
            lines.append("- cells:")
            for tpe, elems in self.cells:
                lines.append("    {}: {}".format(tpe,len(elems)))
        else:
            lines.append("  No cells.")

        if self.point_data:
            lines.append('- point_data:')
            for k,v in self.point_data.items():
                lines.append('  "{}" : {}'.format(k,v.shape))

        if self.cell_data:
            names = ", ".join(self.cell_data.keys())
            lines.append("  Cell data: {}".format(names))

        return "\n".join(lines)


    def toDataFrame(self):
        return None
        


#         Save FlowData Object to vtk
#         """
#         n_points = self.dimensions.x1 * self.dimensions.x2 * self.dimensions.x3
#         vtk_file = Output(filename)
          #self.file = open(self.filename, "w")
          #self.ln = "\n"
#         vtk_file.write_line('# vtk DataFile Version 3.0')
#         vtk_file.write_line('array.mean0D')
#         vtk_file.write_line('ASCII')
#         vtk_file.write_line('DATASET STRUCTURED_POINTS')
#         vtk_file.write_line('DIMENSIONS {}'.format(self.dimensions))
#         vtk_file.write_line('ORIGIN {}'.format(self.origin))
#         vtk_file.write_line('SPACING {}'.format(self.spacing))
#         vtk_file.write_line('POINT_DATA {}'.format(n_points))
#         vtk_file.write_line('FIELD attributes 1')
#         vtk_file.write_line('UAvg 3 {} float'.format(n_points))
#         for u, v, w in zip(self.u, self.v, self.w):
#             vtk_file.write_line('{}'.format(Vec3(u, v, w)))

# --- Paraview
# except: from paraview.simple import *
#    sliceFile = sliceDir + '/' + tStartStr + '/' + sliceName
#    print '     Slice file 1a: ' + sliceFile
#    slice_1a_vtk = LegacyVTKReader( FileNames=[sliceFile] )
#    sliceFile = sliceDir + '/' + tEndStr + '/' + sliceName
#    DataRepresentation3 = GetDisplayProperties(slice_1a_vtk)
#    DataRepresentation3.Visibility = 0
#    SetActiveSource(slice_1a_vtk)

# --- VTK
# import np
# from vtk import vtkStructuredPointsReader
# from vtk.util import np as VN
# 
# reader = vtkStructuredPointsReader()
# reader.SetFileName(filename)
# reader.ReadAllVectorsOn()
# reader.ReadAllScalarsOn()
# reader.Update()
# 
# data = reader.GetOutput()
# 
# dim = data.GetDimensions()
# vec = list(dim)
# vec = [i-1 for i in dim]
# vec.append(3)
# 
# u = VN.vtk_to_np(data.GetCellData().GetArray('velocity'))
# b = VN.vtk_to_numpy(data.GetCellData().GetArray('cell_centered_B'))
# 
# u = u.reshape(vec,order='F')
# b = b.reshape(vec,order='F')
# 
# x = zeros(data.GetNumberOfPoints())
# y = zeros(data.GetNumberOfPoints())
# z = zeros(data.GetNumberOfPoints())
# 
# for i in range(data.GetNumberOfPoints()):
#         x[i],y[i],z[i] = data.GetPoint(i)
# 
# x = x.reshape(dim,order='F')
# y = y.reshape(dim,order='F')
# z = z.reshape(dim,order='F')

# --- vtk
# import vtk
# import numpy
# import vtk_demo.version as version
# 
# 
# def main():
#     """
#     :return: The render window interactor.
#     """
# 
#     chessboard_resolution = 5
#     n_lut_colors = 256
#     data_max_value = 1
# 
#     # Provide some geometry
#     chessboard = vtk.vtkPlaneSource()
#     chessboard.SetXResolution(chessboard_resolution)
#     chessboard.SetYResolution(chessboard_resolution)
#     num_squares = chessboard_resolution * chessboard_resolution
#     #  Force an update so we can set cell data
#     chessboard.Update()
# 
#     # Make some arbitrary data to show on the chessboard geometry
#     data = vtk.vtkFloatArray()
#     for i in range(num_squares):
#         if i == 4:
#             # This square should in principle light up with color given by SetNanColor below
#             data.InsertNextTuple1(numpy.nan)
#         else:
#             thing = (i * data_max_value) / (num_squares - 1)
#             data.InsertNextTuple1(thing)
# 
#     # Make a LookupTable
#     lut = vtk.vtkLookupTable()
#     lut.SetNumberOfColors(n_lut_colors)
#     lut.Build()
#     lut.SetTableRange(0, data_max_value)
#     lut.SetNanColor(.1, .5, .99, 1.0) # <------ This color gets used
#     for i in range(n_lut_colors):
#         # Fill it with arbitrary colors, e.g. grayscale
#         x = data_max_value*i/(n_lut_colors-1)
#         lut.SetTableValue(i, x, x, x, 1.0)
#     lut.SetNanColor(.99, .99, .1, 1.0) # <----- This color gets ignored! ...except by GetNanColor
# 
#     print(lut.GetNanColor())  # <-- Prints the color set by the last SetNanColor call above!
# 
#     chessboard.GetOutput().GetCellData().SetScalars(data)
# 
#     mapper = vtk.vtkPolyDataMapper()
#     mapper.SetInputConnection(chessboard.GetOutputPort())
#     mapper.SetScalarRange(0, data_max_value)
#     mapper.SetLookupTable(lut)
#     mapper.Update()
# 
#     actor = vtk.vtkActor()
#     actor.SetMapper(mapper)
# 
#     renderer = vtk.vtkRenderer()
#     ren_win = vtk.vtkRenderWindow()
#     ren_win.AddRenderer(renderer)
#     renderer.SetBackground(vtk.vtkNamedColors().GetColor3d('MidnightBlue'))
#     renderer.AddActor(actor)
# 
#     iren = vtk.vtkRenderWindowInteractor()

# --- pyvtk
#         """Read vtk-file stored previously with tovtk."""
#         p = pyvtk.VtkData(filename)
#         xn = array(p.structure.points)
#         dims = p.structure.dimensions
#         try:
#             N = eval(p.header.split(" ")[-1])

    # --- Extract
    #     # Convert the center of the turbine coordinate in m to High-Resolution domains left most corner in (i,j)
    #     xe_index = int(origin_at_precusr[0]/10)
    #     ye_index = int(origin_at_precusr[1]/10)
    #     
    #     # Read the full domain from VTK
    #     reader = vtk.vtkStructuredPointsReader()
    #     reader.SetFileName(in_vtk)
    #     reader.Update()
    #     
    #     # Extract the High Resolution domain at same spacial spacing by specifying the (i,i+14),(j,j+14),(k,k+20) tuples
    #     extract = vtk.vtkExtractVOI()
    #     extract.SetInputConnection(reader.GetOutputPort())
    #     extract.SetVOI(xe_index, xe_index+14, ye_index, ye_index+14, 0, 26)
    #     extract.SetSampleRate(1, 1, 1)
    #     extract.Update()
    #     
    #     # Write the extract as VTK
    #     points = extract.GetOutput()
    #     vec = points.GetPointData().GetVectors('Amb')
    # 
    #     with open(out_vtk, 'a') as the_file:
    #         the_file.write('# vtk DataFile Version 3.0\n')
    #         the_file.write('High\n')
    #         the_file.write('ASCII\n')
    #         the_file.write('DATASET STRUCTURED_POINTS\n')
    #         the_file.write('DIMENSIONS %d %d %d\n' % points.GetDimensions())
    #         the_file.write('ORIGIN %f %f %f\n' % origin_at_stitch)
    #         the_file.write('SPACING %f %f %f\n' % points.GetSpacing())
    #         the_file.write('POINT_DATA %d\n' % points.GetNumberOfPoints())
    #         the_file.write('VECTORS Amb float\n')
    #         for i in range(points.GetNumberOfPoints()):
    #             the_file.write('%f %f %f\n' % vec.GetTuple(i) )

    # --- Stitch
    #     reader = vtk.vtkStructuredPointsReader()
    #     reader.SetFileName(in_vtk)
    #     reader.Update()
    #     
    #     hAppend = vtk.vtkImageAppend()
    #     hAppend.SetAppendAxis(0)
    #     for i in range(nx):
    #         hAppend.AddInputData(reader.GetOutput())
    #     hAppend.Update()
    #     
    #     vAppend = vtk.vtkImageAppend()
    #     vAppend.SetAppendAxis(1)
    #     for i in range(ny):
    #         vAppend.AddInputData(hAppend.GetOutput())
    #     vAppend.Update()
    #     
    #     points = vAppend.GetOutput()
    #     vec = points.GetPointData().GetVectors('Amb')
    # 
    #     with open(out_vtk, 'a') as the_file:
    #         the_file.write('# vtk DataFile Version 3.0\n')
    #         the_file.write('Low\n')
    #         the_file.write('ASCII\n')
    #         the_file.write('DATASET STRUCTURED_POINTS\n')
    #         the_file.write('DIMENSIONS %d %d %d\n' % points.GetDimensions())
    #         the_file.write('ORIGIN %f %f %f\n' % points.GetOrigin())
    #         the_file.write('SPACING %f %f %f\n' % points.GetSpacing())
    #         the_file.write('POINT_DATA %d\n' % points.GetNumberOfPoints())
    #         the_file.write('VECTORS Amb float\n')
    #         for i in range(points.GetNumberOfPoints()):
    #             the_file.write('%f %f %f\n' % vec.GetTuple(i) )


    # 


# --------------------------------------------------------------------------------}
# --- The code below is taken from meshio 
#     https://github.com/nschloe/meshio
#     The MIT License (MIT)
#     Copyright (c) 2015-2020 meshio developers
# --------------------------------------------------------------------------------{
ReadError  = BrokenFormatError
WriteError = BrokenFormatError

def _vtk_to_meshio_order(vtk_type, numnodes, dtype=int):
    # meshio uses the same node ordering as VTK for most cell types. However, for the
    # linear wedge, the ordering of the gmsh Prism [1] is adopted since this is found in
    # most codes (Abaqus, Ansys, Nastran,...). In the vtkWedge [2], the normal of the
    # (0,1,2) triangle points outwards, while in gmsh this normal points inwards.
    # [1] http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
    # [2] https://vtk.org/doc/nightly/html/classvtkWedge.html
    if vtk_type == 13:
        return numpy.array([0, 2, 1, 3, 5, 4], dtype=dtype)
    else:
        return numpy.arange(0, numnodes, dtype=dtype)

def _meshio_to_vtk_order(meshio_type, numnodes, dtype=int):
    if meshio_type == "wedge":
        return numpy.array([0, 2, 1, 3, 5, 4], dtype=dtype)
    else:
        return numpy.arange(0, numnodes, dtype=dtype)

vtk_to_meshio_type = {
    0: "empty",
    1: "vertex",
    # 2: 'poly_vertex',
    3: "line",
    # 4: 'poly_line',
    5: "triangle",
    # 6: 'triangle_strip',
    7: "polygon",
    8: "pixel",
    9: "quad",
    10: "tetra",
    # 11: 'voxel',
    12: "hexahedron",
    13: "wedge",
    14: "pyramid",
    15: "penta_prism",
    16: "hexa_prism",
    21: "line3",
    22: "triangle6",
    23: "quad8",
    24: "tetra10",
    25: "hexahedron20",
    26: "wedge15",
    27: "pyramid13",
    28: "quad9",
    29: "hexahedron27",
    30: "quad6",
    31: "wedge12",
    32: "wedge18",
    33: "hexahedron24",
    34: "triangle7",
    35: "line4",
    42: "polyhedron",
    #
    # 60: VTK_HIGHER_ORDER_EDGE,
    # 61: VTK_HIGHER_ORDER_TRIANGLE,
    # 62: VTK_HIGHER_ORDER_QUAD,
    # 63: VTK_HIGHER_ORDER_POLYGON,
    # 64: VTK_HIGHER_ORDER_TETRAHEDRON,
    # 65: VTK_HIGHER_ORDER_WEDGE,
    # 66: VTK_HIGHER_ORDER_PYRAMID,
    # 67: VTK_HIGHER_ORDER_HEXAHEDRON,
    # Arbitrary order Lagrange elements
    68: "VTK_LAGRANGE_CURVE",
    69: "VTK_LAGRANGE_TRIANGLE",
    70: "VTK_LAGRANGE_QUADRILATERAL",
    71: "VTK_LAGRANGE_TETRAHEDRON",
    72: "VTK_LAGRANGE_HEXAHEDRON",
    73: "VTK_LAGRANGE_WEDGE",
    74: "VTK_LAGRANGE_PYRAMID",
    # Arbitrary order Bezier elements
    75: "VTK_BEZIER_CURVE",
    76: "VTK_BEZIER_TRIANGLE",
    77: "VTK_BEZIER_QUADRILATERAL",
    78: "VTK_BEZIER_TETRAHEDRON",
    79: "VTK_BEZIER_HEXAHEDRON",
    80: "VTK_BEZIER_WEDGE",
    81: "VTK_BEZIER_PYRAMID",
}
meshio_to_vtk_type = {v: k for k, v in vtk_to_meshio_type.items()}


# --------------------------------------------------------------------------------}
# --- Mesh 
# --------------------------------------------------------------------------------{
class CellBlock(collections.namedtuple("CellBlock", ["type", "data"])):
    def __repr__(self):
        return "<meshio CellBlock, type: {}, num cells: {}>".format(self.type,len(self.data))
# --------------------------------------------------------------------------------}
# --- File _vtk.py  from meshio
# --------------------------------------------------------------------------------{
vtk_type_to_numnodes = numpy.array(
    [
        0,  # empty
        1,  # vertex
        -1,  # poly_vertex
        2,  # line
        -1,  # poly_line
        3,  # triangle
        -1,  # triangle_strip
        -1,  # polygon
        -1,  # pixel
        4,  # quad
        4,  # tetra
        -1,  # voxel
        8,  # hexahedron
        6,  # wedge
        5,  # pyramid
        10,  # penta_prism
        12,  # hexa_prism
        -1,
        -1,
        -1,
        -1,
        3,  # line3
        6,  # triangle6
        8,  # quad8
        10,  # tetra10
        20,  # hexahedron20
        15,  # wedge15
        13,  # pyramid13
        9,  # quad9
        27,  # hexahedron27
        6,  # quad6
        12,  # wedge12
        18,  # wedge18
        24,  # hexahedron24
        7,  # triangle7
        4,  # line4
    ]
)


# These are all VTK data types.
# One sometimes finds 'vtktypeint64', but this is ill-formed.
vtk_to_numpy_dtype_name = {
    "bit": "bool",
    "unsigned_char": "uint8",
    "char": "int8",
    "unsigned_short": "uint16",
    "short": "int16",
    "unsigned_int": "uint32",
    "int": "int32",
    "unsigned_long": "uint64",
    "long": "int64",
    "float": "float32",
    "double": "float64",
    "vtktypeint32": "int32",  # vtk DataFile Version 5.1
    "vtktypeint64": "int64",  # vtk DataFile Version 5.1
    "vtkidtype": "int32",  # may be either 32-bit or 64-bit (VTK_USE_64BIT_IDS)
}

numpy_to_vtk_dtype = {
    v: k for k, v in vtk_to_numpy_dtype_name.items() if "vtk" not in k
}

# supported vtk dataset types
vtk_dataset_types = [
    "UNSTRUCTURED_GRID",
    "STRUCTURED_POINTS",
    "STRUCTURED_GRID",
    "RECTILINEAR_GRID",
]
# additional infos per dataset type
vtk_dataset_infos = {
    "UNSTRUCTURED_GRID": [],
    "STRUCTURED_POINTS": [
        "DIMENSIONS",
        "ORIGIN",
        "SPACING",
        "ASPECT_RATIO",  # alternative for SPACING in version 1.0 and 2.0
    ],
    "STRUCTURED_GRID": ["DIMENSIONS"],
    "RECTILINEAR_GRID": [
        "DIMENSIONS",
        "X_COORDINATES",
        "Y_COORDINATES",
        "Z_COORDINATES",
    ],
}

# all main sections in vtk
vtk_sections = [
    "METADATA",
    "DATASET",
    "POINTS",
    "CELLS",
    "CELL_TYPES",
    "POINT_DATA",
    "CELL_DATA",
    "LOOKUP_TABLE",
    "COLOR_SCALARS",
]





def read(filename):
    """Reads a VTK vtk file."""
    with open(filename, "rb") as f:
        out = read_buffer(f)
    return out


def read_buffer(f):
    # initialize output data
    info = VTKFile()

    # skip header and title
    f.readline()
    f.readline()

    data_type = f.readline().decode("utf-8").strip().upper()
    if data_type not in ["ASCII", "BINARY"]:
        raise WrongFormatError("Unknown VTK data type ",data_type)
    info.is_ascii = data_type == "ASCII"

    while True:
        line = f.readline().decode("utf-8")
        if not line:
            # EOF
            break

        line = line.strip()
        if len(line) == 0:
            continue

        info.split = line.split()
        info.section = info.split[0].upper()

        if info.section in vtk_sections:
            _read_section(f, info)
        else:
            _read_subsection(f, info)

    _check_mesh(info)
    cells, cell_data = translate_cells(info.c, info.ct, info.cell_data_raw)

    info.cells     = cells
    info.cell_data = cell_data

    return info


def _read_section(f, info):
    if info.section == "METADATA":
        _skip_meta(f)

    elif info.section == "DATASET":
        info.active = "DATASET"
        info.dataset["type"] = info.split[1].upper()
        if info.dataset["type"] not in vtk_dataset_types:
            raise BrokenFormatError(
                "Only VTK '{}' supported (not {}).".format(
                    "', '".join(vtk_dataset_types), info.dataset["type"]
                )
            )

    elif info.section == "POINTS":
        info.active = "POINTS"
        info.num_points = int(info.split[1])
        data_type = info.split[2].lower()
        info.points = _read_points(f, data_type, info.is_ascii, info.num_points)

    elif info.section == "CELLS":
        info.active = "CELLS"
        last_pos = f.tell()
        try:
            line = f.readline().decode("utf-8")
        except UnicodeDecodeError:
            line = ""
        if "OFFSETS" in line:
            # vtk DataFile Version 5.1 - appearing in Paraview 5.8.1 outputs
            # No specification found for this file format.
            # See the question on ParaView Discourse Forum:
            # <https://discourse.paraview.org/t/specification-of-vtk-datafile-version-5-1/5127>.
            info.num_offsets = int(info.split[1])
            info.num_items = int(info.split[2])
            dtype = numpy.dtype(vtk_to_numpy_dtype_name[line.split()[1]])
            offsets = _read_cells(f, info.is_ascii, info.num_offsets, dtype)
            line = f.readline().decode("utf-8")
            assert "CONNECTIVITY" in line
            dtype = numpy.dtype(vtk_to_numpy_dtype_name[line.split()[1]])
            connectivity = _read_cells(f, info.is_ascii, info.num_items, dtype)
            info.c = (offsets, connectivity)
        else:
            f.seek(last_pos)
            info.num_items = int(info.split[2])
            info.c = _read_cells(f, info.is_ascii, info.num_items)

    elif info.section == "CELL_TYPES":
        info.active = "CELL_TYPES"
        info.num_items = int(info.split[1])
        info.ct = _read_cell_types(f, info.is_ascii, info.num_items)

    elif info.section == "POINT_DATA":
        info.active = "POINT_DATA"
        info.num_items = int(info.split[1])

    elif info.section == "CELL_DATA":
        info.active = "CELL_DATA"
        info.num_items = int(info.split[1])

    elif info.section == "LOOKUP_TABLE":
        info.num_items = int(info.split[2])
        numpy.fromfile(f, count=info.num_items * 4, sep=" ", dtype=float)
        # rgba = data.reshape((info.num_items, 4))

    elif info.section == "COLOR_SCALARS":
        nValues = int(info.split[2])
        # re-use num_items from active POINT/CELL_DATA
        num_items = info.num_items
        dtype = numpy.ubyte
        if info.is_ascii:
            dtype = float
        numpy.fromfile(f, count=num_items * nValues, dtype=dtype)


def _read_subsection(f, info):
    if info.active == "POINT_DATA":
        d = info.point_data
    elif info.active == "CELL_DATA":
        d = info.cell_data_raw
    elif info.active == "DATASET":
        d = info.dataset
    else:
        d = info.field_data

    if info.section in vtk_dataset_infos[info.dataset["type"]]:
        if info.section[1:] == "_COORDINATES":
            info.num_points = int(info.split[1])
            data_type = info.split[2].lower()
            d[info.section] = _read_coords(f, data_type, info.is_ascii, info.num_points)
        else:
            if info.section == "DIMENSIONS":
                d[info.section] = list(map(int, info.split[1:]))
            else:
                d[info.section] = list(map(float, info.split[1:]))
            if len(d[info.section]) != 3:
                raise BrokenFormatError(
                    "Wrong number of info in section '{}'. Need 3, got {}.".format(
                        info.section, len(d[info.section])
                    )
                )
    elif info.section == "SCALARS":
        d.update(_read_scalar_field(f, info.num_items, info.split, info.is_ascii))
    elif info.section == "VECTORS":
        d.update(_read_field(f, info.num_items, info.split, [3], info.is_ascii))
    elif info.section == "TENSORS":
        d.update(_read_field(f, info.num_items, info.split, [3, 3], info.is_ascii))
    elif info.section == "FIELD":
        d.update(_read_fields(f, int(info.split[2]), info.is_ascii))
    else:
        raise WrongFormatError("Unknown section ",info.section)


def _check_mesh(info):
    if info.dataset["type"] == "UNSTRUCTURED_GRID":
        if info.c is None:
            raise ReadError("Required section CELLS not found.")
        if info.ct is None:
            raise ReadError("Required section CELL_TYPES not found.")

    elif info.dataset["type"] == "STRUCTURED_POINTS":
        dim = info.dataset["DIMENSIONS"]
        ori = info.dataset["ORIGIN"]
        spa = (
            info.dataset["SPACING"]
            if "SPACING" in info.dataset
            else info.dataset["ASPECT_RATIO"]
        )
        axis = [
            numpy.linspace(ori[i], ori[i] + (dim[i] - 1.0) * spa[i], dim[i])
            for i in range(3)
        ]
        info.xp_grid=axis[0]
        info.yp_grid=axis[1]
        info.zp_grid=axis[2]

        info.points = _generate_points(axis)
        info.c, info.ct = _generate_cells(dim=info.dataset["DIMENSIONS"])

    elif info.dataset["type"] == "RECTILINEAR_GRID":
        axis = [
            info.dataset["X_COORDINATES"],
            info.dataset["Y_COORDINATES"],
            info.dataset["Z_COORDINATES"],
        ]
        info.xp_grid=axis[0]
        info.yp_grid=axis[1]
        info.zp_grid=axis[2]

        info.points = _generate_points(axis)
        info.c, info.ct = _generate_cells(dim=info.dataset["DIMENSIONS"])

    elif info.dataset["type"] == "STRUCTURED_GRID":
        info.c, info.ct = _generate_cells(dim=info.dataset["DIMENSIONS"])
        # TODO x_grid, y_grid, z_grid points



def _generate_cells(dim):
    ele_dim = [d - 1 for d in dim if d > 1]
    ele_no = numpy.prod(ele_dim, dtype=int)
    spatial_dim = len(ele_dim)

    if spatial_dim == 1:
        # cells are lines in 1D
        cells = numpy.empty((ele_no, 3), dtype=int)
        cells[:, 0] = 2
        cells[:, 1] = numpy.arange(ele_no, dtype=int)
        cells[:, 2] = cells[:, 1] + 1
        cell_types = numpy.full(ele_no, 3, dtype=int)

    elif spatial_dim == 2:
        # cells are quad in 2D
        cells = numpy.empty((ele_no, 5), dtype=int)
        cells[:, 0] = 4
        cells[:, 1] = numpy.arange(0, ele_no, dtype=int)
        cells[:, 1] += numpy.arange(0, ele_no, dtype=int) // ele_dim[0]
        cells[:, 2] = cells[:, 1] + 1
        cells[:, 3] = cells[:, 1] + 2 + ele_dim[0]
        cells[:, 4] = cells[:, 3] - 1
        cell_types = numpy.full(ele_no, 9, dtype=int)
    else:
        # cells are hex in 3D
        cells = numpy.empty((ele_no, 9), dtype=int)
        cells[:, 0] = 8
        cells[:, 1] = numpy.arange(ele_no)
        cells[:, 1] += (ele_dim[0] + ele_dim[1] + 1) * (
            numpy.arange(ele_no) // (ele_dim[0] * ele_dim[1])
        )
        cells[:, 1] += (numpy.arange(ele_no) % (ele_dim[0] * ele_dim[1])) // ele_dim[0]
        cells[:, 2] = cells[:, 1] + 1
        cells[:, 3] = cells[:, 1] + 2 + ele_dim[0]
        cells[:, 4] = cells[:, 3] - 1
        cells[:, 5] = cells[:, 1] + (1 + ele_dim[0]) * (1 + ele_dim[1])
        cells[:, 6] = cells[:, 5] + 1
        cells[:, 7] = cells[:, 5] + 2 + ele_dim[0]
        cells[:, 8] = cells[:, 7] - 1
        cell_types = numpy.full(ele_no, 12, dtype=int)

    return cells.reshape(-1), cell_types

def _generate_points(axis):
    x_dim = len(axis[0])
    y_dim = len(axis[1])
    z_dim = len(axis[2])
    pnt_no = x_dim * y_dim * z_dim
    x_id, y_id, z_id = numpy.mgrid[0:x_dim, 0:y_dim, 0:z_dim]
    points = numpy.empty((pnt_no, 3), dtype=axis[0].dtype)
    # VTK sorts points and cells in Fortran order
    points[:, 0] = axis[0][x_id.reshape(-1, order="F")]
    points[:, 1] = axis[1][y_id.reshape(-1, order="F")]
    points[:, 2] = axis[2][z_id.reshape(-1, order="F")]
    return points


def _read_coords(f, data_type, is_ascii, num_points):
    dtype = numpy.dtype(vtk_to_numpy_dtype_name[data_type])
    if is_ascii:
        coords = numpy.fromfile(f, count=num_points, sep=" ", dtype=dtype)
    else:
        # Binary data is big endian, see
        # <https://www.vtk.org/Wiki/VTK/Writing_VTK_files_using_python#.22legacy.22>.
        dtype = dtype.newbyteorder(">")
        coords = numpy.fromfile(f, count=num_points, dtype=dtype)
        line = f.readline().decode("utf-8")
        if line != "\n":
            raise ReadError()
    return coords


def _read_points(f, data_type, is_ascii, num_points):
    dtype = numpy.dtype(vtk_to_numpy_dtype_name[data_type])
    if is_ascii:
        points = numpy.fromfile(f, count=num_points * 3, sep=" ", dtype=dtype)
    else:
        # Binary data is big endian, see
        # <https://www.vtk.org/Wiki/VTK/Writing_VTK_files_using_python#.22legacy.22>.
        dtype = dtype.newbyteorder(">")
        points = numpy.fromfile(f, count=num_points * 3, dtype=dtype)
        line = f.readline().decode("utf-8")
        if line != "\n":
            raise ReadError()
    return points.reshape((num_points, 3))


def _read_cells(f, is_ascii, num_items, dtype=numpy.dtype("int32")):
    if is_ascii:
        c = numpy.fromfile(f, count=num_items, sep=" ", dtype=dtype)
    else:
        dtype = dtype.newbyteorder(">")
        c = numpy.fromfile(f, count=num_items, dtype=dtype)
        line = f.readline().decode("utf-8")
        if line != "\n":
            raise ReadError()
    return c


def _read_cell_types(f, is_ascii, num_items):
    if is_ascii:
        ct = numpy.fromfile(f, count=int(num_items), sep=" ", dtype=int)
    else:
        # binary
        ct = numpy.fromfile(f, count=int(num_items), dtype=">i4")
        line = f.readline().decode("utf-8")
        # Sometimes, there's no newline at the end
        if line.strip() != "":
            raise ReadError()
    return ct


def _read_scalar_field(f, num_data, split, is_ascii):
    data_name = split[1]
    data_type = split[2].lower()
    try:
        num_comp = int(split[3])
    except IndexError:
        num_comp = 1

    # The standard says:
    # > The parameter numComp must range between (1,4) inclusive; [...]
    if not (0 < num_comp < 5):
        raise ReadError("The parameter numComp must range between (1,4) inclusive")

    dtype = numpy.dtype(vtk_to_numpy_dtype_name[data_type])
    lt, _ = f.readline().decode("utf-8").split()
    if lt.upper() != "LOOKUP_TABLE":
        raise ReadError()

    if is_ascii:
        data = numpy.fromfile(f, count=num_data * num_comp, sep=" ", dtype=dtype)
    else:
        # Binary data is big endian, see
        # <https://www.vtk.org/Wiki/VTK/Writing_VTK_files_using_python#.22legacy.22>.
        dtype = dtype.newbyteorder(">")
        data = numpy.fromfile(f, count=num_data * num_comp, dtype=dtype)
        line = f.readline().decode("utf-8")
        if line != "\n":
            raise ReadError()

    data = data.reshape(-1, num_comp)
    return {data_name: data}


def _read_field(f, num_data, split, shape, is_ascii):
    data_name = split[1]
    data_type = split[2].lower()

    dtype = numpy.dtype(vtk_to_numpy_dtype_name[data_type])
    # prod()
    # <https://stackoverflow.com/q/2104782/353337>
    k = reduce((lambda x, y: x * y), shape)

    if is_ascii:
        data = numpy.fromfile(f, count=k * num_data, sep=" ", dtype=dtype)
    else:
        # Binary data is big endian, see
        # <https://www.vtk.org/Wiki/VTK/Writing_VTK_files_using_python#.22legacy.22>.
        dtype = dtype.newbyteorder(">")
        data = numpy.fromfile(f, count=k * num_data, dtype=dtype)
        line = f.readline().decode("utf-8")
        if line != "\n":
            raise ReadError()

    data = data.reshape(-1, *shape)
    return {data_name: data}


def _read_fields(f, num_fields, is_ascii):
    data = {}
    for _ in range(num_fields):
        line = f.readline().decode("utf-8").split()
        if line[0] == "METADATA":
            _skip_meta(f)
            name, shape0, shape1, data_type = f.readline().decode("utf-8").split()
        else:
            name, shape0, shape1, data_type = line

        shape0 = int(shape0)
        shape1 = int(shape1)
        dtype = numpy.dtype(vtk_to_numpy_dtype_name[data_type.lower()])

        if is_ascii:
            dat = numpy.fromfile(f, count=shape0 * shape1, sep=" ", dtype=dtype)
        else:
            # Binary data is big endian, see
            # <https://www.vtk.org/Wiki/VTK/Writing_VTK_files_using_python#.22legacy.22>.
            dtype = dtype.newbyteorder(">")
            dat = numpy.fromfile(f, count=shape0 * shape1, dtype=dtype)
            line = f.readline().decode("utf-8")
            if line != "\n":
                raise ReadError()

        if shape0 != 1:
            dat = dat.reshape((shape1, shape0))

        data[name] = dat

    return data


def _skip_meta(f):
    # skip possible metadata
    # https://vtk.org/doc/nightly/html/IOLegacyInformationFormat.html
    while True:
        line = f.readline().decode("utf-8").strip()
        if not line:
            # end of metadata is a blank line
            break


def translate_cells(data, types, cell_data_raw):
    # https://www.vtk.org/doc/nightly/html/vtkCellType_8h_source.html
    # Translate it into the cells array.
    # `data` is a one-dimensional vector with
    # (num_points0, p0, p1, ... ,pk, numpoints1, p10, p11, ..., p1k, ...
    # or a tuple with (offsets, connectivity)
    has_polygon = numpy.any(types == meshio_to_vtk_type["polygon"])

    cells = []
    cell_data = {}
    if has_polygon:
        numnodes = numpy.empty(len(types), dtype=int)
        # If some polygons are in the VTK file, loop over the cells
        numcells = len(types)
        offsets = numpy.empty(len(types), dtype=int)
        offsets[0] = 0
        for idx in range(numcells - 1):
            numnodes[idx] = data[offsets[idx]]
            offsets[idx + 1] = offsets[idx] + numnodes[idx] + 1

        idx = numcells - 1
        numnodes[idx] = data[offsets[idx]]
        if not numpy.all(numnodes == data[offsets]):
            raise ReadError()

        # TODO: cell_data
        for idx, vtk_cell_type in enumerate(types):
            start = offsets[idx] + 1
            cell_idx = start + _vtk_to_meshio_order(
                vtk_cell_type, numnodes[idx], offsets.dtype
            )
            cell = data[cell_idx]

            cell_type = vtk_to_meshio_type[vtk_cell_type]
            if cell_type == "polygon":
                cell_type += str(data[offsets[idx]])

            if len(cells) > 0 and cells[-1].type == cell_type:
                cells[-1].data.append(cell)
            else:
                cells.append(CellBlock(cell_type, [cell]))

        # convert data to numpy arrays
        for k, c in enumerate(cells):
            cells[k] = CellBlock(c.type, numpy.array(c.data))
    else:
        # Deduct offsets from the cell types. This is much faster than manually going
        # through the data array. Slight disadvantage: This doesn't work for cells with
        # a custom number of points.
        numnodes = vtk_type_to_numnodes[types]
        if not numpy.all(numnodes > 0):
            raise ReadError("File contains cells that meshio cannot handle.")
        if isinstance(data, tuple):
            offsets, conn = data
            if not numpy.all(numnodes == numpy.diff(offsets)):
                raise ReadError()
            idx0 = 0
        else:
            offsets = numpy.cumsum(numnodes + 1) - (numnodes + 1)

            if not numpy.all(numnodes == data[offsets]):
                raise ReadError()
            idx0 = 1
            conn = data

        b = numpy.concatenate(
            [[0], numpy.where(types[:-1] != types[1:])[0] + 1, [len(types)]]
        )
        for start, end in zip(b[:-1], b[1:]):
            meshio_type = vtk_to_meshio_type[types[start]]
            n = numnodes[start]
            cell_idx = idx0 + _vtk_to_meshio_order(types[start], n, dtype=offsets.dtype)
            indices = numpy.add.outer(offsets[start:end], cell_idx)
            cells.append(CellBlock(meshio_type, conn[indices]))
            for name, d in cell_data_raw.items():
                if name not in cell_data:
                    cell_data[name] = []
                cell_data[name].append(d[start:end])

    return cells, cell_data


def write(filename, mesh, binary=True):
    def pad(array):
        return numpy.pad(array, ((0, 0), (0, 1)), "constant")

    if mesh.points.shape[1] == 2:
        points = pad(mesh.points)
    else:
        points = mesh.points

    if mesh.point_data:
        for name, values in mesh.point_data.items():
            if len(values.shape) == 2 and values.shape[1] == 2:
                mesh.point_data[name] = pad(values)

    for name, data in mesh.cell_data.items():
        for k, values in enumerate(data):
            if len(values.shape) == 2 and values.shape[1] == 2:
                data[k] = pad(data[k])

    with open(filename, "wb") as f:
        f.write(b"# vtk DataFile Version 4.2\n")
        f.write("written \n".encode("utf-8"))
        f.write(("BINARY\n" if binary else "ASCII\n").encode("utf-8"))
        f.write(b"DATASET UNSTRUCTURED_GRID\n")

        # write points and cells
        _write_points(f, points, binary)
        _write_cells(f, mesh.cells, binary)

        # write point data
        if mesh.point_data:
            num_points = mesh.points.shape[0]
            f.write("POINT_DATA {}\n".format(num_points).encode("utf-8"))
            _write_field_data(f, mesh.point_data, binary)

        # write cell data
        if mesh.cell_data:
            total_num_cells = sum(len(c.data) for c in mesh.cells)
            f.write("CELL_DATA {}\n".format(total_num_cells).encode("utf-8"))
            _write_field_data(f, mesh.cell_data, binary)


def _write_points(f, points, binary):
    f.write(
        "POINTS {} {}\n".format(
            len(points), numpy_to_vtk_dtype[points.dtype.name]
        ).encode("utf-8")
    )

    if binary:
        # Binary data must be big endian, see
        # <https://www.vtk.org/Wiki/VTK/Writing_VTK_files_using_python#.22legacy.22>.
        # if points.dtype.byteorder == "<" or (
        #     points.dtype.byteorder == "=" and sys.byteorder == "little"
        # ):
        points.astype(points.dtype.newbyteorder(">")).tofile(f, sep="")
    else:
        # ascii
        points.tofile(f, sep=" ")
    f.write(b"\n")


def _write_cells(f, cells, binary):
    total_num_cells = sum([len(c.data) for c in cells])
    total_num_idx = sum([c.data.size for c in cells])
    # For each cell, the number of nodes is stored
    total_num_idx += total_num_cells
    f.write("CELLS {} {}\n".format(total_num_cells,total_num_idx).encode("utf-8"))
    if binary:
        for c in cells:
            n = c.data.shape[1]
            cell_idx = _meshio_to_vtk_order(c.type, n)
            dtype = numpy.dtype(">i4")
            # One must force endianness here:
            # <https://github.com/numpy/numpy/issues/15088>
            numpy.column_stack(
                [
                    numpy.full(c.data.shape[0], n, dtype=dtype),
                    c.data[:, cell_idx].astype(dtype),
                ],
            ).astype(dtype).tofile(f, sep="")
        f.write(b"\n")
    else:
        # ascii
        for c in cells:
            n = c.data.shape[1]
            cell_idx = _meshio_to_vtk_order(c.type, n)
            # prepend a column with the value n
            numpy.column_stack(
                [
                    numpy.full(c.data.shape[0], n, dtype=c.data.dtype),
                    c.data[:, cell_idx],
                ]
            ).tofile(f, sep="\n")
            f.write(b"\n")

    # write cell types
    f.write("CELL_TYPES {}\n".format(total_num_cells).encode("utf-8"))
    if binary:
        for c in cells:
            key_ = c.type[:7] if c.type[:7] == "polygon" else c.type
            vtk_type = meshio_to_vtk_type[key_]
            numpy.full(len(c.data), vtk_type, dtype=numpy.dtype(">i4")).tofile(
                f, sep=""
            )
        f.write(b"\n")
    else:
        # ascii
        for c in cells:
            key_ = c.type[:7] if c.type[:7] == "polygon" else c.type
            numpy.full(len(c.data), meshio_to_vtk_type[key_]).tofile(f, sep="\n")
            f.write(b"\n")


def _write_field_data(f, data, binary):
    f.write(("FIELD FieldData {}\n".format(len(data))).encode("utf-8"))
    for name, values in data.items():
        if isinstance(values, list):
            values = numpy.concatenate(values)
        if len(values.shape) == 1:
            num_tuples = values.shape[0]
            num_components = 1
        else:
            num_tuples = values.shape[0]
            num_components = values.shape[1]

        if " " in name:
            raise WriteError("VTK doesn't support spaces in field names", name)

        f.write(
            (
                "{} {} {} {}\n".format(
                    name,
                    num_components,
                    num_tuples,
                    numpy_to_vtk_dtype[values.dtype.name],
                )
            ).encode("utf-8")
        )
        if binary:
            values.astype(values.dtype.newbyteorder(">")).tofile(f, sep="")
        else:
            # ascii
            values.tofile(f, sep=" ")
            # numpy.savetxt(f, points)
        f.write(b"\n")



if __name__=='__main__':
    #plane=VTKFile('tests/_TODO/FastFarm.Low.DisXY1.t1200.vtk')
    #plane=VTKFile('tests/_TODO/FastFarm.Low.DisXZ1.t1200.vtk')
    plane=VTKFile('tests/_TODO/FastFarm.Low.DisXY1.t0_fake.vtk')
    print(plane.points)
    #plane=VTKFile('tests/_TODO/Main_NM80_OF24_vc.FVW_Hub.AllSeg.000000130.vtk')
    print(plane)
#     print(plane.points)
#     print(plane.cells)
#     print(plane.cell_data_raw)
#     print(plane.cell_data)
#     print('x_grid',plane.x_grid)
#     print('PointData',plane.point_data.keys())
#     print('PointData',plane.point_data_grid.keys())
#     print('PointData',plane.points.shape)
#     print(plane.dataset)
#     if len(plane.z_grid)==1:
#         print('PointData',plane.point_data['DisXY'].shape)
#         D=plane.point_data['DisXY']
#         print(len(plane.x_grid), len(plane.y_grid), len(plane.z_grid),D.shape[1])
# 
#         DD= D.reshape(len(plane.x_grid), len(plane.y_grid), len(plane.z_grid),D.shape[1], order='F')
#         print(DD.shape)
#         import matplotlib.pyplot as plt
#         plt.contourf(plane.x_grid, plane.y_grid, DD[:,:,0,0].T)
#         plt.show()
#     elif len(plane.y_grid)==1:
# 
#         print('PointData',plane.point_data['DisXZ'].shape)
#         D=plane.point_data['DisXZ']
#         print(len(plane.x_grid), len(plane.y_grid), len(plane.z_grid),D.shape[1])
# 
#         DD= D.reshape(len(plane.x_grid), len(plane.y_grid), len(plane.z_grid),D.shape[1], order='F')
#         print(DD.shape)
#         import matplotlib.pyplot as plt
#         #plt.contourf(plane.x_grid, plane.z_grid, DD[:,0,:,1].T, antialiased=False)
#         plt.pcolor(plane.x_grid, plane.z_grid, DD[:,0,:,1].T)
#         plt.show()
