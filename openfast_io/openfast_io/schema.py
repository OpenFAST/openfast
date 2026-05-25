"""Human-written parameter schema for openfast_io modules.

DO NOT auto-generate this file from Fortran Registry files.
Entries are written deliberately — description and units require domain knowledge.

To find what parameters are missing for a new OpenFAST version, run:
    python -m openfast_io.tools.check_registry_drift --openfast-root /path/to/openfast
"""

# Schema dict: version → module → param_name → metadata
_SCHEMA: dict[str, dict[str, dict]] = {
    '5.0.0': {
        'ElastoDyn': {
            'Echo':       {'type': bool,  'desc': 'Echo input to .ech file', 'units': None},
            'Method':     {'type': int,   'desc': 'Integration method {1:RK4, 2:AB4, 3:ABM4}',
                           'units': None, 'enum': [1, 2, 3]},
            'DT':         {'type': float, 'desc': 'ElastoDyn integration time step', 'units': 's'},
            'FlapDOF1':   {'type': bool,  'desc': 'First flapwise blade mode DOF', 'units': None},
            'FlapDOF2':   {'type': bool,  'desc': 'Second flapwise blade mode DOF', 'units': None},
            'EdgeDOF':    {'type': bool,  'desc': 'First edgewise blade mode DOF', 'units': None},
            'PitchDOF':   {'type': bool,  'desc': 'Blade pitch DOF', 'units': None},
            'TeetDOF':    {'type': bool,  'desc': 'Rotor-teeter DOF [2-blade only]', 'units': None},
            'DrTrDOF':    {'type': bool,  'desc': 'Drivetrain rotational-flexibility DOF', 'units': None},
            'GenDOF':     {'type': bool,  'desc': 'Generator DOF', 'units': None},
            'YawDOF':     {'type': bool,  'desc': 'Nacelle-yaw DOF', 'units': None},
            'TwFADOF1':   {'type': bool,  'desc': 'First fore-aft tower bending-mode DOF', 'units': None},
            'TwFADOF2':   {'type': bool,  'desc': 'Second fore-aft tower bending-mode DOF', 'units': None},
            'TwSSDOF1':   {'type': bool,  'desc': 'First side-to-side tower bending-mode DOF', 'units': None},
            'TwSSDOF2':   {'type': bool,  'desc': 'Second side-to-side tower bending-mode DOF', 'units': None},
            'PtfmSgDOF':  {'type': bool,  'desc': 'Platform horizontal surge translation DOF', 'units': None},
            'PtfmSwDOF':  {'type': bool,  'desc': 'Platform horizontal sway translation DOF', 'units': None},
            'PtfmHvDOF':  {'type': bool,  'desc': 'Platform vertical heave translation DOF', 'units': None},
            'PtfmRDOF':   {'type': bool,  'desc': 'Platform roll tilt rotation DOF', 'units': None},
            'PtfmPDOF':   {'type': bool,  'desc': 'Platform pitch tilt rotation DOF', 'units': None},
            'PtfmYDOF':   {'type': bool,  'desc': 'Platform yaw rotation DOF', 'units': None},
            'OoPDefl':    {'type': float, 'desc': 'Initial out-of-plane blade-tip displacement', 'units': 'm'},
            'IPDefl':     {'type': float, 'desc': 'Initial in-plane blade-tip displacement', 'units': 'm'},
            'BlPitch1':   {'type': float, 'desc': 'Blade 1 initial pitch', 'units': 'deg'},
            'BlPitch2':   {'type': float, 'desc': 'Blade 2 initial pitch', 'units': 'deg'},
            'BlPitch3':   {'type': float, 'desc': 'Blade 3 initial pitch', 'units': 'deg'},
            'TeetDefl':   {'type': float, 'desc': 'Initial or fixed teeter angle', 'units': 'deg'},
            'Azimuth':    {'type': float, 'desc': 'Initial azimuth angle for blade 1', 'units': 'deg'},
            'RotSpeed':   {'type': float, 'desc': 'Initial or fixed rotor speed', 'units': 'rpm'},
            'NacYaw':     {'type': float, 'desc': 'Initial or fixed nacelle-yaw angle', 'units': 'deg'},
            'NumBl':      {'type': int,   'desc': 'Number of blades', 'units': None, 'enum': [1, 2, 3]},
            'TipRad':     {'type': float, 'desc': 'Blade tip-to-hub radius (preconed)', 'units': 'm'},
            'HubRad':     {'type': float, 'desc': 'Hub radius (preconed)', 'units': 'm'},
            'TowerHt':    {'type': float, 'desc': 'Height of tower above ground or MSL', 'units': 'm'},
            'TowerBsHt':  {'type': float, 'desc': 'Height of tower base above ground or MSL', 'units': 'm'},
            'GBoxEff':    {'type': float, 'desc': 'Gearbox efficiency', 'units': '%'},
            'GBRatio':    {'type': float, 'desc': 'Gearbox ratio', 'units': None},
            'DTTorSpr':   {'type': float, 'desc': 'Drivetrain torsional spring', 'units': 'N-m/rad'},
            'DTTorDmp':   {'type': float, 'desc': 'Drivetrain torsional damper', 'units': 'N-m/(rad/s)'},
            'BldFile':    {'type': list,  'desc': 'Blade property file paths', 'units': None,
                           'is_array': True, 'item_type': str, 'is_file_ref': True},
            'TwrFile':    {'type': str,   'desc': 'Tower property file path', 'units': None,
                           'is_file_ref': True},
            'BldNodes':   {'type': int,   'desc': 'Number of blade nodes per blade for analysis', 'units': None},
            'TwrNodes':   {'type': int,   'desc': 'Number of tower nodes for analysis', 'units': None},
            'HubMass':    {'type': float, 'desc': 'Hub mass', 'units': 'kg'},
            'NacMass':    {'type': float, 'desc': 'Nacelle mass', 'units': 'kg'},
            'PtfmMass':   {'type': float, 'desc': 'Platform mass', 'units': 'kg'},
            'GenIner':    {'type': float, 'desc': 'Generator inertia about HSS', 'units': 'kg m^2'},
            'ShftTilt':   {'type': float, 'desc': 'Rotor shaft tilt angle', 'units': 'deg'},
            'OverHang':   {'type': float, 'desc': 'Distance from yaw axis to rotor apex', 'units': 'm'},
            'Twr2Shft':   {'type': float, 'desc': 'Vertical distance from tower-top to rotor shaft', 'units': 'm'},
        },
        'Fst': {
            'Echo':       {'type': bool,  'desc': 'Echo input to .ech file', 'units': None},
            'AbortLevel': {'type': str,   'desc': 'Error level for aborting', 'units': None},
            'TMax':       {'type': float, 'desc': 'Total run time', 'units': 's'},
            'DT':         {'type': float, 'desc': 'Global recommended time step', 'units': 's'},
            'ModCoupling': {'type': int,  'desc': 'Module coupling method', 'units': None, 'enum': [1, 2, 3]},
            'InterpOrder': {'type': int,  'desc': 'Interpolation order for input-output extrap', 'units': None},
            'NRotors':    {'type': int,   'desc': 'Number of rotors', 'units': None},
            'CompElast':  {'type': int,   'desc': 'Structural dynamics module {1:ElastoDyn, 2:ED+BD, 3:SimpleED}',
                           'units': None, 'enum': [1, 2, 3]},
            'CompInflow': {'type': int,   'desc': 'Inflow wind module {0:still air, 1:InflowWind, 2:ExtInflow}',
                           'units': None, 'enum': [0, 1, 2]},
            'CompAero':   {'type': int,   'desc': 'Aerodynamics module {0:none, 1:AeroDisk, 2:AeroDyn, 3:ExtLoads}',
                           'units': None, 'enum': [0, 1, 2, 3]},
            'CompServo':  {'type': int,   'desc': 'Controller/drive module {0:none, 1:ServoDyn}',
                           'units': None, 'enum': [0, 1]},
            'CompHydro':  {'type': int,   'desc': 'Hydrodynamics module {0:none, 1:HydroDyn}',
                           'units': None, 'enum': [0, 1]},
            'CompSeaSt':  {'type': int,   'desc': 'Sea state module {0:none, 1:SeaState}',
                           'units': None, 'enum': [0, 1]},
            'CompSub':    {'type': int,   'desc': 'Substructure module {0:none, 1:SubDyn, 2:ExtPtfm}',
                           'units': None, 'enum': [0, 1, 2]},
            'CompMooring': {'type': int,  'desc': 'Mooring module {0:none, 1:MAP, 2:FEAMooring, 3:MoorDyn, 4:OrcaFlex}',
                            'units': None, 'enum': [0, 1, 2, 3, 4]},
            'CompIce':    {'type': int,   'desc': 'Ice module {0:none, 1:IceFloe, 2:IceDyn}',
                           'units': None, 'enum': [0, 1, 2]},
            'CompSoil':   {'type': int,   'desc': 'Soil module {0:none, 1:SoilDyn}',
                           'units': None, 'enum': [0, 1]},
            'MHK':        {'type': int,   'desc': 'MHK turbine type {0:Not MHK, 1:Fixed MHK, 2:Floating MHK}',
                           'units': None, 'enum': [0, 1, 2]},
            'Gravity':    {'type': float, 'desc': 'Gravitational acceleration', 'units': 'm/s^2'},
            'AirDens':    {'type': float, 'desc': 'Air density', 'units': 'kg/m^3'},
            'WtrDens':    {'type': float, 'desc': 'Water density', 'units': 'kg/m^3'},
            'WtrDpth':    {'type': float, 'desc': 'Water depth', 'units': 'm'},
            'EDFile':     {'type': str,   'desc': 'ElastoDyn input file path', 'units': None,
                           'is_file_ref': True},
            'AeroFile':   {'type': str,   'desc': 'AeroDyn input file path', 'units': None,
                           'is_file_ref': True},
            'ServoFile':  {'type': str,   'desc': 'ServoDyn input file path', 'units': None,
                           'is_file_ref': True},
            'InflowFile': {'type': str,   'desc': 'InflowWind input file path', 'units': None,
                           'is_file_ref': True},
            'HydroFile':  {'type': str,   'desc': 'HydroDyn input file path', 'units': None,
                           'is_file_ref': True},
            'SeaStFile':  {'type': str,   'desc': 'SeaState input file path', 'units': None,
                           'is_file_ref': True},
            'SubFile':    {'type': str,   'desc': 'SubDyn/ExtPtfm input file path', 'units': None,
                           'is_file_ref': True},
            'MooringFile': {'type': str,  'desc': 'Mooring input file path', 'units': None,
                            'is_file_ref': True},
            'SoilFile':   {'type': str,   'desc': 'SoilDyn input file path', 'units': None,
                           'is_file_ref': True},
        },
    },
    '4.0.0': {
        # Params present in 4.x but removed in 5.0.0
        'AeroDyn': {
            'Buoyancy': {'type': bool, 'desc': 'Enable buoyancy effects (removed in v5.0.0)', 'units': None},
        },
        'BeamDyn': {
            'UsePitchAct': {'type': bool,  'desc': 'Use a pitch actuator (removed in v5.0.0)', 'units': None},
            'PitchJ':      {'type': float, 'desc': 'Pitch actuator inertia (removed in v5.0.0)', 'units': 'kg*m^2'},
            'PitchK':      {'type': float, 'desc': 'Pitch actuator stiffness (removed in v5.0.0)', 'units': 'N*m/rad'},
            'PitchC':      {'type': float, 'desc': 'Pitch actuator damping (removed in v5.0.0)', 'units': 'N*m*s/rad'},
        },
    },
}

# File-reference parameters by module — used by validation
FILE_REF_PARAMS: dict[str, list[str]] = {
    'Fst':       ['EDFile', 'BDBldFile(1)', 'BDBldFile(2)', 'BDBldFile(3)',
                  'InflowFile', 'AeroFile', 'ServoFile',
                  'SeaStFile', 'HydroFile', 'SubFile', 'MooringFile', 'SoilFile'],
    'ElastoDyn': ['BldFile1', 'BldFile2', 'BldFile3', 'TwrFile'],
    'AeroDyn':   ['ADBlFile(1)', 'ADBlFile(2)', 'ADBlFile(3)'],
    'ServoDyn':  ['DLL_FileName'],
}


def get_schema(module: str, version: str = '5.0.0') -> dict:
    """Get parameter schema for a module at a given version.

    For v5.0.0, returns the v5 schema directly.
    For v4.0.0, returns the v5 schema merged with v4-specific params (those removed in v5).
    """
    base = _SCHEMA.get('5.0.0', {}).get(module, {}).copy()
    if version != '5.0.0' and version in _SCHEMA:
        base.update(_SCHEMA[version].get(module, {}))
    return base


def get_param_info(module: str, param: str, version: str = '5.0.0') -> dict:
    """Get metadata for a single parameter. Returns {} if not in schema."""
    return get_schema(module, version).get(param, {})
