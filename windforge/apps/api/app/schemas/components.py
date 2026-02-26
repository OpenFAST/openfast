"""Component schemas: Tower, Blade, Airfoil, Controller, TurbineModel."""

from datetime import datetime
from uuid import UUID

from pydantic import BaseModel, Field


# ---------------------------------------------------------------------------
# Nested station / polar schemas
# ---------------------------------------------------------------------------
class TowerStationSchema(BaseModel):
    """A single tower distributed-property station."""

    frac: float = Field(..., ge=0.0, le=1.0, description="Fractional height (0=base, 1=top)")
    mass_den: float = Field(..., gt=0, description="Mass density (kg/m)")
    fa_stiff: float = Field(..., gt=0, description="Fore-aft stiffness (N m^2)")
    ss_stiff: float = Field(..., gt=0, description="Side-side stiffness (N m^2)")
    outer_diameter: float = Field(..., gt=0, description="Outer diameter (m)")
    wall_thickness: float = Field(..., gt=0, description="Wall thickness (m)")


class BladeStructuralStationSchema(BaseModel):
    """A single blade structural station."""

    frac: float = Field(..., ge=0.0, le=1.0)
    pitch_axis: float = Field(..., ge=0.0, le=1.0)
    struct_twist: float  # deg
    mass_den: float = Field(..., gt=0)  # kg/m
    flap_stiff: float = Field(..., gt=0)  # N m^2
    edge_stiff: float = Field(..., gt=0)  # N m^2


class BladeAeroStationSchema(BaseModel):
    """A single blade aerodynamic station."""

    frac: float = Field(..., ge=0.0, le=1.0)
    chord: float = Field(..., gt=0)  # m
    aero_twist: float  # deg
    airfoil_id: str  # UUID as string
    aero_center: float = Field(default=0.25, ge=0.0, le=1.0)


class AirfoilPolarSchema(BaseModel):
    """Polar data for a single Reynolds number."""

    re: float = Field(..., gt=0, description="Reynolds number")
    alpha: list[float] = Field(..., description="Angle of attack array (deg)")
    cl: list[float] = Field(..., description="Lift coefficient array")
    cd: list[float] = Field(..., description="Drag coefficient array")
    cm: list[float] = Field(..., description="Moment coefficient array")


# ---------------------------------------------------------------------------
# Tower
# ---------------------------------------------------------------------------
class TowerCreate(BaseModel):
    name: str = Field(..., min_length=1, max_length=255)
    tower_height: float = Field(..., gt=0)
    tower_base_height: float = Field(default=10.0, ge=0)
    tower_fa_damping_1: float = Field(default=1.0, ge=0)
    tower_fa_damping_2: float = Field(default=1.0, ge=0)
    tower_ss_damping_1: float = Field(default=1.0, ge=0)
    tower_ss_damping_2: float = Field(default=1.0, ge=0)
    stations: list[TowerStationSchema] | None = None
    fa_mode_1_coeffs: list[float] | None = None
    fa_mode_2_coeffs: list[float] | None = None
    ss_mode_1_coeffs: list[float] | None = None
    ss_mode_2_coeffs: list[float] | None = None


class TowerUpdate(BaseModel):
    name: str | None = Field(None, min_length=1, max_length=255)
    tower_height: float | None = Field(None, gt=0)
    tower_base_height: float | None = Field(None, ge=0)
    tower_fa_damping_1: float | None = Field(None, ge=0)
    tower_fa_damping_2: float | None = Field(None, ge=0)
    tower_ss_damping_1: float | None = Field(None, ge=0)
    tower_ss_damping_2: float | None = Field(None, ge=0)
    stations: list[TowerStationSchema] | None = None
    fa_mode_1_coeffs: list[float] | None = None
    fa_mode_2_coeffs: list[float] | None = None
    ss_mode_1_coeffs: list[float] | None = None
    ss_mode_2_coeffs: list[float] | None = None


class TowerResponse(BaseModel):
    id: UUID
    project_id: UUID
    name: str
    version: int
    tower_height: float
    tower_base_height: float
    tower_fa_damping_1: float
    tower_fa_damping_2: float
    tower_ss_damping_1: float
    tower_ss_damping_2: float
    stations: list[TowerStationSchema] | None = None
    fa_mode_1_coeffs: list[float] | None = None
    fa_mode_2_coeffs: list[float] | None = None
    ss_mode_1_coeffs: list[float] | None = None
    ss_mode_2_coeffs: list[float] | None = None
    is_active: bool
    created_at: datetime

    model_config = {"from_attributes": True}


# ---------------------------------------------------------------------------
# Blade
# ---------------------------------------------------------------------------
class BladeCreate(BaseModel):
    name: str = Field(..., min_length=1, max_length=255)
    blade_length: float = Field(..., gt=0)
    structural_stations: list[BladeStructuralStationSchema] | None = None
    aero_stations: list[BladeAeroStationSchema] | None = None
    flap_mode_1_coeffs: list[float] | None = None
    flap_mode_2_coeffs: list[float] | None = None
    edge_mode_1_coeffs: list[float] | None = None
    blade_flap_damping: float = Field(default=2.0, ge=0)
    blade_edge_damping: float = Field(default=2.0, ge=0)


class BladeUpdate(BaseModel):
    name: str | None = Field(None, min_length=1, max_length=255)
    blade_length: float | None = Field(None, gt=0)
    structural_stations: list[BladeStructuralStationSchema] | None = None
    aero_stations: list[BladeAeroStationSchema] | None = None
    flap_mode_1_coeffs: list[float] | None = None
    flap_mode_2_coeffs: list[float] | None = None
    edge_mode_1_coeffs: list[float] | None = None
    blade_flap_damping: float | None = Field(None, ge=0)
    blade_edge_damping: float | None = Field(None, ge=0)


class BladeResponse(BaseModel):
    id: UUID
    project_id: UUID
    name: str
    version: int
    blade_length: float
    structural_stations: list[BladeStructuralStationSchema] | None = None
    aero_stations: list[BladeAeroStationSchema] | None = None
    flap_mode_1_coeffs: list[float] | None = None
    flap_mode_2_coeffs: list[float] | None = None
    edge_mode_1_coeffs: list[float] | None = None
    blade_flap_damping: float
    blade_edge_damping: float
    is_active: bool
    created_at: datetime

    model_config = {"from_attributes": True}


# ---------------------------------------------------------------------------
# Airfoil
# ---------------------------------------------------------------------------
class AirfoilCreate(BaseModel):
    name: str = Field(..., min_length=1, max_length=255)
    family: str | None = None
    thickness_ratio: float | None = Field(None, gt=0, le=1.0)
    polars: list[AirfoilPolarSchema] | None = None
    source: str | None = None


class AirfoilUpdate(BaseModel):
    name: str | None = Field(None, min_length=1, max_length=255)
    family: str | None = None
    thickness_ratio: float | None = None
    polars: list[AirfoilPolarSchema] | None = None
    source: str | None = None


class AirfoilResponse(BaseModel):
    id: UUID
    org_id: UUID
    name: str
    family: str | None = None
    thickness_ratio: float | None = None
    polars: list[AirfoilPolarSchema] | None = None
    source: str | None = None
    created_at: datetime

    model_config = {"from_attributes": True}


# ---------------------------------------------------------------------------
# Controller
# ---------------------------------------------------------------------------
class ControllerCreate(BaseModel):
    name: str = Field(..., min_length=1, max_length=255)
    controller_type: str = Field(default="baseline")
    pcmode: int = Field(default=0, ge=0)
    vscontrl: int = Field(default=1, ge=0)
    parameters: dict | None = None
    dll_filename: str | None = None
    dll_procname: str | None = None


class ControllerUpdate(BaseModel):
    name: str | None = Field(None, min_length=1, max_length=255)
    controller_type: str | None = None
    pcmode: int | None = Field(None, ge=0)
    vscontrl: int | None = Field(None, ge=0)
    parameters: dict | None = None
    dll_filename: str | None = None
    dll_procname: str | None = None


class ControllerResponse(BaseModel):
    id: UUID
    project_id: UUID
    name: str
    version: int
    controller_type: str
    pcmode: int
    vscontrl: int
    parameters: dict | None = None
    dll_filename: str | None = None
    dll_procname: str | None = None
    is_active: bool
    created_at: datetime

    model_config = {"from_attributes": True}


# ---------------------------------------------------------------------------
# TurbineModel
# ---------------------------------------------------------------------------
class TurbineModelCreate(BaseModel):
    name: str = Field(..., min_length=1, max_length=255)
    tower_id: UUID | None = None
    blade_id: UUID | None = None
    controller_id: UUID | None = None
    gearbox_ratio: float | None = None
    generator_inertia: float | None = None
    drivetrain_stiffness: float | None = None
    drivetrain_damping: float | None = None
    hub_mass: float | None = None
    hub_inertia: float | None = None
    nacelle_mass: float | None = None
    nacelle_inertia: float | None = None
    overhang: float | None = None
    shaft_tilt: float | None = None
    precone: float | None = None
    rotor_speed_rated: float | None = None
    dof_flags: dict | None = None


class TurbineModelUpdate(BaseModel):
    name: str | None = Field(None, min_length=1, max_length=255)
    tower_id: UUID | None = None
    blade_id: UUID | None = None
    controller_id: UUID | None = None
    gearbox_ratio: float | None = None
    generator_inertia: float | None = None
    drivetrain_stiffness: float | None = None
    drivetrain_damping: float | None = None
    hub_mass: float | None = None
    hub_inertia: float | None = None
    nacelle_mass: float | None = None
    nacelle_inertia: float | None = None
    overhang: float | None = None
    shaft_tilt: float | None = None
    precone: float | None = None
    rotor_speed_rated: float | None = None
    dof_flags: dict | None = None


class TurbineModelResponse(BaseModel):
    id: UUID
    project_id: UUID
    name: str
    version: int
    tower_id: UUID | None = None
    blade_id: UUID | None = None
    controller_id: UUID | None = None
    gearbox_ratio: float | None = None
    generator_inertia: float | None = None
    drivetrain_stiffness: float | None = None
    drivetrain_damping: float | None = None
    hub_mass: float | None = None
    hub_inertia: float | None = None
    nacelle_mass: float | None = None
    nacelle_inertia: float | None = None
    overhang: float | None = None
    shaft_tilt: float | None = None
    precone: float | None = None
    rotor_speed_rated: float | None = None
    dof_flags: dict | None = None
    is_active: bool
    created_at: datetime

    model_config = {"from_attributes": True}
