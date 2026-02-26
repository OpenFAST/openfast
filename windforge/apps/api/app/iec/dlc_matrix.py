"""
IEC 61400-1 Ed.4 Design Load Case (DLC) matrix definition.

Defines all 20 DLCs from IEC 61400-1 Ed.4 Table 4 as structured
dataclass entries, along with wind model and analysis type enumerations.

Each DLC template specifies the wind model, operating condition,
analysis type, partial safety factor, number of seeds, and any
special conditions (faults, yaw misalignment, special wind speeds).

References
----------
IEC 61400-1:2019 (Ed.4) Table 4 - Design load cases
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Optional


class WindModel(Enum):
    """IEC wind condition models (Section 6.3)."""
    NTM = auto()    # Normal Turbulence Model
    ETM = auto()    # Extreme Turbulence Model
    EWM = auto()    # Extreme Wind speed Model (50-yr or 1-yr)
    EOG = auto()    # Extreme Operating Gust
    EDC = auto()    # Extreme Direction Change
    ECD = auto()    # Extreme Coherent gust with Direction change
    EWS = auto()    # Extreme Wind Shear
    NWP = auto()    # Normal Wind Profile (deterministic, for steady)


class AnalysisType(Enum):
    """Load analysis type for each DLC."""
    ULTIMATE = auto()   # Ultimate strength analysis
    FATIGUE = auto()    # Fatigue damage analysis


class OperatingCondition(Enum):
    """Turbine operating condition during the DLC."""
    POWER_PRODUCTION = auto()
    POWER_PRODUCTION_PLUS_FAULT = auto()
    START_UP = auto()
    NORMAL_SHUT_DOWN = auto()
    EMERGENCY_SHUT_DOWN = auto()
    PARKED_STANDING_STILL = auto()
    PARKED_PLUS_FAULT = auto()
    TRANSPORT_ASSEMBLY_MAINTENANCE = auto()


class FaultType(Enum):
    """Type of fault applicable to the DLC."""
    NONE = auto()
    CONTROL_SYSTEM = auto()
    PROTECTION_SYSTEM = auto()
    ELECTRICAL = auto()
    BLADE_SEIZURE = auto()
    PITCH_RUNAWAY = auto()
    GRID_LOSS = auto()
    YAW_SYSTEM = auto()


@dataclass(frozen=True)
class DLCTemplate:
    """Template definition for one IEC Design Load Case.

    Attributes
    ----------
    number : str
        DLC number (e.g., "1.1", "6.2a").
    description : str
        Human-readable description of the load case.
    wind_model : WindModel
        Wind condition model to apply.
    operating_condition : OperatingCondition
        Turbine operating state during this DLC.
    analysis_type : AnalysisType
        Whether this is an ultimate or fatigue analysis.
    default_safety_factor : float
        Default partial safety factor for loads (gamma_f).
    default_num_seeds : int
        Default number of turbulent seeds per wind speed bin.
    needs_fault : bool
        Whether this DLC requires a fault condition.
    fault_type : FaultType
        Type of fault to simulate (if needs_fault is True).
    yaw_misalignment : float
        Default yaw misalignment angle (degrees). 0 means no misalignment
        or use DLC-specific logic.
    special_wind_speed : Optional[str]
        If set, indicates a special wind speed requirement:
        "Vr" (rated), "Vout" (cut-out), "Ve50" (50-yr extreme),
        "Ve1" (1-yr extreme), "Vin" (cut-in), "Vr+/-2", "Vhub".
    wind_speed_range : str
        Description of the wind speed range for this DLC.
        "Vin_to_Vout" for the full operating range,
        or a specific value description.
    """
    number: str
    description: str
    wind_model: WindModel
    operating_condition: OperatingCondition
    analysis_type: AnalysisType
    default_safety_factor: float
    default_num_seeds: int
    needs_fault: bool = False
    fault_type: FaultType = FaultType.NONE
    yaw_misalignment: float = 0.0
    special_wind_speed: Optional[str] = None
    wind_speed_range: str = "Vin_to_Vout"


# ---------------------------------------------------------------------------
# IEC 61400-1 Ed.4 Table 4 -- Complete DLC Matrix
# ---------------------------------------------------------------------------

IEC_DLC_TABLE: list[DLCTemplate] = [
    # ===================================================================
    # 1. Power production
    # ===================================================================
    DLCTemplate(
        number="1.1",
        description="Power production - Normal turbulence, ultimate",
        wind_model=WindModel.NTM,
        operating_condition=OperatingCondition.POWER_PRODUCTION,
        analysis_type=AnalysisType.ULTIMATE,
        default_safety_factor=1.35,
        default_num_seeds=6,
        yaw_misalignment=0.0,
        wind_speed_range="Vin_to_Vout",
    ),
    DLCTemplate(
        number="1.2",
        description="Power production - Normal turbulence, fatigue",
        wind_model=WindModel.NTM,
        operating_condition=OperatingCondition.POWER_PRODUCTION,
        analysis_type=AnalysisType.FATIGUE,
        default_safety_factor=1.0,
        default_num_seeds=6,
        yaw_misalignment=0.0,
        wind_speed_range="Vin_to_Vout",
    ),
    DLCTemplate(
        number="1.3",
        description="Power production - Extreme turbulence model",
        wind_model=WindModel.ETM,
        operating_condition=OperatingCondition.POWER_PRODUCTION,
        analysis_type=AnalysisType.ULTIMATE,
        default_safety_factor=1.35,
        default_num_seeds=6,
        yaw_misalignment=0.0,
        wind_speed_range="Vin_to_Vout",
    ),
    DLCTemplate(
        number="1.4",
        description="Power production - Extreme coherent gust with direction change",
        wind_model=WindModel.ECD,
        operating_condition=OperatingCondition.POWER_PRODUCTION,
        analysis_type=AnalysisType.ULTIMATE,
        default_safety_factor=1.35,
        default_num_seeds=1,
        yaw_misalignment=0.0,
        special_wind_speed="Vr-2_Vr_Vr+2",
        wind_speed_range="Vr-2_Vr_Vr+2",
    ),
    DLCTemplate(
        number="1.5",
        description="Power production - Extreme wind shear",
        wind_model=WindModel.EWS,
        operating_condition=OperatingCondition.POWER_PRODUCTION,
        analysis_type=AnalysisType.ULTIMATE,
        default_safety_factor=1.35,
        default_num_seeds=1,
        yaw_misalignment=0.0,
        wind_speed_range="Vin_to_Vout",
    ),

    # ===================================================================
    # 2. Power production plus occurrence of fault
    # ===================================================================
    DLCTemplate(
        number="2.1",
        description="Power production + fault - Control system fault",
        wind_model=WindModel.NTM,
        operating_condition=OperatingCondition.POWER_PRODUCTION_PLUS_FAULT,
        analysis_type=AnalysisType.ULTIMATE,
        default_safety_factor=1.35,
        default_num_seeds=6,
        needs_fault=True,
        fault_type=FaultType.CONTROL_SYSTEM,
        wind_speed_range="Vin_to_Vout",
    ),
    DLCTemplate(
        number="2.2",
        description="Power production + fault - Protection system fault",
        wind_model=WindModel.NTM,
        operating_condition=OperatingCondition.POWER_PRODUCTION_PLUS_FAULT,
        analysis_type=AnalysisType.ULTIMATE,
        default_safety_factor=1.1,
        default_num_seeds=6,
        needs_fault=True,
        fault_type=FaultType.PROTECTION_SYSTEM,
        wind_speed_range="Vin_to_Vout",
    ),
    DLCTemplate(
        number="2.3",
        description="Power production + fault - Extreme operating gust with fault",
        wind_model=WindModel.EOG,
        operating_condition=OperatingCondition.POWER_PRODUCTION_PLUS_FAULT,
        analysis_type=AnalysisType.ULTIMATE,
        default_safety_factor=1.1,
        default_num_seeds=1,
        needs_fault=True,
        fault_type=FaultType.ELECTRICAL,
        special_wind_speed="Vr-2_Vr_Vr+2_Vout",
        wind_speed_range="Vr-2_Vr_Vr+2_Vout",
    ),
    DLCTemplate(
        number="2.4",
        description="Power production + fault - NTM with fault, fatigue",
        wind_model=WindModel.NTM,
        operating_condition=OperatingCondition.POWER_PRODUCTION_PLUS_FAULT,
        analysis_type=AnalysisType.FATIGUE,
        default_safety_factor=1.0,
        default_num_seeds=6,
        needs_fault=True,
        fault_type=FaultType.CONTROL_SYSTEM,
        wind_speed_range="Vin_to_Vout",
    ),

    # ===================================================================
    # 3. Start up
    # ===================================================================
    DLCTemplate(
        number="3.1",
        description="Start up - Normal wind profile",
        wind_model=WindModel.NWP,
        operating_condition=OperatingCondition.START_UP,
        analysis_type=AnalysisType.FATIGUE,
        default_safety_factor=1.0,
        default_num_seeds=1,
        special_wind_speed="Vin_Vr_Vout",
        wind_speed_range="Vin_Vr_Vout",
    ),
    DLCTemplate(
        number="3.2",
        description="Start up - Extreme operating gust",
        wind_model=WindModel.EOG,
        operating_condition=OperatingCondition.START_UP,
        analysis_type=AnalysisType.ULTIMATE,
        default_safety_factor=1.35,
        default_num_seeds=1,
        special_wind_speed="Vin_Vr-2_Vr_Vr+2",
        wind_speed_range="Vin_Vr-2_Vr_Vr+2",
    ),
    DLCTemplate(
        number="3.3",
        description="Start up - Extreme direction change",
        wind_model=WindModel.EDC,
        operating_condition=OperatingCondition.START_UP,
        analysis_type=AnalysisType.ULTIMATE,
        default_safety_factor=1.35,
        default_num_seeds=1,
        special_wind_speed="Vin_Vr-2_Vr_Vr+2",
        wind_speed_range="Vin_Vr-2_Vr_Vr+2",
    ),

    # ===================================================================
    # 4. Normal shut down
    # ===================================================================
    DLCTemplate(
        number="4.1",
        description="Normal shut down - Normal wind profile",
        wind_model=WindModel.NWP,
        operating_condition=OperatingCondition.NORMAL_SHUT_DOWN,
        analysis_type=AnalysisType.FATIGUE,
        default_safety_factor=1.0,
        default_num_seeds=1,
        special_wind_speed="Vr-2_Vr_Vr+2_Vout",
        wind_speed_range="Vr-2_Vr_Vr+2_Vout",
    ),
    DLCTemplate(
        number="4.2",
        description="Normal shut down - Extreme operating gust",
        wind_model=WindModel.EOG,
        operating_condition=OperatingCondition.NORMAL_SHUT_DOWN,
        analysis_type=AnalysisType.ULTIMATE,
        default_safety_factor=1.35,
        default_num_seeds=1,
        special_wind_speed="Vr-2_Vr_Vr+2_Vout",
        wind_speed_range="Vr-2_Vr_Vr+2_Vout",
    ),

    # ===================================================================
    # 5. Emergency shut down
    # ===================================================================
    DLCTemplate(
        number="5.1",
        description="Emergency shut down",
        wind_model=WindModel.NTM,
        operating_condition=OperatingCondition.EMERGENCY_SHUT_DOWN,
        analysis_type=AnalysisType.ULTIMATE,
        default_safety_factor=1.35,
        default_num_seeds=6,
        special_wind_speed="Vr-2_Vr_Vr+2",
        wind_speed_range="Vr-2_Vr_Vr+2",
    ),

    # ===================================================================
    # 6. Parked (standing still or idling)
    # ===================================================================
    DLCTemplate(
        number="6.1",
        description="Parked - Extreme wind model 50-yr recurrence",
        wind_model=WindModel.EWM,
        operating_condition=OperatingCondition.PARKED_STANDING_STILL,
        analysis_type=AnalysisType.ULTIMATE,
        default_safety_factor=1.35,
        default_num_seeds=6,
        yaw_misalignment=0.0,
        special_wind_speed="Ve50",
        wind_speed_range="Ve50",
    ),
    DLCTemplate(
        number="6.2",
        description="Parked - Extreme wind model 50-yr with loss of electrical network",
        wind_model=WindModel.EWM,
        operating_condition=OperatingCondition.PARKED_PLUS_FAULT,
        analysis_type=AnalysisType.ULTIMATE,
        default_safety_factor=1.1,
        default_num_seeds=6,
        needs_fault=True,
        fault_type=FaultType.GRID_LOSS,
        yaw_misalignment=180.0,
        special_wind_speed="Ve50",
        wind_speed_range="Ve50",
    ),
    DLCTemplate(
        number="6.3",
        description="Parked - Extreme wind model 1-yr recurrence",
        wind_model=WindModel.EWM,
        operating_condition=OperatingCondition.PARKED_STANDING_STILL,
        analysis_type=AnalysisType.ULTIMATE,
        default_safety_factor=1.35,
        default_num_seeds=6,
        yaw_misalignment=20.0,
        special_wind_speed="Ve1",
        wind_speed_range="Ve1",
    ),
    DLCTemplate(
        number="6.4",
        description="Parked - Normal turbulence model, fatigue",
        wind_model=WindModel.NTM,
        operating_condition=OperatingCondition.PARKED_STANDING_STILL,
        analysis_type=AnalysisType.FATIGUE,
        default_safety_factor=1.0,
        default_num_seeds=6,
        wind_speed_range="Vin_to_Vout",
    ),

    # ===================================================================
    # 7. Parked and fault conditions
    # ===================================================================
    DLCTemplate(
        number="7.1",
        description="Parked + fault - Extreme wind model 1-yr recurrence",
        wind_model=WindModel.EWM,
        operating_condition=OperatingCondition.PARKED_PLUS_FAULT,
        analysis_type=AnalysisType.ULTIMATE,
        default_safety_factor=1.1,
        default_num_seeds=6,
        needs_fault=True,
        fault_type=FaultType.YAW_SYSTEM,
        yaw_misalignment=180.0,
        special_wind_speed="Ve1",
        wind_speed_range="Ve1",
    ),

    # ===================================================================
    # 8. Transport, assembly, maintenance, and repair
    # ===================================================================
    DLCTemplate(
        number="8.1",
        description="Transport, assembly, maintenance and repair",
        wind_model=WindModel.EWM,
        operating_condition=OperatingCondition.TRANSPORT_ASSEMBLY_MAINTENANCE,
        analysis_type=AnalysisType.ULTIMATE,
        default_safety_factor=1.5,
        default_num_seeds=1,
        special_wind_speed="Vmaint",
        wind_speed_range="Vmaint",
    ),
]


def get_dlc_by_number(dlc_number: str) -> Optional[DLCTemplate]:
    """Look up a DLC template by its number string.

    Parameters
    ----------
    dlc_number : str
        DLC number, e.g. "1.1", "6.2", "2.3".

    Returns
    -------
    Optional[DLCTemplate]
        The matching DLC template, or None if not found.
    """
    for dlc in IEC_DLC_TABLE:
        if dlc.number == dlc_number:
            return dlc
    return None


def get_dlcs_by_type(analysis_type: AnalysisType) -> list[DLCTemplate]:
    """Get all DLCs of a specific analysis type.

    Parameters
    ----------
    analysis_type : AnalysisType
        ULTIMATE or FATIGUE.

    Returns
    -------
    list[DLCTemplate]
        All DLC templates matching the specified analysis type.
    """
    return [dlc for dlc in IEC_DLC_TABLE if dlc.analysis_type == analysis_type]
