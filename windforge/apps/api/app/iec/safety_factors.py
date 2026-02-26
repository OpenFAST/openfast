"""
IEC 61400-1 Ed.4 partial safety factors.

Implements the partial safety factor system from IEC 61400-1:2019
Section 7.6 and Tables 3 and 5, covering:
  - Load factors (gamma_f) per DLC
  - Material factors (gamma_m) for normal and special conditions
  - Consequence factors (gamma_n) based on component failure class

The design load is computed as:
  F_d = gamma_f * gamma_n * F_k

And the design resistance as:
  R_d = R_k / gamma_m

where F_k is the characteristic load and R_k is the characteristic
resistance.

References
----------
IEC 61400-1:2019 (Ed.4), Tables 3 and 5, Section 7.6
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum, auto
from typing import Optional


class ComponentClass(Enum):
    """Consequence of failure classification (Table 5).

    Component Class 1: Failure does not result in nacelle or blade
                       failure, limited secondary damage.
    Component Class 2: Failure may result in significant structural
                       damage but no collapse.
    Component Class 3: Failure leads to turbine collapse or thrown
                       components posing risk to people.
    """
    CLASS_1 = auto()  # Lower consequence
    CLASS_2 = auto()  # Normal consequence
    CLASS_3 = auto()  # Higher consequence


class MaterialCategory(Enum):
    """Material safety factor category."""
    NORMAL = auto()    # Standard materials with well-characterized properties
    SPECIAL = auto()   # Materials with less-characterized properties


@dataclass(frozen=True)
class PartialSafetyFactors:
    """Complete set of partial safety factors for a design check.

    Attributes
    ----------
    gamma_f : float
        Partial safety factor for loads.
    gamma_m : float
        Partial safety factor for materials.
    gamma_n : float
        Consequence of failure factor.
    combined_load_factor : float
        Product gamma_f * gamma_n applied to characteristic loads.
    """
    gamma_f: float
    gamma_m: float
    gamma_n: float

    @property
    def combined_load_factor(self) -> float:
        """Combined load and consequence factor: gamma_f * gamma_n."""
        return self.gamma_f * self.gamma_n

    @property
    def total_safety_factor(self) -> float:
        """Total safety factor combining load, consequence, and material."""
        return self.gamma_f * self.gamma_n * self.gamma_m


# ---------------------------------------------------------------------------
# Table 3: Partial safety factors for loads (gamma_f)
# ---------------------------------------------------------------------------

# Mapping from DLC number to default gamma_f
# These are the "normal" design situation values.
DLC_LOAD_FACTORS: dict[str, float] = {
    # Power production
    "1.1": 1.35,
    "1.2": 1.00,  # Fatigue
    "1.3": 1.35,
    "1.4": 1.35,
    "1.5": 1.35,
    # Power production + fault
    "2.1": 1.35,
    "2.2": 1.10,
    "2.3": 1.10,
    "2.4": 1.00,  # Fatigue
    # Start up
    "3.1": 1.00,  # Fatigue
    "3.2": 1.35,
    "3.3": 1.35,
    # Normal shut down
    "4.1": 1.00,  # Fatigue
    "4.2": 1.35,
    # Emergency shut down
    "5.1": 1.35,
    # Parked
    "6.1": 1.35,
    "6.2": 1.10,
    "6.3": 1.35,
    "6.4": 1.00,  # Fatigue
    # Parked + fault
    "7.1": 1.10,
    # Transport / maintenance
    "8.1": 1.50,
}


def get_load_factor(dlc_number: str) -> float:
    """Get the partial safety factor for loads for a specific DLC.

    Parameters
    ----------
    dlc_number : str
        DLC number string (e.g., "1.1", "6.2").

    Returns
    -------
    float
        Partial safety factor gamma_f.

    Raises
    ------
    KeyError
        If the DLC number is not recognized.
    """
    if dlc_number not in DLC_LOAD_FACTORS:
        raise KeyError(
            f"Unknown DLC number '{dlc_number}'. "
            f"Valid DLCs: {sorted(DLC_LOAD_FACTORS.keys())}"
        )
    return DLC_LOAD_FACTORS[dlc_number]


# ---------------------------------------------------------------------------
# Table 5: Material partial safety factors (gamma_m)
# ---------------------------------------------------------------------------

def get_material_factor(
    category: MaterialCategory = MaterialCategory.NORMAL,
) -> float:
    """Get the material partial safety factor.

    IEC 61400-1 Ed.4, Table 5:
      - Normal materials:  gamma_m = 1.1
      - Special materials: gamma_m = 1.3

    Additional factors may be applied for specific material types
    (composites, concrete, etc.) per the relevant material standards.

    Parameters
    ----------
    category : MaterialCategory
        Material category (NORMAL or SPECIAL).

    Returns
    -------
    float
        Material partial safety factor gamma_m.
    """
    if category == MaterialCategory.NORMAL:
        return 1.1
    elif category == MaterialCategory.SPECIAL:
        return 1.3
    else:
        return 1.1  # Default to normal


# ---------------------------------------------------------------------------
# Consequence of failure factor (gamma_n)
# ---------------------------------------------------------------------------

def get_consequence_factor(
    component_class: ComponentClass = ComponentClass.CLASS_2,
) -> float:
    """Get the consequence of failure factor.

    IEC 61400-1 Ed.4, Section 7.6.2:
      - Class 1 (lower consequence):  gamma_n = 0.9
      - Class 2 (normal):             gamma_n = 1.0
      - Class 3 (higher consequence): gamma_n = 1.3

    Parameters
    ----------
    component_class : ComponentClass
        Component failure consequence class.

    Returns
    -------
    float
        Consequence factor gamma_n.
    """
    factors = {
        ComponentClass.CLASS_1: 0.9,
        ComponentClass.CLASS_2: 1.0,
        ComponentClass.CLASS_3: 1.3,
    }
    return factors.get(component_class, 1.0)


# ---------------------------------------------------------------------------
# Combined safety factor computation
# ---------------------------------------------------------------------------

def compute_safety_factors(
    dlc_number: str,
    material_category: MaterialCategory = MaterialCategory.NORMAL,
    component_class: ComponentClass = ComponentClass.CLASS_2,
) -> PartialSafetyFactors:
    """Compute the complete set of partial safety factors.

    Combines load, material, and consequence factors for a specific
    design load case and component configuration.

    Parameters
    ----------
    dlc_number : str
        DLC number (e.g., "1.1").
    material_category : MaterialCategory
        Material category for gamma_m.
    component_class : ComponentClass
        Component failure consequence class for gamma_n.

    Returns
    -------
    PartialSafetyFactors
        Complete set of partial safety factors.

    Examples
    --------
    >>> sf = compute_safety_factors("1.1")
    >>> print(f"gamma_f={sf.gamma_f}, gamma_m={sf.gamma_m}, gamma_n={sf.gamma_n}")
    gamma_f=1.35, gamma_m=1.1, gamma_n=1.0
    >>> print(f"Total safety factor: {sf.total_safety_factor:.3f}")
    Total safety factor: 1.485
    """
    gamma_f = get_load_factor(dlc_number)
    gamma_m = get_material_factor(material_category)
    gamma_n = get_consequence_factor(component_class)

    return PartialSafetyFactors(
        gamma_f=gamma_f,
        gamma_m=gamma_m,
        gamma_n=gamma_n,
    )


def fatigue_safety_factor(
    material_category: MaterialCategory = MaterialCategory.NORMAL,
) -> float:
    """Get the combined fatigue safety factor.

    For fatigue DLCs the load factor gamma_f = 1.0, and the
    consequence factor is typically gamma_n = 1.0 for normal
    design. The material factor still applies.

    IEC 61400-1 Ed.4, Section 7.6.3:
      gamma_mf = 1.1 for normal materials
      gamma_mf = 1.3 for special materials

    Parameters
    ----------
    material_category : MaterialCategory
        Material category.

    Returns
    -------
    float
        Fatigue safety factor for materials.
    """
    return get_material_factor(material_category)
