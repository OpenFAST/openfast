#ifndef MODMESH_H
#define MODMESH_H

#include <stdint.h>

#ifdef OPENFAST_DOUBLE_PRECISION
typedef double Real_t; // Default kind for floating-point numbers
#else
typedef float Real_t; // Default kind for floating-point numbers
#endif

typedef enum
{
    COMPONENT_INPUT = 1,  // parameter for "input mesh"
    COMPONENT_OUTPUT = 2, // parameter for "output mesh"
    COMPONENT_STATE = 3,  // parameter for "state mesh" (not recommended to use)
} ComponentType;

typedef enum
{
    MASKID_FORCE = 1,           // parameter for fields holding force
    MASKID_MOMENT = 2,          // parameter for fields holding moment
    MASKID_ORIENTATION = 3,     // parameter for fields holding orientation
    MASKID_TRANSLATIONDISP = 4, // parameter for fields holding translational displacement
    MASKID_TRANSLATIONVEL = 5,  // parameter for fields holding translational velocity
    MASKID_ROTATIONVEL = 6,     // parameter for fields holding rotational velocity
    MASKID_TRANSLATIONACC = 7,  // parameter for fields holding translational acceleration
    MASKID_ROTATIONACC = 8,     // parameter for fields holding rotational acceleration
    MASKID_SCALAR = 9,          // parameter for fields holding scalars
    FIELDMASK_SIZE = 9,         // maximum number of fields in a mesh
} MeshField;

typedef enum
{
    ELEMENT_POINT = 1,    // parameter for elements of point
    ELEMENT_LINE2 = 2,    // parameter for elements of 2-point lines
    ELEMENT_LINE3 = 3,    // parameter for elements of 3-point lines (currently unused)
    ELEMENT_TRI3 = 4,     // parameter for elements (currently unused)
    ELEMENT_TRI6 = 5,     // parameter for elements (currently unused)
    ELEMENT_QUAD4 = 6,    // parameter for elements (currently unused)
    ELEMENT_QUAD8 = 7,    // parameter for elements (currently unused)
    ELEMENT_TET4 = 8,     // parameter for elements (currently unused)
    ELEMENT_TET10 = 9,    // parameter for elements (currently unused)
    ELEMENT_HEX8 = 10,    // parameter for elements (currently unused)
    ELEMENT_HEX20 = 11,   // parameter for elements (currently unused)
    ELEMENT_WEDGE6 = 12,  // parameter for elements (currently unused)
    ELEMENT_WEDGE15 = 13, // parameter for elements (currently unused)
    NELEMKINDS = 13,      // parameter for maximum number of element kinds
} MeshElementType;

typedef struct
{
    Real_t *Position[3];          // XYZ coordinate of node (3,:)
    double *RefOrientation[3][3]; // Original/reference orientation [DCM] (3,3,:)
    Real_t *Force[3];             // Field: Force vectors (3,NNodes)
    Real_t *Moment[3];            // Field: Moment vectors (3,NNodes)
    double *Orientation[3][3];    // Field: Direction Cosine Matrix (DCM) (3,3,NNodes)
    double *TranslationDisp[3];   // Field: Translational displacements (3,NNodes)
    Real_t *RotationVel[3];       // Field: Rotational velocities (3,NNodes)
    Real_t *TranslationVel[3];    // Field: Translational velocities (3,NNodes)
    Real_t *RotationAcc[3];       // Field: Rotational accelerations (3,NNodes)
    Real_t *TranslationAcc[3];    // Field: Translational accelerations (3,NNodes)
    Real_t **Scalars;             // Scalars (nScalars,NNodes)
    int32_t nScalars;             // Stores value of nScalars when created
    int32_t fields;               // Bit packed integer describing the active fields
} Mesh_t;

#endif // MODMESH_H