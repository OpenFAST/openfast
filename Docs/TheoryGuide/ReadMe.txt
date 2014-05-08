The documents make the most sense if studied in the following order:

1) FASTDOFs.xls -- Contains a listing of the DOF indices used by the equations of motion in FAST.
2) FASTCoordinateSystems.doc -- Documents the transformation matrices relating each coordinate system in FAST.  Unfortunately, I do not have pictures in this document that diagram these coordinate systems.  Hopefully, you can visualize them by means of the transformation matrices.
3) FASTKinematics.doc -- Documents the linear position, velocity, and acceleration vectors of each "important" point in the system and documents the angular velocity and acceleration vectors of each "important" reference frame in the system.  Also included is documentation of the partial velocity vectors needed by Kane's dynamics.
4) FASTKinetics.doc -- Documents the derivation of the equations of motion using Kane's dynamics.
5) FASTLoads.doc -- Documents how the output loads are computed using terms from the equations of motion.
6) FASTMotions.doc -- Documents how the output motions are computed using variables from the equations of motion.
7) FASTLogicFlow.doc -- Contains a listing of the subroutine names used by FAST.  The names are listed in the order they are called within the program


bjj: note that these documents were written for FAST v6. As of FAST v8, the term "FAST" in this document applies mostly to ElastoDyn.