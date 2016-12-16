The two NWTC_Library-related types files cannot be generated through the registry.  At the moment it is a manual process.
   
The NWTC registry input file gets split into two sections: one for the mesh mapping and everything else. 
It's not an automatic process since you have to copy the SetErrStat routine into NWTC_Library_Types.f90, 
and you have to copy the mesh-related types/routines into the ModMesh_Types.f90 file. 
Originally, we also had to change some other parts, too, but I've hard-coded some stuff 
in the registry source code for when it is trying to generate types for the NWTC_Library module. 
We could hard-code the registry to generate SetErrStat() at some point, too.