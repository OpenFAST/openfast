.. _sd_appendix_D:

Appendix D. List of Output Channels
===================================

This is a list of all possible output parameters for the SubDyn module.
The names are grouped by meaning, but can be ordered in the OUTPUT
CHANNELS section of the SubDyn input file as the user sees fit. :math:`M \alpha N \beta`,
refers to node :math:`\beta` of member :math:`\alpha`, where :math:`\alpha` is a number in the range [1,99] and
corresponds to row :math:`\alpha` in the MEMBER OUTPUT LIST table (see :numref:`SD_Member_Output`) and
:math:`\beta` is a number in the range [1,9] and corresponds to node :math:`\beta` in the
**NodeCnt** list of that table entry.

Some outputs are in the SS reference coordinate system (global
inertial-frame coordinate system), and end with the suffix `ss`; others
refer to the local (member) reference system and they have suffixes
"Xe", "Ye", or "Ze" (see Section 7).

Table C-1. List of Output Channels.

+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| Channel Name(s)                       | Units                                                        | Description                                                                                                                         |
+=======================================+==============================================================+=====================================================================================================================================+
| *Base and Interface Reaction Loads*                                                                                                                                                                                                        |           
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| ReactFXss, ReactFYss, ReactFZss,      | (N), (N), (N),                                               | Total base reaction forces and moments                                                                                              |
|                                       |                                                              |                                                                                                                                     |
| ReactMXss, ReactMYss, ReactMZss,      | (Nm), (Nm), (Nm)                                             | at the (0.,0.,-**WtrDpth**) location in SS coordinate system                                                                        |
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| Intf?FXss, Intf?FYss, Intf?FZss,      | (N), (N), (N),                                               | Total interface reaction forces and moments at the TP reference points (platform reference points) in SS coordinate system.         |
|                                       |                                                              | ? can be replaced with any number between 1 and 9 to indicate which transition piece to output.                                     |
| Intf?MXss, Intf?MYss, Intf?MZss,      | (Nm), (Nm), (Nm)                                             | Omitting ? defaults to transition piece 1 for backward compatibility.                                                               |
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| *Interface Kinematics*                                                                                                                                                                                                                     |        
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| Intf?TDXss, Intf?TDYss, Intf?TDZss,   | (m), (m), (m),                                               | Displacements and rotations of the TP reference points in SS coordinate system. The rotation angles are Tait-Bryan angles following |
|                                       |                                                              | the convention of intrinsic yaw first, pitch second, and roll last. ? can be replaced with any number between 1 and 9 to indicate   |
| Intf?RDXss, Intf?RDYss, Intf?RDZss    | (rad), (rad), (rad)                                          | which transition piece to output. Omitting ? defaults to transition piece 1 for backward compatibility.                             |
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| Intf?TDXe, Intf?TDYe, Intf?TDZe,      | (m), (m), (m),                                               | Elastic part of the TP reference point displacements and (small angle) rotations in the rigid-body coordinate system.               |
|                                       |                                                              | ? can be replaced with any number between 1 and 9 to indicate which transition piece to output.                                     |
| Intf?RDXe, Intf?RDYe, Intf?RDZe       | (rad), (rad), (rad)                                          | Omitting ? defaults to transition piece 1 for backward compatibility.                                                               |
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| Intf?TAXss, Intf?TAYss, Intf?TAZss,   | (:math:`{m/s^2}`), (:math:`{m/s^2}`), (:math:`{m/s^2}`),     | Translational and rotational accelerations of the TP reference points in SS coordinate system.                                      |
|                                       |                                                              | ? can be replaced with any number between 1 and 9 to indicate which transition piece to output.                                     |
| Intf?RAXss, Intf?RAYss, Intf?RAZss    | (:math:`{rad/s^2}`), (:math:`{rad/s^2}`), (:math:`{rad/s^2}`)| Omitting ? defaults to transition piece 1 for backward compatibility.                                                               |
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| *Rigid-Body Kinematics (floating only)*                                                                                                                                                                                                    |        
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| RBTDXss, RBTDYss, RBTDZss,            | (m), (m), (m),                                               | Displacements and rotations of the rigid-body reference point in SS coordinate system.                                              |
|                                       |                                                              |                                                                                                                                     |
| RBRDXss, RBRDYss, RBRDZss             | (rad), (rad), (rad)                                          | The rotation angles are Tait-Bryan angles following the convention of intrinsic yaw first, pitch second, and roll last.             |
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| RBTVXss, RBTVYss, RBTVZss,            | (m/s), (m/s), (m/s),                                         | Translational and rotational velocities of the rigid-body reference point in SS coordinate system.                                  |
|                                       |                                                              |                                                                                                                                     |
| RBRVXss, RBRVYss, RBRVZss             | (rad/s), (rad/s), (rad/s)                                    |                                                                                                                                     |
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| RBTAXss, RBTAYss, RBTAZss,            | (:math:`{m/s^2}`), (:math:`{m/s^2}`), (:math:`{m/s^2}`),     | Translational and rotational accelerations of the rigid-body reference point in SS coordinate system.                               |
|                                       |                                                              |                                                                                                                                     |
| RBRAXss, RBRAYss, RBRAZss             | (:math:`{rad/s^2}`), (:math:`{rad/s^2}`), (:math:`{rad/s^2}`)|                                                                                                                                     |
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| *Modal Parameters*                                                                                                                                                                                                                         |              
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| SSqm01-SSqm99                         | (-)                                                          | C-B modal variables (up to first 99)                                                                                                |
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| SSqmd01-SSqmd99                       | (1/s)                                                        | First time-derivatives of C-B modal variables (up to first 99)                                                                      |
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| SSqmdd01-SSqmdd99                     | (:math:`{1/s^2}`)                                            | Second time-derivatives of C-B modal variables (up to first 99)                                                                     |
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| *Node Kinematics*                                                                                                                                                                                                                          |           
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| :math:`{M \alpha N \beta}` TDxss,     | (m),                                                         | Nodal translational displacements of :math:`M \alpha N \beta` (total displacement: rigid-body motion + elastic deflection)	     |
|					|							       | 																     |
| :math:`{M \alpha N \beta}` TDyss, 	| (m),                                                         |  									                                                             |
|					|							       | (up to 99 x 9 = 891 designated locations) in SS coordinate system								     |
| :math:`{M \alpha N \beta}` TDzss,	| (m),                                                         |																     |
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| :math:`{M \alpha N \beta}` RDxe,      | (rad),  						       | Nodal rotational displacements of :math:`M \alpha N \beta` (small elastic part only; no rigid-body rotation)                        |
|					|							       |																     |
| :math:`{M \alpha N \beta}` RDye,	| (rad),  						       |																     |
|                                       |                                                              | (up to 99 x 9 = 891 designated locations) in member local coordinate system                                                         |
| :math:`{M \alpha N \beta}` RDze	| (rad)  						       |																     |
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| :math:`{M \alpha N \beta}` TAxe,	| (:math:`{m/s^2}`),					       | Nodal translational accelerations of :math:`M \alpha N \beta`                                                                       |
|					|							       |																     |
| :math:`{M \alpha N \beta}` TAye,	| (:math:`{m/s^2}`),					       |																     |
|					|							       | (up to 99 x 9 = 891 designated locations) in member local coordinate system                                                         |
| :math:`{M \alpha N \beta}` TAze       | (:math:`{m/s^2}`)					       |																     |
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| :math:`{M \alpha N \beta}` RAxe,	| (:math:`{rad/s^2}`),					       | Nodal rotational accelerations of :math:`M \alpha N \beta`                                                                          |
|					|							       |																     |
| :math:`{M \alpha N \beta}` RAye,	| (:math:`{rad/s^2}`),					       |																     |
|					|							       | (up to 99 x 9 = 891 designated locations) in member local coordinate system                                                         |
| :math:`{M \alpha N \beta}` RAze       | (:math:`{rad/s^2}`)					       |																     |
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| *Node Forces and Moments*                                                                                                                                                                                                                  |           
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| :math:`{M \alpha N \beta}` FKxe,	| (N),                                                         |  Static (elastic) component of reaction forces and moments                                               			     |
|					|                                                              |                                                                                                     				     |
| :math:`{M \alpha N \beta}` FKye,	| (N),                                                         |  at :math:`M \alpha N \beta`  along local member coordinate system                                                                  |               
|					|							       |																     |
| :math:`{M \alpha N \beta}` FKze       | (N),							       | 																     |
|					|							       |																     |
| :math:`{M \alpha N \beta}` MKxe,	| (Nm),							       | 																     |
|					|							       |																     |
| :math:`{M \alpha N \beta}` MKye,	| (Nm),							       | 																     |
|					|							       |																     |
| :math:`{M \alpha N \beta}` MKze       | (Nm)							       | 																     |
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| :math:`{M \alpha N \beta}` FMxe,	| (N),                                                         |  Dynamic (inertial) component of reaction forces and moments                                               			     |
|					|                                                              |                                                                                                     				     |
| :math:`{M \alpha N \beta}` FMye,	| (N),                                                         |  at :math:`M \alpha N \beta`  along local member coordinate system                                                                  |               
|					|							       |																     |
| :math:`{M \alpha N \beta}` FMze       | (N),							       | 																     |
|					|							       |																     |
| :math:`{M \alpha N \beta}` MMxe,	| (Nm),							       | 																     |
|					|							       |																     |
| :math:`{M \alpha N \beta}` MMye,	| (Nm),							       | 																     |
|					|							       |																     |
| :math:`{M \alpha N \beta}` MMze       | (Nm)							       | 																     |
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+

