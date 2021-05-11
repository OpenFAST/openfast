.. _sd_appendix_D:

Appendix D. List of Output Channels
===================================

This is a list of all possible output parameters for the SubDyn module.
The names are grouped by meaning, but can be ordered in the OUTPUT
CHANNELS section of the SubDyn input file as the user sees fit. :math:`M \alpha N \beta`,
refers to node :math:`\beta` of member :math:`\alpha`, where :math:`\alpha` is a number in the range [1,9] and
corresponds to row :math:`\alpha` in the MEMBER OUTPUT LIST table (see Section ) and
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
| *Base and Interface Reaction Loads*   |                                                                                                                                                                                                    |           
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| ReactFXss, ReactFYss, ReactFZss,      | (N), (N), (N),                                               | Total base reaction forces and moments                                                                                              |
|                                       |                                                              |                                                                                                                                     |
| ReactMXss, ReactMYss, ReactMZss,      | (Nm), (Nm), (Nm)                                             | at the (0.,0.,-**WtrDpth**) location in SS coordinate system                                                                        |
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| IntfFXss, IntfFYss, IntfFZss,         | (N), (N), (N),                                               | Total interface reaction forces and moments                                                                                         |
|                                       |                                                              |                                                                                                                                     |
| IntfMXss, IntfMYss, IntfMZss,         | (Nm), (Nm), (Nm)                                             | at the TP reference point (platform reference point) location in SS coordinate system                                               |
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| Interface Kinematics                  |                                                                                                                                                                                                    |        
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| IntfTDXss, IntfTDYss, IntfTDZss,      | (m), (m), (m),                                               | Displacements and rotations of the TP reference point                                                                               |
|                                       |                                                              |                                                                                                                                     |
| IntfRDXss, IntfRDYss IntfRDZss        | (rad), (rad), (rad)                                          | (platform reference point) location in SS coordinate system                                                                         |
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| IntfTAXss, IntfTAYss, IntfTAZss,      | (:math:`{m/s^2}`), (:math:`{m/s^2}`), (:math:`{m/s^2}`),     | Translational and rotational accelerations of the TP reference point                                                                |
|                                       |                                                              |                                                                                                                                     |
| IntfRAXss, IntfRAYss IntfRAZss        | (:math:`{rad/s^2}`), (:math:`{rad/s^2}`), (:math:`{rad/s^2}`)| (platform reference point) location in SS coordinate system                                                                         |
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| *Modal Parameters*                    |                                                                                                                                                                                                    |              
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| SSqm01-SSqm99                         | (-)                                                          | C-B modal variables (up to first 99)                                                                                                |
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| SSqmd01-SSqmd99                       | (1/s)                                                        | First time-derivatives of C-B modal variables (up to first 99)                                                                      |
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| SSqmdd01-SSqmdd99                     | (:math:`{1/s^2}`)                                            | Second time-derivatives of C-B modal variables (up to first 99)                                                                     |
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| *Node Kinematics*                     |                                                                                                                                                                                                    |           
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| :math:`{M \alpha N \beta}` TDxss,     | (m)                                                          | Nodal translational displacements of :math:`M \alpha N \beta`							                     |
|					|							       | 																     |
| :math:`M \alpha N \beta` TDyss, 	|							       |  									                                                             |
|					|							       | (up to 81 designated locations) in SS coordinate system									     |
| :math:`M \alpha N \beta` TDzss,	|							       |																     |
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| :math:`{M \alpha N \beta}` RDxe,      | (rad)  						       | Nodal rotational displacements of :math:`M \alpha N \beta`                                                                          |
|					|							       |																     |
| :math:`{M \alpha N \beta}` RDye,	|							       |																     |
|                                       |                                                              | (up to 81 designated locations) in member local coordinate system                                                                   |
| :math:`{M \alpha N \beta}` RDze	|							       |																     |
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| :math:`{M \alpha N \beta}` TAxe,	| (:math:`{m/s^2}`)					       | Nodal translational accelerations  of :math:`M \alpha N \beta`                                                                      |
|					|							       |																     |
| :math:`{M \alpha N \beta}` TAye,	|							       |																     |
|					|							       | (up to 81 designated locations) in member local coordinate system                                                                   |
| :math:`{M \alpha N \beta}` TAze       |                                                              |																     |
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| :math:`{M \alpha N \beta}` RAxe,	| (:math:`{rad/s^2}`)					       | Nodal rotational accelerations  of :math:`M \alpha N \beta`                                                                         |
|					|							       |																     |
| :math:`{M \alpha N \beta}` RAye,	|							       |																     |
|					|							       | (up to 81 designated locations) in member local coordinate system                                                                   |
| :math:`{M \alpha N \beta}` RAze       |                                                              |																     |
+---------------------------------------+--------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------+
| *Node Forces and Moments*             |                                                                                                                                                                                                    |           
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

