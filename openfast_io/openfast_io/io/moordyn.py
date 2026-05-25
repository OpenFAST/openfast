"""
MoorDynIO – read / write MoorDyn input files.

Produces ``{'MoorDyn': md}``

The MoorDyn file format is section-header-based rather than line-sequential,
so parsing loops on header detection.
"""
from __future__ import annotations

import os
import re
from typing import Any, Dict, Optional, Callable

from .base import ModuleIO
from ..parsing import (
    float_read,
    readline_filterComments,
)


class MoorDynIO(ModuleIO):
    """Read / write MoorDyn input files."""

    # ------------------------------------------------------------------
    def read(
        self,
        file_path: str,
        base_dir: str = '',
        *,
        outlist: Optional[dict] = None,
        read_outlist_fn: Optional[Callable] = None,
        **kwargs,
    ) -> dict:
        md: Dict[str, Any] = {}
        md_file = os.path.normpath(os.path.join(base_dir, file_path)) if base_dir else file_path

        # Init optional headers
        md['Rod_Name'] = []
        md['Body_ID'] = []
        md['Rod_ID'] = []

        f = open(md_file)
        data_line = f.readline()

        while data_line:
            tag = ''.join(data_line.strip().split()).lower()

            # ---- LINE TYPES ----
            if 'linetypes' in tag or 'linedictionary' in tag:
                f.readline(); f.readline()
                for k in ['Name', 'Diam', 'MassDen', 'EA', 'NonLinearEA', 'BA_zeta', 'EI',
                          'Cd', 'Ca', 'CdAx', 'CaAx', 'Cl', 'dF', 'cF']:
                    md[k] = []
                dl = readline_filterComments(f).split()
                while dl[0] and dl[0][:3] != '---':
                    md['Name'].append(str(dl[0]))
                    md['Diam'].append(float(dl[1]))
                    md['MassDen'].append(float(dl[2]))
                    md['EA'].append([float_read(x) for x in dl[3].split('|')])
                    md['BA_zeta'].append([float(x) for x in dl[4].split('|')])
                    md['EI'].append(float(dl[5]))
                    md['Cd'].append(float(dl[6]))
                    md['Ca'].append(float(dl[7]))
                    md['CdAx'].append(float(dl[8]))
                    md['CaAx'].append(float(dl[9]))
                    if len(dl) == 10:
                        md['Cl'].append(None); md['dF'].append(None); md['cF'].append(None)
                    elif len(dl) == 11:
                        md['Cl'].append(float(dl[10])); md['dF'].append(None); md['cF'].append(None)
                    elif len(dl) >= 13:
                        md['Cl'].append(float(dl[10])); md['dF'].append(float(dl[11])); md['cF'].append(float(dl[12]))

                    if isinstance(md['EA'][-1], list) and len(md['EA'][-1]) == 1 and isinstance(md['EA'][-1][0], str):
                        ea_file = os.path.normpath(os.path.join(os.path.dirname(md_file), md['EA'][-1][0]))
                        md['NonLinearEA'].append(_read_nonlinear_ea(ea_file))
                    else:
                        md['NonLinearEA'].append(None)

                    dl = readline_filterComments(f).split()
                data_line = ''.join(dl)

            # ---- ROD TYPES ----
            elif 'rodtypes' in tag or 'roddictionary' in tag:
                f.readline(); f.readline()
                for k in ['Rod_Diam', 'Rod_MassDen', 'Rod_Cd', 'Rod_Ca', 'Rod_CdEnd', 'Rod_CaEnd']:
                    md[k] = []
                dl = readline_filterComments(f).split()
                while dl[0] and dl[0][:3] != '---':
                    md['Rod_Name'].append(dl[0])
                    md['Rod_Diam'].append(float(dl[1]))
                    md['Rod_MassDen'].append(float(dl[2]))
                    md['Rod_Cd'].append(float(dl[3]))
                    md['Rod_Ca'].append(float(dl[4]))
                    md['Rod_CdEnd'].append(float(dl[5]))
                    md['Rod_CaEnd'].append(float(dl[6]))
                    dl = readline_filterComments(f).split()
                data_line = ''.join(dl)

            # ---- BODIES ----
            elif 'bodies' in tag or 'bodylist' in tag or 'bodyproperties' in tag:
                f.readline(); f.readline()
                for k in ['Body_Attachment', 'X0', 'Y0', 'Z0', 'r0', 'p0', 'y0',
                          'Body_Mass', 'Body_CG', 'Body_I', 'Body_Volume', 'Body_CdA', 'Body_Ca']:
                    md[k] = []
                md['Body_ID'] = []
                dl = readline_filterComments(f).split()
                while dl[0] and dl[0][:3] != '---':
                    md['Body_ID'].append(int(dl[0]))
                    md['Body_Attachment'].append(dl[1])
                    md['X0'].append(float(dl[2])); md['Y0'].append(float(dl[3])); md['Z0'].append(float(dl[4]))
                    md['r0'].append(float(dl[5])); md['p0'].append(float(dl[6])); md['y0'].append(float(dl[7]))
                    md['Body_Mass'].append(float(dl[8]))
                    md['Body_CG'].append([float(x) for x in dl[9].split('|')])
                    md['Body_I'].append([float(x) for x in dl[10].split('|')])
                    md['Body_Volume'].append(float(dl[11]))
                    md['Body_CdA'].append([float(x) for x in dl[12].split('|')])
                    md['Body_Ca'].append([float(x) for x in dl[13].split('|')])
                    dl = readline_filterComments(f).split()
                data_line = ''.join(dl)

            # ---- RODS ----
            elif 'rods' in tag or 'rodlist' in tag or 'rodproperties' in tag:
                f.readline(); f.readline()
                for k in ['Rod_Type', 'Rod_Attachment', 'Xa', 'Ya', 'Za', 'Xb', 'Yb', 'Zb', 'Rod_NumSegs', 'RodOutputs']:
                    md[k] = []
                md['Rod_ID'] = []
                dl = readline_filterComments(f).split()
                while dl[0] and dl[0][:3] != '---':
                    md['Rod_ID'].append(dl[0])
                    md['Rod_Type'].append(dl[1])
                    md['Rod_Attachment'].append(dl[2])
                    md['Xa'].append(float(dl[3])); md['Ya'].append(float(dl[4])); md['Za'].append(float(dl[5]))
                    md['Xb'].append(float(dl[6])); md['Yb'].append(float(dl[7])); md['Zb'].append(float(dl[8]))
                    md['Rod_NumSegs'].append(int(dl[9]))
                    md['RodOutputs'].append(dl[10])
                    dl = readline_filterComments(f).split()
                data_line = ''.join(dl)

            # ---- POINTS ----
            elif 'points' in tag or 'connectionproperties' in tag or \
                 'nodeproperties' in tag or 'pointproperties' in tag or 'pointlist' in tag:
                f.readline(); f.readline()
                for k in ['Point_ID', 'Attachment', 'X', 'Y', 'Z', 'M', 'V', 'CdA', 'CA']:
                    md[k] = []
                dl = readline_filterComments(f).split()
                while dl[0] and dl[0][:3] != '---':
                    md['Point_ID'].append(int(dl[0]))
                    md['Attachment'].append(str(dl[1]))
                    md['X'].append(float(dl[2])); md['Y'].append(float(dl[3])); md['Z'].append(float(dl[4]))
                    md['M'].append(float(dl[5])); md['V'].append(float(dl[6]))
                    md['CdA'].append(float(dl[7])); md['CA'].append(float(dl[8]))
                    dl = readline_filterComments(f).split()
                data_line = ''.join(dl)

            # ---- LINES ----
            elif 'lines' in tag or 'lineproperties' in tag or 'linelist' in tag:
                f.readline(); f.readline()
                for k in ['Line_ID', 'LineType', 'AttachA', 'AttachB', 'UnstrLen', 'NumSegs', 'Outputs']:
                    md[k] = []
                dl = readline_filterComments(f).split()
                while dl[0] and dl[0][:3] != '---':
                    md['Line_ID'].append(int(dl[0]))
                    md['LineType'].append(str(dl[1]))
                    md['AttachA'].append(str(dl[2])); md['AttachB'].append(str(dl[3]))
                    md['UnstrLen'].append(float(dl[4]))
                    md['NumSegs'].append(int(dl[5]))
                    md['Outputs'].append(str(dl[6]))
                    dl = readline_filterComments(f).split()
                data_line = ''.join(dl)

            # ---- FAILURE ----
            elif 'failure' in tag:
                f.readline(); f.readline()
                for k in ['Failure_ID', 'Failure_Point', 'Failure_Line(s)', 'FailTime', 'FailTen']:
                    md[k] = []
                dl = readline_filterComments(f).split()
                while dl[0] and dl[0][:3] != '---':
                    md['Failure_ID'].append(int(dl[0]))
                    md['Failure_Point'].append(dl[1])
                    md['Failure_Line(s)'].append([int(x) for x in dl[2].split(',')])
                    md['FailTime'].append(float(dl[3]))
                    md['FailTen'].append(float(dl[4]))
                    dl = readline_filterComments(f).split()
                data_line = ''.join(dl)

            # ---- CONTROL ----
            elif 'control' in tag:
                f.readline(); f.readline()
                md['ChannelID'] = []
                md['Lines_Control'] = []
                dl = readline_filterComments(f).split()
                while dl[0] and dl[0][:3] != '---':
                    md['ChannelID'].append(int(dl[0]))
                    control_lines = []
                    for lines in dl[1:]:
                        for line in lines.split(','):
                            control_lines.append(line.strip(','))
                    while '' in control_lines:
                        control_lines.remove('')
                    md['Lines_Control'].append(control_lines)
                    dl = readline_filterComments(f).split()
                data_line = ''.join(dl)

            # ---- EXTERNAL ----
            elif 'external' in tag:
                f.readline(); f.readline()
                for k in ['External_ID', 'Object', 'Fext', 'Blin', 'Bquad', 'CSys']:
                    md[k] = []
                dl = readline_filterComments(f).split()
                while dl[0] and dl[0][:3] != '---':
                    md['External_ID'].append(int(dl[0]))
                    md['Object'].append(dl[1])
                    md['Fext'].append([float(x) for x in dl[2].split('|')])
                    md['Blin'].append([float(x) for x in dl[3].split('|')])
                    md['Bquad'].append([float(x) for x in dl[4].split('|')])
                    md['CSys'].append(dl[5])
                    dl = readline_filterComments(f).split()
                data_line = ''.join(dl)

            # ---- OPTIONS ----
            elif 'options' in tag:
                md['option_values'] = []
                md['option_names'] = []
                md['option_descriptions'] = []
                dl = readline_filterComments(f).split()
                while dl[0] and dl[0][:3] != '---':
                    option_value = dl[0].upper()
                    option_name = dl[1].upper()
                    option_description = ' '.join(dl[2:]) if len(dl) > 2 else '-'
                    if option_name == 'WATERKIN':
                        md['WaterKin'] = option_value.strip('"')
                    md['option_values'].append(float_read(option_value.strip('"')))
                    md['option_names'].append(option_name)
                    md['option_descriptions'].append(option_description)
                    dl = readline_filterComments(f).split()
                data_line = ''.join(dl)

            # ---- OUTPUTS ----
            elif 'outputs' in tag:
                outlist_md = {}
                dl = readline_filterComments(f)
                while (dl and dl[0:3] != '---') and ('END' not in dl):
                    if '"' in dl:
                        # Strip surrounding quotes from each channel token
                        dl = dl.replace('"', '')
                    channels = [c.strip() for c in dl.split(',')]
                    for c in channels:
                        if c:
                            outlist_md[c] = True
                    dl = readline_filterComments(f)
                if outlist is not None:
                    outlist['MoorDyn'] = outlist_md
                md['_outlist'] = outlist_md
                f.close()
                break

            else:
                data_line = f.readline()

        return {'MoorDyn': md}

    # ------------------------------------------------------------------
    def write(
        self,
        data: dict,
        file_path: str,
        base_dir: str = '',
        *,
        outlist: Optional[dict] = None,
        **kwargs,
    ) -> None:
        md = data['MoorDyn']

        with open(file_path, 'w') as f:
            f.write('--------------------- MoorDyn Input File ------------------------------------\n')
            f.write('Generated with OpenFAST_IO\n')

            # Line Types
            if md.get('Name'):
                f.write('----------------------- LINE TYPES ------------------------------------------\n')
                f.write('Name    Diam    MassDen    EA    BA/-zeta    EI    Cd    Ca    CdAx    CaAx\n')
                f.write('(-)     (m)     (kg/m)     (N)   (N-s/-)    (-)   (-)   (-)   (-)     (-)\n')
                for i in range(len(md['Name'])):
                    ea_str = '|'.join([str(v) for v in md['EA'][i]])
                    ba_str = '|'.join([str(v) for v in md['BA_zeta'][i]])
                    ln = f"{md['Name'][i]}   {md['Diam'][i]}   {md['MassDen'][i]}   {ea_str}   {ba_str}   {md['EI'][i]}   {md['Cd'][i]}   {md['Ca'][i]}   {md['CdAx'][i]}   {md['CaAx'][i]}"
                    if md.get('Cl') and md['Cl'][i] is not None:
                        ln += f"   {md['Cl'][i]}"
                    if md.get('dF') and md['dF'][i] is not None:
                        ln += f"   {md['dF'][i]}   {md['cF'][i]}"
                    f.write(ln + '\n')

            # Rod Types
            if md.get('Rod_Name'):
                f.write('----------------------- ROD TYPES -------------------------------------------\n')
                f.write('TypeName  Diam  MassDenInAir  Cd  Ca  CdEnd  CaEnd\n')
                f.write('(name)    (m)   (kg/m)        (-) (-) (-)    (-)\n')
                for i in range(len(md['Rod_Name'])):
                    f.write(f"{md['Rod_Name'][i]}   {md['Rod_Diam'][i]}   {md['Rod_MassDen'][i]}   {md['Rod_Cd'][i]}   {md['Rod_Ca'][i]}   {md['Rod_CdEnd'][i]}   {md['Rod_CaEnd'][i]}\n")

            # Bodies
            if md.get('Body_ID'):
                f.write('----------------------- BODIES ----------------------------------------------\n')
                f.write('ID  Attachment  X0  Y0  Z0  r0  p0  y0  Mass  CG  I  Volume  CdA  Ca\n')
                f.write('(-) (-)         (m) (m) (m) (deg)(deg)(deg)(kg) (m) (kg-m^2) (m^3) (m^2) (-)\n')
                for i in range(len(md['Body_ID'])):
                    cg = '|'.join([str(v) for v in md['Body_CG'][i]])
                    bi = '|'.join([str(v) for v in md['Body_I'][i]])
                    cda = '|'.join([str(v) for v in md['Body_CdA'][i]])
                    ca = '|'.join([str(v) for v in md['Body_Ca'][i]])
                    f.write(f"{md['Body_ID'][i]}  {md['Body_Attachment'][i]}  {md['X0'][i]}  {md['Y0'][i]}  {md['Z0'][i]}  {md['r0'][i]}  {md['p0'][i]}  {md['y0'][i]}  {md['Body_Mass'][i]}  {cg}  {bi}  {md['Body_Volume'][i]}  {cda}  {ca}\n")

            # Rods
            if md.get('Rod_ID'):
                f.write('----------------------- RODS ------------------------------------------------\n')
                f.write('ID  RodType  Attachment  Xa  Ya  Za  Xb  Yb  Zb  NumSegs  RodOutputs\n')
                f.write('(-) (-)      (-)         (m) (m) (m) (m) (m) (m) (-)      (-)\n')
                for i in range(len(md['Rod_ID'])):
                    f.write(f"{md['Rod_ID'][i]}  {md['Rod_Type'][i]}  {md['Rod_Attachment'][i]}  {md['Xa'][i]}  {md['Ya'][i]}  {md['Za'][i]}  {md['Xb'][i]}  {md['Yb'][i]}  {md['Zb'][i]}  {md['Rod_NumSegs'][i]}  {md['RodOutputs'][i]}\n")

            # Points
            if md.get('Point_ID'):
                f.write('----------------------- POINTS -----------------------------------------------\n')
                f.write('ID  Attachment  X     Y     Z     M     V    CdA    CA\n')
                f.write('(-) (-)         (m)   (m)   (m)   (kg)  (m^3) (m^2) (-)\n')
                for i in range(len(md['Point_ID'])):
                    f.write(f"{md['Point_ID'][i]}   {md['Attachment'][i]}   {md['X'][i]}   {md['Y'][i]}   {md['Z'][i]}   {md['M'][i]}   {md['V'][i]}   {md['CdA'][i]}   {md['CA'][i]}\n")

            # Lines
            if md.get('Line_ID'):
                f.write('----------------------- LINES -----------------------------------------------\n')
                f.write('ID  LineType  AttachA  AttachB  UnstrLen  NumSegs  LineOutputs\n')
                f.write('(-) (-)       (-)      (-)      (m)       (-)      (-)\n')
                for i in range(len(md['Line_ID'])):
                    f.write(f"{md['Line_ID'][i]}   {md['LineType'][i]}   {md['AttachA'][i]}   {md['AttachB'][i]}   {md['UnstrLen'][i]}   {md['NumSegs'][i]}   {md['Outputs'][i]}\n")

            # Failure
            if md.get('Failure_ID'):
                f.write('----------------------- FAILURE ----------------------------------------------\n')
                f.write('ID  Point  Line(s)  FailTime  FailTen\n')
                f.write('(-) (-)    (-)      (s)       (N)\n')
                for i in range(len(md['Failure_ID'])):
                    lines_str = ','.join([str(x) for x in md['Failure_Line(s)'][i]])
                    f.write(f"{md['Failure_ID'][i]}   {md['Failure_Point'][i]}   {lines_str}   {md['FailTime'][i]}   {md['FailTen'][i]}\n")

            # Control
            if md.get('ChannelID'):
                f.write('----------------------- CONTROL -----------------------------------------------\n')
                f.write('ChannelID  Line(s)\n')
                f.write('(-)        (-)\n')
                for i in range(len(md['ChannelID'])):
                    cl = ','.join(md['Lines_Control'][i])
                    f.write(f"{md['ChannelID'][i]}   {cl}\n")

            # External
            if md.get('External_ID'):
                f.write('----------------------- EXTERNAL -----------------------------------------------\n')
                f.write('ID  Object  Fext  Blin  Bquad  CSys\n')
                f.write('(-) (-)     (N)   (N/(m/s))  (N/(m/s)^2)  (-)\n')
                for i in range(len(md['External_ID'])):
                    fext = '|'.join([str(v) for v in md['Fext'][i]])
                    blin = '|'.join([str(v) for v in md['Blin'][i]])
                    bquad = '|'.join([str(v) for v in md['Bquad'][i]])
                    f.write(f"{md['External_ID'][i]}   {md['Object'][i]}   {fext}   {blin}   {bquad}   {md['CSys'][i]}\n")

            # Options
            if md.get('option_names'):
                f.write('----------------------- OPTIONS -----------------------------------------------\n')
                for i in range(len(md['option_names'])):
                    val = md['option_values'][i]
                    name = md['option_names'][i]
                    desc = md['option_descriptions'][i]
                    f.write(f"{val}   {name}   {desc}\n")

            # Outputs
            f.write('----------------------- OUTPUTS ------------------------------------------------\n')
            out = md.get('_outlist', {})
            if outlist and 'MoorDyn' in outlist:
                out = outlist['MoorDyn']
            for ch in out:
                f.write(f'"{ch}"\n')
            f.write('END\n')
            f.write('----------------------- need this line ------------------\n')


def _read_nonlinear_ea(filepath):
    """Read a non-linear EA file (strain-EA curves)."""
    try:
        import numpy as np
        data = np.loadtxt(filepath)
        return data
    except Exception:
        return None
