


def dfstructure2stfile(dfStructure, H2_stfile):
    """ 
    Write a HAWC2 st file from a pandas dataframe of structural data.

    For fully populated matrices, the columns are:
        ['r_[m]','m_[kg/m]','x_cg_[m]','y_cg_[m]','ri_x_[m]','ri_y_[m]','pitch_[deg]','x_e_[m]','y_e_[m]','K11','K12','K13','K14','K15','K16','K22','K23','K24','K25','K26','K33','K34','K35','K36','K44','K45','K46','K55','K56','K66']
    For Timoshenko:
        ['r_[m]','m_[kg/m]','x_cg_[m]','y_cg_[m]','ri_x_[m]','ri_y_[m]','x_sh_[m]','y_sh_[m]','E_[N/m^2]','G_[N/m^2]','I_x_[m^4]','I_y_[m^4]','I_p_[m^4]','k_x_[-]','k_y_[-]','A_[m^2]','pitch_[deg]','x_e_[m]','y_e_[m]']


    TODO: 
      - make the input independent of column order
      - use hawc2_st_file.py instead

    """

    FPM = 'K11' in dfStructure
    if FPM:
        cols=['r_[m]','m_[kg/m]','x_cg_[m]','y_cg_[m]','ri_x_[m]','ri_y_[m]','pitch_[deg]','x_e_[m]','y_e_[m]','K11','K12','K13','K14','K15','K16','K22','K23','K24','K25','K26','K33','K34','K35','K36','K44','K45','K46','K55','K56','K66']
    else:
        cols=['r_[m]','m_[kg/m]','x_cg_[m]','y_cg_[m]','ri_x_[m]','ri_y_[m]', 'x_sh_[m]','y_sh_[m]','E_[N/m^2]','G_[N/m^2]','I_x_[m^4]','I_y_[m^4]','I_p_[m^4]','k_x_[-]','k_y_[-]','A_[m^2]','pitch_[deg]','x_e_[m]','y_e_[m]']

    with open(H2_stfile, 'w') as f:
        f.write('%i ; number of sets, Nset\n' % 1)
        f.write('-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n')
        f.write('#%i ; set number\n' % 1)
        f.write('\t'.join(['{:20s}'.format(s) for s in cols])+'\n')
        f.write('$%i %i\n' % (1, dfStructure.shape[0])) 
        f.write('\n'.join('\t'.join('%19.13e' %x for x in y) for y in dfStructure.values)) # TODO: this assumes a column order
