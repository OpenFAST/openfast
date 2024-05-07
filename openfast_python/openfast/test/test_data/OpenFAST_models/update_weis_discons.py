'''
Update the DISCON.IN examples in the WEIS repository using the Tune_Case/ .yaml files

'''
import os
from rosco.toolbox.ofTools.fast_io.update_discons import update_discons

weis_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
iea15_dir = os.path.join(weis_dir,'examples/01_aeroelasticse/OpenFAST_models/IEA-15-240-RWT/')

if __name__=="__main__":

    # {tune_yaml:discon}
    discon_map = {
        os.path.join(iea15_dir,'IEA-15-240-RWT-UMaineSemi/IEA15MW-UMaineSemi.yaml'): os.path.join(iea15_dir,'IEA-15-240-RWT-UMaineSemi/DISCON-UMaineSemi.IN'),
        os.path.join(iea15_dir,'IEA-15-240-RWT-Monopile/IEA15MW-Monopile.yaml'): os.path.join(iea15_dir,'IEA-15-240-RWT-Monopile/DISCON-Monopile.IN'),
        # 'IEA15MW.yaml': 'IEA-15-240-RWT-UMaineSemi/DISCON-UMaineSemi.IN',
        # 'BAR.yaml': 'BAR_10/BAR_10_DISCON.IN'
    }

    # # Directories
    # test_dir = os.path.dirname(os.path.abspath(__file__))
    # tune_dir = os.path.realpath(os.path.join(test_dir,'../Tune_Cases'))

    # # Make paths absolute
    # map_abs = {}
    # for tune, test in map_rel.items():
    #     tune = os.path.join(tune_dir,tune)
    #     map_abs[tune] = os.path.join(test_dir,test)

    # Make discons    
    update_discons(discon_map)
