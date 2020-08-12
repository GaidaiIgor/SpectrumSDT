# Author:  SLEPc Team
# Contact: slepc-maint@upv.es

def get_slepc_dir():
    import os
    return os.path.dirname(__file__)

def get_config():
    conf = {}
    conf['SLEPC_DIR'] = get_slepc_dir()
    return conf
