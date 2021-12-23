from cx_Freeze import setup, Executable

base = None    

executables = [Executable("Long_track_CA_NEW_GUI.py", base=base)]

packages = ["os", "json", "idna", "tkinter", "subprocess", "psutil"]
options = {
    'build_exe': {    
        'packages':packages,
    },    
}

setup(
    name = "SLM Microstructure Simulation",
    options = options,
    version = "1.0",
    description = 'CA model by Andro, IMT, May 2021',
    executables = executables
)
