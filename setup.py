import cx_Freeze
import sys
import os.path

PYTHON_INSTALL_DIR = os.path.dirname(os.path.dirname(os.__file__))
os.environ['TCL_LIBRARY'] = os.path.join(PYTHON_INSTALL_DIR, 'tcl', 'tcl8.6')
os.environ['TK_LIBRARY'] = os.path.join(PYTHON_INSTALL_DIR, 'tcl', 'tk8.6')

additional_mods = ['numpy.core._methods', 'numpy.lib.format','idna.idnadata']
base = None

if sys.platform == 'win32':
    base = "Win32GUI"

executables = [cx_Freeze.Executable("swiss2go.py", base=base)]

cx_Freeze.setup(
    name = "Swiss2GO App ",
    options = {"build_exe": {"packages":["tkinter","Bio","Bio.Blast","csv", "requests", "time"],'includes':additional_mods,'include_files':[
            os.path.join(PYTHON_INSTALL_DIR, 'DLLs', 'tk86t.dll'),
            os.path.join(PYTHON_INSTALL_DIR, 'DLLs', 'tcl86t.dll'),
         ]}},
    version = "0.01",
    description = "Swiss2GO application",
    executables = executables
    )
