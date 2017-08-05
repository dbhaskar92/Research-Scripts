import string
import sys,time
from os import environ
from os import getcwd

sys.path.append(environ["PYTHON_MODULE_PATH"])
sys.path.append(environ["SWIG_LIB_INSTALL_DIR"])

import CompuCellSetup

sim,simthread = CompuCellSetup.getCoreSimulationObjects()     

CompuCellSetup.initializeSimulationObjects(sim,simthread)

steppableRegistry = CompuCellSetup.getSteppableRegistry()

from SmallGTPaseSteppables import ConstraintInitializerSteppable
ConstraintInitializerSteppableInstance = ConstraintInitializerSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(ConstraintInitializerSteppableInstance)

from SmallGTPaseSteppables import RacRhoSignallingClass
RacRhoSignallingInstance = RacRhoSignallingClass(sim,_frequency=1)
steppableRegistry.registerSteppable(RacRhoSignallingInstance)

from SmallGTPaseSteppables import GrowthSteppable
GrowthSteppableInstance = GrowthSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(GrowthSteppableInstance)

from SmallGTPaseSteppables import ExtraFields
extraFieldsInstance=ExtraFields(sim,_frequency=1)
steppableRegistry.registerSteppable(extraFieldsInstance)

from SmallGTPaseSteppables import ExtraPlotSteppable
ExtraPlotSteppableInstance = ExtraPlotSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(ExtraPlotSteppableInstance)

from SmallGTPaseSteppables import XMLDataSteppable
XMLDataSteppableInstance = XMLDataSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(XMLDataSteppableInstance)
        
CompuCellSetup.mainLoop(sim,simthread,steppableRegistry)