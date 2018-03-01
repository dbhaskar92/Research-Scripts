import sys
from os import environ
from os import getcwd
import string

sys.path.append(environ["PYTHON_MODULE_PATH"])

import CompuCellSetup

sim,simthread = CompuCellSetup.getCoreSimulationObjects()      

CompuCellSetup.initializeSimulationObjects(sim,simthread)

steppableRegistry = CompuCellSetup.getSteppableRegistry()

from TissueInvasionSteppables import ConstraintInitializerSteppable
ConstraintInitializerSteppableInstance = ConstraintInitializerSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(ConstraintInitializerSteppableInstance)

from TissueInvasionSteppables import GrowthSteppable
GrowthSteppableInstance = GrowthSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(GrowthSteppableInstance)

from TissueInvasionSteppables import MitosisSteppable
MitosisSteppableInstance = MitosisSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(MitosisSteppableInstance)

from TissueInvasionSteppables import DeathSteppable
DeathSteppableInstance = DeathSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(DeathSteppableInstance)

from TissueInvasionSteppables import ExtraPlotSteppable
ExtraPlotSteppableInstance = ExtraPlotSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(ExtraPlotSteppableInstance)

from TissueInvasionSteppables import XMLDataSteppable
XMLDataSteppableInstance = XMLDataSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(XMLDataSteppableInstance)
        
CompuCellSetup.mainLoop(sim,simthread,steppableRegistry)
