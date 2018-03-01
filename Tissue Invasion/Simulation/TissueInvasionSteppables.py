from PySteppables import *
import CompuCell
import random
import sys
import math
import string
import lxml.etree as ET
from PySteppablesExamples import MitosisSteppableBase
            

class ConstraintInitializerSteppable(SteppableBasePy):
    
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        
    def start(self):
        for cell in self.cellList:
            cell.targetVolume = random.randint(25, 45)
            cell.lambdaVolume = 2.0
            cell.lambdaVecX = 0.0
            cell.lambdaVecY = 0.0
        
    
class GrowthSteppable(SteppableBasePy):
    
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)            
    
    def step(self,mcs):
    
        for cell in self.cellList:
            
            totalMediumSurfacearea = 0.0
            totalcommonSurfacearea = 0.0
            
            for neighbor , commonSurfaceArea in self.getCellNeighborDataList(cell):                
                if neighbor:
                    totalcommonSurfacearea += commonSurfaceArea
                else:
                    totalcommonSurfacearea += commonSurfaceArea
                    totalMediumSurfacearea += commonSurfaceArea
            
            # cell growth is modulated by fraction of perimeter in contact with ECM (contact inhibition)
            sarearatio = totalMediumSurfacearea/totalcommonSurfacearea
            growthfactor = 1.0/(1.0 + math.exp(-2.0*sarearatio))
            
            # apply force if cell is in contact with ECM
            if sarearatio > 0.0 :
                cell.targetVolume += growthfactor*0.01
                # force vector values: -10, -20, -30
                cell.lambdaVecY = -5.0                   
            else:
                cell.targetVolume += growthfactor*0.01
                cell.lambdaVecY = 0.0
            # print "Cell: ",cell.id," Growth: ",growthfactor," TVolume: ",cell.targetVolume," Force: ",cell.lambdaVecY


class MitosisSteppable(MitosisSteppableBase):
    
    def __init__(self,_simulator,_frequency=1):
        MitosisSteppableBase.__init__(self,_simulator, _frequency)
    
    def step(self,mcs):
        cells_to_divide = []
        
        for cell in self.cellList:
            if cell.volume > 50: 
                cells_to_divide.append(cell)
                
        for cell in cells_to_divide:
            # to switch mitosis on/off, leave one of the below lines uncommented
            self.divideCellRandomOrientation(cell)
            # alternative division options:
            # self.divideCellOrientationVectorBased(cell,1,0,0)              
            # self.divideCellAlongMajorAxis(cell)                            
            # self.divideCellAlongMinorAxis(cell)                           

    def updateAttributes(self):
        parentCell=self.mitosisSteppable.parentCell
        childCell=self.mitosisSteppable.childCell
        
        parentCell.targetVolume = random.randint(25,45)
        childCell.targetVolume = random.randint(25,45)
        childCell.lambdaVolume = parentCell.lambdaVolume
        childCell.type = parentCell.type


class DeathSteppable(SteppableBasePy):
    
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        self.scalarCLField = self.createScalarFieldCellLevelPy("cellPressure")
    
    def step(self,mcs):
        # kill cell if it is under too much pressure
        self.scalarCLField.clear()
        
        for cell in self.cellList:
            pressure = cell.targetVolume - cell.volume
            self.scalarCLField[cell] = pressure
            if pressure > 25:
                print "Warning Cell: ",cell.id," Type: ",cell.type," Pressure: ",pressure," TVolume: ",cell.targetVolume


class XMLDataSteppable(SteppableBasePy):
 
    sim = ET.Element("simulation")
    
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        self.simulator=_simulator
        self.inventory=self.simulator.getPotts().getCellInventory()
        self.cellList=CellList(self.inventory)

    def start(self):
        XMLDataSteppable.sim.set("simulation_type", "CPM_Tissue_Invasion")
        XMLDataSteppable.sim.set("type_list", "MCFpIRES (green),MCFmPodo (red)")
        XMLDataSteppable.sim.set("axis_division", "random")
        XMLDataSteppable.sim.set("cell_cycle_model", "contact_inhibited")
        XMLDataSteppable.sim.set("time_step", "1")
        XMLDataSteppable.sim.set("extra_sim_info", "none")
        
    def step(self,mcs):        
        
        time = ET.SubElement(XMLDataSteppable.sim, "time")
        time.set("t", str(mcs))
        
        heterotypicSurfacearea = 0.0
        
        for cell in self.cellList:
            ce = ET.SubElement(time, "cell")
            ce.set("cell_id", str(cell.id))
            ce.set("type", str(cell.type))
            ce.set("x", str(cell.xCOM))
            ce.set("y", str(cell.yCOM))
            ce.set("area", str(cell.volume))
            ce.set("perimeter", str(cell.surface))
            ce.set("pressure", str(cell.targetVolume-cell.volume))
            ce.set("ext_force", str(cell.lambdaVecX)+" "+str(cell.lambdaVecY))
            
            totalMediumSurfacearea = 0.0
            totalcommonSurfacearea = 0.0
            nlist = []
            
            for neighbor , commonSurfaceArea in self.getCellNeighborDataList(cell):                
                if neighbor:
                    nlist.append(neighbor.id)
                    totalcommonSurfacearea += commonSurfaceArea
                    if neighbor.type != cell.type:
                            heterotypicSurfacearea += commonSurfaceArea
                else:
                    totalcommonSurfacearea += commonSurfaceArea
                    totalMediumSurfacearea += commonSurfaceArea
            
            ce.set("contact_perimeter", str(totalMediumSurfacearea))
            ce.set("neighbors", " ".join(map(str, nlist)))
            
            # print "Cell: ",cell.id," Vol: ",cell.volume," SArea: ",cell.surface," COM x: ",xCOM," y: ",yCOM
            time.set("heterotypic_boundary_length", str(heterotypicSurfacearea/2.0))
            
    def finish(self):
        # this function is evaluated after the last MCS
        tree = ET.ElementTree(XMLDataSteppable.sim)
        tree.write("/home/dbhaskar92/Downloads/simulation_ti_results.xml", pretty_print=True)        


class ExtraPlotSteppable(SteppablePy):   
    
    def __init__(self,_simulator,_frequency=1):
        SteppablePy.__init__(self,_frequency)
        self.simulator=_simulator
        self.inventory=self.simulator.getPotts().getCellInventory()
        self.cellList=CellList(self.inventory)

    def start(self):
        import CompuCellSetup  
        
        self.pW=CompuCellSetup.viewManager.plotManager.getNewPlotWindow()
        
        if not self.pW:
            return
        
        #Plot Title - properties           
        self.pW.setTitle("Number of Cells")
        self.pW.setTitleSize(11)
        self.pW.setTitleColor("black")
        
        #plot background
        self.pW.setPlotBackgroundColor("white")
        
        # properties of x axis
        self.pW.setXAxisTitle("Monte Carlo Step (MCS)")
        self.pW.setXAxisTitleSize(10)      
        self.pW.setXAxisTitleColor("black")              
        
        # properties of y axis
        self.pW.setYAxisTitle("Number of Cells")        
        # self.pW.setYAxisLogScale()
        # If you use logscale, it will diverge when one of the two types disappears.
        self.pW.setYAxisTitleSize(10)        
        self.pW.setYAxisTitleColor("black")                      
        
        # choices for style are NoCurve,Lines,Sticks,Steps,Dots
        self.pW.addPlot('MCFpIRES',_style='Lines',_color='green',_size=2)
        self.pW.addPlot('MCFmPodo',_style='Lines',_color='red',_size=2)
        self.pW.addPlot('Total',_style='Lines',_color='blue',_size=2)
     
        self.pW.addGrid()
        self.pW.addAutoLegend("top")
        self.clearFlag=False
    
    def step(self,mcs):
    
        if not self.pW:
            print "To get scientific plots working you need extra packages installed:"
            print "Windows/OSX Users: Make sure you have numpy installed. For instructions please visit www.compucell3d.org/Downloads"
            print "Linux Users: Make sure you have numpy and PyQwt installed. Consult manual pages on how to install those packages"
            return        
        
        numpIRES = 0    # epithelial cells (green)
        numPodo = 0     # mesenchymal cells (red)
        numTotal = 0
        
        for cell in self.cellList: 
            if cell.type == 1:
                numpIRES = numpIRES + 1
            elif cell.type == 2:
                numPodo = numPodo + 1
          
        numTotal = len(self.cellList)
        
        self.pW.addDataPoint("MCFpIRES", mcs, numpIRES)
        self.pW.addDataPoint("MCFmPodo", mcs, numPodo)
        self.pW.addDataPoint("Total", mcs, numTotal)
        
        self.pW.showAllPlots()
