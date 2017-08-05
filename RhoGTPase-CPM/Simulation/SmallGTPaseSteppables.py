import os
import sys
import math
import string
import random
import CompuCell
import CompuCellSetup
import numpy as NP
import lxml.etree as ET
from PySteppables import *
from PlayerPython import *
from XMLUtils import dictionaryToMapStrStr as d2mss
from PySteppablesExamples import MitosisSteppableBase

class ConstraintInitializerSteppable(SteppableBasePy):
    
    def __init__(self,_simulator,_frequency=1):
        
        SteppableBasePy.__init__(self,_simulator,_frequency)
        
    def start(self):
        
        for cell in self.cellList:
            
            cell.targetVolume = cell.volume
            cell.targetSurface = 0.0
            cell.lambdaVolume = 1.0
            cell.lambdaSurface = 0.0


class RacRhoSignallingClass(SteppableBasePy):
    
    def __init__(self,_simulator,_frequency):
        
        SteppableBasePy.__init__(self,_simulator,_frequency)
    
    def start(self):
             
        options = {'steps':10, 'stiff':True}
        self.setSBMLGlobalOptions(options)
        modelFile = 'Simulation/SmallGTPase.sbml'
        self.addSBMLToCellTypes(_modelFile=modelFile,_modelName='RR',_types=[self.RED,self.GRAY],_stepSize=0.001)
    
        state={}
    
        for cell in self.cellList:
            
            state['G'] = random.uniform(0,1)  
            state['Gavg'] = state['G']
            state['At'] = cell.targetVolume
            state['beta'] = random.uniform(0.05,0.3)
            state['A'] = cell.volume
            
            self.setSBMLState(_modelName='RR',_cell=cell,_state=state)
            
            cellDict=self.getDictionaryAttribute(cell)
            cellDict['GTPase']=state['G']
        
    def step(self,mcs):
        
        for cell in self.cellList:
            
            G=0.0; nn=0        
            for neighbor , commonSurfaceArea in self.getCellNeighborDataList(cell):
                if neighbor:
                    nn+=1
                    state=self.getSBMLState(_modelName='RR',_cell=neighbor)
                    G+=state['G']
                    
            state=self.getSBMLState(_modelName='RR',_cell=cell)
            gtpase = state['G']
            target_vol = state['At']
            cell_beta = state['beta']
            G+=state['G']
            
            nn+=1
            if (nn > 0):
                G=G/nn
                
            state={}  
            state['Gavg'] = G
            state['G'] = gtpase
            state['At'] = target_vol
            state['beta'] = cell_beta
            state['A'] = cell.volume
            
            self.setSBMLState(_modelName='RR',_cell=cell,_state=state)
            state = self.getSBMLState(_modelName='RR',_cell=cell)
            
            cellDict = self.getDictionaryAttribute(cell)
            cellDict['GTPase'] = state['G']
            
            t_vol = state['At']
            
            if cellDict['GTPase'] < 1.0:
                cell.type = self.GRAY
            else:
                cell.type = self.RED
            
            cell.targetVolume = t_vol
            cell.targetSurface = 0.0
            cell.lambdaVolume = 1.0
            cell.lambdaSurface = 0.0

        self.timestepSBML()
        
        
class ExtraFields(SteppableBasePy):
    
  def __init__(self,_simulator,_frequency=1):
      
    SteppableBasePy.__init__(self,_simulator,_frequency)
    self.scalarFieldG=CompuCellSetup.createScalarFieldCellLevelPy("CellularGTPase")
   
  def step(self,mcs):
      
    self.scalarFieldG.clear()
    
    for cell in self.cellList:
      if cell:
        cellDict=CompuCell.getPyAttrib(cell)
        self.scalarFieldG[cell] = cellDict['GTPase']

    
class GrowthSteppable(SteppableBasePy):
    
    def __init__(self,_simulator,_frequency=1):
        
        SteppableBasePy.__init__(self,_simulator,_frequency)            
    
    def step(self,mcs):
    
        for cell in self.cellList:
            
            state=self.getSBMLState(_modelName='RR',_cell=cell)
            GTPase=state['G']
            cell.lambdaVolume = 1.0
            cell.lambdaSurface = 0.0
            cell.targetVolume = state['At']
            cell.targetSurface = 0.0
            
            DEBUG = 0
            
            if DEBUG == 1:
                if mcs % 10 == 0:
                    print "Cell: ",cell.id," TVolume: ",cell.targetVolume," Volume: ",cell.volume," GTPase: ",GTPase


class XMLDataSteppable(SteppableBasePy):
 
    sim = ET.Element("simulation")
    
    def __init__(self,_simulator,_frequency=1):
        
        SteppableBasePy.__init__(self,_simulator,_frequency)
        self.simulator=_simulator
        self.inventory=self.simulator.getPotts().getCellInventory()
        self.cellList=CellList(self.inventory)

    def start(self):
        
        XMLDataSteppable.sim.set("simulation_type", "CPM_Small_GTPase")
        XMLDataSteppable.sim.set("type_list", "Low_GTP (Gray),High_GTP (Red)")
        XMLDataSteppable.sim.set("axis_division", "none")
        XMLDataSteppable.sim.set("cell_cycle_model", "none")
        XMLDataSteppable.sim.set("time_step", "1")
        
    def step(self,mcs):
        
        time = ET.SubElement(XMLDataSteppable.sim, "time")
        time.set("t", str(mcs))
        
        for cell in self.cellList:
            
            ce = ET.SubElement(time, "cell")
            state=self.getSBMLState(_modelName='RR',_cell=cell)
            
            ce.set("cell_id", str(cell.id))
            ce.set("type", str(cell.type))
            ce.set("x", str(cell.xCOM))
            ce.set("y", str(cell.yCOM))
            ce.set("area", str(cell.volume))
            ce.set("perimeter", str(cell.surface))
            ce.set("pressure", str(state['At']-cell.volume))
            ce.set("target_volume", str(state['At']))
            ce.set("gtpase", str(state['G']))
            
    def finish(self):

        tree = ET.ElementTree(XMLDataSteppable.sim)
        tree.write("/home/dbhaskar92/Downloads/sim_sml_gtpase.xml", pretty_print=True)
        
        
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
        
        # plot title          
        self.pW.setTitle("Rho GTPase Activity")
        self.pW.setTitleSize(11)
        self.pW.setTitleColor("black")
        
        # plot background
        self.pW.setPlotBackgroundColor("white")
        
        # properties of x axis
        self.pW.setXAxisTitle("Monte Carlo Step (MCS)")
        self.pW.setXAxisTitleSize(10)      
        self.pW.setXAxisTitleColor("black")              
        
        # properties of y axis
        self.pW.setYAxisTitle("G")        
        # self.pW.setYAxisLogScale()
        self.pW.setYAxisTitleSize(10)        
        self.pW.setYAxisTitleColor("black")                      
        
        # choices for style are NoCurve, Lines, Sticks, Steps, Dots
        self.pW.addPlot('C1', _style='Lines', _color='green', _size=2)
        self.pW.addPlot('C2', _style='Lines', _color='red', _size=2)
        self.pW.addPlot('C3', _style='Lines', _color='blue', _size=2)
        self.pW.addPlot('C4', _style='Lines', _color='purple', _size=2)
        self.pW.addPlot('C5', _style='Lines', _color='orange', _size=2)
        self.pW.addPlot('C6', _style='Lines', _color='magenta', _size=2)
        self.pW.addPlot('C7', _style='Lines', _color='cyan', _size=2)
        self.pW.addPlot('C8', _style='Lines', _color='brown', _size=2)
        self.pW.addPlot('C9', _style='Lines', _color='gold', _size=2)
     
        self.pW.addGrid()
        self.pW.addAutoLegend("top")
        self.clearFlag=False
    
    def step(self,mcs):
        
        if not self.pW:
            
            print "To get scientific plots working you need extra packages installed:"
            print "Windows/OSX Users: Make sure you have numpy installed. For instructions please visit www.compucell3d.org/Downloads"
            print "Linux Users: Make sure you have numpy and PyQwt installed. Consult manual pages on how to install those packages"
            return
            
        area = 0
        tarea = 0
        gtpase = 0
        
        for cell in self.cellList:
            
            cellDict=CompuCell.getPyAttrib(cell)
            gtpase = cellDict['GTPase']
            
            if cell.id == 1:
                self.pW.addDataPoint("C1", mcs, gtpase)
            elif cell.id ==2:
                self.pW.addDataPoint("C2", mcs, gtpase)
            elif cell.id ==3:
                self.pW.addDataPoint("C3", mcs, gtpase)
            elif cell.id ==4:
                self.pW.addDataPoint("C4", mcs, gtpase)
            elif cell.id ==5:
                self.pW.addDataPoint("C5", mcs, gtpase)
            elif cell.id ==6:
                self.pW.addDataPoint("C6", mcs, gtpase)
            elif cell.id ==7:
                self.pW.addDataPoint("C7", mcs, gtpase)
            elif cell.id ==8:
                self.pW.addDataPoint("C8", mcs, gtpase)
            elif cell.id ==9:
                self.pW.addDataPoint("C9", mcs, gtpase)
        
        self.pW.showAllPlots()
        