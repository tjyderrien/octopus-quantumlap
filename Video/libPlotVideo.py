##################################################################################
## Copyright (C) 2016 M. Compostella, M. Rampp, H. Appel
## Modified by TJY Derrien, July 2016
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA.
##
##################################################################################
## General Notes:
## This is a python script for the visualization sotfware VisIt.
## Submit the script with "visit -cli -s script_name.py"
## or "visit -cli -nowin -s script.py" for a version without plotting windows
##
## Last update: 29.06.2016
## Tested on: VisIt 2.10 (patch 2)
## at: Max Planck Computing and Data Facility, Garching
##################################################################################


# This script creates a multi-windows 3D plot of a laser pulse propagating towards Azulene.
# TODO: add vertical red bar on the laser graph (difficult! not included into Visit yet.) 
# TODO: convert units for input the laser intensity to W/cm2 ! V/m is awful to manage. 

# Variables that you can set before running the script:
#
# inputdata     - File with the FULL name (ending with "/") of the directories containing the data that you
#                 want to plot. Each line is an input step for the plot.
#
# outdir        - Directory where the .png images are saved.
#
# elfvtk        - Name of the .vtk file containing the Octopus results. The script assumes a variable called 'Re' is
#                 present in this file.
#
# geometry      - Name of the file containing the geometry of the molecule (e.g. geometry.xyz).
#
# maxe_x        - Name of the .vtk file containing the values for the electric field along the x axis.
#                 The script assumes the variable named 'Re' is present in this file.
#
# maxe_y        - Name of the .vtk file containing the values for the electric field along the y axis.
#                 The script assumes the variable named 'Re' is present in this file.
#
# maxe_z        - Name of the .vtk file containing the values for the electric field along the z axis.
#                 The script assumes the variable named 'Re' is present in this file.
#
# bonds         - Set to True if you want to plot the bonds of the molecule. Note that bonds are created according to the
#                 distance between the elements of the molecule.
#


# Set the variable below as DEFAULT values
# First row of <inputdata> should be equal to the ground state data. 
# All other time steps will be normalized to this one. 
# import scipy.constants as phys

# We define some functions usefull for all windows

""" ============ FUNCTIONS ORIGINALLY DEFINED BY MPSD / MPG =================
"""

BohrRadius=5.2917721067e-1 #Angstroms
# TODO: replace with link to the Python library maybe.
# BohrRadius = phys.physical_constants['Bohr radius']
# print "Bohr radius"
# print BohrRadius

def extension(dataset, var):     # this function creates a temporary window, calculates the extensions of a dataset and destroys the window
    AddWindow()
    ReOpenDatabase(dataset)
    AddPlot("Pseudocolor", var)
    DrawPlots()
    Query("SpatialExtents")
    ext = GetQueryOutputValue()
    DeleteWindow()
    return ext

def rescmol(reszoom):
    AddOperator("Transform" , 0)
    res = TransformAttributes()
    res.doScale = 1
    res.scaleOrigin = (0, 0, 0)
    res.scaleX = reszoom
    res.scaleY = reszoom
    res.scaleZ = reszoom
    SetOperatorOptions(res)

def createbonds():
    AddOperator("CreateBonds", 0)
    bonds = CreateBondsAttributes()
    bonds.minDist = (0.4, 0.4)
    bonds.maxDist = (3, 4)
    SetOperatorOptions(bonds)
    
    # Defining an arrow if we are going to zoom on the molecule
#if (molzoom>1):
    #ar1 = CreateAnnotationObject("Line2D")
    #ar1.visible = 1
    #ar1.position  = (0.543, 0.52)
    #ar1.position2 = (0.80, 0.88)
    #ar1.width = 2
    #ar1.beginArrow = 2
    #ar1.useForegroundForLineColor = 0
    #ar1.color = (160, 160, 160, 255)

def winatt1():
    w1 = View3DAttributes()
    w1.viewNormal = (-0.56, 0.36, -0.75)
    w1.viewUp = (0.26, 0.93, 0.25)
    w1.viewAngle = 30
    w1.imageZoom = 1.1
    # The following lines should work for all molzoom, fixing the view to the scale of the electric field
    w1.focus = (0, 0, (ext1[2]+ext1[4])/10.)
    w1.parallelScale = ext1[5]
    w1.nearPlane = 2.*ext1[4]
    w1.farPlane = 2.*ext1[5]
    # The following lines work if you choose a molzoom that is consistent with the scale of the electric field
    #w1.focus = (0, 0, -20*molzoom)
    #w1.parallelScale = 170.*molzoom
    #w1.nearPlane = -340.*molzoom
    #w1.farPlane = 340.*molzoom
    SetView3D(w1)
    
    
    
 ## == ADDINGS BY THIBAULT JY DERRIEN
 
 #def plotVerticalBar(varX=step, colorcode=(0, 0, 0, 255)): #{{{
  #OpenDatabase(filepath)
  #AddPlot("Scatter", varY, 1, 1)
  #ScatterAtts = ScatterAttributes()
  #ScatterAtts.var1 = varX
  #ScatterAtts.var1Role = ScatterAtts.Coordinate0
  #ScatterAtts.var2 = varY
  #ScatterAtts.var2Role = ScatterAtts.Coordinate1
  #ScatterAtts.colorType = ScatterAtts.ColorBySingleColor
  #ScatterAtts.singleColor = colorcode
  #ScatterAtts.pointSize = 0.05
  #ScatterAtts.pointSizePixels = 3
  #ScatterAtts.pointType = ScatterAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
  #SetDefaultPlotOptions(ScatterAtts)
  #SetPlotOptions(ScatterAtts)
  #DrawPlots()
  ##}}}
from visit import *
from visit_utils import *
#from visit_utils.qplot import *


def AddTracerObject(filepath, time): #{{{
  # Plots a line on a plot to mark the current time
  # /!\ This uses PropertyTree as the simple Python interface is not yet developed in VisIt 2.10.3. 
  # See: http://www.visitusers.org/index.php?title=Visit-tutorial-Advanced-movie-making#Creating_the_curve_images_over_time
  # and: http://fossies.org/dox/visit2.10.3/plots_8py_source.html
  # NOTE: Under development. 
  # filepath = "laser.txt"
  
  n_curves  = 1

  curves = [PropertyTree() for i in xrange(n_curves)]
  curves[0].file = filepath
  curves[0].index = 0

  p=PropertyTree()
  p.plots = [PropertyTree() for i in xrange(n_curves)]
  p.plots[0].type = "line"
  p.plots[0].curve = curves[0]
  tracer_line = PropertyTree()
  tracer_line.curve = [1]
  tracer_line.type = "tracer_line"
  tracer_line.color = (0, 185, 50, 255)
  tracer_line.tracer_x = time #put here your time
  p.plots.append(tracer_line)
  scene = CurveScene(p)
#}}}

def plot1D(filepath, varX, varY, colorcode=(255, 0, 0, 255)): #{{{
  OpenDatabase(filepath)
  AddPlot("Scatter", varY, 1, 1)
  ScatterAtts = ScatterAttributes()
  ScatterAtts.var1 = varX
  ScatterAtts.var1Role = ScatterAtts.Coordinate0
  ScatterAtts.var2 = varY
  ScatterAtts.var2Role = ScatterAtts.Coordinate1
  ScatterAtts.colorType = ScatterAtts.ColorBySingleColor
  ScatterAtts.singleColor = colorcode
  ScatterAtts.pointSize = 0.05
  ScatterAtts.pointSizePixels = 3
  ScatterAtts.pointType = ScatterAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
  SetDefaultPlotOptions(ScatterAtts)
  SetPlotOptions(ScatterAtts)
  DrawPlots()
  #}}}

def plot2Ddifference(inxyz, indensityvtk, step, caption=1, OriginX=0., OriginY=0., OriginZ=0., axis="X", extreme=0.03, molvar="element"):#{{{
  SetViewExtentsType(1)	
  ReOpenDatabase(indensityvtk)
  AddPlot("Pseudocolor", "Difference", 1, 1)
      
  #Defining the color scale
  PseudocolorAtts = PseudocolorAttributes()
  PseudocolorAtts.colorTableName = "difference" #"PuBu"
  PseudocolorAtts.minFlag = 1
  PseudocolorAtts.maxFlag = 1
  PseudocolorAtts.legendFlag = caption
  PseudocolorAtts.min = -extreme
  PseudocolorAtts.max = extreme
  SetPlotOptions(PseudocolorAtts)
  #DrawPlots()
  #SetActivePlots(1)
      
  ## Setting box size
  #AddOperator("Box", 1)
  #BoxAtts = BoxAttributes()
  #BoxAtts.amount = BoxAtts.Some  # Some, All
  #BoxAtts.minx = -10
  #BoxAtts.maxx = 50
  #BoxAtts.miny = -40
  #BoxAtts.maxy = 40
  #BoxAtts.minz = -10
  #BoxAtts.maxz = 50
  #BoxAtts.inverse = 0
  #SetOperatorOptions(BoxAtts, 1)

  # Trying to restrict viewing to the box
  SetViewExtentsType(1)

  ## Adding molecules
  #SetActiveWindow(1)
  #SetActivePlots((0,1))
  ReOpenDatabase(inxyz)
  AddPlot("Molecule", molvar, 1, 1)
  mol1 = MoleculeAttributes()
  mol1.drawAtomsAs = mol1.SphereAtoms  # NoAtoms, SphereAtoms, ImposterAtoms
  mol1.colorBonds = mol1.ColorByAtom  # ColorByAtom, SingleColor
  mol1.bondSingleColor = (128, 128, 128, 255)
  mol1.scaleRadiusBy = mol1.Fixed # Fixed, Covalent, Atomic, Variable
  mol1.radiusFixed = 0.3 #0.3 ?
  mol1.radiusVariable = "default"
  mol1.radiusScaleFactor = 1.0
  mol1.drawBondsAs = mol1.CylinderBonds
  mol1.bondRadius = 0.12
  mol1.bondLineWidth = 0
  mol1.legendFlag = caption
  mol1.elementColorTable = "cpk_jmol"
  mol1.residueTypeColorTable = "amino_shapely"
  mol1.atomSphereQuality = mol1.Super
  mol1.bondLineStyle = mol1.SOLID
  mol1.residueSequenceColorTable = "Default"
  mol1.continuousColorTable = "Default"
  mol1.bondCylinderQuality = mol1.Super
  #MoleculeAtts.legendFlag = 1
  #MoleculeAtts.minFlag = 0
  #MoleculeAtts.scalarMin = 0
  #MoleculeAtts.maxFlag = 0
  #MoleculeAtts.scalarMax = 1
  SetPlotOptions(mol1)

  # Build a 2D slide from 3D system
  SetActivePlots((0, 1))
  #SetActivePlots(0)
  AddOperator("Slice", 1)
  SliceAtts = SliceAttributes()
  if(axis == "X"): #Aligned with H bonds :) 
    SliceAtts.normal = (1, 0, 0)
    SliceAtts.axisType = SliceAtts.XAxis  # XAxis, YAxis, ZAxis, Arbitrary, ThetaPhi
    SliceAtts.upAxis = (0, 0, 1)
    SliceAtts.flip = 1
    SliceAtts.theta = 90
    SliceAtts.phi = 0
    SliceAtts.project2d = 1
  
  if(axis == "Y"): #Not aligned with Si chemical bonds :-(
    SliceAtts.normal = (0, -1, 0)
    SliceAtts.axisType = SliceAtts.YAxis  # XAxis, YAxis, ZAxis, Arbitrary, ThetaPhi
    SliceAtts.upAxis = (0, 0, 1)
    SliceAtts.flip = 0
    SliceAtts.theta = 180
    SliceAtts.phi = 0
    SliceAtts.project2d = 1
    
  if(axis == "SiBond"): #Aligned with Si chemical bonds :)
    SliceAtts.normal = (1, 1, 1)
    SliceAtts.originType = SliceAtts.Point  # Point, Intercept, Percent, Zone, Node
    SliceAtts.originPoint = (OriginX, OriginY, OriginZ)
    SliceAtts.axisType = SliceAtts.Arbitrary  # XAxis, YAxis, ZAxis, Arbitrary, ThetaPhi
    SliceAtts.upAxis = (0, 0, 1)
    SliceAtts.flip = 0
    SliceAtts.theta = 0
    SliceAtts.phi = 0
    SliceAtts.project2d = 0
    
  if(axis == "Z"):
    SliceAtts.normal = (0, 0, -1)
    SliceAtts.axisType = SliceAtts.ZAxis  # XAxis, YAxis, ZAxis, Arbitrary, ThetaPhi
    SliceAtts.upAxis = (1, 0, 0)
    SliceAtts.flip = 1
    SliceAtts.theta = 0
    SliceAtts.phi = -90
    SliceAtts.project2d = 1
  
  SliceAtts.interactive = 1
  SetOperatorOptions(SliceAtts, 1)

  #Adding a label for the time step    
  l1s = CreateAnnotationObject("Text2D")
  l1s.text = "step # "+step
  l1s.position = (0.05, 0.05)
  
  # Drawing plot to run query
  DrawPlots()
  #}}}
  
def plot3Ddifference(inxyz, indensityvtk, step, caption=1, OriginX=0., OriginY=0., OriginZ=0., axis="X", extreme=0.01, ZoomValue=1.0, Enable_Perspective=0, molvar="element"):#{{{
  SetViewExtentsType(1)
  ## Adding molecules
  #SetActiveWindow(1)
  #SetActivePlots((0,1))
  ReOpenDatabase(inxyz)
  AddPlot("Molecule", molvar, 1, 1)
  mol1 = MoleculeAttributes()
  mol1.drawAtomsAs = mol1.SphereAtoms  # NoAtoms, SphereAtoms, ImposterAtoms
  mol1.colorBonds = mol1.ColorByAtom  # ColorByAtom, SingleColor
  mol1.bondSingleColor = (128, 128, 128, 255)
  mol1.scaleRadiusBy = mol1.Fixed # Fixed, Covalent, Atomic, Variable
  mol1.radiusFixed = 0.3 #0.3 ?
  mol1.radiusVariable = "default"
  mol1.radiusScaleFactor = 1.0
  mol1.drawBondsAs = mol1.CylinderBonds
  mol1.bondRadius = 0.12
  mol1.bondLineWidth = 0
  mol1.elementColorTable = "cpk_jmol"
  mol1.residueTypeColorTable = "amino_shapely"
  mol1.atomSphereQuality = mol1.Super
  mol1.bondLineStyle = mol1.SOLID
  mol1.residueSequenceColorTable = "Default"
  mol1.continuousColorTable = "Default"
  mol1.bondCylinderQuality = mol1.Super
  #MoleculeAtts.legendFlag = 1
  #MoleculeAtts.minFlag = 0
  #MoleculeAtts.scalarMin = 0
  #MoleculeAtts.maxFlag = 0
  #MoleculeAtts.scalarMax = 1
  SetPlotOptions(mol1)
  
  ReOpenDatabase(indensityvtk)
  AddPlot("Contour", "Difference", 1, 1)
  #SetActivePlots((0, 1))
  ContourAtts = ContourAttributes()
  ContourAtts.contourNLevels = 6
  ContourAtts.contourMethod = ContourAtts.Level  # Level, Value, Percent
  ContourAtts.minFlag = 1
  ContourAtts.maxFlag = 1
  ContourAtts.min = -extreme
  ContourAtts.max = extreme
  # 3 color map (depletion: blue, overdensity: red)
  #ContourAtts.SetMultiColor(0, (25, 94, 255, 255))
  #ContourAtts.SetMultiColor(1, (255, 255, 255, 0))
  #ContourAtts.SetMultiColor(2, (255, 0, 0, 255))
  # 5 color map (depletion: blue, overdensity: red)
  #ContourAtts.SetMultiColor(0, (51, 102, 255, 89))
  #ContourAtts.SetMultiColor(1, (153, 204, 255, 57))
  #ContourAtts.SetMultiColor(2, (255, 255, 255, 0))
  #ContourAtts.SetMultiColor(3, (255, 153, 0, 58))
  #ContourAtts.SetMultiColor(4, (255, 0, 0, 124))
  # 6 color map
  ContourAtts.SetMultiColor(0, (51, 102, 255, 89))
  ContourAtts.SetMultiColor(1, (153, 204, 255, 57))
  ContourAtts.SetMultiColor(2, (204, 255, 255, 10))
  ContourAtts.SetMultiColor(3, (255, 204, 0, 5))
  ContourAtts.SetMultiColor(4, (255, 153, 0, 58))
  ContourAtts.SetMultiColor(5, (255, 0, 0, 121))  
  SetPlotOptions(ContourAtts)
  # Magnification on one atomic site below the surface
  SetWindowMode("zoom")
  View3DAtts = View3DAttributes()
  View3DAtts.viewAngle = 30
  View3DAtts.focus = (0., 0., 0.) # (0.108315, 0.108315, 0.)
  View3DAtts.viewUp = (0., 0., 1.)
  View3DAtts.viewAngle = 30
  View3DAtts.parallelScale = 57.1578
  View3DAtts.nearPlane = -114.316
  View3DAtts.farPlane = 114.316
  View3DAtts.imageZoom = ZoomValue
  View3DAtts.perspective = Enable_Perspective
  View3DAtts.eyeAngle = 2
  View3DAtts.centerOfRotationSet = 0
  View3DAtts.centerOfRotation = (0.0, 0.0, 0.0)
  View3DAtts.axis3DScaleFlag = 0
  View3DAtts.axis3DScales = (1, 1, 1)
  View3DAtts.shear = (0, 0, 1)
  View3DAtts.windowValid = 1
  
  if (axis=="X"):
    View3DAtts.viewNormal = (0, -1, 0)
    View3DAtts.viewUp = (0, 0, 1)
    View3DAtts.imagePan = (OriginY, OriginZ)
    
  elif (axis=="Y"):
    View3DAtts.viewNormal = (1, 0, 0)
    View3DAtts.viewUp = (0, 0, 1)
    View3DAtts.imagePan = (OriginX, OriginZ)

  elif (axis=="3D"):
    View3DAtts.viewNormal = (0.95, -0.04, 0.18)
    View3DAtts.focus = (-0.21, -0.21, -0.21)
    View3DAtts.viewUp = (-0.18, -0.011, 0.98)
    View3DAtts.imagePan = (OriginX, OriginZ)
  
#  if(zoomPicture != 1.0):
#    View3DAtts.imagePan = (0.0211978, -0.240702)
#  View3DAtts.imageZoom = 9.13208
#  View3DAtts.perspective = 0
  SetView3D(View3DAtts)
  
  #Adding a label for the time step    
  l1s = CreateAnnotationObject("Text2D")
  l1s.text = "step # "+step
  l1s.position = (0.05, 0.05)
  
  DrawPlots()
  #}}}
