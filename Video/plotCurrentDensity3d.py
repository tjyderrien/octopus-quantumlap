# To use this script: 
# visit -nocli -nowin -s /scratch/tjyderri/octopus/Video/plotDensity.py

##################################################################################
## Copyright (C) 2016 M. Compostella, M. Rampp, H. Appel
## Modified by TJY Derrien, 21st July 2016
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
# TODO: Plot both polarization direections on the same plot !!! Laser direction oscillates, not the other one! 
# TODO: add vertical red bar on the laser graph
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

inputdata = 'listoffiles.txt'
outdir    = './'

densityvtk= 'current.vtk'
#DensityGsVtk= '/../td.0000000/density.vtk' #This should work for any run in our frame. 
#densityvtk= 'elf_rs.vtk'

#elfvtk    = 'elf_rs.vtk'
geometry  = 'geometry.vtk'

bonds     = False

BohrRadius=5.2917721067e-1 #Angstroms

OriginX=0.
OriginY=0.
OriginZ=0. #0.6788749999999979/BohrRadius

MaxValue=0.01
ZoomValue=1.0
EnablePerspective=0

# BohrRadius = phys.physical_constants['Bohr radius']
# print "Bohr radius"
# print BohrRadius

###############################################
# Setting up Python interpreter, if necessary #
###############################################
import sys

# Detect if running in CLI ('visit -cli') and if not check Python version
isCLI = True
if 'VisItException' not in dir():
    if sys.version_info[1] < 5:
        print "\nRequires Python version 2.5 or greater\n"
        sys.exit()
    isCLI=False

# Import visit into a standard Python interpreter, if required (get error
# if run Launch() from CLI). Launch() is required to start the Viewer and
# import the rest of the API.
if not isCLI:
    print
    print "Running VisIt's script from Python interpreter"
    print
    print "Insert path to the 'site-packages' installation directory"
    visitpath = raw_input("E.g: /home/user/VisIt/current/linux-x86_64/lib/site-packages \n")
    winornot = raw_input("Do you want to use a plotting window? (y/n) \n")
    sys.path.append(visitpath)
    from visit import *
    if (winornot=="y"):
        Launch()
    else:
        LaunchNowin()


################################
# Setting up python and inputs #
################################
# Some useful modules
import math, time
import glob, os, shutil
import subprocess
from visit_utils import *

start_script = time.time()

# Setting initial variables or using default input
subprocess.call(["clear", ""])

#print "      *** Welcome to visualization script for VisIt ***      "
#print ""
#print "Description: This script creates a set of 3D plots of a laser"
#print "             pulse propagating towards Si sample."
#print ""
#print "Enter your choice:"
#print "1- Use default value for the setup variables"
#print "2- Specify the location of the file containing the list of input"
#print "3- Specify single inputs"
#answer = raw_input("\n")

#if (answer=="1"):
fileex = os.path.exists(inputdata)
if (fileex==False):
    print "File do not exists!"
    quit()
with open(inputdata, 'r') as infile:
    indata = infile.read().split('\n')
outex = os.path.exists(outdir)
if (outex==False):
    print "Directory does not exists!"
    quit()
molzoom = 1
savename  = 'CurrentDensity'
moviename = "CurrentDensity.mpg"
moviespeed = '2'
subprocess.call(["clear", ""])

#elif (answer=="2"):
    #print "Enter full path to the file containing the list of input:"
    #inputdata = raw_input("(e.g: /home/user/inputdata.txt)\n")
    #fileex = os.path.exists(inputdata)
    #if (fileex==False):
        #print "File does not exists!"
        #quit()
    #with open(inputdata, 'r') as infile:
        #indata = infile.read().split('\n')
    #outdir = raw_input("Enter output directory: (press Enter for current directory)\n")
    #if (len(outdir)==0):
        #outdir = subprocess.check_output(["pwd"]).rstrip("\n")+"/"
    #elif (outdir[len(outdir)-1]!="/"):
        #outdir=outdir+"/"
    #outex = os.path.exists(outdir)
    #if (outex==False):
        #print "Directory does not exists!"
        #quit()
    #molzoom = raw_input("Zoom for the molecule (default is 1):\n")
    #if (len(molzoom)==0):
        #molzoom = 1.
    #else:
        #molzoom = float(molzoom)
    #savename = raw_input("Enter name for output images:\n")
    #if (len(savename)==0):
        #savename = 'density'
    #print "Enter name for output movie with extension: (press Enter to skip movie creation)"
    #moviename = raw_input("You can choose between .wmv, .mpg and .mov file formats\n")
    #if (len(moviename)!=0):
        #moviespeed = raw_input("Set movie speed (1 fast - 100 slow)\n")
    #subprocess.call(["clear", ""])
#elif (answer=="3"):
    #inputdir = raw_input("Enter full path to the directory containing the input timesteps:\n")
    #if (inputdir[len(inputdir)-1]!="/"):
        #inputdir=inputdir+"/"
    #inputstart = raw_input("Enter the name of the first input (e.g. td.0000000/):\n")
    #if (inputstart[len(inputstart)-1]!="/"):
        #inputstart=inputstart+"/"
    #inputend = raw_input("Enter the name of the last input (e.g. td.0030000/):\n")
    #if (inputend[len(inputend)-1]!="/"):
        #inputend=inputend+"/"
    #print "Note that we assume that each input directory contains"
    #print "a elf_rs.vtk, a geometry.xyz and a maxwell_e_field-[x,y,z].vtk file"
    #inputstride = raw_input("Enter the timestep distance between inputs (e.g. 15):\n")
    #if (int(inputstride)==0):
        #totinput=1
    #else:
        #totinput=(int( inputend[len(inputend)-8:len(inputend)-1] )-int( inputstart[len(inputstart)-8:len(inputstart)-1] ))/int(inputstride) +1
    #indata = [None]*(totinput)
    #outdir = raw_input("Enter output directory: (press Enter for current directory)\n")
    #if (len(outdir)==0):
        #outdir = subprocess.check_output(["pwd"]).rstrip("\n")+"/"
    #elif (outdir[len(outdir)-1]!="/"):
        #outdir=outdir+"/"
    #outex = os.path.exists(outdir)
    #if (outex==False):
        #print "Directory does not exists!"
        #quit()
    #molzoom = raw_input("Zoom for the molecule (default is 1):\n")
    #if (len(molzoom)==0):
        #molzoom = 1.
    #else:
        #molzoom = float(molzoom)
    #savename = raw_input("Enter name for output images:\n")
    #if (len(savename)==0):
        #savename = 'azul'
    #print "Enter name for output movie with extension: (press Enter to skip movie creation)"
    #print "Note that you need an encoder like ffmpeg on this machine in order to be able to create a movie"
    #moviename = raw_input("You can choose between .wmv, .mpg and .mov file formats\n")
    #if (len(moviename)!=0):
        #moviespeed = raw_input("Set movie speed (1 fast - 100 slow)\n")
    #subprocess.call(["clear", ""])
    #for i in range(0,totinput):
        #step=int( inputstart[len(inputstart)-8:len(inputstart)-1] )+int(inputstride)*i
        #indata[i]=inputdir + inputstart[0:len(inputend)-8] +"{0:07d}".format(step)+"/"
#else:
    #print "Value not recognized."
    #quit()

# Checking if current directory coincides with output directory
changedir=False
currentpath = subprocess.check_output(["pwd"]).rstrip("\n")+"/"
if (currentpath!=outdir):
    changedir=True

# Checking if input files exist
indata = filter(None, indata)
for i in range(0,len(indata)):
    direx = os.path.exists(indata[i])
    if (direx==False):
        print "Warning! Input directory ",indata[i]," DOES NOT exist! Removing it from input list...\n"
        indata[i]=None

# Removing empty lines
indata = filter(None, indata)

# Check if there are still input files available
if (len(indata)==0):
    print "Total number of input files is 0!"
    print "Please check the paths and the name of the inputs"
    print ""
    print "Note that the script will continue in case you want to produce a video with old images"
    press = raw_input("Press enter to continue or type \"q\" to quit\n")
    if (press=='q'):
        quit()

# Summary of setup variables
print "Summary of the setup variables:\n"
print "-Total number of input: ", len(indata)
if (len(indata)!=0):
    print "-First input: ", indata[0]
if (len(indata)!=0):
    print "-Last input:  ", indata[len(indata)-1]
print "-Output directory: ", outdir
print "-Molecule zoom: ", molzoom
print "-Images name: ", savename+"####.png"
if (len(moviename)==0):
    print "-Movie file will not be created"
else:
    print "-Movie file "+moviename+" will be created"
    print "-Movie speed: ", moviespeed
print ""
press = raw_input("Press enter to continue or type \"q\" to quit\n")
if (press=='q'):
    quit()
subprocess.call(["clear", ""])

#########################
# We now start plotting #
#########################
#
# Setting saving options
s = SaveWindowAttributes()
s.format = s.PNG
s.fileName = savename
s.width, s.height = 2048,1536
s.resConstraint = 0
if (changedir==True):
   s.outputToCurrentDirectory = 0
   s.outputDirectory = outdir
else:
   s.outputToCurrentDirectory = 1
s.advancedMultiWindowSave = 1

# Window structure
# #====#====#====#
# #    #    # 3  #
# # 1  # 2  #====#
# #    #    # 4  #
# #====#====#====#

s.subWindowAtts.win1.size = ( int(s.width/2.), int(s.height)/1.)
s.subWindowAtts.win1.position = (0,0)
s.subWindowAtts.win2.size = ( int(s.width/2.), int(s.height)/1.)
s.subWindowAtts.win2.position = (int(s.width/2.), 0)
#s.subWindowAtts.win3.size = ( int(s.width/3.), int(s.height)/2.)
#s.subWindowAtts.win3.position = ( int(2.*s.width/3.), int(s.height)/2. )
#s.subWindowAtts.win4.size = ( int(s.width/3.), int(s.height)/2.)
#s.subWindowAtts.win4.position = ( int(2.*s.width/3.), 0. )
SetSaveWindowAttributes(s)

# Window 1 is available by default
AddWindow() # Create window 2
#AddWindow() # Create window 3
#AddWindow() # Create window 4

# We define some functions usefull for all windows
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


#### Window 1 ####
SetActiveWindow(1)
#SetActivePlots(1)
SetViewExtentsType(1)

#=================================================
#> End spontaneous state
#=================================================

#InvertBackgroundColor()

#Creating a label for the time step
l1s = CreateAnnotationObject("Text2D")
l1s.text = "Time step"
l1s.position = (0.05, 0.05)
#l1s.height = 0.025
#l1s.fontBold = 1
#l1s.textColor = (255, 255, 255, 255)
l1s.useForegroundForTextColor = 1

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

#from visit_utils import *
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

def plot2Ddifference(axis="X", extreme=0.03, caption=1):#{{{
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
  
  if(axis == "Y"): #Not aligned with Si chemical bonds :-(
    SliceAtts.normal = (0, -1, 0)
    SliceAtts.axisType = SliceAtts.YAxis  # XAxis, YAxis, ZAxis, Arbitrary, ThetaPhi
    SliceAtts.upAxis = (0, 0, 1)
    SliceAtts.flip = 0
    SliceAtts.theta = 180
    SliceAtts.phi = 0
    
  if(axis == "SiBond"): #Aligned with Si chemical bonds :)
    SliceAtts.normal = (0.5, 0, 0.5)
    SliceAtts.originType = SliceAtts.Point  # Point, Intercept, Percent, Zone, Node
    SliceAtts.originPoint = (OriginX, OriginY, OriginZ)
    SliceAtts.axisType = SliceAtts.Arbitrary  # XAxis, YAxis, ZAxis, Arbitrary, ThetaPhi
    SliceAtts.upAxis = (0, 0, 1)
    SliceAtts.flip = 0
    SliceAtts.theta = 0
    SliceAtts.phi = 0
    
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
  #l1s = GetAnnotationObject("Text2D")
  l1s.text = "step # "+step
  l1s.position = (0.05, 0.05)
  
  # Drawing plot to run query
  DrawPlots()
  #}}}
  
def plot3Ddifference(axis="X", extreme=MaxValue, zoomPicture=ZoomValue, EnablePerspective=0):#{{{
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
  ContourAtts.contourNLevels = 5
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
  ContourAtts.SetMultiColor(0, (51, 102, 255, 89))
  ContourAtts.SetMultiColor(1, (153, 204, 255, 57))
  ContourAtts.SetMultiColor(2, (255, 255, 255, 0))
  ContourAtts.SetMultiColor(3, (255, 153, 0, 58))
  ContourAtts.SetMultiColor(4, (255, 0, 0, 124))
  SetPlotOptions(ContourAtts)
  # Magnification on one atomic site below the surface
  #SetWindowMode("zoom")
  
  View3DAtts = View3DAttributes()
  View3DAtts.viewAngle = 30
  View3DAtts.parallelScale = 57.1578
  View3DAtts.nearPlane = -114.316
  View3DAtts.farPlane = 114.316
  View3DAtts.imagePan = (0, 0)
  View3DAtts.imageZoom = zoomPicture
  View3DAtts.perspective = EnablePerspective
  View3DAtts.eyeAngle = 2
  View3DAtts.centerOfRotationSet = 0
  View3DAtts.centerOfRotation = (0.108315, 0.108315, 4.76837e-05)
  View3DAtts.axis3DScaleFlag = 0
  View3DAtts.axis3DScales = (1, 1, 1)
  
  if (axis=="X"):
    View3DAtts.viewNormal = (0, -1, 0)
    View3DAtts.viewUp = (0, 0, 1)
    
  elif (axis=="Y"):
    View3DAtts.viewNormal = (1, 0, 0)
    View3DAtts.viewUp = (0, 0, 1)
  
#  View3DAtts.imagePan = (0.0211978, -0.240702)
#  View3DAtts.imageZoom = 9.13208
#  View3DAtts.perspective = 0
  SetView3D(View3DAtts)
  
  #Adding a label for the time step    
  #l1s = GetAnnotationObject("Text2D")
  l1s.text = "step # "+step
  l1s.position = (0.05, 0.05)
  
  DrawPlots()
  #}}}

def plot3Dcurrent(axis="X", extreme=MaxValue, zoomPicture=ZoomValue, EnablePerspective=0):#{{{
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
  ContourAtts.SetMultiColor(0, (0, 0, 255, 22))
  ContourAtts.SetMultiColor(1, (0, 204, 255, 56))
  ContourAtts.SetMultiColor(2, (153, 204, 255, 55))
  ContourAtts.SetMultiColor(3, (255, 255, 153, 58))
  ContourAtts.SetMultiColor(4, (255, 153, 0, 65))
  ContourAtts.SetMultiColor(5, (255, 0, 0, 66))
  SetPlotOptions(ContourAtts)
  # Magnification on one atomic site below the surface
  #SetWindowMode("zoom")
  
  View3DAtts = View3DAttributes()
  View3DAtts.viewAngle = 30
  View3DAtts.parallelScale = 57.1578
  View3DAtts.nearPlane = -114.316
  View3DAtts.farPlane = 114.316
  View3DAtts.imagePan = (0, 0)
  View3DAtts.imageZoom = zoomPicture
  View3DAtts.perspective = EnablePerspective
  View3DAtts.eyeAngle = 2
  View3DAtts.centerOfRotationSet = 0
  View3DAtts.centerOfRotation = (0.108315, 0.108315, 4.76837e-05)
  View3DAtts.axis3DScaleFlag = 0
  View3DAtts.axis3DScales = (1, 1, 1)
  
  if (axis=="X"):
    View3DAtts.viewNormal = (0, -1, 0)
    View3DAtts.viewUp = (0, 0, 1)
    
  elif (axis=="Y"):
    View3DAtts.viewNormal = (1, 0, 0)
    View3DAtts.viewUp = (0, 0, 1)
  
#  View3DAtts.imagePan = (0.0211978, -0.240702)
#  View3DAtts.imageZoom = 9.13208
#  View3DAtts.perspective = 0
  SetView3D(View3DAtts)
  
  #Adding a label for the time step    
  #l1s = GetAnnotationObject("Text2D")
  l1s.text = "step # "+step
  l1s.position = (0.05, 0.05)
  
  DrawPlots()
  #}}}

##############################
# Cycling on all input files #
##############################
imgtstart = time.time()
#for i in range(0,1):     # Use only the first input for tests
for i in range(0,len(indata)):
    activepath = indata[i]
    if (activepath[len(activepath)-1]!="/"):
        activepath = activepath+"/"
    #indensityGsVtk = activepath+DensityGsVtk #TJYD: we will plot the difference with ground state
    indensityvtk =   activepath+densityvtk
    #invtk = activepath+elfvtk
    inxyz = activepath+geometry
    #inefx = activepath+maxe_x
    #inefy = activepath+maxe_y
    #inefz = activepath+maxe_z
    step = activepath[len(activepath)-8:len(activepath)-1]

    # We open all databases for the plots
    #OpenDatabase(indensityGsVtk)
    elfvar="Re"
    #DefineScalarExpression("GroundStateData", "conn_cmfe(<"+indensityGsVtk+":Re>, <mesh>)") #TODO: script the definition of this mesh
    OpenDatabase(indensityvtk)
    #DefineScalarExpression("ExcitedStateData", "conn_cmfe(<"+indensityvtk+":Re>, <mesh>)")
    # DefineScalarExpression("Difference", "Re") #-GroundStateData") #DEBUG
    #DefineScalarExpression("Difference", "ExcitedStateData - GroundStateData") #Re: TD-data, GroundStateData: GS-density. This is the performed operation on the scalar data. 
    
    #OpenDatabase(invtk)
    #elfvar="Re"
    OpenDatabase(inxyz)
    molvar="element"
    #OpenDatabase(inefx)
    #fldx="Re"
    # We specify the variables of the y and z components of the E field, but it is not necessary to open these databases
    #fldy="Re"
    #fldz="Re"
    
    ## TODO: Building an electric field object - useful later 
    ## We define a vector field from the [x,y,z] scalar components of the electric field
    #expression = "{conn_cmfe(<"+inefx+":"+fldx+">,mesh), conn_cmfe(<"+inefy+":"+fldy+">,mesh), conn_cmfe(<"+inefz+":"+fldz+">,mesh)}"
    #DefineVectorExpression("Evec", expression)

    #expression2D = "{conn_cmfe(<"+inefx+":"+fldx+">,mesh), conn_cmfe(<"+inefy+":"+fldy+">,mesh), 0.}"
    #DefineVectorExpression("Evec2D", expression2D)

    ## Retriving extention of the E-field along x dataset
    #ext = extension(inefx, fldx)
    
    #### Window 1 ####
    SetActiveWindow(1)
    #SetActivePlots(0)
    plot3Dcurrent(axis="X", extreme=MaxValue, zoomPicture=ZoomValue)
    
    #### Window 2 ####
    #CloneWindow()
    #SetActivePlots(1)
    SetActiveWindow(2)
    #plot2Ddifference(axis="SiBond", extreme=MaxValue, caption=0)
    plot3Dcurrent(axis="Y", extreme=MaxValue, zoomPicture=ZoomValue)
    
    #### Window 3 ####
    #CloneWindow()
    #SetActivePlots(0)
    #SetActiveWindow(3)
    #plot2Ddifference(axis="Z", extreme=MaxValue, caption=0)
    
    #### Window 4 ####
    ##SetActiveWindow(4)
    
    # Adding the plot of the laser pulse
    
    # NOTE: filtering for the td.general files is done by ./Video/buildFileList.sh, 
    #plot1D(filepath="../td.general/laser.txt", varX="var01", varY="var02", colorcode=(255, 0, 0, 255))
    #AddTracerObject(filepath="../td.general/laser.txt", time=step)
    #TODO: put to the scale
    #plot1D(filepath="../td.general/laser2.txt", varX="var01", varY="var02", colorcode=(0, 0, 255, 255))
    #plot1D(filepath="../td.general/current.txt", varX="var01", varY="var02", colorcode=(0, 0, 0, 255))
    #plotVerticalBar(varX=step, colorcode=(0, 0, 0, 255)) #black vertical line
    
    #TODO: add a vertical bar for the current time step
    # the same index than <step> variable should be found in the laser file
    
    #### All Windows ####
    if (i==0):
        zoomtstr = 0.
        zoomtend = 0.

    # Saving the image
    #SetViewExtentsType(1)
    SaveWindow()

    # Clear and close
    for winind in range(1,5):
        SetActiveWindow(winind)
        DeleteAllPlots()
        RemoveAllOperators()
   
    CloseDatabase(indensityvtk)
    CloseDatabase(indensityGsVtk)
    CloseDatabase(inxyz)

    # Estimating remaining time
    print "Time left (sec): ",((time.time()-imgtstart-(zoomtend-zoomtstr))/(i+1))*(len(indata)-(i+1))

# End of the loop over the input files

imgtend = time.time()
subprocess.call(["clear", ""])
print "Creation of the images completed."

# Movie creation
if (len(moviename)!=0):
    vidtstart= time.time() 
    input_pattern = outdir+savename+"%04d.png"
    output_pattern = outdir+moviename
    encoding.encode(input_pattern,output_pattern,fdup=int(moviespeed))
    vidtend = time.time()


# Reporting working time
print "\n"
print " *** Time statistics ***"
if (len(indata)!=0):
    print "Images creation time:"
    print("--- %s seconds ---" % (imgtend - imgtstart-(zoomtend-zoomtstr)))
    print( "On average %s seconds per image" % ((imgtend - imgtstart-(zoomtend-zoomtstr))/len(indata)) )
    if (molzoom>1.):
        print ""
        print "Additional time for creating the zooming:"
        print("--- %s seconds ---" % (zoomtend-zoomtstr))
    print ""
print "Video creation time:"
if (len(moviename)!=0):
    print("--- %s seconds ---" % (vidtend - vidtstart))
else:
    print"--- Not Available ---"
print ""
print "Script total execution time:"
print("--- %s seconds ---" % (time.time() - start_script))
print ""

quit()

exit()

##################################################################
