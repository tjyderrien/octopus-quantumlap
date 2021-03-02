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

inputdata = 'listoffiles.txt'
outdir    = './'

# densityvtk= 'density.vtk'
densityvtk= 'elf_rs.vtk'

#elfvtk    = 'elf_rs.vtk'
#geometry  = 'geometry.xyz'

bonds     = False


###############################################
# Setting up Python interpreter, if necessary #
###############################################
import sys

# Detect if running in CLI ('visit -cli') and if not check Python version
isCLI = True
if 'VisItException' not in dir():
    if sys.version_info[1] < 5:
        print("\nRequires Python version 2.5 or greater\n")
        sys.exit()
    isCLI=False

# Import visit into a standard Python interpreter, if required (get error
# if run Launch() from CLI). Launch() is required to start the Viewer and
# import the rest of the API.
if not isCLI:
    print()
    print("Running VisIt's script from Python interpreter")
    print()
    print("Insert path to the 'site-packages' installation directory")
    visitpath = input("E.g: /home/user/VisIt/current/linux-x86_64/lib/site-packages \n")
    winornot = input("Do you want to use a plotting window? (y/n) \n")
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
import time
import os
import subprocess
from visit_utils import *

start_script = time.time()

# Setting initial variables or using default input
subprocess.call(["clear", ""])

print(" *** Welcome to Azulene3D.py, a visualization script for VisIt ***")
print("")
print("Description: This script creates a set of 3D plots of a laser")
print("             pulse propagating towards Azulene C10H8.")
print("")
print("Enter your choice:")
print("1- Use default value for the setup variables")
print("2- Specify the location of the file containing the list of input")
print("3- Specify single inputs")
answer = input("\n")

if (answer=="1"):
    fileex = os.path.exists(inputdata)
    if (fileex==False):
        print("File do not exists!")
        quit()
    with open(inputdata, 'r') as infile:
        indata = infile.read().split('\n')
    outex = os.path.exists(outdir)
    if (outex==False):
        print("Directory does not exists!")
        quit()
    molzoom = 1.
    savename  = 'density'
    moviename = "density.mpg"
    moviespeed = '50'
    subprocess.call(["clear", ""])
elif (answer=="2"):
    print("Enter full path to the file containing the list of input:")
    inputdata = input("(e.g: /home/user/inputdata.txt)\n")
    fileex = os.path.exists(inputdata)
    if (fileex==False):
        print("File does not exists!")
        quit()
    with open(inputdata, 'r') as infile:
        indata = infile.read().split('\n')
    outdir = input("Enter output directory: (press Enter for current directory)\n")
    if (len(outdir)==0):
        outdir = subprocess.check_output(["pwd"]).rstrip("\n")+"/"
    elif (outdir[len(outdir)-1]!="/"):
        outdir=outdir+"/"
    outex = os.path.exists(outdir)
    if (outex==False):
        print("Directory does not exists!")
        quit()
    molzoom = input("Zoom for the molecule (default is 1):\n")
    if (len(molzoom)==0):
        molzoom = 1.
    else:
        molzoom = float(molzoom)
    savename = input("Enter name for output images:\n")
    if (len(savename)==0):
        savename = 'density'
    print("Enter name for output movie with extension: (press Enter to skip movie creation)")
    moviename = input("You can choose between .wmv, .mpg and .mov file formats\n")
    if (len(moviename)!=0):
        moviespeed = input("Set movie speed (1 fast - 100 slow)\n")
    subprocess.call(["clear", ""])
elif (answer=="3"):
    inputdir = input("Enter full path to the directory containing the input timesteps:\n")
    if (inputdir[len(inputdir)-1]!="/"):
        inputdir=inputdir+"/"
    inputstart = input("Enter the name of the first input (e.g. td.0000000/):\n")
    if (inputstart[len(inputstart)-1]!="/"):
        inputstart=inputstart+"/"
    inputend = input("Enter the name of the last input (e.g. td.0030000/):\n")
    if (inputend[len(inputend)-1]!="/"):
        inputend=inputend+"/"
    print("Note that we assume that each input directory contains")
    print("a elf_rs.vtk, a geometry.xyz and a maxwell_e_field-[x,y,z].vtk file")
    inputstride = input("Enter the timestep distance between inputs (e.g. 15):\n")
    if (int(inputstride)==0):
        totinput=1
    else:
        totinput=(int( inputend[len(inputend)-8:len(inputend)-1] )-int( inputstart[len(inputstart)-8:len(inputstart)-1] ))/int(inputstride) +1
    indata = [None]*(totinput)
    outdir = input("Enter output directory: (press Enter for current directory)\n")
    if (len(outdir)==0):
        outdir = subprocess.check_output(["pwd"]).rstrip("\n")+"/"
    elif (outdir[len(outdir)-1]!="/"):
        outdir=outdir+"/"
    outex = os.path.exists(outdir)
    if (outex==False):
        print("Directory does not exists!")
        quit()
    molzoom = input("Zoom for the molecule (default is 1):\n")
    if (len(molzoom)==0):
        molzoom = 1.
    else:
        molzoom = float(molzoom)
    savename = input("Enter name for output images:\n")
    if (len(savename)==0):
        savename = 'azul'
    print("Enter name for output movie with extension: (press Enter to skip movie creation)")
    print("Note that you need an encoder like ffmpeg on this machine in order to be able to create a movie")
    moviename = input("You can choose between .wmv, .mpg and .mov file formats\n")
    if (len(moviename)!=0):
        moviespeed = input("Set movie speed (1 fast - 100 slow)\n")
    subprocess.call(["clear", ""])
    for i in range(0,totinput):
        step=int( inputstart[len(inputstart)-8:len(inputstart)-1] )+int(inputstride)*i
        indata[i]=inputdir + inputstart[0:len(inputend)-8] +"{0:07d}".format(step)+"/"
else:
    print("Value not recognized.")
    quit()

# Checking if current directory coincides with output directory
changedir=False
currentpath = subprocess.check_output(["pwd"]).rstrip("\n")+"/"
if (currentpath!=outdir):
    changedir=True

# Checking if input files exist
indata = [_f for _f in indata if _f]
for i in range(0,len(indata)):
    direx = os.path.exists(indata[i])
    if (direx==False):
        print("Warning! Input directory ",indata[i]," DOES NOT exist! Removing it from input list...\n")
        indata[i]=None

# Removing empty lines
indata = [_f for _f in indata if _f]

# Check if there are still input files available
if (len(indata)==0):
    print("Total number of input files is 0!")
    print("Please check the paths and the name of the inputs")
    print("")
    print("Note that the script will continue in case you want to produce a video with old images")
    press = input("Press enter to continue or type \"q\" to quit\n")
    if (press=='q'):
        quit()

# Summary of setup variables
print("Summary of the setup variables:\n")
print("-Total number of input: ", len(indata))
if (len(indata)!=0):
    print("-First input: ", indata[0])
if (len(indata)!=0):
    print("-Last input:  ", indata[len(indata)-1])
print("-Output directory: ", outdir)
print("-Molecule zoom: ", molzoom)
print("-Images name: ", savename+"####.png")
if (len(moviename)==0):
    print("-Movie file will not be created")
else:
    print("-Movie file "+moviename+" will be created")
    print("-Movie speed: ", moviespeed)
print("")
press = input("Press enter to continue or type \"q\" to quit\n")
if (press=='q'):
    quit()
subprocess.call(["clear", ""])

#########################
# We now start plotting #
#########################
# This script will create 4 different windows (plus a 5th window used for calculations, but then destroyed). 
# Therefore, we start by defining the saving options for these images.
# We want to obtain something like:
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

## Multiwindows support
#s.advancedMultiWindowSave = 1
#s.subWindowAtts.win1.size = (1024,1024)
#s.subWindowAtts.win1.position = (0,0)
#s.subWindowAtts.win2.size = (1024,1024)
#s.subWindowAtts.win2.position = (1024,0)
#s.subWindowAtts.win3.size = (2048,256)
#s.subWindowAtts.win3.position = (0,1024)
#s.subWindowAtts.win4.size = (2048,256)
#s.subWindowAtts.win4.position = (0,1280)
SetSaveWindowAttributes(s)

# First, we prepare the 4 plotting windows and set their attributes
# Window 1 is available by default
#AddWindow() # Create window 2
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
SetViewExtentsType(1)

# Setting annotations
a1 = AnnotationAttributes()
a1.userInfoFlag = 0
a1.databaseInfoFlag = 0
a1.axes3D.visible = 0
a1.axes3D.triadFlag = 0
a1.axes3D.bboxFlag = 1
a1.axes3D.lineWidth = 1
SetAnnotationAttributes(a1)

# Setting box size
BoxAtts = BoxAttributes()
BoxAtts.amount = BoxAtts.All  # Some, All
BoxAtts.minx = -10
BoxAtts.maxx = 50
BoxAtts.miny = -5
BoxAtts.maxy = 60
BoxAtts.minz = -10
BoxAtts.maxz = 30
BoxAtts.inverse = 0
SetOperatorOptions(BoxAtts, 0)

# Trying to restrict viewing to it
SetViewExtentsType(1)
# Begin spontaneous state
View2DAtts = View2DAttributes()
View2DAtts.windowCoords = (-5.19914, 28.5951, -2.17319, 60.4145)
View2DAtts.viewportCoords = (0.2, 0.95, 0.15, 0.95)
View2DAtts.fullFrameActivationMode = View2DAtts.Auto  # On, Off, Auto
View2DAtts.fullFrameAutoThreshold = 100
View2DAtts.xScale = View2DAtts.LINEAR  # LINEAR, LOG
View2DAtts.yScale = View2DAtts.LINEAR  # LINEAR, LOG
View2DAtts.windowValid = 1
SetView2D(View2DAtts)
# End spontaneous state

ViewCurveAtts = ViewCurveAttributes()
ViewCurveAtts.domainCoords = (0, 1)
ViewCurveAtts.rangeCoords = (0, 1)
ViewCurveAtts.viewportCoords = (0.2, 0.95, 0.15, 0.95)
ViewCurveAtts.domainScale = ViewCurveAtts.LINEAR  # LINEAR, LOG
ViewCurveAtts.rangeScale = ViewCurveAtts.LINEAR  # LINEAR, LOG
SetViewCurve(ViewCurveAtts)
View2DAtts = View2DAttributes()
View2DAtts.windowCoords = (-5.19914, 28.5951, -2.17319, 60.4145)
View2DAtts.viewportCoords = (0.2, 0.95, 0.15, 0.95)
View2DAtts.fullFrameActivationMode = View2DAtts.Auto  # On, Off, Auto
View2DAtts.fullFrameAutoThreshold = 100
View2DAtts.xScale = View2DAtts.LINEAR  # LINEAR, LOG
View2DAtts.yScale = View2DAtts.LINEAR  # LINEAR, LOG
View2DAtts.windowValid = 1
SetView2D(View2DAtts)
View3DAtts = View3DAttributes()
View3DAtts.viewNormal = (0, 0, 1)
View3DAtts.focus = (-0.216648, -0.216648, 0)
View3DAtts.viewUp = (0, 1, 0)
View3DAtts.viewAngle = 30
View3DAtts.parallelScale = 90.1424
View3DAtts.nearPlane = -180.285
View3DAtts.farPlane = 180.285
View3DAtts.imagePan = (0, 0)
View3DAtts.imageZoom = 1
View3DAtts.perspective = 1
View3DAtts.eyeAngle = 2
View3DAtts.centerOfRotationSet = 0
View3DAtts.centerOfRotation = (-0.216648, -0.216648, 0)
View3DAtts.axis3DScaleFlag = 0
View3DAtts.axis3DScales = (1, 1, 1)
View3DAtts.shear = (0, 0, 1)
View3DAtts.windowValid = 1
SetView3D(View3DAtts)
ViewAxisArrayAtts = ViewAxisArrayAttributes()
ViewAxisArrayAtts.domainCoords = (0, 1)
ViewAxisArrayAtts.rangeCoords = (0, 1)
ViewAxisArrayAtts.viewportCoords = (0.15, 0.9, 0.1, 0.85)
SetViewAxisArray(ViewAxisArrayAtts)


#InvertBackgroundColor()

# Creating a label for the description
#l1 = CreateAnnotationObject("Text2D")
#l1.text = "Pulse propagation"
#l1.position = (0.68, 0.05)
#l1.height = 0.025
#l1.fontBold = 1
#l1.textColor = (255, 255, 255, 255)
#l1.useForegroundForTextColor = 0

# Creating a label for the time step
#l1s = CreateAnnotationObject("Text2D")
#l1s.text = "Time step"
#l1s.position = (0.05, 0.05)
#l1s.height = 0.025
#l1s.fontBold = 1
##l1s.textColor = (255, 255, 255, 255)
#l1s.useForegroundForTextColor = 1

## Creating a label for the molecule zoom
#l1m = CreateAnnotationObject("Text2D")
#l1m.text = "Molecule"
#l1m.position = (0.68, 0.90)
#l1m.height = 0.025
#l1m.fontBold = 1
#l1m.textColor = (255, 255, 255, 255)
#l1m.useForegroundForTextColor = 0

# Defining an arrow if we are going to zoom on the molecule
if (molzoom>1):
    ar1 = CreateAnnotationObject("Line2D")
    ar1.visible = 1
    ar1.position  = (0.543, 0.52)
    ar1.position2 = (0.80, 0.88)
    ar1.width = 2
    ar1.beginArrow = 2
    ar1.useForegroundForLineColor = 0
    ar1.color = (160, 160, 160, 255)

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

#### ADDED MY PERSONAL SCRIPT HERE

##OpenDatabase("localhost:/eos/scratch/tjyderri/octopus/02-Slab2D/Si-4x4x6-k-1x1x1-L10-MoveIons/lambda800e-9-tau10e-15-field1E12-CEP0.0/lambda400e-9-tau10e-15-field5E8-CEP0.0-delay0.0/lambda400e-9-tau10e-15-field5E8-CEP0.0-delay0.0/output_iter/td.0000000/density.vtk", 0)
#AddPlot("Pseudocolor", "Re", 1, 1)
#AddOperator("Slice", 1)
##SetActivePlots(0)
##SetActivePlots(0)
#SliceAtts = SliceAttributes()
##SliceAtts.originType = SliceAtts.Intercept  # Point, Intercept, Percent, Zone, Node
##SliceAtts.originPoint = (0, 0, 0)
##SliceAtts.originIntercept = 0
##SliceAtts.originPercent = 0
##SliceAtts.originZone = 0
##SliceAtts.originNode = 0
#SliceAtts.normal = (0, -1, 0)
#SliceAtts.axisType = SliceAtts.YAxis  # XAxis, YAxis, ZAxis, Arbitrary, ThetaPhi
#SliceAtts.upAxis = (0, 0, 1)
#SliceAtts.project2d = 1
#SliceAtts.interactive = 1
##SliceAtts.flip = 0
##SliceAtts.originZoneDomain = 0
##SliceAtts.originNodeDomain = 0
##SliceAtts.meshName = "mesh"
##SliceAtts.theta = 0
##SliceAtts.phi = 0
#SetOperatorOptions(SliceAtts, 1)
#DrawPlots()
# The SaveView RPC is not supported in the VisIt module so it will not be logged.
#SaveWindowAtts = SaveWindowAttributes()
#SaveWindowAtts.outputToCurrentDirectory = 1
#SaveWindowAtts.outputDirectory = "."
#SaveWindowAtts.fileName = "visit"
#SaveWindowAtts.family = 1
#SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
#SaveWindowAtts.width = 1024
#SaveWindowAtts.height = 1024
#SaveWindowAtts.screenCapture = 0
#SaveWindowAtts.saveTiled = 0
#SaveWindowAtts.quality = 80
#SaveWindowAtts.progressive = 0
#SaveWindowAtts.binary = 0
#SaveWindowAtts.stereo = 0
#SaveWindowAtts.compression = SaveWindowAtts.PackBits  # None, PackBits, Jpeg, Deflate
#SaveWindowAtts.forceMerge = 0
#SaveWindowAtts.resConstraint = SaveWindowAtts.ScreenProportions  # NoConstraint, EqualWidthHeight, ScreenProportions
#SaveWindowAtts.advancedMultiWindowSave = 0
#SaveWindow()



##### Window 2 ####
#SetActiveWindow(2)

## Setting annotations
#a2 = AnnotationAttributes()
#a2.userInfoFlag = 0
#a2.databaseInfoFlag = 0
#a2.axes3D.visible = 0
#a2.axes3D.triadFlag = 0
#a2.axes3D.bboxFlag = 0
#a2.axes2D.visible = 0
#SetAnnotationAttributes(a2)

#InvertBackgroundColor()

## Creating a label for the description
#l2 = CreateAnnotationObject("Text2D")
#l2.text = "ELF and E-field at Z=0"
#l2.position = (0.58, 0.05)
#l2.height = 0.025
#l2.fontBold = 1
#l2.textColor = (255, 255, 255, 255)
#l2.useForegroundForTextColor = 0

#def slice2():
    #AddOperator("Slice", 0)
    #s2 = SliceAttributes()
    #s2.originType = s2.Percent
    #s2.originPercent = 50
    #s2.axisType = s2.ZAxis
    #SetOperatorOptions(s2)

#def winatt2():
    #w2 = View2DAttributes()
    #w2.windowCoords = (-17.*molzoom, 17.*molzoom, -17.*molzoom, 17.*molzoom)
    #w2.viewportCoords = (0, 1, 0, 1)
    #w2.fullFrameActivationMode = w2.On
    #SetView2D(w2)


##### Window 3 ####
#SetActiveWindow(3)

## Setting annotations
#a3 = AnnotationAttributes()
#a3.userInfoFlag = 0
#a3.databaseInfoFlag = 0
#a3.axes2D.visible = 0
#SetAnnotationAttributes(a3)

#InvertBackgroundColor()

## Creating a label for the description
#l3 = CreateAnnotationObject("Text2D")
#l3.text = "E-field x slice along Y-Z"
#l3.position = (0.77, 0.05)
#l3.height = 0.1
#l3.fontBold = 1
#l3.textColor = (0, 0, 0, 255)  # Black
##l3.textColor = (255, 255, 255, 255)  # White
#l3.useForegroundForTextColor = 0

#def slice3():
    #AddOperator("Slice", 0)
    #s3 = SliceAttributes()
    #s3.originType = s3.Percent
    #s3.originPercent = 50
    #s3.axisType = s3.XAxis
    #SetOperatorOptions(s3)

#def winatt3():
    #w3 = View2DAttributes()
    #w3.windowCoords = (ext[4], ext[5], ext[2], ext[3])
    #w3.viewportCoords = (0, 1, 0, 1)
    #w3.fullFrameActivationMode = w3.On
    #SetView2D(w3)


##### Window 4 ####
#SetActiveWindow(4)

## Setting annotations
#a4 = AnnotationAttributes()
#a4.userInfoFlag = 0
#a4.databaseInfoFlag = 0
#a4.axes2D.visible = 0
#SetAnnotationAttributes(a2)

#InvertBackgroundColor()

## Setting Lineout attributes
#lo4 = GetGlobalLineoutAttributes()
#lo4.windowId = 4
#lo4.createWindow = 0
#SetGlobalLineoutAttributes(lo4)

## Creating a label for the description
#l4 = CreateAnnotationObject("Text2D")
#l4.text = "Pulse intensity along Z"
#l4.position = (0.77, 0.05)
#l4.height = 0.1
#l4.fontBold = 1
#l4.textColor = (255, 255, 255, 255)
#l4.useForegroundForTextColor = 0

#def winatt4():
    #length = ext[5]-ext[4]
    #w4 = ViewCurveAttributes()
    #w4.domainCoords = (0, length)
    #w4.rangeCoords = (-0.05, 0.05)
    #w4.viewportCoords = (0, 1, 0, 1)
    #SetViewCurve(w4)



##############################
# Cycling on all input files #
##############################
imgtstart = time.time()
#for i in range(0,1):     # Use only the first input for tests
for i in range(0,len(indata)):
    activepath = indata[i]
    if (activepath[len(activepath)-1]!="/"):
        activepath = activepath+"/"
    indensityvtk = activepath+densityvtk
    #invtk = activepath+elfvtk
    #inxyz = activepath+geometry
    #inefx = activepath+maxe_x
    #inefy = activepath+maxe_y
    #inefz = activepath+maxe_z
    step = activepath[len(activepath)-8:len(activepath)-1]

    # We open all databases for the plots
    OpenDatabase(indensityvtk)
    elfvar="Re"
    #OpenDatabase(invtk)
    #elfvar="Re"
    #OpenDatabase(inxyz)
    #molvar="element"
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
    SetViewExtentsType(1)
    ## TJYD: Adding pseudo-color plot on slices
    ReOpenDatabase(indensityvtk)
    AddPlot("Pseudocolor", "Re", 1, 1)
    AddOperator("Slice", 1)
    SliceAtts = SliceAttributes()
    SliceAtts.normal = (0, -1, 0)
    SliceAtts.axisType = SliceAtts.YAxis  # XAxis, YAxis, ZAxis, Arbitrary, ThetaPhi
    SliceAtts.upAxis = (0, 0, 1)
    SliceAtts.project2d = 1
    SliceAtts.interactive = 1
    SetOperatorOptions(SliceAtts, 1)
    
    #Defining the color scale
    PseudocolorAtts = PseudocolorAttributes()
    PseudocolorAtts.colorTableName = "PuBu"
    SetPlotOptions(PseudocolorAtts)
    
    SetActivePlots(0)
    
    
    # Setting box size
    AddOperator("Box", 1)
    BoxAtts = BoxAttributes()
    BoxAtts.amount = BoxAtts.Some  # Some, All
    BoxAtts.minx = -10
    BoxAtts.maxx = 50
    BoxAtts.miny = -5
    BoxAtts.maxy = 60
    BoxAtts.minz = -10
    BoxAtts.maxz = 30
    BoxAtts.inverse = 0
    SetOperatorOptions(BoxAtts, 1)

    # Trying to restrict viewing to the box
    SetViewExtentsType(1)
    # Begin spontaneous state
    View2DAtts = View2DAttributes()
    View2DAtts.windowCoords = (-5.19914, 28.5951, -2.17319, 60.4145)
    View2DAtts.viewportCoords = (0.2, 0.95, 0.15, 0.95)
    View2DAtts.fullFrameActivationMode = View2DAtts.Auto  # On, Off, Auto
    View2DAtts.fullFrameAutoThreshold = 100
    View2DAtts.xScale = View2DAtts.LINEAR  # LINEAR, LOG
    View2DAtts.yScale = View2DAtts.LINEAR  # LINEAR, LOG
    View2DAtts.windowValid = 1
    SetView2D(View2DAtts)
    # End spontaneous state

    ViewCurveAtts = ViewCurveAttributes()
    ViewCurveAtts.domainCoords = (0, 1)
    ViewCurveAtts.rangeCoords = (0, 1)
    ViewCurveAtts.viewportCoords = (0.2, 0.95, 0.15, 0.95)
    ViewCurveAtts.domainScale = ViewCurveAtts.LINEAR  # LINEAR, LOG
    ViewCurveAtts.rangeScale = ViewCurveAtts.LINEAR  # LINEAR, LOG
    SetViewCurve(ViewCurveAtts)
    View2DAtts = View2DAttributes()
    View2DAtts.windowCoords = (-5.19914, 28.5951, -2.17319, 60.4145)
    View2DAtts.viewportCoords = (0.2, 0.95, 0.15, 0.95)
    View2DAtts.fullFrameActivationMode = View2DAtts.Auto  # On, Off, Auto
    View2DAtts.fullFrameAutoThreshold = 100
    View2DAtts.xScale = View2DAtts.LINEAR  # LINEAR, LOG
    View2DAtts.yScale = View2DAtts.LINEAR  # LINEAR, LOG
    View2DAtts.windowValid = 1
    SetView2D(View2DAtts)
    View3DAtts = View3DAttributes()
    View3DAtts.viewNormal = (0, 0, 1)
    View3DAtts.focus = (-0.216648, -0.216648, 0)
    View3DAtts.viewUp = (0, 1, 0)
    View3DAtts.viewAngle = 30
    View3DAtts.parallelScale = 90.1424
    View3DAtts.nearPlane = -180.285
    View3DAtts.farPlane = 180.285
    View3DAtts.imagePan = (0, 0)
    View3DAtts.imageZoom = 1
    View3DAtts.perspective = 1
    View3DAtts.eyeAngle = 2
    View3DAtts.centerOfRotationSet = 0
    View3DAtts.centerOfRotation = (-0.216648, -0.216648, 0)
    View3DAtts.axis3DScaleFlag = 0
    View3DAtts.axis3DScales = (1, 1, 1)
    View3DAtts.shear = (0, 0, 1)
    View3DAtts.windowValid = 1
    SetView3D(View3DAtts)
    ViewAxisArrayAtts = ViewAxisArrayAttributes()
    ViewAxisArrayAtts.domainCoords = (0, 1)
    ViewAxisArrayAtts.rangeCoords = (0, 1)
    ViewAxisArrayAtts.viewportCoords = (0.15, 0.9, 0.1, 0.85)
    SetViewAxisArray(ViewAxisArrayAtts)

    
    ## Plotting E-field contours
    #ReOpenDatabase(inefx)
    #AddPlot("Contour", fldx)
    #c1 = ContourAttributes()
    #c1.colorType = c1.ColorByMultipleColors
    #c1.contourMethod = c1.Value
    #c1.contourValue = (-0.02, 0.02)
    ##c1.contourMethod = c1.Percent
    ##c1.contourPercent = (30, 70)
    #c1.SetMultiColor(0,(0,0,255,127))
    #c1.SetMultiColor(1,(255,0,0,127))
    #c1.lineWidth= 1
    #c1.legendFlag = 0
    #SetPlotOptions(c1)

    # Drawing plot to run query
    #SetViewExtentsType(1)
    DrawPlots()
    #SetViewExtentsType(1)
    ## Obtaining the spatial extent of the E-field
    #Query("SpatialExtents")
    #ext1 = GetQueryOutputValue()

    # Plotting the molecule
    #ReOpenDatabase(inxyz)
    #AddPlot("Molecule", molvar)
    #mol1 = MoleculeAttributes()
    #mol1.drawAtomsAs = mol1.SphereAtoms
    #mol1.scaleRadiusBy = mol1.Atomic
    #mol1.radiusScaleFactor = molzoom/2.
    #mol1.drawBondsAs = mol1.CylinderBonds
    #mol1.bondRadius = 0.04
    #mol1.legendFlag = 0
    #mol1.atomSphereQuality = mol1.Super
    #mol1.bondCylinderQuality = mol1.Super
    #SetPlotOptions(mol1)

    # Rescaling the molecule
    #rescmol(molzoom)

    ## Plotting some ELF contours around the molecule
    #ReOpenDatabase(invtk)
    #AddPlot("Contour", elfvar)
    #c1 = ContourAttributes()
    #c1.colorType = c1.ColorByMultipleColors
    #c1.contourMethod = c1.Value
    #c1.contourValue = (0.01, 0.1, 0.26, 0.65, 0.85)
    #c1.SetMultiColor(0,(60,60,255,20))
    #c1.SetMultiColor(1,(111,255,79,20))
    #c1.SetMultiColor(2,(255,255,60,20))
    #c1.SetMultiColor(3,(255,208,74,20))
    #c1.SetMultiColor(4,(255,187,146,20))
    #c1.lineWidth= 1
    #c1.legendFlag = 0
    #SetPlotOptions(c1)

    #rescmol(molzoom)
   
    # Adding a label for the time step    
    #l1s = GetAnnotationObject("Text2D")
    #l1s.text = "step # "+step

    # Adding a label for the molecule zoom
    #lm = GetAnnotationObject("Text2D3")
    #lm.text = "Molecule x"+str(molzoom)

    # Drawing
    #DrawPlots()

    # Setting window attribute
    #winatt1()
    
    ##### Window 2 ####
    #SetActiveWindow(2)
    
    ## Adding the plot of the Azulene molecule
    #ReOpenDatabase(inxyz)
    #AddPlot("Pseudocolor", "element")
    #mol2 = PseudocolorAttributes()
    #mol2.colorTableName = "cpk_jmol"
    #mol2.minFlag,mol2.min = 1, -1.
    #mol2.maxFlag,mol2.max = 1, 118.
    #mol2.legendFlag = 0
    #mol2.pointSizePixels = 25
    #mol2.pointType = mol2.Sphere
    #mol2.lineType = mol2.Line
    #mol2.lineWidth = 3
    #mol2.renderSurfaces = 1
    #mol2.renderPoints = 1
    #mol2.smoothingLevel = 2
    #SetPlotOptions(mol2)

    #AddOperator("Project", 0)

    ## Creating bonds
    #if (bonds==True):
        #createbonds()

    ## Rescaling the molecule
    #rescmol(molzoom)
   
    ## Plotting ELF contours
    #ReOpenDatabase(invtk)
    #AddPlot("Contour", elfvar)
    #c2 = ContourAttributes()
    #c2.colorType = c2.ColorByMultipleColors
    #c2.contourMethod = c2.Value
    #c2.contourValue = (0.01, 0.1, 0.26, 0.65, 0.85)
    #c2.SetMultiColor(0,(60,60,255,255))
    #c2.SetMultiColor(1,(111,255,79,255))
    #c2.SetMultiColor(2,(255,255,60,255))
    #c2.SetMultiColor(3,(255,208,74,255))
    #c2.SetMultiColor(4,(255,187,146,255))
    #c2.lineWidth= 1
    #c2.legendFlag = 0
    #SetPlotOptions(c2)

    #slice2()

    ## Rescaling the contours
    #rescmol(molzoom)

    ## Plotting Electric vector field
    #ReOpenDatabase(inefx)
    #AddPlot("Vector", "Evec2D")
    #v2 = VectorAttributes()
    #v2.useLegend = 0
    #v2.glyphLocation = v2.UniformInSpace
    #v2.nVectors = 500
    #v2.lineWidth = 1
    #v2.scale = 0.15
    #v2.scaleByMagnitude = 0
    #v2.autoScale = 1
    #v2.headSize = 0.15
    #v2.lineStem = v2.Cylinder
    #v2.stemWidth = 0.03
    #v2.vectorOrigin = v2.Middle
    #v2.geometryQuality = v2.High
    ##v2.colorByMag = 1
    ##v2.colorTableName = "hot_and_cold"
    #v2.colorByMag = 0
    #v2.vectorColor = (0,255,255,100)
    #SetPlotOptions(v2)

    #slice2()
    
    ## Drawing
    #DrawPlots()

    ## Setting window attributes
    #winatt2()
    
    ##### Window 3 #### (plus, we create Window 5)
    #SetActiveWindow(3)

    ## Plotting E-field along x
    #ReOpenDatabase(inefx)
    #AddPlot("Pseudocolor", fldx)
    #p3 = PseudocolorAttributes()
    #p3.colorTableName = "hot_desaturated"
    #p3.legendFlag = 0
    #p3.minFlag,p3.min = 1, -0.05
    #p3.maxFlag,p3.max = 1, 0.05
    #SetPlotOptions(p3)

    #slice3()

    ## We clone this window to use it for the laser pulse intensity in Window 4
    #CloneWindow()
    #SetActiveWindow(3)

    ## Plotting Azulene molecule with a volume plot
    ##ReOpenDatabase(invtk)
    ##AddPlot("Volume", elfvar)
    ##v3 = VolumeAttributes()
    ##v3.legendFlag = 0
    ##v3.smoothData = 1
    ##v3.opacityAttenuation = 1.
    ##v3.SetColorControlPoints(GetColorTable("gray"))
    ##SetPlotOptions(v3)

    ## Plotting Azulene molecule with a multi-contour plot
    #ReOpenDatabase(invtk)
    #AddPlot("Contour", elfvar)
    #c3 = ContourAttributes()
    #c3.colorType = c3.ColorByMultipleColors
    #c3.contourMethod = c3.Value
    #c3.contourValue = (0.01, 0.1, 0.26, 0.65, 0.85)
    #c3.SetMultiColor(0,(60,60,255,255))
    #c3.SetMultiColor(1,(111,255,79,255))
    #c3.SetMultiColor(2,(255,255,60,255))
    #c3.SetMultiColor(3,(255,208,74,255))
    #c3.SetMultiColor(4,(255,187,146,255))
    #c3.lineWidth= 1
    #c3.legendFlag = 0
    #SetPlotOptions(c3)
    
    #slice3()
    #rescmol(molzoom)

    ## Drawing
    #DrawPlots()

    ## Setting window attributes
    #winatt3()
    
    ##### Window 4 #### (plus, we destroy Window 5)
    #SetActiveWindow(5)
    
    ## We plot the Pseudocolor slice of the E-field and evaluate the value along a line 
    #DrawPlots()  

    ## This command created a curve plot in window 4
    #Lineout((ext[4],0,0),(ext[5],0,0))    

    ## Curve attributes
    #SetActiveWindow(4)
    #c4 = CurveAttributes()
    #c4.curveColor = (255, 255, 255, 255)
    #c4.lineWidth = 2
    #c4.showLegend = 0
    #c4.showLabels = 0
    #SetPlotOptions(c4)

    ## Setting window attributes
    #winatt4()
    
    ## We now delete window 5
    #SetActiveWindow(5)
    #DeleteWindow()
    
    #### All Windows ####
    if (i==0):
        zoomtstr = 0.
        zoomtend = 0.

    # Zooming on the molecule
    #if (i==0 and molzoom>1.):
        #zoomtstr = time.time() 
        #zoomsteps = 100
        #for zoomidx in range(0,zoomsteps+1):
            #if (float(molzoom)/float(zoomsteps)>1.):
                #thiszoom = 1.+( float(zoomidx)*(float(molzoom)/float(zoomsteps))-1. )
                #if (thiszoom==0):
                    #thiszoom = 1.
            #else:
                #thiszoom = 1.+  float(zoomidx)*(float(molzoom)-1.)/float(zoomsteps)
            ## We change window 1
            #SetActiveWindow(1)
            ## Changing the molecule
            #SetActivePlots(1)
            #mol = GetPlotOptions()
            #mol.radiusScaleFactor = thiszoom/2.
            #SetPlotOptions(mol)
            #RemoveLastOperator()
            #rescmol(thiszoom)
            ## Changing the contours of the molecule
            #SetActivePlots(2)
            #RemoveLastOperator()
            #rescmol(thiszoom)

            #lm = GetAnnotationObject("Text2D3")
            #lm.text = "Molecule x"+str(thiszoom)

            #DrawPlots()
            
            #if (zoomidx==zoomsteps):
                #ar1 = GetAnnotationObject("Line2D1")
                #ar1.visible = 0
            
            ## We keep window 2 the same
            
            ## We change window 3
            #SetActiveWindow(3)
            #SetActivePlots(1)
            #RemoveLastOperator()
            #rescmol(thiszoom)

            #DrawPlots()
            
            ## We keep window 4 the same
            
            ## Saving the image
            #if (zoomidx!=zoomsteps):
                #SaveWindow()
        #zoomtend = time.time() 

    # Saving the image
    #SetViewExtentsType(1)
    SaveWindow()

    # Clear and close
    for winind in range(1,5):
        SetActiveWindow(winind)
        DeleteAllPlots()
        RemoveAllOperators()
    #CloseDatabase(inefx)
    CloseDatabase(indensityvtk)
    #CloseDatabase(inxyz)
    #CloseDatabase(invtk)

    # Estimating remaining time
    print("Time left (sec): ",((time.time()-imgtstart-(zoomtend-zoomtstr))/(i+1))*(len(indata)-(i+1)))

# End of the loop over the input files

imgtend = time.time()
subprocess.call(["clear", ""])
print("Creation of the images completed.")

# Movie creation
if (len(moviename)!=0):
    vidtstart= time.time() 
    input_pattern = outdir+savename+"%04d.png"
    output_pattern = outdir+moviename
    encoding.encode(input_pattern,output_pattern,fdup=int(moviespeed))
    vidtend = time.time()


# Reporting working time
print("\n")
print(" *** Time statistics ***")
if (len(indata)!=0):
    print("Images creation time:")
    print(("--- %s seconds ---" % (imgtend - imgtstart-(zoomtend-zoomtstr))))
    print(( "On average %s seconds per image" % ((imgtend - imgtstart-(zoomtend-zoomtstr))/len(indata)) ))
    if (molzoom>1.):
        print("")
        print("Additional time for creating the zooming:")
        print(("--- %s seconds ---" % (zoomtend-zoomtstr)))
    print("")
print("Video creation time:")
if (len(moviename)!=0):
    print(("--- %s seconds ---" % (vidtend - vidtstart)))
else:
    print("--- Not Available ---")
print("")
print("Script total execution time:")
print(("--- %s seconds ---" % (time.time() - start_script)))
print("")

quit()

exit()

##################################################################
