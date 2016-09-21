##################################################################################
## Copyright (C) 2016 M. Compostella, M. Rampp, H. Appel
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


# This script creates a 3D plot of the scatter of a proton off a molecule of C2H4.


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
# elfcube       - Name of the .cube file containing the Octopus results and the position of the elements of the molecule. 
#                 The script assumes the variables named 'electron_density' and 'elements' are present in this file.
#
# bonds         - Set to True if you want to plot the bonds of the molecule. Note that bonds are created according to the
#                 distance between the elements of the molecule.
#
# usecube       - You can choose to use the .vtk+.xyz input files (usecube=False) or the .cube file (usecube=True) by setting 
#                 this variable.
#


# Set the variable below as DEFAULT values

inputdata = '/home/user/inputdata.txt'
outdir    = '/home/user/output/'

elfvtk    = 'elf_rs.vtk'
geometry  = 'geometry.xyz'
elfcube   = 'elf_rs.cube'

bonds     = False
usecube   = False



#################################################################################
# !!!!!!!!! You are not supposed to modify what is below this point !!!!!!!!!!! #
#################################################################################

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

print " *** Welcome to C2H4_3D.py, a visualization script for VisIt ***"
print ""
print "Description: This script creates a set of 3D plots of the scatter"
print "             of a proton off a molecule of C2H4."
print ""
print "Enter your choice:"
print "1- Use default value for the setup variables"
print "2- Specify the location of the file containing the list of input"
print "3- Specify single inputs"
answer = raw_input("\n")

if (answer=="1"):
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
    zoom      = 2.5
    plottype  = '3'
    backgrid  = False
    rotation  = False
    keeprot   = False
    savename  = 'mol3D'
    moviename = ""
    subprocess.call(["clear", ""])
elif (answer=="2"):
    print "Enter full path to the file containing the list of input:"
    inputdata = raw_input("(e.g: /home/user/inputdata.txt)\n")
    fileex = os.path.exists(inputdata)
    if (fileex==False):
        print "File does not exists!"
        quit()
    with open(inputdata, 'r') as infile:
        indata = infile.read().split('\n')
    outdir = raw_input("Enter output directory: (press Enter for current directory)\n")
    if (len(outdir)==0):
        outdir = subprocess.check_output(["pwd"]).rstrip("\n")+"/"
    elif (outdir[len(outdir)-1]!="/"):
        outdir=outdir+"/"
    outex = os.path.exists(outdir)
    if (outex==False):
        print "Directory does not exists!"
        quit()
    print "Enter zoom value:"
    zoom = raw_input("Suggested value for this plot is 2.5\n")
    if (len(zoom)==0):
        zoom = 2.5
    print "Choose the plot type:"
    print "1 - Volume plot"
    print "2 - Volume plot with ray casting"
    plottype = raw_input("3 - Multi-contour plot\n")
    if (plottype!='1' and plottype!='2'):
        plottype = 3
    plotgrid = raw_input("Do you want to plot a background grid? (y/n)\n")
    if (plotgrid=="y"):
        backgrid=True
    else:
        backgrid=False
    plotrot = raw_input("Do you want to plot an initial rotation around the molecule? (y/n)\n")
    if (plotrot=="y"):
        rotation=True
    else:
        rotation=False
    plotkeeprot = raw_input("Do you want to keep the camera rotating during the scattering? (y/n)\n")
    if (plotkeeprot=="y"):
        keeprot=True
    else:
        keeprot=False
    savename = raw_input("Enter name for output images:\n")
    if (len(savename)==0):
        savename = 'mol3D'
    print "Enter name for output movie with extension: (press Enter to skip movie creation)"
    moviename = raw_input("You can choose between .wmv, .mpg and .mov file formats\n")
    if (len(moviename)!=0):
        moviespeed = raw_input("Set movie speed (1 fast - 100 slow)\n")
    subprocess.call(["clear", ""])
elif (answer=="3"):
    inputdir = raw_input("Enter full path to the directory containing the input timesteps:\n")
    if (inputdir[len(inputdir)-1]!="/"):
        inputdir=inputdir+"/"
    inputstart = raw_input("Enter the name of the first input (e.g. td.0000000/):\n")
    if (inputstart[len(inputstart)-1]!="/"):
        inputstart=inputstart+"/"
    inputend = raw_input("Enter the name of the last input (e.g. td.0030000/):\n")
    if (inputend[len(inputend)-1]!="/"):
        inputend=inputend+"/"
    print "Note that we assume that each input directory contains"
    print "a elf_rs.vtk, a geometry.xyz and a elf_rs.cube file"
    inputstride = raw_input("Enter the timestep distance between inputs (e.g. 15):\n")
    if (int(inputstride)==0):
        totinput=1
    else:
        totinput=(int( inputend[len(inputend)-8:len(inputend)-1] )-int( inputstart[len(inputstart)-8:len(inputstart)-1] ))/int(inputstride) +1
    indata = [None]*(totinput)
    outdir = raw_input("Enter output directory: (press Enter for current directory)\n")
    if (len(outdir)==0):
        outdir = subprocess.check_output(["pwd"]).rstrip("\n")+"/"
    elif (outdir[len(outdir)-1]!="/"):
        outdir=outdir+"/"
    outex = os.path.exists(outdir)
    if (outex==False):
        print "Directory does not exists!"
        quit()
    print "Enter zoom value:"
    zoom = raw_input("Suggested value for this plot is 2.5\n")
    if (len(zoom)==0):
        zoom = 2.5
    print "Choose the plot type:"
    print "1 - Volume plot"
    print "2 - Volume plot with ray casting"
    plottype = raw_input("3 - Multi-contour plot\n")
    if (plottype!='1' and plottype!='2'):
        plottype = 3
    plotgrid = raw_input("Do you want to plot a background grid? (y/n)\n")
    if (plotgrid=="y"):
        backgrid=True
    else:
        backgrid=False
    plotrot = raw_input("Do you want to plot an initial rotation around the molecule? (y/n)\n")
    if (plotrot=="y"):
        rotation=True
    else:
        rotation=False
    plotkeeprot = raw_input("Do you want to keep the camera rotating during the scattering? (y/n)\n")
    if (plotkeeprot=="y"):
        keeprot=True
    else:
        keeprot=False
    savename = raw_input("Enter name for output images:\n")
    if (len(savename)==0):
        savename = 'mol3D'
    print "Enter name for output movie with extension: (press Enter to skip movie creation)"
    print "Note that you need an encoder like ffmpeg on this machine in order to be able to create a movie"
    moviename = raw_input("You can choose between .wmv, .mpg and .mov file formats\n")
    if (len(moviename)!=0):
        moviespeed = raw_input("Set movie speed (1 fast - 100 slow)\n")
    subprocess.call(["clear", ""])
    for i in range(0,totinput):
        step=int( inputstart[len(inputstart)-8:len(inputstart)-1] )+int(inputstride)*i
        indata[i]=inputdir + inputstart[0:len(inputend)-8] +"{0:07d}".format(step)+"/"
else:
    print "Value not recognized."
    quit()

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
print "-Zoom value: ", zoom
if (plottype=='1'):
    print "-Volume plot will be created"
elif (plottype=='2'):
    print "-Volume plot with ray casting will be created"
else:
    print "-Contour plot will be created"
if (backgrid==True):
    print "-Plot will have a background grid"
else:
    print "-Plot will not have a background grid"
if (rotation==True):
    print "-Initial rotation around the molecule will be displayed"
else:
    print "-Initial rotation around the molecule will not be displayed"
if (keeprot==True):
    print "-Camera will keep rotating around the molecule during scattering"
else:
    print "-Camera will be static during the scattering"
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
# Setting annotations
a = AnnotationAttributes()
a.userInfoFlag = 0
a.databaseInfoFlag = 0
a.axes3D.visible = 0
a.axes3D.triadFlag = 0
a.axes3D.bboxFlag = 0
SetAnnotationAttributes(a)

# Setting some visual and rendering options
if (plottype=='1'):
    InvertBackgroundColor()
    r = GetRenderingAttributes()
    r.scalableActivationMode = r.Always
    r.specularFlag = 1
    r.specularCoeff = 0.3 #Strength
    r.specularPower = 40  #Sharpness
    r.specularColor = (255, 255, 255, 255)
    SetRenderingAttributes(r)
    cmapname=GetActiveContinuousColorTable()
    if (answer!="1"):
        print "Current color table is "+cmapname
        cmans = raw_input("Do you want to change it? (y/n)\n")
        if (cmans=="y"):
            print "The available color tables are:"
            for ctnames in ColorTableNames():
                print "- "+ctnames
            print ""
            print "Note: if you want to define a new color table, please refer to"
            print "      the VisIt manual."
            cmapname = raw_input("Enter the new color table:\n")
elif (plottype=='2'):
    InvertBackgroundColor()
    r = GetRenderingAttributes()
    r.scalableActivationMode = r.Always
    r.specularFlag = 1
    r.specularCoeff = 0.3 #Strength
    r.specularPower = 40  #Sharpness
    r.specularColor = (255, 255, 255, 255)
    SetRenderingAttributes(r)
    cmapname=GetActiveContinuousColorTable()
    if (answer!="1"):
        print "Current color table can be change only in the script."
        cmans = raw_input("Press enter to continue.\n")
else:
    InvertBackgroundColor()
    r = GetRenderingAttributes()
    r.scalableActivationMode = r.Always
    r.doShadowing = 1
    r.shadowStrength = 0.3
    r.specularFlag = 1
    r.specularCoeff = 0.3 #Strength
    r.specularPower = 40  #Sharpness
    r.specularColor = (255, 255, 255, 255)
    SetRenderingAttributes(r)
    if (answer!="1"):
        print "This multi-contour plot has a predefined color table"
        cmans = raw_input("Do you want to change it? (y/n)\n")
        if (cmans=="y"):
            print "The available color tables are:"
            for ctnames in ColorTableNames():
                print "- "+ctnames
            print ""
            print "Note: if you want to define a new color table, please refer to"
            print "      the VisIt manual."
            cmapname = raw_input("Enter the new color table:\n")
            lcontval = [None]*GetColorTable(cmapname).GetNumControlPoints()
            colval   = [None]*GetColorTable(cmapname).GetNumControlPoints()
            for numindex in range(0,len(lcontval)):
                lcontval[numindex] = GetColorTable(cmapname).GetControlPoints(numindex).position
                colval[numindex]   = ( GetColorTable(cmapname).GetControlPoints(numindex).colors[0], GetColorTable(cmapname).GetControlPoints(numindex).colors[1], GetColorTable(cmapname).GetControlPoints(numindex).colors[2], float(20.+float(numindex)*(80-20)/float(len(lcontval)-1)) )
            lcontval[0]=lcontval[0]+0.01
            lcontval[len(lcontval)-1]=lcontval[len(lcontval)-1]-0.01
            # Converting list in tuple
            contval = tuple(lcontval)
        else:
            contval = (0.99, 0.9, 0.75, 0.5, 0.1)
            colval  = [None]*len(contval)
            colval[0] = (102,   0,   0, 89) #dark red 35%
            colval[1] = (204, 255, 229, 51) #light green 20%
            colval[2] = (255, 128,   0, 38) #orange 15%
            colval[3] = (192, 192, 192, 25) #gray 10%
            colval[4] = (255, 255, 255, 25) #white 10%

# Setting saving options
s = SaveWindowAttributes()
s.format = s.PNG
s.fileName = savename
s.width, s.height = 1024,1024
s.resConstraint = 0
if (changedir==True):
   s.outputToCurrentDirectory = 0
   s.outputDirectory = outdir
else:
   s.outputToCurrentDirectory = 1
SetSaveWindowAttributes(s)

# Creating a label for future use
l = CreateAnnotationObject("Text2D")
l.text = "step # 0000000"
l.position = (0.05, 0.05)
l.height = 0.03
l.fontBold = 1

# Camera movements

# Defining a box
def AddBox():
    AddOperator("Box", 0)
    b = BoxAttributes()
    b.amount = 1
    b.minx = -7
    b.maxx = 7
    b.miny = -7
    b.maxy = 7
    b.minz = -7
    b.maxz = 7
    SetOperatorOptions(b)

# Background grid
def grid():
    AddPlot("Mesh", "mesh")
    m = MeshAttributes()
    m.legendFlag = 0
    SetPlotOptions(m)
    AddOperator("Resample", 0)
    res = ResampleAttributes()
    res.samplesX = 10
    res.samplesY = 10
    res.is3D = 0
    SetOperatorOptions(res)

def rotatecamera():
    # Create the control points for the views.
    c0 = View3DAttributes()
    c0.viewNormal = (0, 0, 1)
    c0.focus = (0, 0, 0)
    c0.viewUp = (0, 1, 0)
    c0.viewAngle = 30
    c0.parallelScale = 12.
    c0.nearPlane = -24.
    c0.farPlane = 24.
    c0.imageZoom = float(zoom)
    c0.perspective = 1
 
    c1 = View3DAttributes()
    c1.viewNormal = (-1, 0, 0)
    c1.focus = (0, 0, 0)
    c1.viewUp = (0, 1, 0)
    c1.viewAngle = 30
    c1.parallelScale = 12.
    c1.nearPlane = -24.
    c1.farPlane = 24.
    c1.imageZoom = float(zoom)
    c1.perspective = 1
 
    c2 = View3DAttributes()
    c2.viewNormal = (0, 0, -1)
    c2.focus = (0, 0, 0)
    c2.viewUp = (0, 1, 0)
    c2.viewAngle = 30
    c2.parallelScale = 12.
    c2.nearPlane = -24.
    c2.farPlane = 24.
    c2.imageZoom = float(zoom)
    c2.perspective = 1
 
    c3 = View3DAttributes()
    c3.viewNormal = (1, 0, 0)
    c3.focus = (0, 0, 0)
    c3.viewUp = (0, 1, 0)
    c3.viewAngle = 30
    c3.parallelScale = 12.
    c3.nearPlane = -24.
    c3.farPlane = 24.
    c3.imageZoom = float(zoom)
    c3.perspective = 1
 
    c4 = c0
 
    # Create a tuple of camera values and x values. The x values are weights
    # that help to determine where in [0,1] the control points occur.
    cpts = (c0, c1, c2, c3, c4)
    x=[]
    for i in range(5):
        x = x + [float(i) / float(4.)]
 
    # Animate the camera. Note that we use the new built-in EvalCubicSpline
    # function which takes a t value from [0,1] a tuple of t values and a tuple
    # of control points. In this case, the control points are View3DAttributes
    # objects that we are using to animate the camera but they can be any object
    # that supports +, * operators.
    nsteps = 100
    for i in range(nsteps):
        t = float(i) / float(nsteps - 1)
        c = EvalCubicSpline(t, x, cpts)
        SetView3D(c)
        if (i!=nsteps-1):
            SaveWindow()
    
# Cycling on all input files
imgtstart = time.time()
#for i in range(0,1):     # Use only the first input for tests
for i in range(0,len(indata)):
    activepath = indata[i]
    if (activepath[len(activepath)-1]!="/"):
        activepath = activepath+"/"
    invtk = activepath+elfvtk
    inxyz = activepath+geometry
    incube = activepath+elfcube
    step = activepath[len(activepath)-8:len(activepath)-1]

# Adding a label for the time step    
    l = GetAnnotationObject("Text2D1")
    l.text = "step # "+step

# Opening Database
    if (usecube==False):
        OpenDatabase(invtk)
        elfvar="Re"
    else:
        OpenMDServer()
        DBopt = GetDefaultFileOpenOptions("Cube")
        DBopt["ExtendVolumeByOneCell"] = 0
        SetDefaultFileOpenOptions("Cube", DBopt)
        OpenDatabase(incube, 0, "Cube")
        elfvar="electron_density"
    
    if (plottype=='1'):
        # Adding Volume plot of the elf
        AddPlot("Volume", elfvar)
        vol = VolumeAttributes()
        vol.legendFlag = 0
        vol.smoothData = 1
        vol.opacityAttenuation = 0.07
        vol.num3DSlices = 500
        vol.resampleTarget = 1000000
        vol.rendererType = vol.Texture3D
        vol.lowGradientLightingReduction = vol.Lowest
        vol.SetColorControlPoints(GetColorTable(cmapname))
        SetPlotOptions(vol)
    elif (plottype=='2'):
        AddPlot("Volume", elfvar)
        ray = VolumeAttributes()
        ray.legendFlag = 0
        ray.lightingFlag = 1
        ray.colorControlPoints.GetControlPoints(0).colors = (0, 255, 255, 255)
        ray.colorControlPoints.GetControlPoints(0).position = 0.436155
        ray.colorControlPoints.GetControlPoints(1).colors = (255, 102, 0, 255)
        ray.colorControlPoints.GetControlPoints(1).position = 0.536585
        ray.colorControlPoints.GetControlPoints(2).colors = (255, 255, 0, 255)
        ray.colorControlPoints.GetControlPoints(2).position = 0.75
        ray.colorControlPoints.GetControlPoints(3).colors = (255, 0, 0, 255)
        ray.colorControlPoints.GetControlPoints(3).position = 1
        ray.opacityAttenuation = 0.18
        ray.opacityMode = ray.FreeformMode
        ray.resampleTarget = 50000
        ray.freeformOpacity = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 4, 5, 7, 8, 11, 14, 19, 22, 23, 26, 28, 28, 29, 31, 32, 34, 34, 35, 37, 38, 40, 41, 43, 43, 43, 43, 44, 46, 47, 49, 49, 49, 49, 49, 49, 49, 49, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 49, 46, 43, 41, 40, 38, 37, 35, 35, 34, 31, 29, 29, 28, 25, 20, 16, 11, 8, 5, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 9, 14, 20, 25, 28, 31, 34, 37, 40, 42, 46, 49, 52, 55, 56, 59, 61, 61, 62, 64, 65, 67, 68, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 67, 65, 65, 64, 61, 59, 58, 55, 52, 49, 46, 43, 40, 35, 32, 29, 25, 19, 14, 11, 8, 5, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 8, 20, 31, 42, 49, 53, 55, 60, 65, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70)
        ray.smoothData = 1
        ray.samplesPerRay = 200
        ray.rendererType = ray.RayCasting
        SetPlotOptions(ray)
    else:
        # Adding 3D contour plot
        AddPlot("Contour", elfvar)
        # Contour plot option
        c = ContourAttributes()
        c.colorType = c.ColorByMultipleColors
        c.contourMethod = c.Value
        c.contourValue = contval
        for numpoint in range(0,len(contval)):
            c.SetMultiColor(numpoint, colval[numpoint])
        c.lineWidth= 1
        c.legendFlag = 0
        SetPlotOptions(c)
    
# Adding the background grid, if required
    if (backgrid == True):
        grid()
    
# Choosing the Database and adding the molecule
    if (usecube==False):
        # Opening .xyz Database
        OpenDatabase(inxyz)

# Adding molecule
    AddPlot("Molecule", "element")

# Molecule options
    mol = MoleculeAttributes()
    mol.drawAtomsAs = mol.SphereAtoms
    mol.scaleRadiusBy = mol.Atomic
    mol.radiusScaleFactor = 0.15
    mol.drawBondsAs = mol.CylinderBonds
    mol.bondRadius = 0.04
    mol.legendFlag = 0
    mol.atomSphereQuality = mol.Super
    mol.bondCylinderQuality = mol.Super
    SetPlotOptions(mol)

# Adding box operator
    AddBox()

# Creating bonds
    if (bonds==True):
        AddOperator("CreateBonds", 0)    
    
# Drawing
    DrawPlots()

# Setting the view and rotation around the molecule
    if (i==0):
        # Setting View
        wa = View3DAttributes()
        wa.viewNormal = (0, 0, 1)
        wa.focus = (0, 0, 0)
        wa.viewUp = (0, 1, 0)
        wa.viewAngle = 30
        wa.parallelScale = 12.
        wa.nearPlane = -24.
        wa.farPlane = 24.
        wa.imageZoom = float(zoom)
        wa.perspective = 1
        SetView3D(wa)
        if (rotation==True):
            rotstart = time.time()
            rotatecamera()
            rotend = time.time()
        else:
            rotstart = 0.
            rotend = 0.

# Keep rotating the camera during scattering
    if (keeprot==True):
        kr = View3DAttributes()
        x = math.sin(2.*math.pi*(float(i)/float(len(indata)-1)))
        z = math.cos(2.*math.pi*(float(i)/float(len(indata)-1)))
        kr.viewNormal = (-x, 0, z)
        kr.focus = (0, 0, 0)
        kr.viewUp = (0, 1, 0)
        kr.viewAngle = 30
        kr.parallelScale = 12.
        kr.nearPlane = -24.
        kr.farPlane = 24.
        kr.imageZoom = float(zoom)
        kr.perspective = 1
        kr.centerOfRotation = (0., 0., 0.)
        SetView3D(kr)

# Saving the image
    SaveWindow()
    
# Delay before cleaning
#    time.sleep(1)

# Clear and close
    DeleteAllPlots()
    if (usecube==False):
        CloseDatabase(inxyz)
        CloseDatabase(invtk)
    else:
        CloseDatabase(incube)

# Estimating remaining time
    print "Time left (sec): ",((time.time()-imgtstart-(rotend-rotstart))/(i+1))*(len(indata)-(i+1))

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
    print("--- %s seconds ---" % (imgtend - imgtstart-(rotend-rotstart)))
    print( "On average %s seconds per image" % ((imgtend - imgtstart-(rotend-rotstart))/len(indata)) )
    if (rotation==True):
        print ""
        print "Additional time for creating the rotation:"
        print("--- %s seconds ---" % (rotend-rotstart))
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

