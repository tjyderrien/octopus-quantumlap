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


# This script creates a 2D plot of the scatter of a proton off a molecule of C2H4.


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

inputdata = '/scratch/tjyderri/octopus/VisualizationWorkshop/instructions/C2H4_2D/inputdata_eos.txt'
outdir    = '/scratch/tjyderri/octopus/VisualizationWorkshop/instructions/C2H4_2D/output/'

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

print " *** Welcome to C2H4_2D.py, a visualization script for VisIt ***"
print ""
print "Description: This script creates a set of 2D plots of the scatter"
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
    limits    = [-5,5,-5,5]
    savename  = 'mol2D'
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
    print "Enter plot limits separated by spaces: (xmin xmax ymin ymax)"
    limits = raw_input("Suggested values for this plot are -5 5 -5 5\n").split(" ")
    if (len(limits)!=4):
        limits = [-5,5,-5,5]
    savename = raw_input("Enter name for output images:\n")
    if (len(savename)==0):
        savename = 'mol2D'
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
    print "Enter plot limits separated by spaces: (xmin xmax ymin ymax)"
    limits = raw_input("Suggested values for this plot are -5 5 -5 5\n").split(" ")
    if (len(limits)!=4):
        limits = [-5,5,-5,5]
    savename = raw_input("Enter name for output images:\n")
    if (len(savename)==0):
        savename = 'mol2D'
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
print "-Image limits: ",limits[0],"< X <", limits[1]
print "               ",limits[2],"< Y <", limits[3]
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
# Creating color table
ct = ((60,60,255,0.), (0,255,255,0.1), (65,255,255,0.16), (111,255,79,0.19), (255,255,60,0.26), (255,208,74,0.65), (255,187,146,0.85), (255,255,255,1))
ccpl = ColorControlPointList()
for pt in ct:
    p = ColorControlPoint()
    p.colors = (pt[0], pt[1], pt[2], 255)
    p.position = pt[3]
    ccpl.AddControlPoints(p)
AddColorTable("elf", ccpl)
cmapname="elf"
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

# Setting annotations
a = AnnotationAttributes()
a.userInfoFlag = 0
a.databaseInfoFlag = 0
a.axes2D.visible = 0
SetAnnotationAttributes(a)

# Setting View
wa = GetView2D()
wa.windowCoords = (float(limits[0]),float(limits[1]),float(limits[2]),float(limits[3]))
wa.viewportCoords = (0, 1, 0, 1)
wa.fullFrameActivationMode = wa.On
SetView2D(wa)
SetViewExtentsType("actual")
ToggleLockViewMode()
ToggleMaintainViewMode()

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
    
# Adding Psedocolor of the elf
    AddPlot("Pseudocolor", elfvar)

# Adding a slice operator to this plot
    AddOperator("Slice", 0)
    sli = SliceAttributes()
    sli.originType = sli.Point
    sli.originPoint = (0, 0, 0.075)
    sli.axisType = sli.ZAxis
    SetOperatorOptions(sli)
    
# Pseudocolor plot option
    p = PseudocolorAttributes()
    p.colorTableName = cmapname
    p.legendFlag = 0
    p.minFlag,p.min = 1, 0.
    p.maxFlag,p.max = 1, 1.
    p.scaling = 0
    SetPlotOptions(p)

# Adding contour at 0.8
    AddPlot("Contour", elfvar)

# Contour plot option
    c = ContourAttributes()
    c.colorType = c.ColorBySingleColor
    c.singleColor = (118, 42, 0, 255)
    c.contourMethod = 1
    c.contourValue = (0.8)
    c.lineWidth= 1
    c.legendFlag = 0
    SetPlotOptions(c)
    
# Adding a slice operator to this plot
    AddOperator("Slice", 0)
    sli = SliceAttributes()
    sli.originType = sli.Point
    sli.originPoint = (0, 0, 0.075)
    sli.axisType = sli.ZAxis
    SetOperatorOptions(sli)

# Choosing the Database and adding the plot
    if (usecube==False):
        # Opening .xyz Database
        OpenDatabase(inxyz)

# Adding molecule and projecting it to the plane
    AddPlot("Pseudocolor", "element")
    mol = PseudocolorAttributes()
    mol.colorTableName = "cpk_jmol"
    mol.minFlag,mol.min = 1, -1.
    mol.maxFlag,mol.max = 1, 118.
    mol.legendFlag = 0
    mol.pointSizePixels = 25
    mol.pointType = mol.Sphere
    mol.lineType = mol.Line
    mol.lineWidth = 1
    mol.renderSurfaces = 1
    mol.renderPoints = 1
    mol.smoothingLevel = 2
    SetPlotOptions(mol)

    AddOperator("Project", 0)
    pro = ProjectAttributes()
    SetOperatorOptions(pro)

# Creating bonds
    if (bonds==True):
        AddOperator("CreateBonds", 0)    

# Drawing
    DrawPlots()

#  Saving the image
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
    print "Time left (sec): ",((time.time()-imgtstart)/(i+1))*(len(indata)-(i+1))

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
    print("--- %s seconds ---" % (imgtend - imgtstart))
    print( "On average %s seconds per image" % ((imgtend - imgtstart)/len(indata)) )
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



