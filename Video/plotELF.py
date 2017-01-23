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

#densityvtk= 'density.vtk'
densityvtk= 'elf_rs.vtk'
DensityGsVtk= '/../td.0000000/elf_rs.vtk' #This should work for any run in our frame. 

#elfvtk    = 'elf_rs.vtk'
geometry  = 'geometry.vtk' #'geometry.xyz'

ShowCaption = 1 
bonds     = False

OriginX=0.
OriginY=0.
OriginZ=0. #0.6788749999999979/BohrRadius

MaxValue=0.5
ZoomValue=0.7 #1:Bulk #3:Slabs
EnablePerspective=0

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
from libPlotVideo import *

start_script = time.time()

# Setting initial variables or using default input
subprocess.call(["clear", ""])

print " *** Welcome to Azulene3D.py, a visualization script for VisIt ***"
print ""
print "Description: This script creates a set of 3D plots of a laser"
print "             pulse propagating towards Azulene C10H8."
print ""
print "Enter your choice:"
print "1- Use default value for the setup variables"
print "2- Specify the location of the file containing the list of input"
print "3- Specify single inputs"
answer = raw_input("\n")

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
molzoom = 1.
savename  = 'density'
moviename = "density.mpg"
moviespeed = '50'
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




#### Window 1 ####
SetActiveWindow(1)
#SetActivePlots(1)
SetViewExtentsType(1)

#=================================================
#> End spontaneous state
#=================================================

#InvertBackgroundColor()

#Creating a label for the time step
#l1s = CreateAnnotationObject("Text2D")
#l1s.text = "Time step"
#l1s.position = (0.05, 0.05)
#l1s.height = 0.025
#l1s.fontBold = 1
#l1s.textColor = (255, 255, 255, 255)
#l1s.useForegroundForTextColor = 1

##############################
# Cycling on all input files #
##############################
imgtstart = time.time()
#for i in range(0,1):     # Use only the first input for tests
for i in range(0,len(indata)):
    activepath = indata[i]
    if (activepath[len(activepath)-1]!="/"):
        activepath = activepath+"/"
    indensityGsVtk = activepath+DensityGsVtk #TJYD: we will plot the difference with ground state
    indensityvtk = activepath+densityvtk
    #invtk = activepath+elfvtk
    inxyz = activepath+geometry
    #inefx = activepath+maxe_x
    #inefy = activepath+maxe_y
    #inefz = activepath+maxe_z
    step = activepath[len(activepath)-8:len(activepath)-1]

    # We open all databases for the plots
    OpenDatabase(indensityvtk)
    elfvar="elf_rs" #Old simulations: Re. "Re"
    DefineScalarExpression("GroundStateData", "conn_cmfe(<"+indensityGsVtk+":"+elfvar+">, <mesh>)") #TODO: script the definition of this mesh
    OpenDatabase(indensityvtk)
    DefineScalarExpression("ExcitedStateData", "conn_cmfe(<"+indensityvtk+":"+elfvar+">, <mesh>)")
    # DefineScalarExpression("Difference", "Re") #-GroundStateData") #DEBUG
    DefineScalarExpression("Difference", "ExcitedStateData - GroundStateData") #Re: TD-data, GroundStateData: td_value - GS. This is the performed operation on the scalar data. 
    
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
    plot3Ddifference(inxyz, indensityvtk, step, ShowCaption, OriginX, OriginY, OriginZ, "X", MaxValue, ZoomValue, EnablePerspective, molvar)
    
    #### Window 2 ####
    #CloneWindow()
    #SetActivePlots(1)
    SetActiveWindow(2)
    #plot2Ddifference(inxyz, indensityvtk, step, ShowCaption, OriginX, OriginY, OriginZ, axis="SiBond", extreme=MaxValue, ZoomValue, EnablePerspective, molvar)
    plot3Ddifference(inxyz, indensityvtk, step, ShowCaption, OriginX, OriginY, OriginZ, "Y", MaxValue, ZoomValue, EnablePerspective, molvar)
    
    #### Window 3 ####
    #CloneWindow()
    #SetActivePlots(0)
    #SetActiveWindow(3)
    #plot2Ddifference(inxyz, indensityvtk, axis="Z", extreme=MaxValue, caption=0)
    
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
    CloseDatabase(indensityvtk)

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
