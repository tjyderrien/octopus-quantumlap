# Now adding plotting script remade from scratch. 
AddPlot("Pseudocolor", "Re", 1, 1)
AddOperator("Slice", 1)
SetActivePlots(0)
SetActivePlots(0)
AddOperator("Box", 1)
BoxAtts = BoxAttributes()
BoxAtts.amount = BoxAtts.Some  # Some, All
BoxAtts.minx = -20
BoxAtts.maxx = 70
BoxAtts.miny = -20
BoxAtts.maxy = 70
BoxAtts.minz = -20
BoxAtts.maxz = 70
BoxAtts.inverse = 0
SetOperatorOptions(BoxAtts, 1)
DrawPlots()


# THis is not the right one
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
gg

:q

########
SetViewExtentsType(1)
# Begin spontaneous state
View2DAtts = View2DAttributes()
View2DAtts.windowCoords = (-29.0284, 28.5951, -80.4078, 80.4078)
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
View2DAtts.windowCoords = (-29.0284, 28.5951, -80.4078, 80.4078)
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
DrawPlots()


#####3333


SetViewExtentsType(1)
# Begin spontaneous state
View2DAtts = View2DAttributes()
View2DAtts.windowCoords = (-29.0284, 28.5951, -80.4078, 80.4078)
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
View2DAtts.windowCoords = (-29.0284, 28.5951, -80.4078, 80.4078)
View2DAtts.viewportCoords = (0.2, 0.95, 0.15, 0.95)
View2DAtts.fullFrameActivationMode = View2DAtts.Auto  # On, Off, Auto
View2DAtts.fullFrameAutoThreshold = 100
View2DAtts.xScale = View2DAtts.LINEAR  # LINEAR, LOG
View2DAtts.yScale = View2DAtts.LINEAR  # LINEAR, LOG
View2DAtts.windowValid = 1
SetView2D(View2DAtts)
View3DAtts = View3DAttributes()
View3DAtts.viewNormal = (-0.193224, -0.873456, 0.446922)
View3DAtts.focus = (-0.216648, -0.216648, 0)
View3DAtts.viewUp = (-0.0920526, 0.469636, 0.878048)
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


PseudocolorAtts = PseudocolorAttributes()
PseudocolorAtts.scaling = PseudocolorAtts.Linear  # Linear, Log, Skew
PseudocolorAtts.skewFactor = 1
PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData  # OriginalData, CurrentPlot
PseudocolorAtts.minFlag = 0
PseudocolorAtts.min = 0
PseudocolorAtts.maxFlag = 0
PseudocolorAtts.max = 0.2
PseudocolorAtts.centering = PseudocolorAtts.Natural  # Natural, Nodal, Zonal
PseudocolorAtts.colorTableName = "PuBu"
PseudocolorAtts.invertColorTable = 0
PseudocolorAtts.opacityType = PseudocolorAtts.FullyOpaque  # ColorTable, FullyOpaque, Constant, Ramp, VariableRange
PseudocolorAtts.opacityVariable = ""
PseudocolorAtts.opacity = 1
PseudocolorAtts.opacityVarMin = 0
PseudocolorAtts.opacityVarMax = 1
PseudocolorAtts.opacityVarMinFlag = 0
PseudocolorAtts.opacityVarMaxFlag = 0
PseudocolorAtts.pointSize = 0.05
PseudocolorAtts.pointType = PseudocolorAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
PseudocolorAtts.pointSizeVarEnabled = 0
PseudocolorAtts.pointSizeVar = "default"
PseudocolorAtts.pointSizePixels = 2
PseudocolorAtts.lineType = PseudocolorAtts.Line  # Line, Tube, Ribbon
PseudocolorAtts.lineStyle = PseudocolorAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
PseudocolorAtts.lineWidth = 0
PseudocolorAtts.tubeDisplayDensity = 10
PseudocolorAtts.tubeRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
PseudocolorAtts.tubeRadiusAbsolute = 0.125
PseudocolorAtts.tubeRadiusBBox = 0.005
PseudocolorAtts.varyTubeRadius = 0
PseudocolorAtts.varyTubeRadiusVariable = ""
PseudocolorAtts.varyTubeRadiusFactor = 10
PseudocolorAtts.endPointType = PseudocolorAtts.None  # None, Tails, Heads, Both
PseudocolorAtts.endPointStyle = PseudocolorAtts.Spheres  # Spheres, Cones
PseudocolorAtts.endPointRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
PseudocolorAtts.endPointRadiusAbsolute = 1
PseudocolorAtts.endPointRadiusBBox = 0.005
PseudocolorAtts.endPointRatio = 2
PseudocolorAtts.renderSurfaces = 1
PseudocolorAtts.renderWireframe = 0
PseudocolorAtts.renderPoints = 0
PseudocolorAtts.smoothingLevel = 0
PseudocolorAtts.legendFlag = 1
PseudocolorAtts.lightingFlag = 1
SetPlotOptions(PseudocolorAtts)

