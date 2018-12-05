# -*- coding: mbcs -*-
# Do not delete the following import lines
from abaqus import *
from abaqusConstants import *
import __main__
import section
import regionToolset
import displayGroupMdbToolset as dgm
import mesh
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import part
import material
import assembly
import step
import interaction
import load
import optimization
import job
import sketch
import connectorBehavior
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    UNDEFORMED, ORIENT_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(
    renderStyle=WIREFRAME, )
session.viewports['Viewport: 1'].odbDisplay.materialOrientationOptions.setValues(
    axis1Color='#FF0000', axis2Color='#00FF00', axis3Color='#0000FF')
session.viewports['Viewport: 1'].odbDisplay.materialOrientationOptions.setValues(
    lineThickness=MEDIUM)
session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(
    edgeColorWireHide='#000000')
session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(
    edgeLineThickness=THIN)
session.viewports['Viewport: 1'].odbDisplay.superimposeOptions.setValues(
    renderStyle=WIREFRAME)
session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(
    renderStyle=SHADED, fillColor='#FFFFFF', translucency=ON)	
session.viewports['Viewport: 1'].view.setProjection(projection=PARALLEL)
session.viewports['Viewport: 1'].view.setValues(nearPlane=2.83742, 
    farPlane=5.75095, cameraPosition=(3.32453, 3.29065, 2.57238), 
    cameraUpVector=(-0.505782, -0.499829, 0.703104), cameraTarget=(0.6, 
    0.55, 0.7))
session.viewports['Viewport: 1'].view.fitView()		


