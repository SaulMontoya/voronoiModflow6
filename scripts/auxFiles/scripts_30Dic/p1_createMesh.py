from vorModflow.geoVoronoi import createVoronoi
import matplotlib.pyplot as plt
import numpy as np
import time

#Create mesh object
vorMesh = createVoronoi()

#Define base refinement and refinement levels
vorMesh.defineParameters(maxRef = 500, refLevel=3, refTimes = 1)

#Open limit layers and refinement definition layers
vorMesh.addLimitLayer('basin','../shps/Angascancha_Basin_Extension.shp')
vorMesh.addDiscretizationLayer('river','../shps/rios.shp')

#Generate point pair array
vorMesh.domainPairArray(txtFile='../txt/pairArray.txt')

#Relaxation of voronoi polygons
vorMesh.relaxVertices()
#vorMesh.relaxVertices()
#vorMesh.relaxVertices()

#Create polygon object and clip to the limit layer
vorGraph = vorMesh.createClipVoronoi(shapePath='../shps/voronoiGridRelaxed.shp')