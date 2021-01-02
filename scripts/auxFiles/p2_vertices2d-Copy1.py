import os, json
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
from tqdm import tqdm

# Open vector and raster files
vorMesh = gpd.read_file('../shps/voronoiGridRelaxed.shp')

# Get grid index
intervalNumber = 10
meshBounds = vorMesh.total_bounds
gridXarray = np.linspace(meshBounds[0],meshBounds[2],intervalNumber+1)
gridYarray = np.linspace(meshBounds[1],meshBounds[3],intervalNumber+1)

totalVerticesList = []
cell2dArrays = []
polygonCentroidList = []
gridIndexList =[]

#defining function
def findIndex(var, coordArray):
    for interval in range(intervalNumber):
#        if var > coordArray[interval] and var < coordArray[interval+1]:
        if var >= coordArray[interval] and var < coordArray[interval+1]:
            return interval
            break

print('\nCreating a unique list of vertices [[x1,y1],[x2,y2],...]')
for index,row in tqdm(vorMesh.iterrows(), total= vorMesh.shape[0]):
    #vertices xy
    coords = row.geometry.exterior.coords.xy
    totalVerticesList += list(zip(coords[0],coords[1]))

uniqueVerticesArray = np.unique(np.array(totalVerticesList),axis=0)
uniqueVerticesList = uniqueVerticesArray.tolist()

vertexIndexDict = {}
for index, vertex in enumerate(uniqueVerticesList):
    srtVertex = str(vertex)
    vertexIndexDict[srtVertex]=index

print('\nExtracting cell2d data and grid index')
for index,row in tqdm(vorMesh.iterrows(), total= vorMesh.shape[0]):
    rowGeometry = row.geometry
    coords = rowGeometry.exterior.coords
    #cell2d array
    cellArray = []
    #add index
    cellArray.append(index)
    #add centroid
    cellArray += list(rowGeometry.centroid.coords[0])
    #working with vertices number and vertex
    vertexIndexList = []
    for vertex in coords:
        srtVertex = str(list(vertex))
        vertexIndexList.append(vertexIndexDict[srtVertex])
    cellArray.append(len(vertexIndexList))
    cellArray += vertexIndexList
    cell2dArrays.append(cellArray)
    #get grid index
    xmin = min(coords.xy[0])
    xmax = max(coords.xy[0])
    ymin = min(coords.xy[1])
    ymax = max(coords.xy[1])
    xInterBeg = findIndex(xmin,gridXarray)
    xInterEnd = findIndex(xmax,gridXarray)
    yInterBeg = findIndex(ymin,gridYarray)
    yInterEnd = findIndex(ymax,gridYarray)
    gridIndexList.append([[xInterBeg,xInterEnd],[yInterBeg,yInterEnd]])

uniqueVerticesArray = np.unique(np.array(totalVerticesList),axis=0)
uniqueVerticesList = uniqueVerticesArray.tolist()
indexedVerticesList = [[index, row[0], row[1]] for index, row in enumerate(uniqueVerticesList)]

disvDict = {}
disvDict['NCPL'] = len(vorMesh.index)
disvDict['NVERT'] = len(uniqueVerticesList)
disvDict['uniqueVerticesList']=uniqueVerticesList
disvDict['indexedVerticesList']=indexedVerticesList
disvDict['cell2dArrays'] = cell2dArrays

spatialIndexDict = {}
spatialIndexDict['intervalNumber'] = intervalNumber
spatialIndexDict['gridXarray'] = list(gridXarray)
spatialIndexDict['gridYarray'] = list(gridYarray)
spatialIndexDict['gridIndexList'] = gridIndexList

with open('../txt/disvDict.json', 'w') as outf:
    json.dump(disvDict, outf)

with open('../txt/spatialIndexDict.json', 'w') as outf:
    json.dump(spatialIndexDict, outf)
