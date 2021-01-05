import os, json
import numpy as np
#import geopandas as gpd
import fiona
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
from tqdm import tqdm

vorMesh = fiona.open('../modelData/shps/voronoiGridRelaxed.shp')

# Get grid index
intervalNumber = 10
meshBounds = vorMesh.bounds
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
for index, row in enumerate(tqdm(vorMesh, total= len(vorMesh))):
    #vertices xy
    if len(row['geometry']['coordinates']) == 1:
        #print(row['geometry']['coordinates'])
        xyList = [[i[0],i[1]] for i in row['geometry']['coordinates'][0]]
        totalVerticesList += xyList
    elif len(row['geometry']['coordinates']) > 1:
        print(row['geometry']['coordinates'])
        print(index)
        for vertexList in row['geometry']['coordinates']:
            #print(vertexList)
            xyList = [[i[0],i[1]] for i in vertexList[0]]
            #print(xyList)
            totalVerticesList += xyList
    else:
        pass
    
uniqueVerticesArray = np.unique(np.array(totalVerticesList), axis=0)
uniqueVerticesList = uniqueVerticesArray.tolist()

vertexIndexDict = {}
for index, vertex in enumerate(uniqueVerticesList):
    strVertex = str(vertex)
    vertexIndexDict[strVertex]=index
    
    
print('\nExtracting cell2d data and grid index')
for index,row in enumerate(tqdm(vorMesh)):#.iterrows(), total= vorMesh.shape[0]):
    rowCoords = row['geometry']['coordinates'][0]
    rowPoly = Polygon(rowCoords)
    #print(index)
    #print(rowPoly.bounds)
    #coords = rowGeometry.exterior.coords
    #cell2d array
    cellArray = []
    #add index
    cellArray.append(index)
    #add centroid
    cellArray += list(rowPoly.centroid.coords[0])
    #working with vertices number and vertex
    vertexIndexList = []
    for vertex in rowCoords:
        #print(vertex)
        strVertex = str(list(vertex))
        #print(vertexIndexDict[strVertex])
        vertexIndexList.append(vertexIndexDict[strVertex])
    cellArray.append(len(vertexIndexList))
    cellArray += vertexIndexList
    cell2dArrays.append(cellArray)
    #get grid index
    xmin = rowPoly.bounds[0] #min(coords.xy[0])
    xmax = rowPoly.bounds[2] #max(coords.xy[0])
    ymin = rowPoly.bounds[1] #min(coords.xy[1])
    ymax = rowPoly.bounds[3] #max(coords.xy[1])
    xInterBeg = findIndex(xmin,gridXarray)
    xInterEnd = findIndex(xmax,gridXarray)
    yInterBeg = findIndex(ymin,gridYarray)
    yInterEnd = findIndex(ymax,gridYarray)
    gridIndexList.append([[xInterBeg,xInterEnd],[yInterBeg,yInterEnd]])

uniqueVerticesArray = np.unique(np.array(totalVerticesList),axis=0)
uniqueVerticesList = uniqueVerticesArray.tolist()
indexedVerticesList = [[index, row[0], row[1]] for index, row in enumerate(uniqueVerticesList)]

disvDict = {}
disvDict['NCPL'] = len(vorMesh)
disvDict['NVERT'] = len(uniqueVerticesList)
disvDict['uniqueVerticesList']=uniqueVerticesList
disvDict['indexedVerticesList']=indexedVerticesList
disvDict['cell2dArrays'] = cell2dArrays

spatialIndexDict = {}
spatialIndexDict['intervalNumber'] = intervalNumber
spatialIndexDict['gridXarray'] = list(gridXarray)
spatialIndexDict['gridYarray'] = list(gridYarray)
spatialIndexDict['gridIndexList'] = gridIndexList

with open('../modelData/txt/disvDict.json', 'w') as outf:
    json.dump(disvDict, outf)

with open('../modelData/txt/spatialIndexDict.json', 'w') as outf:
    json.dump(spatialIndexDict, outf)