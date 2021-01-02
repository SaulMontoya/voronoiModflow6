import geopandas as gpd
from workFunctions import timeit
from shapely.ops import linemerge
import matplotlib.pyplot as plt

riversDf =  gpd.read_file('../shps/rios.shp')
meshDf = gpd.read_file('../shps/voronoiGrid.shp')
riversDf = riversDf.dissolve(by='TIPO')
indexList = []

#RCH package

def copyBaseFile(basePath,workPath,key,value):
    baseFile = open(basePath,'r')
    baseRead = baseFile.read()
    workFile = open(workPath,'w')
    baseRead = baseRead.replace(key,value)
    workFile.write(baseRead)
    workFile.close()

copyBaseFile('../model/baseModel/flow.rch','../model/flow.rch','%%MAXBOUND%%',str(len(meshDf.index)))
copyBaseFile('../model/baseModel/flow.evt','../model/flow.evt','%%MAXBOUND%%',str(len(meshDf.index)))
@timeit
def getPolygons(**kwargs):
    indexList = []
    for polyIndex,polyRow in meshDf.iterrows():
        for riverIndex, riverRow in riversDf.iterrows():
            if riverRow.geometry.crosses(polyRow.geometry):
                indexList.append(polyIndex)
    return indexList

#print(riversDf.head())
#riversDf = riversDf.dissolve(by='TIPO')
#riversDf.plot()
#plt.show()

crossedPolygons = getPolygons()

baseDrnFile = open('../model/baseModel/flow.drn','r')
baseDrn = baseDrnFile.read()
workDrnFile = open('../model/flow.drn','w')
baseDrn = baseDrn.replace('%%MAXBOUND%%',str(len(crossedPolygons)))
workDrnFile.write(baseDrn)
workDrnFile.close()

cellTopFile = open('../model/disv_mtop.disv','r').read()
cellTopList = cellTopFile.split(' ')

drnCells = open('../model/drnCells.drn','w')
for index, cellTop in enumerate(cellTopList):
    if index in crossedPolygons:
        drnCells.write('%d %d 0.001 \n'%(index+1, int(cellTop)-1))
#workDrnFile = open('../model/drnCells.drn','w')
drnCells.close()
