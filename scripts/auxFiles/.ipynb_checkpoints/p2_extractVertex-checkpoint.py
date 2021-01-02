import geopandas as gpd
import numpy as np
import json, rasterio

vorMesh = gpd.read_file('../shps/voronoiGrid.shp')
demRaster = rasterio.open('../rst/ASTGTM2_18S_xyz2.asc')

#vorMesh.set_index('FID',inplace=True)
print(vorMesh.index)

#print(vorMesh.head())
totalVerticesList = []

for index,row in vorMesh.iterrows():
    coords = row.geometry.exterior.coords.xy
    totalVerticesList += list(zip(coords[0],coords[1]))

uniqueVerticesArray = np.unique(np.array(totalVerticesList),axis=0)

uniqueVerticesList = uniqueVerticesArray.tolist()

vorMesh['verticesList'] = np.empty((len(vorMesh), 0)).tolist()
polygonVerticesList = []
polygonCentroidList = []

for index,row in vorMesh.iterrows():
    coords = row.geometry.exterior.coords.xy
    polygonCentroidList.append(list(row.geometry.centroid.coords[0]))
    vertexIndexList = []
    for vertex in list(zip(coords[0],coords[1])):
        vertexIndex = uniqueVerticesList.index(list(vertex))
        vertexIndexList.append(vertexIndex)
    polygonVerticesList.append(vertexIndexList)

disvDict = {}
disvDict['NCPL'] = len(vorMesh.index)
disvDict['NLAY'] = 6
disvDict['modelBottom'] = 3200
disvDict['thickRatio'] = [0.05,0.1,0.2,0.4,0.7,1]
disvDict['NVERT'] = len(uniqueVerticesList)
disvDict['uniqueVerticesList']=uniqueVerticesList
disvDict['polyCentroidList'] = polygonCentroidList
disvDict['polygonVerticesSequence'] = polygonVerticesList

disvDict['polyCentroidElevation'] = []
for centroid in disvDict['polyCentroidList']:
    valuesXY = demRaster.sample([centroid])
    disvDict['polyCentroidElevation'].append(int(list(valuesXY)[0][0]))

with open('../txt/disvDict.json', 'w') as outf:
    json.dump(disvDict, outf)


baseDisvFile = open('../model/baseModel/flow.disv','r')
baseDisv = baseDisvFile.read()
workDisvFile = open('../model/flow.disv','w')

baseDisv = baseDisv.replace('%%NCPL%%',str(disvDict['NCPL']))
baseDisv = baseDisv.replace('%%NLAY%%',str(disvDict['NLAY']))
baseDisv = baseDisv.replace('%%NVERT%%',str(disvDict['NVERT']))
workDisvFile.write(baseDisv)
workDisvFile.close()
