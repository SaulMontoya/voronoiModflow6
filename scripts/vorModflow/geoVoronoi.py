import numpy as np
import tqdm, time
from datetime import datetime
from vorModflow.lloydRelax import Field
#import geopandas as gpd
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
from scipy.spatial import Voronoi,cKDTree
#import geospatial libraries
import fiona
from tqdm import tqdm
from shapely.ops import split
from shapely.geometry import LineString, Polygon, mapping
from collections import OrderedDict

class createVoronoi():
    def __init__(self):
        self.discArrays = {}
        self.modelDis = {}
        self.pairArray = None

    def addLimitLayer(self, name, shapePath):
        #shape = gpd.read_file(shapePath)
        limitShape = fiona.open(shapePath)
        limitGeom = Polygon(limitShape[0]['geometry']['coordinates'][0])
        limitBounds = limitGeom.bounds
        self.modelDis['xMin'], self.modelDis['xMax'] = [limitBounds[i] for i in [0,2]]
        self.modelDis['yMin'], self.modelDis['yMax'] = [limitBounds[i] for i in [1,3]]
        self.modelDis['xDim'] = limitBounds[2] - limitBounds[0]
        self.modelDis['yDim'] = limitBounds[3] - limitBounds[1]
        #self.modelDis['refDim'] = np.sqrt([self.modelDis['xDim']**2+self.modelDis['yDim']**2])
        self.modelDis['limitShape'] = limitShape
        self.modelDis['limitGeometry'] = limitGeom
        self.modelDis['crs'] = limitShape.crs

    def defineParameters(self, maxRef = 500, refLevel=5, refTimes = 2):
        self.modelDis['refLevel'] = refLevel
        self.modelDis['maxRef'] = maxRef
        totalRefSizeList = [maxRef/refLevel**i for i in range(refTimes+1)]
        self.modelDis['refSizeList'] = totalRefSizeList[:-1]
        self.modelDis['minRef'] = totalRefSizeList[-1]
        print('\n/--------Sumary of cell discretization-------/')
        print('Maximun refinement: %.2f m.'%self.modelDis['maxRef'])
        print('Minimum refinement: %.2f m.'%self.modelDis['minRef'])
        print('/--------------------------------------------/\n',flush=True)
        #refSizeList = [self.modelDis['maxRef'],self.modelDis['maxRef']/5]


    def addDiscretizationLayer(self, name, shapePath):
        #shape = gpd.read_file(shapePath)
        shape = fiona.open(shapePath)
        self.discArrays[name] = shape

    def vertexAsList(self,shape):
        def getEquidistantPoints(p1, p2, parts):
            return zip(np.linspace(p1[0], p2[0], parts+1),
                       np.linspace(p1[1], p2[1], parts+1))
        vertexList = []
        #for line in shape.iterrows():
        for line in shape:
            partialList = []

            if line['geometry']['type'] == 'Polygon':
                lineGeom = Polygon(line['geometry']['coordinates'][0])
                pointObject = lineGeom.exterior.coords.xy
                pointList = list(zip(pointObject[0],pointObject[1]))
                for index, point in enumerate(pointList):
                    partialList.append(point)
            elif line['geometry']['type'] == 'LineString':
                lineGeom = LineString(line['geometry']['coordinates'])
                pointObject = lineGeom.coords.xy
                pointList = list(zip(pointObject[0],pointObject[1]))
                distDiag = cdist(pointList,pointList,'euclidean').diagonal(1)
                for index, point in enumerate(pointList):
                    lastIndex=len(pointList)-1
                    if  index != lastIndex:
                        if distDiag[index] < self.modelDis['minRef']:
                            partialList.append(point)
                        elif distDiag[index] > self.modelDis['minRef']:
                            parts = int(distDiag[index]//self.modelDis['minRef'])
                            distribPoints = list(getEquidistantPoints(pointList[index], pointList[index+1], parts))
                            partialList+=distribPoints
                    elif index == lastIndex:
                        partialList.append(point)
            else:
                pass
            vertexList += partialList
        return vertexList

    def pairArrayGenerator(self):
        #print("\n Start pair array generation \n")
        pairArrayList = []
        nPointX = int(self.modelDis['xDim']/self.modelDis['maxRef'])
        nPointY = int(self.modelDis['yDim']/self.modelDis['maxRef'])
        xMax = self.modelDis['xMax']+2*self.modelDis['maxRef']
        yMax = self.modelDis['yMax']+2*self.modelDis['maxRef']
        xArray = np.arange(self.modelDis['xMin'], xMax, step=self.modelDis['maxRef'])
        yArray = np.arange(self.modelDis['yMin'], yMax, step=self.modelDis['maxRef'])
        #refSizeList = [self.modelDis['maxRef'],self.modelDis['maxRef']/5]

        def griddedPoint(xArray,yArray,refSize, arrayType):
            tempPointList = []
            xOutList = []
            yOutList = []
            refInterval = refSize*(self.modelDis['refLevel']//2)/self.modelDis['refLevel']
            refArray = np.linspace(-refInterval, refInterval, num=self.modelDis['refLevel'])

            if arrayType=='1d':
                for x, y in tqdm(list(zip(xArray, yArray)), desc='Filling points with refinement of %.2f m'%refSize):
                        distance = self.tree.query([x,y])[0]

                        #if distance < refSize:
                        if distance < 2*refSize:
                            for xRef in refArray:
                                for yRef in refArray:
                                    if [xRef,yRef] == [0,0]: #Adds the same point
                                        xOutList.append(x)
                                        yOutList.append(y)
                                    elif [xRef,yRef] != [0,0]:
                                        refPoint = [x+xRef,y+yRef]
                                        distance = self.tree.query(refPoint)[0]
                                        #if distance < refSize*2:
                                        if distance < refSize:
                                            tempPointList.append(refPoint)
                                            xOutList.append(x+xRef)
                                            yOutList.append(y+yRef)
                                        else:
                                            pass


            elif arrayType=='2d':
                with tqdm(total=xArray.shape[0] * yArray.shape[0], desc='Filling points with refinement of %.2f m'%refSize) as pbar:
                    for x in xArray:
                        for y in yArray:
                            tempPointList.append([x,y])
                            distance = self.tree.query([x,y])[0]

                            for xRef in refArray:
                                for yRef in refArray:
                                    refPoint = [x+xRef,y+yRef]
                                    distance = self.tree.query(refPoint)[0]
                                    if [xRef,yRef] == [0,0]: #Adds the same point
                                        xOutList.append(x)
                                        yOutList.append(y)
                                    elif distance < refSize and [xRef,yRef] != [0,0]:
                                        tempPointList.append(refPoint)
                                        xOutList.append(x+xRef)
                                        yOutList.append(y+yRef)
                                    else:
                                        pass
                            pbar.update(1)
            else:
                pass


            xOutArray = np.array(xOutList)
            yOutArray = np.array(yOutList)
            return tempPointList, xOutArray, yOutArray

        for index, refSize in enumerate(self.modelDis['refSizeList']):
            #print("Start generation \n")
            if index == 0:
                tempPointList, xOutArray, yOutArray = griddedPoint(xArray,yArray,refSize,'2d')
            else:
                tempPointList, xOutArray, yOutArray = griddedPoint(xArray,yArray,refSize,'1d')
            #to debug point in the loop
            #tempPointArray = np.array(tempPointList)
            #np.savetxt('../txt/'+str(datetime.now().time().second)+'.txt',tempPointArray)

            pairArrayList += tempPointList
            xArray = np.copy(xOutArray)
            yArray = np.copy(yOutArray)

        self.pairArray = np.array(pairArrayList+self.modelDis['limitPairList'])



    def domainPairArray(self, txtFile='', probIndex=0.01):
        start = time.time()
        partialPairList = []
        #self.modelDis['refDist']=refDist
        #self.modelDis['maxRef']=maxRef
        #self.modelDis['minRef']=minRef
        #Add all shapes to the tuple list
        for key, shape in self.discArrays.items():
            #print(self.vertexAsList(shape))
            partialPairList += self.vertexAsList(shape)
        #adding limit vertex
        limitPairList = self.vertexAsList(self.modelDis['limitShape'])
        self.modelDis['limitPairList'] = limitPairList
        self.modelDis['partialPairList'] = partialPairList


        #Definitio of distance tree
        #print(self.modelDis['partialPairList'])
        self.tree = cKDTree(self.modelDis['partialPairList'])

        #Geneation of total pair array
        self.pairArrayGenerator()

        if txtFile != '':
            np.savetxt(txtFile,self.pairArray)

        print('\n/----Sumary of points for voronoi meshing----/')
        print('Points from limit layers: %d'%len(limitPairList))
        print('Points from definition layers: %d'%len(partialPairList))
        print('Total number of generated points: %d'%self.pairArray.shape[0])
        print('/--------------------------------------------/')
        end = time.time()
        print('\nTime required for point generation: %.2f seconds'%(end - start), flush=True)

    def relaxVertices(self, relaxLevel=3):
        start = time.time()
        for i in range(relaxLevel):
            field = Field(self.pairArray)
            field.relax()
            self.pairArray = field.get_points()
        end = time.time()
        print('\nTime required for Voronoi relaxation: %.2f seconds'%(end - start), flush=True)

    def createClipVoronoi(self,shapePath=''):
        start = time.time()
        vorNew = Voronoi(self.pairArray)
        polygonList = []
        for ptGroup in vorNew.regions:
            ptGroupList = []
            for pt in ptGroup:
                if pt != -1:
                    ptGroupList.append(tuple(vorNew.vertices[pt]))
            if len(ptGroupList) >= 3:
                polygonList.append(Polygon(ptGroupList))

        clipPolyList = []
        for poly in polygonList:
            limitGeometry = self.modelDis['limitGeometry']
            result = limitGeometry.intersection(poly)
            if result.area:
                clipPolyList.append(result)

        schema_props = OrderedDict([("cell", "int")])
        schema={"geometry": "Polygon", "properties": schema_props}
        if shapePath != '':
            outFile = fiona.open(shapePath,mode = 'w',driver = 'ESRI Shapefile',
                                crs = self.modelDis['crs'], schema=schema)
            #for index, poly in enumerate(polygonList):
            for index, poly in enumerate(clipPolyList):
                feature = {
                    "geometry": mapping(poly),
                    "properties": OrderedDict([("cell",index)]),
                }
                outFile.write(feature)
            outFile.close()
            end = time.time()

        clipPolyAreaList =[poly.area for poly in clipPolyList]
        clipPolyAreaArray = np.array(clipPolyAreaList)
        print('\n/---------Sumary of voronoi meshing----------/')
        print('Total number of cells: %d'%len(clipPolyAreaList))
        print('Area of 95 percentile cell: %.2f m2'%np.percentile(clipPolyAreaArray,95))
        print('Area of 50 percent (median) cell: %.2f m2'%np.percentile(clipPolyAreaArray,50))
        print('Area of 5 percentile cell: %.2f m2'%np.percentile(clipPolyAreaArray,5))
        print('/---------------------------------------------/', flush=True)
        print('\nTime required for Voronoi clipping: %.2f seconds'%(end - start))

        fig = plt.figure(figsize=(12,6))
        for poly in clipPolyList:
            #print(poly.geom_type)
            if poly.geom_type=='Polygon':
                plt.plot(*poly.exterior.xy)
            elif poly.geom_type=='MultiPolygon':
                    tempPolyList = list(poly)
                    for tempPoly in tempPolyList:
                        plt.plot(*tempPoly.exterior.xy)
            else: pass
        return fig
