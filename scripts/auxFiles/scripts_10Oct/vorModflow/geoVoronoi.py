import numpy as np
from vorModflow.lloydRelax import Field
import geopandas as gpd
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
from scipy.spatial import Voronoi,cKDTree
from shapely.geometry import Polygon, Point

class createVoronoi():
    def __init__(self):
        self.arrays = {}
        self.geom = {}
        self.pairArray = None

    def addLimitLayer(self, name, shapePath):
        shape = gpd.read_file(shapePath)
        shapeBounds = shape.bounds.values[0]
        self.geom['xMin'], self.geom['xMax'] = [shapeBounds[i] for i in [0,2]]
        self.geom['yMin'], self.geom['yMax'] = [shapeBounds[i] for i in [1,3]]
        self.geom['xDim'] = shapeBounds[2] - shapeBounds[0]
        self.geom['yDim'] = shapeBounds[3] - shapeBounds[1]
        self.geom['refDim'] = np.sqrt([self.geom['xDim']**2+self.geom['yDim']**2])
        self.geom['limitGeometry'] = shape.iloc[0].geometry

    def addDefinitionLayer(self, name, shapePath):
        shape = gpd.read_file(shapePath)
        self.arrays[name] = shape

    def vertexAsList(self,shape):
        def getEquidistantPoints(p1, p2, parts):
            return zip(np.linspace(p1[0], p2[0], parts+1),
                       np.linspace(p1[1], p2[1], parts+1))
        vertexList = []
        for line in shape.iterrows():
            partialList = []
            if line[1].geometry.geom_type == 'Polygon':
                pointList = line[1].geometry.exterior.coords.xy
            elif line[1].geometry.geom_type == 'LineString':
                pointObject = line[1].geometry.coords.xy
                pointList = list(zip(pointObject[0],pointObject[1]))
                distDiag = cdist(pointList,pointList,'euclidean').diagonal(1)
                for index, point in enumerate(pointList):
                    lastIndex=len(pointList)-1
                    if  index != lastIndex:
                        if distDiag[index] < self.geom['minRef']:
                            partialList.append(point)
                        elif distDiag[index] > self.geom['minRef']:
                            parts = int(distDiag[index]//self.geom['minRef'])
                            distribPoints = list(getEquidistantPoints(pointList[index], pointList[index+1], parts))
                            partialList+=distribPoints
                    elif index == lastIndex:
                        partialList.append(point)
            else:
                pass
            vertexList += partialList
        return vertexList

    def pairArrayGenerator(self):
        print("\n Start pair array generation \n")
        pairArrayList = []
        nPointX = int(self.geom['xDim']/self.geom['maxRef'])
        nPointY = int(self.geom['yDim']/self.geom['maxRef'])
        xArray = np.arange(self.geom['xMin'], self.geom['xMax']+2*self.geom['maxRef'],
                            step=self.geom['maxRef'])
        yArray = np.arange(self.geom['yMin'], self.geom['yMax']+2*self.geom['maxRef'],
                            step=self.geom['maxRef'])
        #sizeArray = np.array([self.geom['maxRef']])#,self.geom['maxRef']/5])
        sizeArray = np.array([self.geom['maxRef'],self.geom['maxRef']/5])

        def griddedPoint(xArray,yArray,refSize, arrayType):
            tempPointList = []
            xOutList = []
            yOutList = []
            refArray = np.linspace(-refSize*0.4, refSize*0.4, num=5)
            i = 1

            if arrayType=='1d':
                print('Processing regular 1d array')
                for x, y in list(zip(xArray, yArray)):
                    distance = self.tree.query([x,y])[0]
                    if distance < refSize:
                        for xRef in refArray:
                            for yRef in refArray:
                                if [xRef,yRef] != [0,0]:
                                    refPoint = [x+xRef,y+yRef]
                                    distance = self.tree.query(refPoint)[0]
                                    if distance < refSize/2:
                                        tempPointList.append(refPoint)
                                        xOutList.append(x+xRef)
                                        yOutList.append(y+yRef)
                                        i += 1
                    else:
                        pass
                    if i%10000 == 0:
                        print('Processing %d points'%i)
            elif arrayType=='2d':
                print('Processing regular 2d array')
                for x in xArray:
                    for y in yArray:
                        tempPointList.append([x,y])
                        distance = self.tree.query([x,y])[0]
                        for xRef in refArray:
                            for yRef in refArray:
                                refPoint = [x+xRef,y+yRef]
                                distance = self.tree.query(refPoint)[0]
                                if [xRef,yRef] == [0,0]:
                                    xOutList.append(x+xRef)
                                    yOutList.append(y+yRef)
                                if distance < refSize/2 and [xRef,yRef] != [0,0]:
                                    tempPointList.append(refPoint)
                                    xOutList.append(x+xRef)
                                    yOutList.append(y+yRef)
                                    i += 1
                                else:
                                    pass
                        if i%10000 == 0:
                            print('Processing %d points'%i)

            xOutArray = np.array(xOutList)
            yOutArray = np.array(yOutList)
            return tempPointList, xOutArray, yOutArray

        i = 0
        for refSize in sizeArray:
            print("Start generation \n")
            if i == 0:
                tempPointList, xOutArray, yOutArray = griddedPoint(xArray,yArray,refSize,'2d')
            else:
                tempPointList, xOutArray, yOutArray = griddedPoint(xArray,yArray,refSize,'1d')
            pairArrayList += tempPointList
            xArray = np.copy(xOutArray)
            yArray = np.copy(yOutArray)
            print('Closing first ref number %i of size %d with total points %d'%(i,refSize,len(pairArrayList)))
            i+=1
        self.pairArray = np.array(pairArrayList)

    def domainPairArray(self, minRef = 20, maxRef=500, probIndex=0.01):
        partialPairList = []
        #self.geom['refDist']=refDist
        self.geom['maxRef']=maxRef
        self.geom['minRef']=minRef
        #Add all shapes to the tuple list
        for key, shape in self.arrays.items():
            partialPairList += self.vertexAsList(shape)
        self.geom['partialPairList'] = partialPairList
        print('Original vextex  number %d'%len(partialPairList))

        #Definitio of distance tree
        self.tree = cKDTree(self.geom['partialPairList'])

        #Geneation of total pair array
        self.pairArrayGenerator()
        print('Total number of vertex %d'%self.pairArray.shape[0])

    def relaxVertices(self, relaxLevel=3):
        for i in range(relaxLevel):
            field = Field(self.pairArray)
            field.relax()
            self.pairArray = field.get_points()

    def createClipVoronoi(self):
        vorNew = Voronoi(self.pairArray)
        polygonList = []
        for ptGroup in vorNew.regions:
            ptGroupList = []
            for pt in ptGroup:
                if pt != -1:
                    ptGroupList.append(tuple(vorNew.vertices[pt]))
            if len(ptGroupList) >= 3:
                polygonList.append(Polygon(ptGroupList))
        vorList = gpd.GeoSeries(polygonList)
        self.vorClip = gpd.clip(vorList, self.geom['limitGeometry'])
