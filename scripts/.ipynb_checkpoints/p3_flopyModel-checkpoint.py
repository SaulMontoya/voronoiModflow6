import flopy
import flopy.utils.binaryfile as bf
import json, rasterio
import geopandas as gpd
import numpy as np
from tqdm import tqdm
import time

#open input files
jsonFile = open('../txt/disvDict.json', 'r')
disvDict = json.load(jsonFile)

jsonFile = open('../txt/spatialIndexDict.json', 'r')
spatialIndexDict = json.load(jsonFile)

demRaster = rasterio.open('../rst/ASTGTM2_18S_xyz2.asc')
riversDf =  gpd.read_file('../shps/rios.shp')
meshDf = gpd.read_file('../shps/voronoiGrid.shp')

### MODFLOW Model ###
modelName = 'Model'
workSpace = '../model'
exeName = '../exe/mf6.exe'

sim = flopy.mf6.MFSimulation(sim_name=modelName, version='mf6',
                             exe_name=exeName,
                             sim_ws=workSpace)

tdis = flopy.mf6.ModflowTdis(sim, time_units='SECONDS',
                             perioddata=[[1.0, 1, 1.]])

gwf = flopy.mf6.ModflowGwf(sim, modelname=modelName, save_flows=True)

ims = flopy.mf6.ModflowIms(sim, print_option='SUMMARY', complexity='complex',
                           outer_hclose=1.e-2, inner_hclose=1.e-2)

#data from dictionary
cell2d = disvDict['cell2dArrays']
vertices = disvDict['indexedVerticesList']
ncpl = disvDict['NCPL']
nvert = disvDict['NVERT']
#another disv data
thickRatio = [0.05,0.1,0.2,0.4,0.7,1]
nlay = 6
bot = 3200

top = []

layBotmList = []
for cell in cell2d:
    cellBotmList = []
    cellCentroid = cell[1:3]
    valuesXY = demRaster.sample([cellCentroid])
    topCell = float(list(valuesXY)[0][0])
    top.append(topCell)
    for lay in range(nlay):
        botm = bot + (topCell - bot)*(1-thickRatio[lay])
        cellBotmList.append(botm)
    layBotmList.append(cellBotmList)

topArray = np.array(top)
botm = np.array(layBotmList).T

dis = flopy.mf6.ModflowGwfdisv(gwf, nlay=nlay, ncpl=ncpl, nvert=nvert,
                               top=topArray, botm=botm,
                               vertices=vertices, cell2d=cell2d)

# Flujo de propiedad
kh = [1E-4, 1E-4, 5E-5, 1E-6, 5E-7, 5E-7]
icelltype = [1, 1, 1, 1, 0, 0]
npf = flopy.mf6.ModflowGwfnpf(gwf,
                              save_specific_discharge=True,
                              icelltype=icelltype,
                              k=kh,
                              k33=kh)

ic = flopy.mf6.ModflowGwfic(gwf, strt=topArray.max())

rch = flopy.mf6.ModflowGwfrcha(gwf, recharge=0.15/86400/365)

evt = flopy.mf6.ModflowGwfevta(gwf, surface=top, rate=1.2/86400/365)


riverSpd = []

# Get grid index
intervalNumber = spatialIndexDict['intervalNumber']
gridXarray = spatialIndexDict['gridXarray']
gridYarray = spatialIndexDict['gridYarray']
gridIndexList = spatialIndexDict['gridIndexList']

#defining function
def findIndex(var, coordArray):
    for interval in range(intervalNumber):
#        if var > coordArray[interval] and var < coordArray[interval+1]:
        if var > coordArray[interval] and var < coordArray[interval+1]:
            return interval
            break

i = 0

print('\nWorking with the spatial index of rivers')
riverCellsList = []
riverSpd = []

riversDf[["xInterBeg","xInterEnd","yInterBeg","yInterEnd"]] = 0
for riverIndex, riverRow in tqdm(riversDf.iterrows(), total= riversDf.shape[0]):
    coords = riverRow.geometry.bounds
    xmin = coords[0]
    xmax = coords[2]
    ymin = coords[1]
    ymax = coords[3]
    xInterBeg = findIndex(xmin,gridXarray)
    xInterEnd = findIndex(xmax,gridXarray)
    yInterBeg = findIndex(ymin,gridYarray)
    yInterEnd = findIndex(ymax,gridYarray)
    riversDf.loc[riverIndex,"xInterBeg"] = xInterBeg
    riversDf.loc[riverIndex,"xInterEnd"] = xInterEnd
    riversDf.loc[riverIndex,"yInterBeg"] = yInterBeg
    riversDf.loc[riverIndex,"yInterEnd"] = yInterEnd
#print(riversDf.head())
riversDf.to_csv('../txt/riverDf')

print('\nSpatial location of rivers cells')
riverCellsList = []
riverSpd = []
timeList = []
start = time.time()
for xInt in tqdm(range(intervalNumber)):
    for yInt in range(intervalNumber):
        #print('Working for cell %d, %d'%(xInt, yInt))
        #filter river cells in the gridrectangle
        #f1 = riversDf[(xInt >= riversDf["xInterBeg"])]
        #print(f1['CODIGO'])

        filterRiverDf = riversDf[(xInt >= riversDf["xInterBeg"]) & ((xInt) <= riversDf["xInterEnd"])
                            & ((yInt) >= riversDf["yInterBeg"]) & ((yInt) <= riversDf["yInterEnd"])]# & (riversDf["yInterEnd"]<=yInt+1)]


        cellIndexList=[]
        #if filterRiverDf.count()[0]>0:
        if filterRiverDf.shape[0]>0:
            #print('Number of rivers on grid: %d'%filterRiverDf.count()[0])
            for cellIndex, cellRow in enumerate(gridIndexList):
                #if cellIndex not in riverCellsList:
                if xInt in cellRow[0] and yInt in cellRow[1]:
                    cellIndexList.append(cellIndex)
        
        if xInt == 1 and yInt == 5:
            print(filterRiverDf)
            print('La longitud de celdas es de %d'%len(cellIndexList))
            print('Estas son las celdas:'+str(cellIndexList))
            print('El count es %d'%filterRiverDf.shape[0])
        #    print(xInt, yInt)
        #iterating though river reaches
        #start = time.time()
        if filterRiverDf.shape[0]>0:
            for riverIndex, riverRow in filterRiverDf.iterrows():
                #iterating through cells
                #print(riverRow)
                for cellIndex in cellIndexList:
                        cellPolyRow = meshDf.iloc[cellIndex]
                        #print(riverRow.geometry)
                        if riverRow.geometry.crosses(cellPolyRow.geometry):
                            if xInt == 1 and yInt == 5:
                                print('Esta celdad esta añadida: %d'%cellIndex)
                            #print(cellPolyRow)
                            riverCellsList.append(cellIndex)


riverCellsList = set(riverCellsList)
riverCellsList = list(riverCellsList)
if 48987 in riverCellsList:
    print('++++++++++++++++++++++++++++++++++++++++++++++++++++ok++++++++++++++++++++++++')
riverCellsList.sort()
for riverCell in riverCellsList:
    riverSpd.append([(0,riverCell),top[riverCell],3e-4])
end = time.time()
print(end - start)
        #timeList.append(end-start)

#print(max(timeList))
#print(riverSpd)
        #filterDf = filterDf[(filterDf["xInterEnd"]<=(xInt+1)) & (filterDf["yInterEnd"]<=(yInt+1))]
        #print(filterDf.describe())
        #if filterDf.count()[0]>0:
        #    print(filterDf.describe())
        #    print(xInt,yInt)

        #get river in this interval



#for polyIndex,polyRow in tqdm(meshDf.iterrows(), total= meshDf.shape[0]):
#    for riverIndex, riverRow in riversDf.iterrows():
#        if riverRow.geometry.crosses(polyRow.geometry):
#            riverSpd.append([(0,polyIndex),top[polyIndex],1e5,top[polyIndex]-2])
riv_spd = {0: riverSpd}
#riv = flopy.mf6.ModflowGwfriv(gwf, stress_period_data=riv_spd)
riv = flopy.mf6.ModflowGwfdrn(gwf, stress_period_data=riv_spd)


# wel
#wel_spd = {0: [[(2, 183), -45000],
#                [(2, 18), -45000]
#                ]}
#wel = flopy.mf6.ModflowGwfwel(gwf,stress_period_data=wel_spd)

# oc
hname = '{}.hds'.format(name)
cname = '{}.cbc'.format(name)
oc = flopy.mf6.ModflowGwfoc(gwf, budget_filerecord=cname,
                            head_filerecord=hname,
                            saverecord=[('HEAD', 'ALL'), ('BUDGET', 'ALL')])

# Escribir archivos y correr simulación
print('\nWriting model final and model simulation')
sim.write_simulation()
sim.run_simulation()
