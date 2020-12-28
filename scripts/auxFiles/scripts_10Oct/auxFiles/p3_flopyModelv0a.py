import flopy
import flopy.utils.binaryfile as bf
import json, rasterio
import geopandas as gpd
import numpy as np
from tqdm import tqdm

#open input files
jsonFile = open('../txt/disvDict.json', 'r')
disvDict = json.load(jsonFile)

demRaster = rasterio.open('../rst/ASTGTM2_18S_xyz2.asc')
riversDf =  gpd.read_file('../shps/rios.shp')
meshDf = gpd.read_file('../shps/voronoiGrid.shp')

### MODFLOW Model ###
name = 'Model'
workspace = '../model'
mf_exe_name = '../exe/mf6'

sim = flopy.mf6.MFSimulation(sim_name=name, version='mf6',
                             exe_name=mf_exe_name,
                             sim_ws=workspace)

tdis = flopy.mf6.ModflowTdis(sim, time_units='DAYS',
                             perioddata=[[1.0, 1, 1.]])

gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)

ims = flopy.mf6.ModflowIms(sim, print_option='SUMMARY', complexity='simple',
                           outer_hclose=1.e-2, inner_hclose=1.e-2)

#print(help(flopy.mf6.ModflowGwfdisv))


cell2d = disvDict['cell2dArrays']
vertices = disvDict['indexedVerticesList']
#nlay = disvDict['NLAY']
ncpl = disvDict['NCPL']
nvert = disvDict['NVERT']
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
kh = [50, 50, 50, 50, 50, 50]
icelltype = [1, 1, 0, 0, 0, 0]
npf = flopy.mf6.ModflowGwfnpf(gwf,
                              save_specific_discharge=True,
                              icelltype=icelltype,
                              k=kh,
                              k33=kh)

ic = flopy.mf6.ModflowGwfic(gwf, strt=topArray.max())

rch = flopy.mf6.ModflowGwfrcha(gwf, recharge=0.001)

evt = flopy.mf6.ModflowGwfevta(gwf, surface=top, rate=0.003)


#ibd = tri.get_boundary_marker_array()
#ibd = np.ma.masked_not_equal(ibd, 2)
#riv_node = []
#for i in range(ibd.shape[0]):
#    try:
#        int(ibd[i])
#        riv_node.append(i)
#    except:
#        pass
#riv_spd = {0: [[(0, i), 320, 1e5, 318] for i in riv_node]}
#riv = flopy.mf6.ModflowGwfriv(gwf, stress_period_data=riv_spd)
riverSpd = []
for polyIndex,polyRow in tqdm(meshDf.iterrows(), total= meshDf.shape[0]):
    for riverIndex, riverRow in riversDf.iterrows():
        if riverRow.geometry.crosses(polyRow.geometry):
            riverSpd.append([(0,polyIndex),top[polyIndex],1e5,top[polyIndex]-2])
riv_spd = {0: riverSpd}
riv = flopy.mf6.ModflowGwfriv(gwf, stress_period_data=riv_spd)


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
sim.write_simulation()
sim.run_simulation()
