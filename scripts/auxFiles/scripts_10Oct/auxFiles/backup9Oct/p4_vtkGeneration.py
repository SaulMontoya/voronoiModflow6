import pyvista as pv
import numpy as np
import os, flopy
import matplotlib.pyplot as plt
name = 'Model'
workspace = '../model'
mf_exe_name = '../exe/mf6'

sim = flopy.mf6.MFSimulation.load(sim_name=name, exe_name=mf_exe_name, sim_ws=workspace)

mfmodel = sim.get_model(model_name='model')

fname = os.path.join(workspace, name + '.hds')
hdobj = flopy.utils.HeadFile(fname, precision='double')
head = hdobj.get_data()

mtop = np.ones([mfmodel.disv.nvert.array ])*mfmodel.disv.top[0]
verticesXYArray = np.dstack((mfmodel.disv.vertices.array.xv,mfmodel.disv.vertices.array.yv))[0]

modelz = np.vstack((mfmodel.disv.top.array,mfmodel.disv.botm.array))
nlay = mfmodel.disv.nlay.array

cellTubes = {}
for index, cell in enumerate(mfmodel.disv.cell2d.array):
    vertexIndexList = [x for x in list(cell)[4:] if x is not None]
    filterVertexIndexList = [len(vertexIndexList)] + list(range(len(vertexIndexList)))
    filterVertexIndexArray = np.array(filterVertexIndexList)
    layerTubes = {}
    for lay in range(nlay):
        filterVertexArray = np.array([list(verticesXYArray[vertex]) + [modelz[lay,index]] for vertex in vertexIndexList])
        cellSurf = pv.PolyData(filterVertexArray, filterVertexIndexArray)
        cellSurf.cell_arrays['head'] = head[lay][0][index]
        cellZ = modelz[:,index]
        cellMesh = cellSurf.extrude([0, 0, cellZ[lay+1]-cellZ[lay]])
        layerTubes[str(lay)] = cellMesh
    blocks = pv.MultiBlock(layerTubes)
    cellTubes[str(index)] = blocks.combine()
totalModelCells = pv.MultiBlock(cellTubes)
totalModelMerged = totalModelCells.combine()
totalModelMerged.save('../vtk/modelMerge.vtk')
