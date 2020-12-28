from vorModflow.geoVoronoi import createVoronoi
import matplotlib.pyplot as plt
import geopandas as gpd
import numpy as np
import time

#Create mesh object
vorMesh = createVoronoi()

#Open limit layers and refinement definition layers
vorMesh.addLimitLayer('basin','../shps/Angascancha_Basin_Extension.shp')
vorMesh.addDiscretizationLayer('river','../shps/rios.shp')

#Generate point pair array
start = time.time()
vorMesh.domainPairArray()
end = time.time()
print('Total time required for: Generation of total points %.2f seconds \n'%(end - start))

#Relaxation of voronoi polygons
start = time.time()
vorMesh.relaxVertices()
end = time.time()
print('Total time required for: Relaxation of total points %.2f seconds \n'%(end - start))

#Create polygon object and clip to the limit layer
start = time.time()
vorMesh.createClipVoronoi()
end = time.time()
print('Total time required for: Voronoi gridding %.2f seconds \n'%(end - start))

#Save clipped polygons to shapefile o geojson
vorGdf = gpd.GeoDataFrame(vorMesh.vorClip)
vorGdf.to_file('../shps/voronoiGrid.shp')

#Generate graph with geo features, refinement vertex, relaxed points and clipped voronois
partialPairArray=np.array(vorMesh.geom['partialPairList'])
fig, ax =  plt.subplots(1, 1, figsize=(20,20), sharex=True, sharey=True)
vorMesh.vorClip.plot(ax=ax, facecolor="none", edgecolor='black')
ax.scatter(vorMesh.pairArray[:,0],vorMesh.pairArray[:,1],s=10,alpha=0.5,c='orangered')
ax.scatter(partialPairArray[:,0],partialPairArray[:,1],marker='*',s=10,c='royalblue')
for key, array in vorMesh.arrays.items():
    array.plot(ax=ax, facecolor="none",edgecolor='dodgerblue', alpha=0.8)
plt.show()
