import os, json

jsonFile = open('../txt/disvDict.json')
disvDict = json.load(jsonFile)

topFile = open('../model/disv_mtop.disv','w')
botmFile = open('../model/disv_botm.disv','w')

botDict = {}
for lay in range(disvDict['NLAY']):
    botDict['botmLay'+str(lay)]=[]

for elev in disvDict['polyCentroidElevation']:
    topFile.write(str(elev)+' ')
    totThick = elev - disvDict['modelBottom']
    for lay in range(disvDict['NLAY']):
        botLay = elev - disvDict['thickRatio'][lay]*totThick
        botDict['botmLay'+str(lay)].append(botLay)

for lay in range(disvDict['NLAY']):
    for item in botDict['botmLay'+str(lay)]:
        botmFile.write(str(item)+' ')
    botmFile.write('\n')

topFile.close()
botmFile.close()


verticeFile = open('../model/disv_vertice.disv','w')
for index, vertice in enumerate(disvDict['uniqueVerticesList']):
    verticeFile.write(str(index+1)+' '+str(vertice[0])+' '+str(vertice[1]) + '\n')
verticeFile.close()

cell2dFile = open('../model/disv_cell2d.disv','w')
for index, vertice in enumerate(disvDict['polyCentroidList']):
    seqList = disvDict['polygonVerticesSequence'][index]
    seqString = ''
    for item in seqList[:-1]:
        seqString += str(item+1)+' '
    cell2dFile.write(str(index+1)+' '+str(vertice[0])+' '+str(vertice[1]) +' '+str(len(seqList)-1)+' '+seqString+'\n')
cell2dFile.close()
