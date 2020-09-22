from pymnet import *
from itertools import islice
import sys
import matplotlib

LigRecFile = sys.argv[1]
RecTFFile = sys.argv[2]
TFTarGeneFile = sys.argv[3]
outPic = sys.argv[4]

net = MultilayerNetwork(aspects = 1,fullyInterconnected = False)

def AddELe(node,alist):
    if node not in alist:
        alist.append(node)
    return alist

def StaNode(node,nodeCountDict):
    if node in nodeCountDict:
        nodeCountDict[node] += 1
    else:
        nodeCountDict[node] = 1
    return nodeCountDict

#main
LigList = []
RecList = []
TFList = []
TarList = []

edgeColorDict0 = {}
nodeColorDict0 = {}
nodeLabelColorDict0 = {}
nodeCountDict = {}
LayNodeListDict = {}
edgeWidthDict0 = {}

#read Net Files
LigRec = open(LigRecFile,'r')
for line in islice(LigRec,0,None):
   line = line.strip("\n")
   Lig = line.split('_')[0]
   Rec = line.split('_')[1]
   net[Lig,"Ligand"][Rec,"Receptor"] = 1
   LigList = AddELe(Lig,LigList)
   RecList = AddELe(Rec,RecList)
   edgeColorDict0[((Lig,"Ligand"),(Rec,"Receptor"))] = '#0000FF'
   nodeCountDict = StaNode(Lig,nodeCountDict)
   nodeCountDict = StaNode(Rec,nodeCountDict)
   edgeWidthDict0[((Lig,"Ligand"),(Rec,"Receptor"))] = 1.5

RecTF = open(RecTFFile,'r')
for line in islice(RecTF,0,None):
   line = line.strip("\n")
   Rec = line.split('_')[0]
   TF = line.split('_')[1]
   net[Rec,"Receptor"][TF,"TF"] = 1
   RecList = AddELe(Rec,RecList)
   TFList = AddELe(TF,TFList)
   edgeColorDict0[((Rec,"Receptor"),(TF,"TF"))] = '#009933'
   nodeCountDict = StaNode(Rec,nodeCountDict)
   nodeCountDict = StaNode(TF,nodeCountDict)
   edgeWidthDict0[((Rec,"Receptor"),(TF,"TF"))] = 0.5

TFTar = open(TFTarGeneFile,'r')
for line in islice(TFTar,0,None):
    line = line.strip("\n")
    TF = line.split('_')[0]
    Tar = line.split('_')[1]
    net[TF,"TF"][Tar,"Target gene"] = 1
    TFList = AddELe(TF,TFList)
    TarList = AddELe(Tar,TarList)
    edgeColorDict0[((TF,"TF"),(Tar,"Target gene"))] = '#CC66FF'
    nodeCountDict = StaNode(TF,nodeCountDict)
    nodeCountDict = StaNode(Tar,nodeCountDict)
    edgeWidthDict0[((TF,"TF"),(Tar,"Target gene"))] = 0.5

#read end.

#get node color
def getColor(nodeList,rawDict,lay,color):
    for i in nodeList:
        rawDict[(i,lay)] = color
    return rawDict

nodeColorDict0 = getColor(LigList,nodeColorDict0,"Ligand","#708090")
nodeColorDict0 = getColor(RecList,nodeColorDict0,"Receptor","#F4A460")
nodeColorDict0 = getColor(TFList,nodeColorDict0,"TF","#90EE90")
nodeColorDict0 = getColor(TarList,nodeColorDict0,"Target gene","#6495ED")

nodeLabelColorDict0 = getColor(LigList,nodeLabelColorDict0 ,"Ligand","#000000")
nodeLabelColorDict0 = getColor(RecList,nodeLabelColorDict0 ,"Receptor","#FF6600")
nodeLabelColorDict0 = getColor(TFList,nodeLabelColorDict0 ,"TF","#006600")
nodeLabelColorDict0 = getColor(TarList,nodeLabelColorDict0 ,"Target gene","#0000FF")

#get node size and node label size
LayNodeListDict = {"Ligand":LigList,"Receptor":RecList,"TF":TFList,"Target gene":TarList}
nodeSizeDict = {}
nodeLabelSizeDict = {}
for l in LayNodeListDict:
    for n in LayNodeListDict[l]:
        if l == "Ligand":
            nodeSizeDict[(n,l)] = 0.03
            nodeLabelSizeDict[(n,l)] = 5
        elif l == "Receptor":
            nodeSizeDict[(n,l)] = 0.03
            nodeLabelSizeDict[(n,l)] = 5
        elif l == "TF":
            nodeSizeDict[(n,l)] = 0.01
            nodeLabelSizeDict[(n,l)] = 3
        elif l == "Target gene":
            nodeSizeDict[(n,l)] = 0.01
            nodeLabelSizeDict[(n,l)] = 2

#get node coords
def Coords(gmp_unique,lay,nodeCoorDict):
	Count = len(gmp_unique)
	dis = 1.6/4 #1.6/int(math.sqrt(Count))
	dis2 = 1.6/8
	for i in range(Count):
		zheng = i//9 #(int(math.sqrt(Count))+1)
		yu = i%9 #(int(math.sqrt(Count))+1)
		nodeCoorDict[(gmp_unique[i],lay)] = (-0.8+yu*dis, 0.8-zheng*dis)
	return nodeCoorDict

def Coords2(gmp_unique,lay,nodeCoorDict):
        Count = len(gmp_unique)
        dis = 0.8
        dis2 = 1.2
        for i in range(Count):
                zheng = i// 4#(int(math.sqrt(Count))+1)
                yu = i % 4 #(int(math.sqrt(Count))+1)
                nodeCoorDict[(gmp_unique[i],lay)] = (-1 + yu * dis,-1 - zheng* dis2)
        return nodeCoorDict

def Coords3(gmp_unique,lay,nodeCoorDict):
        Count = len(gmp_unique)
        for i in range(Count):
                zheng = i// 7 #(int(math.sqrt(Count))+1)
                yu = i % 7  #(int(math.sqrt(Count))+1)
                nodeCoorDict[(gmp_unique[i],lay)] = (-1+yu*0.4, -1-zheng*0.5)
        return nodeCoorDict


nodeCoorDict = {}
nodeCoorDict = Coords2(LigList,"Ligand",nodeCoorDict)
nodeCoorDict = Coords2(RecList,"Receptor",nodeCoorDict)
nodeCoorDict = Coords3(TFList,"TF",nodeCoorDict)
nodeCoorDict = Coords(TarList,"Target gene",nodeCoorDict)


#draw net

fig = draw(net, show=False, azim=-30, elev=25, layergap=0.80,
						layout='shell',
						autoscale=True,
						defaultNodeLabelAlpha=0.95,
						defaultLayerAlpha=0.5,
                        defaultEdgeAlpha = 0.5,
                                                defaultEdgeWidth= 0.5,

                                                figsize = [10,8],

						# color
						nodeColorDict = nodeColorDict0,
						nodeLabelColorDict = nodeLabelColorDict0,
                                                edgeColorDict = edgeColorDict0,

						# size
						nodeSizeDict = nodeSizeDict,
						nodeSizeRule={'scalecoeff': 0.2, 'rule': 'scaled'},
						nodeLabelSizeDict = nodeLabelSizeDict,


                                                #width
                                                edgeWidthDict = edgeWidthDict0,


						# coords
						#nodeCoords = nodeCoorDict,
						nodelayerCoords = nodeCoorDict,

						layerColorDict={'Ligand':'#FF99CC', 'Receptor':'#33FFFF','TF':'#FFFF66','Target gene':'#66CCFF'},
						layerOrderDict={'Ligand':1, 'Receptor':2, 'TF':3, 'Target gene':4},
						# node label style
						defaultNodeLabelStyle='oblique', # normal, italic or oblique
                                                                                                
                                                edgeStyleRule={'inter': '-', 'intra': '-', 'rule': 'edgetype'}, 
                                                defaultEdgeStyle='-',

						layerPadding=0.25)

fig.savefig(outPic,dpi = 200)

print ("output net:",outPic)
