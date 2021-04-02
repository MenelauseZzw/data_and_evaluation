import h5py
import numpy as np
from xml.etree import ElementTree
import collections

def readRawFile(filename, shape, thresholdBelow=0.05):
    rawData      = np.fromfile(filename, dtype=np.float32)
    rawData      = np.reshape(rawData, (-1, 5)) # each row consists of rsp(1), dir(3) and rad(1)

    responses    = rawData[:, 0]

    ignor        = responses < thresholdBelow # ignor is 1-D array
    responses    = responses[~ignor]
    tangentLines = rawData[:, 1:4]
    tangentLines = tangentLines[~ignor]
    radiuses     = rawData[:, 4]
    radiuses     = radiuses[~ignor]

    ignor        = np.reshape(ignor, shape) # ignor is 3-D array
    zs,ys,xs     = np.where(~ignor)
    measurements = np.column_stack((xs,ys,zs))

    tangentLinesPoints1 = np.empty_like(measurements, dtype=np.double)
    tangentLinesPoints2 = np.empty_like(measurements, dtype=np.double)

    for i in xrange(len(measurements)):
        p  = measurements[i]
        lp = tangentLines[i]
        tangentLinesPoints1[i] = p + lp
        tangentLinesPoints2[i] = p - lp
    
    dataset = dict()

    dataset['measurements']        = measurements
    dataset['tangentLinesPoints1'] = tangentLinesPoints1
    dataset['tangentLinesPoints2'] = tangentLinesPoints2
    dataset['radiuses']            = radiuses
    dataset['responses']           = responses

    return dataset

def readGxlFile(filename):
    elemTree = ElementTree.parse(filename)

    indices = dict()

    positions = []
    nodeTypes = []

    for idx,node in enumerate(elemTree.getroot().findall('./graph/node')):
        nodeID          = node.attrib['id']
        nodeType        = node.findtext("./attr[@name=' nodeType']/string")
        position        = tuple(float(item.text) for item in node.findall("./attr[@name=' position']/tup/float"))
        indices[nodeID] = idx

        if nodeType == ' root node ':
            nodeTypes.append('r')
        elif nodeType == ' bifurication ':
            nodeTypes.append('b')
        else:
            assert nodeType == ' terminal node '
            nodeTypes.append('t')

        positions.extend(position)

    n = len(indices)
    positions = np.reshape(positions, newshape=(-1, 3))

    indices1      = []
    indices2      = []
    radiusesPrime = []
    direction     = []
    radiusesDict  = collections.defaultdict(list)

    for edge in elemTree.getroot().findall('./graph/edge'):
        edgeID      = edge.attrib['id']
        sourceID    = edge.attrib['from']
        targetID    = edge.attrib['to']
        radiusPrime = float(edge.findtext("./attr[@name=' radius']/float"))
        flow        = float(edge.findtext("./attr[@name=' flow']/float"))
        
        indices1.append(indices[sourceID])
        indices2.append(indices[targetID])
        radiusesPrime.append(radiusPrime)
        radiusesDict[indices[sourceID]].append(radiusPrime)
        radiusesDict[indices[targetID]].append(radiusPrime)

        dir = positions[indices[targetID]] - positions[indices[sourceID]]
        direction.extend(dir)
   
    temp = sorted(radiusesDict.items(),key=lambda item:item[0])
    maxRange = [2.0 * np.max(v[1]) for v in temp]#2.0 is just a predefined number
    dataset = dict()


    dataset['positions']   = np.reshape(positions, newshape=(-1, 3))
    dataset['nodeTypes']   = np.array(nodeTypes, dtype=object)
    dataset['indices1']    = np.array(indices1, dtype=np.int)
    dataset['indices2']    = np.array(indices2, dtype=np.int)
    dataset['radiusPrime'] = np.array(radiusesPrime, dtype=np.double)
    dataset['direction']   = np.reshape(direction,newshape=(-1,3))
    dataset['MaxRange']    = np.array(maxRange, dtype=np.double)
    
    return dataset

def readH5File(filename, n_dims=3):
    dataset = dict()

    with h5py.File(filename, mode='r') as f:
        for name in ['rootIndex','gap','costValue','indices1', 'indices2', 'radiuses', 'weights', 'responses', 'arcRadiuses', 'arcRadiusesMean', 'arcRadiusesStdDev', 'sourceIndices', 'targetIndices', 'connectedComponentsIndices','objectnessMeasure']:
            if name in f:
                dataset[name] = f[name][()]

        for name in ['positions', 'measurements', 'tangentLinesPoints1', 'tangentLinesPoints2']:
            if name in f:
                item = f[name][()]
                item = np.reshape(item, newshape=(-1,n_dims))
                dataset[name] = item

    return dataset

def writeH5File(filename, dataset):
    with h5py.File(filename, mode='w') as f:
        for name in dataset:
            item = dataset[name]
            f.create_dataset(name, data=item.flatten())
