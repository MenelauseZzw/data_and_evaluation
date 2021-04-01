import IO2
import argparse
import pandas as pd
import os.path
import numpy as np
import numpy.linalg as linalg
from sklearn.neighbors import NearestNeighbors


def resamplePoints(points, indices1, indices2, samplingStep, flag=True):
    resampledPoints = []
    resampledIndices = []

    for k, (index1, index2) in enumerate(zip(indices1, indices2)):
        point1 = points[index1]
        point2 = points[index2]
        num = int(np.ceil(linalg.norm(point1 - point2) / samplingStep) + 1)
        if not flag:
            num = 1
        for i in np.arange(num + 1):
            lambd = i / float(num)
            interp = (1 - lambd) * point1 + lambd * point2
            resampledPoints.append(interp)
            resampledIndices.append(k)

    return np.array(resampledPoints), np.array(resampledIndices)


def resamplePointsOnGT(points, nodetypes, maxrange, indices1, indices2, samplingStep):
    resampledPoints = []
    resampledIndices = []
    isBifurcationPoints = []

    for k, (index1, index2) in enumerate(zip(indices1, indices2)):
        point1 = points[index1]
        point2 = points[index2]
        num = int(np.ceil(linalg.norm(point1 - point2) / samplingStep) + 1)
        for i in np.arange(num + 1):
            lambd = i / float(num)
            interp = (1 - lambd) * point1 + lambd * point2
            resampledPoints.append(interp)
            resampledIndices.append(k)
            if nodetypes[index1] == 'b' and linalg.norm(interp - point1) <= maxrange[index1]:
                isBifurcationPoints.append(True)
            elif nodetypes[index2] == 'b' and linalg.norm(interp - point2) <= maxrange[index2]:
                isBifurcationPoints.append(True)
            else:
                isBifurcationPoints.append(False)

    return np.array(resampledPoints), np.array(resampledIndices), np.array(isBifurcationPoints)


def getBifurcationPoints(points, indices1, indices2, nodeTypes=None, maxRange=None):
    if nodeTypes is not None:
        adjacentIndices = [[] for i in np.arange(points.shape[0])]
        anglesGT = []
        anglesGT2 = []
        anglesGT3 = []
        for i, (index1, index2) in enumerate(zip(indices1, indices2)):
            adjacentIndices[index1].append(index2)
            adjacentIndices[index2].append(-1 * index1)
        for i in np.arange(points.shape[0]):
            if nodeTypes[i] == 'b':
                p0 = points[i]
                p1p2 = [points[adjacentIndices[i][j]] for j in np.arange(3) if adjacentIndices[i][j] > 0]
                vector1 = p1p2[0] - p0
                vector2 = p1p2[1] - p0
                vector1 = vector1 / linalg.norm(vector1)
                vector2 = vector2 / linalg.norm(vector2)
                angle = np.arccos(np.dot(vector1, vector2))
                anglesGT.append(angle)

                p3 = [points[-adjacentIndices[i][j]] for j in np.arange(3) if adjacentIndices[i][j] <= 0]
                vector3 = p0 - p3[0]
                vector3 = vector3 / linalg.norm(vector3)
                tmpres = np.clip(np.dot(vector1, vector3), -1, 1)
                angle = np.arccos(tmpres)
                anglesGT2.append(angle)

                tmpres = np.clip(np.dot(vector2, vector3), -1, 1)
                angle = np.arccos(tmpres)
                anglesGT3.append(angle)

        return points[np.where(nodeTypes == 'b')], maxRange[np.where(nodeTypes == 'b')], np.array(anglesGT), np.array(
            anglesGT2), np.array(anglesGT3)
    else:
        numOfPoints = points.shape[0]
        pointsDegree = np.zeros(numOfPoints, dtype=int)
        adjacentPointsAtBifurcation = [[] for i in np.arange(numOfPoints)]

        for i, (index1, index2) in enumerate(zip(indices1, indices2)):
            pointsDegree[index1] += 1
            pointsDegree[index2] += 1
            adjacentPointsAtBifurcation[index1].append(index2)
            adjacentPointsAtBifurcation[index2].append(index1)

        # branching points only include bifurcation points and trifurcation points
        x = pointsDegree >= 3
        y = pointsDegree <= 4
        z = np.vstack((x, y))
        ind = np.all(z, axis=0)

        IndList = np.arange(numOfPoints)
        revIndices = IndList[ind]

        return points[ind], ind, adjacentPointsAtBifurcation, pointsDegree, revIndices


def computeMetricsForROC(GTpoints, RTpoints, upperBound, lowerBound, radii, maxrad, indicesOrig = None):
    learner = NearestNeighbors(n_neighbors=1)
    learner.fit(RTpoints)

    distOrig, closestIndicesOrig = learner.kneighbors(GTpoints)
    distOrig = np.ravel(distOrig)
    closestIndicesOrig = np.ravel(closestIndicesOrig)

    if indicesOrig is None:
        radiiOrig = radii
    else:
        radiiOrig = radii[indicesOrig]

    arrUpOrig = radiiOrig < upperBound
    arrLoOrig = radiiOrig > lowerBound
    numLo_UpOrig = np.count_nonzero(np.logical_and(arrUpOrig, arrLoOrig))
    if indicesOrig is None:
        closerThanRadiusOrig = np.logical_and(np.logical_and(arrUpOrig, arrLoOrig), distOrig < maxrad)
    else:
        closerThanRadiusOrig = np.logical_and(np.logical_and(arrUpOrig, arrLoOrig), distOrig < maxrad[indicesOrig])
    RTclosestIndicesFlag = np.logical_and(arrUpOrig, arrLoOrig)

    numCloserThanRadiusOrig = np.count_nonzero(closerThanRadiusOrig)  # number of points {p_i} in GroundTruthTree such that {q_j} is closest in ReconstructedTree and ||p_i - q_j|| < maxRad
    numOrig = numLo_UpOrig

    learner = NearestNeighbors(n_neighbors=1)
    learner.fit(GTpoints)

    dist, closestIndices = learner.kneighbors(RTpoints)
    dist = np.ravel(dist)
    closestIndices = np.ravel(closestIndices)

    if indicesOrig is None:
        radii = radii[closestIndices]
    else:
        radii = radii[indicesOrig[closestIndices]]

    arrUp = radii < upperBound
    arrLo = radii > lowerBound
    numLo_Up = np.count_nonzero(np.logical_and(arrUp, arrLo))
    if indicesOrig is None:
        closerThanRadius = np.logical_and(np.logical_and(arrUp, arrLo), dist < maxrad[closestIndices])
    else:
        closerThanRadius = np.logical_and(np.logical_and(arrUp, arrLo), dist < maxrad[indicesOrig[closestIndices]])

    numCloserThanRadius = np.count_nonzero(closerThanRadius)  # number of points {q_j} in ReconstructedTree such that {p_i} is closest in GroundTruthTree and ||p_i - q_j|| < maxRad
    num = numLo_Up

    if indicesOrig is None:
        return numCloserThanRadiusOrig, numOrig, numCloserThanRadius, num, closestIndicesOrig, RTclosestIndicesFlag
    else:
        return numCloserThanRadiusOrig, numOrig, numCloserThanRadius, num


def doComputeOverlapMeasure(args):
    dirname_originalVolume = args.dirname_originalVolume
    dirname_noisyVolume = args.dirname_noisyVolume
    dirname_volume = args.dirname_volume
    dirname_reconstructed = args.dirname_reconstructed

    voxelSize = args.voxelWidth
    samplingStep = args.samplingStep
    pointsArrName = args.points
    flagSquare = args.flagSquare
    upperBound = args.upperBound
    lowerBound = args.lowerBound
    Kth = args.K

    doOutputHeader = args.doOutputHeader
    prependHeaderStr = args.prependHeaderStr
    prependRowStr = args.prependRowStr

    debugMode = args.debugMode

    filename = os.path.join(dirname_originalVolume,dirname_volume, 'tree_structure.xml')
    dataset = IO2.readGxlFile(filename)

    positions = dataset['positions']
    positions = voxelSize * positions

    indicesOrig1 = dataset['indices1']
    indicesOrig2 = dataset['indices2']
    radiuses = dataset['radiusPrime']
    directions = dataset['direction']
    maxRange = dataset['MaxRange']
    nodeTypes = dataset['nodeTypes']

    # root of the groundtruth tree
    rind = np.nonzero(nodeTypes == 'r')[0][0]  # root index
    rpos = positions[rind]  # root position

    filename = os.path.join(dirname_noisyVolume,dirname_volume,dirname_reconstructed)
    dataset = IO2.readH5File(filename)

    Points = dataset[pointsArrName]

    indices1 = dataset['indices1']
    indices2 = dataset['indices2']
    tangentlinepoint1 = dataset['tangentLinesPoints1']
    tangentlinepoint2 = dataset['tangentLinesPoints2']

    # find the reconstructed root point
    learner = NearestNeighbors(n_neighbors=1)
    learner.fit(Points)

    dst, clinx = learner.kneighbors(np.array([rpos]))
    resInxRoot = clinx[0][0]  # reconstructed root index

    ############# all points matching ################
    # resample groundtruth and reconstructed tree
    pointsOrig, indicesOrig, bifurOrig = resamplePointsOnGT(positions, nodeTypes, maxRange, indicesOrig1, indicesOrig2,samplingStep)
    points, indices = resamplePoints(Points, indices1, indices2, samplingStep)
    maxrad = np.maximum(radiuses, np.sqrt(2) * voxelSize / 2)

    NumberOfPointsCloserThanRadiusOrig, NumberOfPointsOrig, NumberOfPointsCloserThanRadius, NumberOfPoints = computeMetricsForROC(pointsOrig, points, upperBound, lowerBound, radiuses, maxrad, indicesOrig)

    PercentageOfPointsCloserThanRadiusOrig = NumberOfPointsCloserThanRadiusOrig / float(NumberOfPointsOrig)
    PercentageOfPointsCloserThanRadius = NumberOfPointsCloserThanRadius / float(NumberOfPoints)
    OverlapMeasure = (NumberOfPointsCloserThanRadiusOrig + NumberOfPointsCloserThanRadius) / (float(NumberOfPointsOrig) + float(NumberOfPoints))

    ############# branching/bifurcation points matching ################
    bifurcationPointsOrig, maxRadiusAtBifurcation, anglesGT, anglesGT2, anglesGT3 = getBifurcationPoints(positions,indicesOrig1,indicesOrig2,nodeTypes,maxRange)  # groundtruth bifurcation points
    bifurcationPoints, isBifurcationPoints, adjacentPointsAtBifurcation, pointsDegree, revIndices = getBifurcationPoints(Points, indices1, indices2)  # reconstructed braching points; we still use bifurcation points as the name

    maxradBifur = np.maximum(maxRadiusAtBifurcation, np.sqrt(3) * voxelSize)

    NumberOfBifurPointsCloserThanRadiusOrig, NumberOfBifurPointsOrig, NumberOfBifurPointsCloserThanRadius, NumberOfBifurPoints, RTclosestIndices, RTclosestIndicesFlag  = computeMetricsForROC(bifurcationPointsOrig, bifurcationPoints, upperBound, lowerBound, maxRadiusAtBifurcation, maxradBifur)

    PercentageOfBifurPointsCloserThanRadiusOrig = float(NumberOfBifurPointsCloserThanRadiusOrig) / NumberOfBifurPointsOrig
    PercentageOfBifurPointsCloserThanRadius = float(NumberOfBifurPointsCloserThanRadius) / NumberOfBifurPoints

    ############# angle error at branching/bifurcation points ################
    # BFS to determine the reconstructed branching point angle
    newAdjacentList = [[] for i in np.arange(Points.shape[0])]
    visited = [False] * Points.shape[0]
    startind = resInxRoot
    bfsQueue = []

    visited[startind] = True
    bfsQueue.extend(adjacentPointsAtBifurcation[startind])
    if isBifurcationPoints[startind]:
        for i in np.arange(pointsDegree[startind]):
            if i == (pointsDegree[startind] - 1):
                newAdjacentList[startind].append(-adjacentPointsAtBifurcation[startind][i])
            else:
                newAdjacentList[startind].append(adjacentPointsAtBifurcation[startind][i])

    while len(bfsQueue) != 0:
        currentNode = bfsQueue[0]
        if visited[currentNode] == True:
            bfsQueue = bfsQueue[1:]
            continue
        else:
            visited[currentNode] = True
            adlist = adjacentPointsAtBifurcation[currentNode]
            bfsQueue.extend(adlist)
            bfsQueue = bfsQueue[1:]

            if pointsDegree[currentNode] >= 3:
                for i in np.arange(pointsDegree[currentNode]):
                    if visited[adlist[i]] == False:
                        newAdjacentList[currentNode].append(adlist[i])
                    else:
                        newAdjacentList[currentNode].append(-adlist[i])
            else:
                newAdjacentList[currentNode].extend(adlist)

    # compute angle error
    angleErrors = []
    # print("mean",np.mean(anglesGT)/np.pi*180, "median",np.median(anglesGT)/np.pi*180,"std",np.std(anglesGT/np.pi*180),"bound",np.mean(np.abs(np.array(anglesGT/np.pi*180) - np.median(anglesGT/np.pi*180))))
    for n,bp in enumerate(RTclosestIndices):
        if not RTclosestIndicesFlag[n]:
            continue
        i = revIndices[bp]
        adjacentIndices = newAdjacentList[i]
        vector12 = []
        vector33 = []
        for j in adjacentIndices:
            if j >= 0:
                tmpvec = Points[j] - Points[i]
                vector12.append(tmpvec)
            else:
                tmpvec = Points[i] - Points[-j]
                vector33.append(tmpvec)
        best_error = 10000000
        for i1 in range(len(vector12)):
            for i2 in range(i1 + 1, len(vector12)):
                vector1 = vector12[i1]
                vector2 = vector12[i2]
                vector1 = vector1 / linalg.norm(vector1)
                vector2 = vector2 / linalg.norm(vector2)
                angleRT = np.clip(np.dot(vector1, vector2), -1, 1)
                angleRT = np.arccos(angleRT)
                cur_error = np.abs(angleRT - anglesGT[n]) / np.pi * 180
                best_error = min(cur_error, best_error)
        assert best_error < 10000000, "%f" % cur_error
        if flagSquare:
            angleErrors.append(best_error ** 2)
        else:
            angleErrors.append(best_error)

    angleErrors = np.array(angleErrors)

    # Hausdorff distance, not average actually
    averageAngleAtBifur = angleErrors[np.argsort(angleErrors)[int(Kth / 100.0 * len(angleErrors)) - 1]]
    # averageAngleAtBifurStd = np.std(angleErrors)

    keyValPairs = [(name, eval(name)) for name in (
    'NumberOfPointsCloserThanRadiusOrig', 'NumberOfPointsOrig', 'NumberOfPointsCloserThanRadius', 'NumberOfPoints',
    'PercentageOfPointsCloserThanRadiusOrig', 'PercentageOfPointsCloserThanRadius', 'OverlapMeasure',
    'NumberOfBifurPointsCloserThanRadiusOrig', 'NumberOfBifurPointsOrig', 'NumberOfBifurPointsCloserThanRadius', 'NumberOfBifurPoints',
    'PercentageOfBifurPointsCloserThanRadiusOrig', 'PercentageOfBifurPointsCloserThanRadius',
    'averageAngleAtBifur', 'Kth')]

    if (doOutputHeader):
        print(prependHeaderStr + (",".join(kvp[0] for kvp in keyValPairs)))

    print(prependRowStr + (",".join(str(kvp[1]) for kvp in keyValPairs)))

    if debugMode:
        pass


def analyzeOverlapMeasure(df, voxelSize):
    outResult = pd.DataFrame(columns=['ThresholdValue', 'PercentageOfPointsCloserThanRadiusOrig',
                                      'PercentageOfPointsCloserThanRadiusOrigStd',
                                      'PercentageOfPointsCloserThanRadius', 'PercentageOfPointsCloserThanRadiusStd',
                                      'OverlapMeasure', 'OverlapMeasureStd',
                                      'PercentageOfBifurPointsCloserThanRadiusOrig','PercentageOfBifurPointsCloserThanRadiusOrigStd',
                                      'PercentageOfBifurPointsCloserThanRadius','PercentageOfBifurPointsCloserThanRadiusStd',
                                      'averageAngleAtBifur','averageAngleAtBifurStd'])

    for i, (ThresholdValue, g) in enumerate(df.groupby(['ThresholdValue'])):
        PercentageOfPointsCloserThanRadiusOrig = g['PercentageOfPointsCloserThanRadiusOrig']
        PercentageOfPointsCloserThanRadius = g['PercentageOfPointsCloserThanRadius']
        OverlapMeasure = g['OverlapMeasure']

        PercentageOfBifurPointsCloserThanRadiusOrig = g['PercentageOfBifurPointsCloserThanRadiusOrig']
        PercentageOfBifurPointsCloserThanRadius = g['PercentageOfBifurPointsCloserThanRadius']

        averageAngleAtBifur = g['averageAngleAtBifur']

        outResult.loc[i] = (ThresholdValue,
                            np.mean(PercentageOfPointsCloserThanRadiusOrig),
                            np.std(PercentageOfPointsCloserThanRadiusOrig),
                            np.mean(PercentageOfPointsCloserThanRadius), np.std(PercentageOfPointsCloserThanRadius),
                            np.mean(OverlapMeasure), np.std(OverlapMeasure),
                            np.mean(PercentageOfBifurPointsCloserThanRadiusOrig), np.std(PercentageOfBifurPointsCloserThanRadiusOrig),
                            np.mean(PercentageOfBifurPointsCloserThanRadius),np.std(PercentageOfBifurPointsCloserThanRadius),
                            np.mean(averageAngleAtBifur),np.std(averageAngleAtBifur))
    return outResult


def doAnalyzeOurOverlapMeasureCsv(args):
    dirname = args.dirname
    basename = args.basename
    voxelSize = args.voxelSize

    df = pd.read_csv(os.path.join(dirname, basename))

    overlapMeasure = pd.DataFrame()

    for i, (ParValue, g) in enumerate(df.groupby(['ParValue'])):
        g = analyzeOverlapMeasure(g, voxelSize)
        g.insert(1, 'ParValue', ParValue)

        overlapMeasure = pd.concat((overlapMeasure, g))

    basename, _ = os.path.splitext(basename)
    filename = os.path.join(dirname, basename + 'Average.csv')
    overlapMeasure.sort_values(['ThresholdValue', 'ParValue'], inplace=True)
    overlapMeasure.to_csv(filename, index=False)


def doComputeConnectivityMeasure(args):
    dirname_originalVolume = args.dirname_originalVolume
    dirname_noisyVolume = args.dirname_noisyVolume
    dirname_volume = args.dirname_volume
    dirname_reconstructed = args.dirname_reconstructed

    voxelSize = args.voxelWidth
    samplingStep = args.samplingStep
    pointsArrName = args.points

    doOutputHeader = args.doOutputHeader
    prependHeaderStr = args.prependHeaderStr
    prependRowStr = args.prependRowStr

    filename = os.path.join(dirname_originalVolume, dirname_volume, 'tree_structure.xml')
    dataset = IO2.readGxlFile(filename)

    positions = dataset['positions']
    positions = voxelSize * positions

    indicesOrig1 = dataset['indices1']
    indicesOrig2 = dataset['indices2']
    radiuses = dataset['radiusPrime']
    maxRange = dataset['MaxRange']
    nodeTypes = dataset['nodeTypes']

    filename = os.path.join(dirname_noisyVolume, dirname_volume, dirname_reconstructed)
    dataset = IO2.readH5File(filename)

    Points = dataset[pointsArrName]

    indices1 = dataset['indices1']
    indices2 = dataset['indices2']

    maxrad = np.maximum(radiuses, np.sqrt(2) * voxelSize)

    adjacentIndicesOrig = [[] for i in np.arange(positions.shape[0])]
    for i, (index1, index2) in enumerate(zip(indicesOrig1, indicesOrig2)):
        adjacentIndicesOrig[index1].append(index2)
        adjacentIndicesOrig[index2].append(-1 * index1)

    pointsOrig, indicesOrig, bifurOrig = resamplePointsOnGT(positions, nodeTypes, maxRange, indicesOrig1, indicesOrig2,
                                                         samplingStep)

    learner = NearestNeighbors(n_neighbors=1)
    learner.fit(pointsOrig)

    dist, closestIndices = learner.kneighbors(Points)
    dist = np.ravel(dist)
    closestIndices = np.ravel(closestIndices)
    closerThanRadius = dist < maxrad[indicesOrig[closestIndices]]

    adjacentIndices = [[] for i in np.arange(Points.shape[0])]
    newIndices1_ = []
    newIndices2_ = []
    for i, (index1, index2) in enumerate(zip(indices1, indices2)):
        if closerThanRadius[index1] and closerThanRadius[index2]:
            adjacentIndices[index1].append(index2)
            adjacentIndices[index2].append(index1)
            newIndices1_.append(index1)
            newIndices2_.append(index2)

    # connectivity construction
    connectivity = [0 for i in np.arange(indicesOrig1.shape[0])]
    FalseCount = 0
    TotalFalse = 0
    TotalConnectivity = len(newIndices1_)
    newIndices1 = []
    newIndices2 = []
    for i, (index1, index2) in enumerate(
            zip(newIndices1_, newIndices2_)):  # limited to correctly detected points forming the tubular graph
        k1 = indicesOrig[closestIndices[index1]]
        k2 = indicesOrig[closestIndices[index2]]

        k1Index1 = indicesOrig1[k1]
        k1Index2 = indicesOrig2[k1]
        k2Index1 = indicesOrig1[k2]
        k2Index2 = indicesOrig2[k2]

        if k1 == k2:
            newIndices1.append(index1)
            newIndices2.append(index2)
        elif k1Index2 == k2Index1:  # k1->k2
            connectivity[k2] = 1
            newIndices1.append(index1)
            newIndices2.append(index2)
        elif k2Index2 == k1Index1:  # k2->k1
            connectivity[k1] = 1
            newIndices1.append(index1)
            newIndices2.append(index2)
        elif k1Index1 == k2Index1:  # after bifurcation
            FalseCount += 1
            TotalFalse += 1
        else:
            TotalFalse += 1
    connectivity = np.array(connectivity)
    # newIndices1 and newIndices2 give the correct connectivities
    newIndices1 = np.array(newIndices1)
    newIndices2 = np.array(newIndices2)

    newPoints, newIndices = resamplePoints(Points, newIndices1, newIndices2, samplingStep)
    learner = NearestNeighbors(n_neighbors=1)
    learner.fit(newPoints)

    distOrig, closestIndicesOrig = learner.kneighbors(pointsOrig)
    distOrig = np.ravel(distOrig)

    closerThanRadiusOrig = distOrig < maxrad[indicesOrig]
    ConnectivityRatio = np.count_nonzero(closerThanRadiusOrig) * 1.0 / (len(pointsOrig))

    FalseRatio = FalseCount * 1.0 / TotalConnectivity

    keyValPairs = [(name, eval(name)) for name in ('FalseRatio', 'ConnectivityRatio', 'TotalConnectivity')]

    if (doOutputHeader):
        print(prependHeaderStr + (",".join(kvp[0] for kvp in keyValPairs)))

    print(prependRowStr + (",".join(str(kvp[1]) for kvp in keyValPairs)))


if __name__ == '__main__':
    # create the top-level parser
    argparser = argparse.ArgumentParser()
    subparsers = argparser.add_subparsers()

    # create the parser for the "doComputeOverlapMeasure" command
    subparser = subparsers.add_parser('doComputeOverlapMeasure')
    subparser.add_argument('dirname_originalVolume')
    subparser.add_argument('dirname_noisyVolume')
    subparser.add_argument('dirname_volume')
    subparser.add_argument('dirname_reconstructed')
    subparser.add_argument('voxelWidth', type=float, help="voxel width of the original volume")
    subparser.add_argument('samplingStep', type=float, help="absolute step size of the sampling")
    subparser.add_argument('upperBound', type=float, default=1e10, help="only vessels with radius smaller than upperBound are used to compute the measure")
    subparser.add_argument('lowerBound', type=float, default=0, help="only vessels with radius larger than lowerBound are used to compute the measure")
    subparser.add_argument('K', type=float, help="Kth percentile in Hausdorff distance")
    subparser.add_argument('--points', default='positions')
    subparser.add_argument('--doOutputHeader', default=False, action='store_true')
    subparser.add_argument('--prependHeaderStr', default="")
    subparser.add_argument('--prependRowStr', default="")
    subparser.add_argument('--debugMode', default=False, action='store_true')
    subparser.add_argument('--flagSquare', default=False, action='store_true')
    subparser.set_defaults(func=doComputeOverlapMeasure)

    # create the parser for the "doAnalyzeOurOverlapMeasureCsv" command
    subparser = subparsers.add_parser('doAnalyzeOurOverlapMeasureCsv')
    subparser.add_argument('dirname')
    subparser.add_argument('basename')
    subparser.add_argument('voxelSize', type=float)
    subparser.set_defaults(func=doAnalyzeOurOverlapMeasureCsv)

    # create the parser for the "ComputeConnectivityMeasure" command
    subparser = subparsers.add_parser('doComputeConnectivityMeasure')
    subparser.add_argument('dirname_originalVolume')
    subparser.add_argument('dirname_noisyVolume')
    subparser.add_argument('dirname_volume')
    subparser.add_argument('dirname_reconstructed')
    subparser.add_argument('voxelWidth', type=float)
    subparser.add_argument('samplingStep', type=float)
    subparser.add_argument('--points', default='positions')
    subparser.add_argument('--doOutputHeader', default=False, action='store_true')
    subparser.add_argument('--prependHeaderStr', default="")
    subparser.add_argument('--prependRowStr', default="")
    subparser.set_defaults(func=doComputeConnectivityMeasure)

    # parse the args and call whatever function was selected
    args = argparser.parse_args()
    args.func(args)
