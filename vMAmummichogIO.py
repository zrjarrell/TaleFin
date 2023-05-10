import pandas as pd
import math
import pathwayObjects

#################################################################
#
#         Object preparation for mummichogMA
#
#################################################################

#reads annotations table from MA output to dataframe
def readAnnotations(path):
    annotations = pd.read_table(path, sep=",")
    return annotations

#reads pathway results from MA output to dataframe, filters based on minimum overlap and maximum p-value supplied by user
def readPathways(path, overlapMin, pvalueMax, pvalueType):
    pathwayResults = pd.read_table(path, sep=",")
    sigPathways = pathwayResults[(pathwayResults['Hits.sig'] >= overlapMin) & (pathwayResults[pvalueType] < pvalueMax)].reset_index(drop=True)
    return sigPathways

#rounds up on 5 at deciding decimal position, regardless of outcome (decimal < 0 works on left of decimal position)
def properRound(n, decimals=0):
    expoN = n * 10 ** decimals
    if abs(expoN) - abs(math.floor(expoN)) < 0.5:
        result =  math.floor(expoN) / 10 ** decimals
    else:
        result = math.ceil(expoN) / 10 ** decimals
    if decimals <= 0:
        return int(result)
    return result

#rounds down at specified decimal position (decimal < 0 works on left of decimal position)
def floorRound(n, decimals=0):
    expoN = n * 10 ** decimals
    result = math.floor(expoN) / 10 ** decimals
    if decimals <= 0:
        return int(result)
    return result

#rounds up at specified decimal position (decimal < 0 works on left of decimal position)
def ceilRound(n, decimals=0):
    expoN = n * 10 ** decimals
    result = math.ceil(expoN) / 10 ** decimals
    if decimals <= 0:
        return int(result)
    return result

#adds three new sets of mz and time to table used for mummichog input: traditional rounding, the floor and the ceiling.
#this is essential for capturing all data when aligning between the mummichog input, the original feature table, and the mummichog annotations
#mummichog v1 does not account for issues with python's inbuilt round() function [floating point representation]
#additionally, any rounding of feature mz or time that may have occurred during upstream analysis in R could introduce discrepancies between tables
#R uses a fourth rounding method, statistical rounding (or "bankers rounding").
#building this set of three references prior to alignment with the feature table and mummichog annotations will ensure all data gets
#captured at the correct precision of mz and time measurements
def floorCeilRoundInputs(inputList, MZdecimal=4, timeDecimal=0):
    floorMZs = []
    ceilMZs = []
    roundMZs = []
    floorTimes = []
    ceilTimes = []
    roundTimes = []
    for i in inputList.index:
        floorMZs += [floorRound(inputList.loc[i, 'mz'], MZdecimal)]
        ceilMZs += [ceilRound(inputList.loc[i, 'mz'], MZdecimal)]
        roundMZs += [properRound(inputList.loc[i, 'mz'], MZdecimal)]
        floorTimes += [floorRound(inputList.loc[i, 'time'], timeDecimal)]
        ceilTimes += [ceilRound(inputList.loc[i, 'time'], timeDecimal)]
        roundTimes += [properRound(inputList.loc[i, 'time'], timeDecimal)]
    inputList['floorMZ'] = floorMZs
    inputList['ceilMZ'] = ceilMZs
    inputList['roundMZ'] = roundMZs
    inputList['floorTime'] = floorTimes
    inputList['ceilTime'] = ceilTimes
    inputList['roundTime'] = roundTimes
    return inputList

#aligns annotations with mummichog input. Will replace mz and time columns in mummichog annotations with the correct mz and time
#which were found to match in the mummichog input
def fixAnnotationRounding(annotations, inputList):
    newAnnotations = pd.DataFrame(columns=annotations.columns)
    newAnnotations['mz'] = ''
    newAnnotations['time'] = ''
    for row in annotations.index:
        #searches two possible results of rounding (floor and ceiling) and finally assigns traditionally rounded values
        matches = inputList[(inputList['floorMZ'] == annotations.iloc[row, 0]) | (inputList['ceilMZ'] == annotations.iloc[row, 0])]
        matches = matches.drop_duplicates()
        for match in matches.index: 
            rowData = annotations.iloc[row, :]
            rowData['mz'] = matches.loc[match,'roundMZ']
            rowData['time'] =  matches.loc[match, 'roundTime']
            newAnnotations.loc[len(newAnnotations)] = rowData
    return newAnnotations

#similar to fixAnnotationRounding, this function creates an mz and time pair for the feature table that will correctly align with that in the mummichog input
def fixFeatureTableRounding(featureTable, inputList, annotations, MZdecimal=4, timeDecimal=0):
    selectedFeatureTable = pd.DataFrame(columns=featureTable.columns)
    selectedFeatureTable['roundMZ'] = ''
    selectedFeatureTable['roundTime'] = ''
    #using annotations table to reduce inputList to only those which will be included in the output, cuts time needed significantly
    annotatedFeatures = annotations.iloc[:, -2:].drop_duplicates(ignore_index=True)
    annotatedFeatures.columns = ['roundMZ', 'roundTime']
    mappedInputs = inputList.merge(annotatedFeatures, how='inner', on=['roundMZ','roundTime'])
    #searches for round-proof matches between feature table and and annotated input mzs
    for feature in featureTable.index:
        floorMZ = floorRound(featureTable.loc[feature, :]['mz'], MZdecimal)
        floorTime = floorRound(featureTable.loc[feature, :]['time'], timeDecimal)
        ceilMZ = ceilRound(featureTable.loc[feature, :]['mz'], MZdecimal)
        ceilTime = ceilRound(featureTable.loc[feature, :]['time'], timeDecimal)
        matches = mappedInputs[((mappedInputs['roundMZ']==floorMZ) & (mappedInputs['roundTime']==floorTime)) |
                            ((mappedInputs['roundMZ']==floorMZ) & (mappedInputs['roundTime']==ceilTime)) |
                            ((mappedInputs['roundMZ']==ceilMZ) & (mappedInputs['roundTime']==floorTime)) |
                            ((mappedInputs['roundMZ']==ceilMZ) & (mappedInputs['roundTime']==ceilTime))]
        matches = matches.drop_duplicates()
        for match in matches.index:
            rowData = featureTable.iloc[feature, :]
            rowData['roundMZ'] = matches.loc[match, 'roundMZ']
            rowData['roundTime'] = matches.loc[match, 'roundTime']
            selectedFeatureTable.loc[len(selectedFeatureTable)] = rowData
    if timeDecimal <= 0:
        selectedFeatureTable['roundTime'] = selectedFeatureTable['roundTime'].astype(int)
    #returns reduced featureTable (only annotated features) with a roundMZ and roundTime matching that in the inputList
    return selectedFeatureTable













#need function that does following:
#   builds a list of mzs in pathway table to resemble v1
#   formats annotation table to resemble v1
#   formats pathway table to resemble v1






#iterates through pathway table to build pathway objects
def buildPathways(sigPathways):
    pathways = []
    for path in sigPathways.index:
        pathways += [pathwayObjects.Pathway(sigPathways, path)]
    return pathways

#iterates through pathway objects in pathway list to build KEGGNode objects, makes sure no redundant Nodes are made
#stores in library, KEGG ID is the key
def buildKEGGNodes(pathways, annotations, inputList):
    allKEGGIDs = {}
    for path in pathways:
        for keggID in path.keggIDs:
            if keggID in allKEGGIDs:
                allKEGGIDs[keggID].addPathway(path)
            else:
                allKEGGIDs[keggID] = pathwayObjects.KEGGNode(keggID, path, annotations)
    for keggID in allKEGGIDs:
        allKEGGIDs[keggID].addStats(inputList)
    return allKEGGIDs

#iterates through list of features to build feature objects, stored in dictionary
def buildUniqueFeatures(fullList):
    uniqueFeatures = {}
    uniqueList = fullList.drop_duplicates(subset=['mz', 'time']).reset_index(drop=True)
    for unique in uniqueList.index:
        uniqueFeatures[str(uniqueList.iloc[unique, 1]) + '_' + str(uniqueList.iloc[unique, 2])] = pathwayObjects.featureNode(unique, uniqueList)
    #goes back through list with duplicates to add matched KEGG IDs to each unique feature
    for row in fullList.index:
        uniqueFeatures[str(fullList.iloc[row, 1]) + '_' + str(fullList.iloc[row, 2])].keggMatches += [fullList.iloc[row, 0]]
    return uniqueFeatures, uniqueList

#generates feature objects
def compileFeatureData(features, nodes, threshold, classes):
    #iterates through each node to build a list of all features, includes duplicates
    df = pd.DataFrame(columns= ['id', 'mz', 'time'])
    d = 0
    for node in nodes:
        c = 0
        for pvalue in nodes[node].pvalues:
            if pvalue < threshold:
                df.loc[d] = [nodes[node].id, nodes[node].mzs[c], nodes[node].times[c]]
                d += 1
            c += 1
    #uses list to create unique feature objects
    uniqueFeatures, uniqueList = buildUniqueFeatures(df)
    features = features.iloc[:, 2:].rename(columns={'roundMZ': 'mz', 'roundTime': 'time'})
    df = uniqueList.merge(features, how='inner', on=['mz','time'])
    #generates an intensity table with independent variable and feature intensity and stores it in each feature object
    for feature in uniqueFeatures:
        featureData = df[(df['mz']== uniqueFeatures[feature].mz) & (df['time']== uniqueFeatures[feature].time)]
        featureData = featureData.iloc[:, 3:].T
        featureData = classes.join(featureData)
        featureData.columns = ['class', 'intensity']
        uniqueFeatures[feature].intensityTable = featureData
    return uniqueFeatures

#reads connected KEGG IDs from each feature object to create a link from KEGG ID to object
def mapFeaturesToKEGG(featureDict, keggDict):
    for feature in featureDict:
        for keggID in featureDict[feature].keggMatches:
            keggDict[keggID].sigFeatures += [feature]
