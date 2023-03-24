import re

#builds pathway objects. Intended to store individual pathways with related KEGG IDs, metabolites and statistical outputs from mummichog output
class Pathway:
    def __init__(self, sigPathways, pathwayIndex):
        self.name = sigPathways.iloc[pathwayIndex, 0]
        self.overlapSize = sigPathways.iloc[pathwayIndex, 1]
        self.pathwaySize = sigPathways.iloc[pathwayIndex, 2]
        self.pValue = sigPathways.iloc[pathwayIndex, 4]
        self.metabolites = list(set(sigPathways.iloc[pathwayIndex, 5].split('$')))
        self.keggIDs = list(set(sigPathways.iloc[pathwayIndex, 6].split(';')))
        self.mzs = list(set(re.split(';|,', sigPathways.iloc[pathwayIndex, 7])))
    
    def printInfo(self):
        print(self.name + " - (" + str(self.overlapSize) + "/" + str(self.pathwaySize) + ")")
        print("p-value: " + str(self.pValue))
        print("KEGG IDs: " + str(self.keggIDs))

#builds KEGG objects which can be referenced from the pathway.keggIDs list. IDs are keys to objects stored in a dictionary
class KEGGNode:
    def __init__(self, keggID, path, annotations):
        features = annotations[annotations['id']==keggID]
        self.id = keggID
        self.mzs = features.iloc[:,-2].tolist()
        self.times = features.iloc[:,-1].tolist()
        self.adducts = features.iloc[:,2].tolist()
        self.pvalues = []
        self.tstatistics = []
        self.pathways = [path]
        self.pathwayNames = [path.name]
        self.sigFeatures = []
    
    #adds additional pathways to pathway list of object. Used to handle multiple pathways for one KEGGNode and limit object redundancy
    def addPathway(self, newPathway):
        self.pathways = self.pathways + [newPathway]
        self.pathwayNames += [newPathway.name]
    
    #uses self.mzs and self.times to find matching pvalue and tstatistic from mummichog input table
    def addStats(self, inputList):
        for i in range(0, len(self.mzs)):
            matches = inputList[(inputList['roundMZ']==self.mzs[i]) & (inputList['roundTime']==self.times[i])]
            matches.drop_duplicates()
            if len(matches) == 1:
                self.pvalues += [matches.iloc[0, 2]]
                self.tstatistics += [matches.iloc[0, 3]]
            else:
                self.pvalues += matches.iloc[:, 2].tolist()
                self.tstatistics += matches.iloc[:, 3].tolist()
    
    def printInfo(self):
        print(self.id + ": ")
        print(', '.join(self.pathwayNames))
        print('\nm/z: ' + str(self.mzs))
        print('times: ' + str(self.times))
        print('adducts: ' + str(self.adducts))
        print('p-values: ' + str(self.pvalues))
        print('t-statistics: ' + str(self.tstatistics))

#stores mz feature specific information
class featureNode:
    def __init__(self, unique, uniqueList):
        self.label = str(uniqueList.iloc[unique, 1]) + '_' + str(uniqueList.iloc[unique, 2])
        self.mz = uniqueList.iloc[unique, 1]
        self.time = uniqueList.iloc[unique, 2]
        self.keggMatches = []