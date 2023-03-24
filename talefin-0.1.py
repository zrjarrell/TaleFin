import pandas as pd
import os
import graphicsReport
import json

#testdata
pathToMummichogOutput = "/Users/zrj/Documents/research/scripts/TaleFin/test-data/cd-xi-f_adef_lung_c18neg_limma2way-pls/1677524736.87.2way-outer-join-fdr20/tsv/mcg_pathwayanalysis_2way-outer-join-fdr20.xlsx"

pathToMummichogInput = "/Users/zrj/Documents/research/scripts/TaleFin/test-data/cd-xi-f_adef_lung_c18neg_limma2way-pls/formummichog-outer-join-2way-sigresults-fdr20.txt"

pathToFeatureTable = "/Users/zrj/Documents/research/scripts/TaleFin/test-data/cd-xi-f_adef_lung_c18neg_featuretable.txt"

pathToClasslist = "/Users/zrj/Documents/research/scripts/TaleFin/test-data/cd-xi-f_adef_lung_classlist.txt"

outputPath = "./output1"

#testdata2
pathToMummichogOutput = "/Users/zrj/Documents/research/scripts/TaleFin/test-data-2/cd-xi_males_testes_hilicpos_limma-pls/1648495652.23.hilicpos_p05vip2/tsv/mcg_pathwayanalysis_hilicpos_p05vip2.xlsx"

pathToMummichogInput = "/Users/zrj/Documents/research/scripts/TaleFin/test-data-2/cd-xi_males_testes_hilicpos_limma-pls/hilicpos_formummichog.txt"

pathToFeatureTable = "/Users/zrj/Documents/research/scripts/TaleFin/test-data-2/cd-xi_males_testes_hilicpos_featuretable.txt"

pathToClasslist = "/Users/zrj/Documents/research/scripts/TaleFin/test-data-2/cd-xi_males_testes_classlist.txt"

outputPath = "./output2"

#testdata3
pathToMummichogAnnotations = "/Users/zrj/Documents/research/scripts/TaleFin/test-data-3-ma/mummichog_matched_compound_all.csv"

pathToMummichogPathwayResults = "/Users/zrj/Documents/research/scripts/TaleFin/test-data-3-ma/mummichog_pathway_enrichment.csv"

pathToMummichogInput = "/Users/zrj/Documents/research/scripts/TaleFin/test-data-3-ma/liver_bothcolumns_formummichog_fdr20vip2.txt"

pathToFeatureTableNeg = "/Users/zrj/Documents/research/scripts/TaleFin/test-data-3-ma/davis-se_liver_c18neg_featuretable.txt"

pathToFeatureTablePos = "/Users/zrj/Documents/research/scripts/TaleFin/test-data-3-ma/davis-se_liver_hilicpos_featuretable.txt"

pathToClasslistNeg = "/Users/zrj/Documents/research/scripts/TaleFin/test-data-3-ma/davis-se_liver_classlist.txt"

pathToClasslistPos = "/Users/zrj/Documents/research/scripts/TaleFin/test-data-3-ma/davis-se_liver_classlist.txt"

outputPath = "./output4"

mummichogSelectedFeatureThreshold = 0.05
mummichogSignificanceThreshold = 0.05
mummichogMinimumOverlap = 4

expSetup = "classification" # 'classification' or 'regression'

organism = "hsa" #human: "hsa", mouse: "mmu", chickenk: "gga"

mummichogVersion = 1

with open('keggMapsIDs.json') as json_file:
    keggPathwayMaps = json.load(json_file)

os.mkdir(outputPath)
os.chdir(outputPath)
mummichogInput = pd.read_table(pathToMummichogInput, sep="\t")

if mummichogVersion == 1:
    import v1mummichogIO
    classes = pd.read_table(pathToClasslist, sep="\t", index_col = 0)
    features = pd.read_table(pathToFeatureTable, sep="\t")
    mummichogInput = v1mummichogIO.floorCeilRoundInputs(mummichogInput, 4)
    sigPathwayTable, annotationsTable = v1mummichogIO.splitMummichogOutput(pathToMummichogOutput, mummichogMinimumOverlap, mummichogSignificanceThreshold)
    pathways = v1mummichogIO.buildPathways(sigPathwayTable)
    annotationsTable = v1mummichogIO.fixAnnotationRounding(annotationsTable, mummichogInput)
    features = v1mummichogIO.fixFeatureTableRounding(features, mummichogInput, annotationsTable)
    allIDs = v1mummichogIO.buildKEGGNodes(pathways, annotationsTable, mummichogInput)
    featureDict = v1mummichogIO.compileFeatureData(features, allIDs, mummichogSelectedFeatureThreshold, classes)
    v1mummichogIO.mapFeaturesToKEGG(featureDict, allIDs)
#elif mummichogVersion == 'MA':
    #import vMAmummichogIO

graphicsReport.createGraphs(featureDict, expSetup)

graphicsReport.reducePathwayCards(pathways, allIDs)

graphicsReport.generateHTMLreport(pathways, pathToMummichogOutput, mummichogSignificanceThreshold, mummichogMinimumOverlap, allIDs, featureDict, mummichogSelectedFeatureThreshold, keggPathwayMaps, "hsa")
graphicsReport.generateStylesheet()

