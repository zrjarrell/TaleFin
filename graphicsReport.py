import seaborn as sns
import matplotlib.pyplot as plt
from datetime import datetime
import os

#################################################################
#
# For graph generation and writing the html & css for the report
#
#################################################################

#makes graphs for each feature object. Makes scatter plot with linear reg line for expSetup = 'regression' or
#a boxplot for expSetup = 'classification'. saves file path in object
def createGraphs(featureDict, expSetup):
    os.mkdir(os.path.join('./graphs'))
    pdfDir = os.path.join('./graphs/pdfs')
    svgDir = os.path.join('./graphs/svgs')
    os.mkdir(pdfDir)
    os.mkdir(svgDir)
    if expSetup != 'classification' and expSetup != 'regression':
        print("Please specify your type of analysis, 'regression' or 'categorical'.")
        return
    else:
        for feature in featureDict:
            if expSetup == 'classification':
                graph = sns.boxplot(x='class', y='intensity', data=featureDict[feature].intensityTable)
                graph = sns.stripplot(x='class', y='intensity', data=featureDict[feature].intensityTable, color='orange', jitter=0.2, size=2.5)
            elif expSetup == 'regression':
                graph = sns.regplot(x='class', y='intensity', data=featureDict[feature].intensityTable, line_kws={"color":"r","alpha":0.7,"lw":5})
            plt.title(featureDict[feature].label)
            pdfPath = os.path.join(pdfDir, featureDict[feature].label + '.pdf')
            svgPath = os.path.join(svgDir, featureDict[feature].label + '.svg')
            plt.savefig(pdfPath)
            plt.savefig(svgPath)
            plt.clf()
            featureDict[feature].pdfPath = pdfPath
            featureDict[feature].svgPath = svgPath

#handles multiple KEGG IDs represented by the same set of annotated features. For duplicates, summarizes duplications in a label attribute
def reducePathwayCards(pathways, keggNodes):
    for pathway in pathways:
        pathway.reducedNodesLabels = []
        pathway.reducedNodesKeys = []
        for newID in pathway.keggIDs:
            repeated = False
            c = 0
            for oldID in pathway.reducedNodesKeys:
                if sorted(keggNodes[oldID].mzs) == sorted(keggNodes[newID].mzs) and sorted(keggNodes[oldID].times) == sorted(keggNodes[newID].times) and sorted(keggNodes[oldID].adducts) == sorted(keggNodes[newID].adducts):
                    repeated = True
                    pathway.reducedNodesLabels[c] += (', ' + newID)
                c += 1
            if not repeated:
                pathway.reducedNodesLabels += [newID]
                pathway.reducedNodesKeys += [newID]

#writes card display to report.html
def writeKEGGcard(keggID, keggDict, featureDict, pathway, c, threshold):
    keggCard = """<div class="keggCard">
        <p><b>%s</b></p>
        <div class="graphContainer">
    """%(pathway.reducedNodesLabels[c])
    for feature in keggDict[keggID].sigFeatures:
        keggCard += """<img src="%s" width="450">"""%(featureDict[feature].svgPath)
    keggCard += """</div>
    <table>
        <tr>
            <th><i>m/z</i></th>
            <th>Time</th>
            <th>Adduct</th>
            <th>p-value</th>
            <th>t-statistic</th>
        </tr>"""
    d = 0
    for feature in keggDict[keggID].mzs:
        if keggDict[keggID].pvalues[d] < threshold:
            keggCard += """<tr class="significant">"""
        else:
            keggCard += """<tr>"""
        keggCard += """<td>%s</td>
                <td>%s</td>
                <td>%s</td>
                <td>%s</td>
                <td>%s</td>
            </tr>"""%(str(round(feature, 4)), str(int(round(keggDict[keggID].times[d], 0))), keggDict[keggID].adducts[d], str(round(keggDict[keggID].pvalues[d], 6)), keggDict[keggID].tstatistics[d])
        d += 1
    keggCard += "</table>"
    keggCard += "</div>"
    return keggCard

#writes pathway container for report.html
def writePathwayBlock(pathway, keggDict, featureDict, threshold, keggMaps, organism="map"):
    keggLink = "https://www.kegg.jp/pathway/" + organism + keggMaps[pathway.name]
    for keggID in pathway.keggIDs:
        keggLink += '+' + keggID
    pathwayTitleButton = """
    <button type="button" class="pathwayButton">
        <div class="pathwayLabel">
            <div class="pathwayName"><b>%s</b></div>
            <div class="pathwaySignificance"><b>%s</b></div>
        </div>
        <div><a href="%s" class="keggLink" target="_blank">View Kegg Map</a></div>
    </button>"""%(pathway.name + " - (" + str(pathway.overlapSize) + "/" + str(pathway.pathwaySize) + ")", "p-value: " + str(round(pathway.pValue, 8)), keggLink)
    pathwayContainer = """<div class="pathwayContainer">"""
    c = 0
    for keggID in pathway.reducedNodesKeys:
        pathwayContainer += writeKEGGcard(keggID, keggDict, featureDict, pathway, c, threshold)
        c += 1
    pathwayContainer += """</div>"""
    return pathwayTitleButton + pathwayContainer

#writes report.html for user friendly mummichog summary
def generateHTMLreport(pathways, mummichogOutput, pathwayThreshold, featureThreshold, keggDict, featureDict, selectionThreshold, keggMaps, organism="map"):
    report = open('report.html', 'w')
    globalLink = "https://www.kegg.jp/pathway/" + organism + "01100"
    for keggID in keggDict:
        globalLink += '+' + keggID
    heading = """<!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="utf-8">
        <title>TaleFin Report</title>
        <link href="style.css" rel="stylesheet">
    </head>
    <body>
    <div class="heading">
        <h2><b>TaleFin Report</b></h2>
        <div>
            <p><b>Report for mummichog output:</b> %s</p>
            <p><a href="%s" class="keggLink" target="_blank">View Global Metabolism</a></p>
        </div>
        <p><b>Date:</b> %s</p>
        <p><b>Pathway significance threshold:</b> %s</p>
        <p><b>Required selected features:</b> %s</p>
    </div>
    """%(os.path.basename(mummichogOutput), globalLink, str(datetime.now()), pathwayThreshold, featureThreshold)
    report.write(heading)
    for pathway in pathways:
        pathwayDiv = writePathwayBlock(pathway, keggDict, featureDict, selectionThreshold, keggMaps)
        report.write(pathwayDiv)
    footing = """<script>
    var pathwayButton = document.getElementsByClassName("pathwayButton");
    var i;
    
    for (i = 0; i < pathwayButton.length; i++) {
        pathwayButton[i].addEventListener("click", function() {
            if (event.target.classList.contains("keggLink")) return;
            this.classList.toggle("active");
            var content = this.nextElementSibling;
            if (content.style.display === "flex") {
                content.style.display = "none";
            } else {
                content.style.display = "flex";
            }
        });
    }
    </script>
    </body>
    </html>"""
    report.write(footing)
    report.close()

#writes a style.css document for styling of report.html
def generateStylesheet():
    stylesheet = open('style.css', 'w')
    styles = """
    * {
    font-family: Arial;
    }

    .heading {
        margin: 40px 10px 0;
    }

    .heading > div {
        display: flex;
        justify-content: space-between;
    }

    .heading > div > p {
        margin: 0;
    }

    h2, p {
        padding: 4px;
    }

    .pathwayButton {
        padding: 16px;
        width: 100%;
        margin: 2px 0 0 0;
        border: none;
        text-align: left;
        outline: none;
        font-size: 15px;
        cursor: pointer;
        background-color: rgb(208, 182, 254);
        color: #111;
        display: flex;
        justify-content: space-between;
        border-radius: 8px;
    }

    .active, .pathwayButton:hover {
        background-color: rgb(171, 122, 255);
    }

    .pathwayLabel {
        display: flex;
        justify-content: space-between;
        width: 90%;
    }

    .keggLink,
    .keggLink:visited {
        color: black;
        font-weight: bold;
        width: 10%;
        min-width: 9em;
        text-align: right;
    }

    .pathwayContainer {
        display: none;
        padding: 16px;
        font-size: 90%;
        gap: 10px;
        flex-wrap: wrap;
    }

    .keggCard {
        min-width: 475px;
        flex-grow: 1;
        background-color: rgb(181, 226, 254);
        display: flex;
        flex-direction: column;
        align-items: center;
        border-radius: 8px;
    }

    .graphContainer {
        width: 450px;
        height: 345px;
        overflow-y: scroll;
        overflow-x: hidden;
    }

    table {
        border-collapse: collapse;
        margin: 8px;
    }

    th, td {
        border: 1px solid rgb(117, 201, 254);
        padding: 0.5em;
        text-align: left;
    }

    th {
        background-color: rgb(197, 233, 255);
    }

    td {
        background-color: rgb(211, 238, 255);
    }

    .significant {
        color: rgb(160, 0, 0);
        font-weight: bold;
    }
    """
    stylesheet.write(styles)
    stylesheet.close()