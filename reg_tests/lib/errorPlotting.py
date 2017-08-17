#
# Copyright 2017 National Renewable Energy Laboratory
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

"""
    This program plots the output vectors (vs. time) of a given solution attribute
    for two OpenFAST solutions, with the second solution assumed to be the baseline for
    comparison. It reads two OpenFAST binary output files (.outb), and
    generates three plots of the given attribute (1) comparing the two tests' respective
    values, (2) showing the difference in values, and (3) showing relative difference,
    as compared to the baseline solution.

    Usage: python plotOpenfastOut.py solution1 solution2 attribute
    Example: python plotOpenfastOut.py output-local/Test01.outb output-baseline/Test01.outb Wind1VelX
"""

import sys
import os
import numpy as np
from fast_io import load_output
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import rtestlib as rtl

def validateAndExpandInputs(argv):
    rtl.validateInputOrExit(argv, 3, "solution1 solution2 attribute")
    testSolution = argv[0]
    baselineSolution = argv[1]
    attribute = argv[2]
    rtl.validateFileOrExit(testSolution)
    rtl.validateFileOrExit(baselineSolution)
    return (testSolution, baselineSolution, attribute)

def parseSolution(solution):
    try:
        data, info, _ = load_output(solution)
        return (data, info)
    except Exception as e:
        rtl.exitWithError("Error: {}".format(e))

def plotError(xseries, y1series, y2series, xlabel, title1, title2):
    plt.figure()
    
    ax = plt.subplot(211)
    plt.title(title1)
    plt.grid(True)
    plt.plot(xseries, y1series, "g", linestyle="solid", linewidth=3, label = "Baseline")
    plt.plot(xseries, y2series, "r", linestyle="solid", linewidth=1, label = "Local")
    plt.legend()
    
    ax = plt.subplot(212)
    plt.title(title2)
    plt.grid(True)
    plt.plot(xseries, y1series - y2series)
    plt.xlabel(xlabel)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
    
    plt.tight_layout()
    
    return plt

def savePlot(plt, path, foutname):
    plt.savefig(os.path.join(path, foutname+".png"))

def plotOpenfastError(testSolution, baselineSolution, attribute):
    testSolution, baselineSolution, attribute = validateAndExpandInputs([testSolution, baselineSolution, attribute])
    dict1, info1 = parseSolution(testSolution)
    dict2, info2 = parseSolution(baselineSolution)

    try:
        channel = info1['attribute_names'].index(attribute)
    except Exception as e:
        rtl.exitWithError("Error: Invalid channel name--{}".format(e))

    title1 = attribute + " (" + info1["attribute_units"][channel] + ")"
    title2 = "Baseline - Local"
    xlabel = 'Time (s)'

    timevec = dict1[:, 0]
    y1series = np.array(dict1[:, channel], dtype = np.float)
    y2series = np.array(dict2[:, channel], dtype = np.float)
    plt = plotError(timevec, y1series, y2series, xlabel, title1, title2)

    basePath = os.path.sep.join(testSolution.split(os.path.sep)[:-1])
    plotPath = os.path.join(basePath, "plots")
    rtl.validateDirOrMkdir(plotPath)
    savePlot(plt, plotPath, attribute)
    
    plt.close()
    
def initializePlotDirectory(testSolution, plotList, relativeNorm, maxNorm):
    basePath = os.path.sep.join(testSolution.split(os.path.sep)[:-1])
    plotPath = os.path.join(basePath, "plots")
    caseName = basePath.split(os.path.sep)[-1]
    rtl.validateDirOrMkdir(plotPath)
    with open(os.path.join(plotPath, "plots.html"), "w") as plotshtml:
        plotshtml.write('<!DOCTYPE html>' + "\n")
        plotshtml.write('<html>' + "\n")
        plotshtml.write('<head>' + "\n")
        plotshtml.write('  <title>{}</title>'.format(caseName) + "\n")
        plotshtml.write('  <link href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css" rel="stylesheet" integrity="sha384-BVYiiSIFeK1dGmJRAkycuHAHRg32OmUcww7on3RYdg4Va+PmSTsz/K68vbdEjh4u" crossorigin="anonymous">' + "\n")
        plotshtml.write('  <script src="https://ajax.googleapis.com/ajax/libs/jquery/1.12.4/jquery.min.js"></script>' + "\n")
        plotshtml.write('  <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js" integrity="sha384-Tc5IQib027qvyjSMfHjOMaLkfuWVxZxUPnCJA7l2mCWNIpG9mGCD8wGNIcPD7Txa" crossorigin="anonymous"></script>' + "\n")
        plotshtml.write('</head>' + "\n")
        plotshtml.write('<body>' + "\n")
        plotshtml.write('  <h2 class="text-center">{}</h2>'.format(caseName) + "\n")
        plotshtml.write('  <div class="container">' + "\n")
        plotshtml.write('    <table class="table table-bordered table-hover table-sm" style="margin: auto; width: 50%">' + "\n")
        plotshtml.write('      <thead>' + "\n")
        plotshtml.write('        <tr>' + "\n")
        plotshtml.write('          <th>#</th>' + "\n")
        plotshtml.write('          <th>Channel</th>' + "\n")
        plotshtml.write('          <th>Relative Norm</th>' + "\n")
        plotshtml.write('          <th>Max Norm</th>' + "\n")
        plotshtml.write('        </tr>' + "\n")
        plotshtml.write('      </thead>' + "\n")
        plotshtml.write('      <tbody>' + "\n")
        for i,plot in enumerate(plotList):
            plotshtml.write('        <tr>' + "\n")
            plotshtml.write('          <th scope="row">{}</th>'.format(i) + "\n")
            plotshtml.write('          <td><a href="#{}">{}</a></td>'.format(plot, plot) + "\n")
            plotshtml.write('          <td>{0:0.4e}</td>'.format(relativeNorm[i]) + "\n")
            plotshtml.write('          <td>{0:0.4e}</td>'.format(maxNorm[i]) + "\n")
            plotshtml.write('        </tr>' + "\n")
        plotshtml.write('      </tbody>' + "\n")
        plotshtml.write('    </table>' + "\n")
        plotshtml.write('    <br>' + "\n")
        plotshtml.write('    <div class="row">' + "\n")
        for i,plot in enumerate(plotList):
            plotshtml.write('      <div id={} class="col-sm-12 col-md-6 col-lg-6">'.format(plot) + "\n")
            plotshtml.write('        <img src="{}" class="center-block img-responsive thumbnail">'.format(plot+".png") + "\n")
            plotshtml.write('      </div>' + "\n")
        plotshtml.write('    </div>' + "\n")
        plotshtml.write('  </div>' + "\n")
        plotshtml.write('</body>' + "\n")
        plotshtml.write('</html>' + "\n")

    plotshtml.close()
