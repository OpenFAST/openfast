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
import matplotlib.pyplot as plt
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
        data, info = load_output(solution)
        return (data, info)
    except Exception as e:
        rtl.exitWithError("Error: {}".format(e))

def plotError(xseries, y1series, y2series, title, xlabel, y1label, y2label):
    diff = y1series - y2series
    plt.figure()
    plt.subplot(211)
    plt.title(title)
    plt.grid(True)
    plt.ylabel(y1label)
    plt.plot(xseries, y1series, "g", linestyle="solid", linewidth=3, label = "Baseline")
    plt.plot(xseries, y2series, "r", linestyle="solid", linewidth=1, label = "Local")
    plt.legend()
    plt.subplot(212)
    plt.grid(True)
    plt.plot(xseries, diff)
    plt.ylabel(y2label)
    plt.xlabel(xlabel)
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

    # get test name -- this could break if .outb file is not used
    testname = testSolution.split("/")[-1]
    testname = testname.split(".")[-2]
    xlabel = 'Time (s)'
    y1label = attribute + " (" + info1["attribute_units"][channel] + ")"
    y2label = "Baseline - Local (" + info1['attribute_units'][channel] + ")"

    timevec = dict1[:, 0]
    y1series = np.array(dict1[:, channel], dtype = np.float)
    y2series = np.array(dict2[:, channel], dtype = np.float)
    plt = plotError(timevec, y1series, y2series, testname, xlabel, y1label, y2label)

    basePath = os.path.sep.join(testSolution.split(os.path.sep)[:-1])
    plotPath = os.path.join(basePath, "plots")
    rtl.validateDirOrMkdir(plotPath)
    savePlot(plt, plotPath, attribute)

def initializePlotDirectory(testSolution, plotList):
    basePath = os.path.sep.join(testSolution.split(os.path.sep)[:-1])
    plotPath = os.path.join(basePath, "plots")
    rtl.validateDirOrMkdir(plotPath)
    with open(os.path.join(plotPath, "plots.html"), "w") as plotshtml:
        plotshtml.write('<!DOCTYPE html>' + "\n")
        plotshtml.write('<html>' + "\n")
        plotshtml.write('<head>' + "\n")
        plotshtml.write('<title>YEEEAA</title>' + "\n")
        plotshtml.write('</head>' + "\n")
        plotshtml.write('<body>' + "\n")
        plotshtml.write('<link href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css" rel="stylesheet" integrity="sha384-BVYiiSIFeK1dGmJRAkycuHAHRg32OmUcww7on3RYdg4Va+PmSTsz/K68vbdEjh4u" crossorigin="anonymous">')
        plotshtml.write('<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/font-awesome/4.5.0/css/font-awesome.min.css">' + "\n")
        plotshtml.write('<section>' + "\n")
        plotshtml.write('<div class="container gal-container">' + "\n")

        for plot in plotList:
            plotshtml.write('<div class="col-md-8 col-sm-12 co-xs-12 gal-item">' + "\n")
            plotshtml.write('<div class="box">' + "\n")
            plotshtml.write('<a href="#" data-toggle="modal" data-target="#1">' + "\n")
            plotshtml.write('<img src="{}">'.format(plot+".png") + "\n")
            plotshtml.write('</a>' + "\n")
            plotshtml.write('<div class="modal fade" id="1" tabindex="-1" role="dialog">' + "\n")
            plotshtml.write('<div class="modal-dialog" role="document">' + "\n")
            plotshtml.write('<div class="modal-content">' + "\n")
            plotshtml.write('<button type="button" class="close" data-dismiss="modal" aria-label="Close"><span aria-hidden="true">×</span></button>' + "\n")
            plotshtml.write('<div class="modal-body">' + "\n")
            plotshtml.write('<img src="{}">'.format(plot+".png") + "\n")
            plotshtml.write('</div>' + "\n")
            plotshtml.write('<div class="col-md-12 description">' + "\n")
            plotshtml.write('<h4>This is the first one on my Gallery</h4>' + "\n")
            plotshtml.write('</div>' + "\n")
            plotshtml.write('</div>' + "\n")
            plotshtml.write('</div>' + "\n")
            plotshtml.write('</div>' + "\n")
            plotshtml.write('</div>' + "\n")
            plotshtml.write('</div>' + "\n")

        plotshtml.write('</body>' + "\n")
        plotshtml.write('</html>' + "\n")

    plotshtml.close()
