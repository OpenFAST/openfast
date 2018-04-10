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
    This library provides tools for plotting the output channels over time of a 
    given solution attribute for two OpenFAST solutions, with the second solution
    assumed to be the baseline for comparison. There are functions for solution
    file I/O, plot creation, and html creation for navigating the plots.
"""

import sys
import os
import numpy as np
from fast_io import load_output
import rtestlib as rtl

def _validateAndExpandInputs(argv):
    rtl.validateInputOrExit(argv, 3, "solution1 solution2 attribute")
    testSolution = argv[0]
    baselineSolution = argv[1]
    attribute = argv[2]
    rtl.validateFileOrExit(testSolution)
    rtl.validateFileOrExit(baselineSolution)
    return (testSolution, baselineSolution, attribute)

def _parseSolution(solution):
    try:
        data, info, _ = load_output(solution)
        return (data, info)
    except Exception as e:
        rtl.exitWithError("Error: {}".format(e))

def _plotError(xseries, y1series, y2series, xlabel, title1, title2):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FormatStrFormatter

    plt.figure()

    ax = plt.subplot(211)
    plt.title(title1)
    plt.grid(True)
    plt.plot(xseries, y2series, "g", linestyle="solid", linewidth=3, label = "Baseline")
    plt.plot(xseries, y1series, "r", linestyle="solid", linewidth=1, label = "Local")
    plt.legend()
    
    ax = plt.subplot(212)
    plt.title(title2)
    plt.grid(True)
    plt.plot(xseries, abs(y2series - y1series))
    plt.xlabel(xlabel)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
    
    plt.tight_layout()
    
    return plt

def _savePlot(plt, path, foutname):
    plt.savefig(os.path.join(path, foutname+".png"))

def plotOpenfastError(testSolution, baselineSolution, attribute):
    testSolution, baselineSolution, attribute = _validateAndExpandInputs([testSolution, baselineSolution, attribute])
    dict1, info1 = _parseSolution(testSolution)
    dict2, info2 = _parseSolution(baselineSolution)

    try:
        channel = info1['attribute_names'].index(attribute)
    except Exception as e:
        rtl.exitWithError("Error: Invalid channel name--{}".format(e))

    title1 = attribute + " (" + info1["attribute_units"][channel] + ")"
    title2 = "Max norm"
    xlabel = 'Time (s)'

    timevec = dict1[:, 0]
    y1series = np.array(dict1[:, channel], dtype = np.float)
    y2series = np.array(dict2[:, channel], dtype = np.float)
    plt = _plotError(timevec, y1series, y2series, xlabel, title1, title2)

    basePath = os.path.sep.join(testSolution.split(os.path.sep)[:-1])
    plotPath = os.path.join(basePath, "plots")
    rtl.validateDirOrMkdir(plotPath)
    _savePlot(plt, plotPath, attribute)
    
    plt.close()
    
def _htmlHead(title):
    head  = '<!DOCTYPE html>' + '\n'
    head += '<html>' + '\n'
    head += '<head>' + '\n'
    head += '  <title>{}</title>'.format(title) + '\n'
    head += '  <link href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css" rel="stylesheet" integrity="sha384-BVYiiSIFeK1dGmJRAkycuHAHRg32OmUcww7on3RYdg4Va+PmSTsz/K68vbdEjh4u" crossorigin="anonymous">' + '\n'
    head += '  <script src="https://ajax.googleapis.com/ajax/libs/jquery/1.12.4/jquery.min.js"></script>' + '\n'
    head += '  <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js" integrity="sha384-Tc5IQib027qvyjSMfHjOMaLkfuWVxZxUPnCJA7l2mCWNIpG9mGCD8wGNIcPD7Txa" crossorigin="anonymous"></script>' + '\n'
    head += '  <style media="screen" type="text/css">'
    head += '    .cell-warning {'
    head += '      background-color: #FF6666;'
    head += '    }'
    head += '    .cell-highlight {'
    head += '      background-color: #E5E589;'
    head += '    }'
    head += '  </style>'
    head += '</head>' + '\n'
    return head

def _htmlTail():
    tail = '</html>' + '\n'
    return tail

def _tableHead(columns):
    head  = '    <table class="table table-bordered table-hover table-sm" style="margin: auto; width: 50%">' + '\n'
    head += '      <thead>' + '\n'
    head += '        <tr>' + '\n'
    head += '          <th>#</th>' + '\n'
    for column in columns:
        head += '          <th>{}</th>'.format(column) + '\n'
    head += '        </tr>' + '\n'
    head += '      </thead>' + '\n'
    return head
    
def initializePlotDirectory(testSolution, plotList, relativeNorm, maxNorm):
    basePath = os.path.sep.join(testSolution.split(os.path.sep)[:-1])
    plotPath = os.path.join(basePath, "plots")
    caseName = basePath.split(os.path.sep)[-1]
    rtl.validateDirOrMkdir(plotPath)
    
    with open(os.path.join(plotPath, "plots.html"), "w") as html:
        
        html.write( _htmlHead(caseName) )
        
        html.write('<body>' + '\n')
        html.write('  <h2 class="text-center">{}</h2>'.format(caseName) + '\n')
        html.write('  <div class="container">' + '\n')
        html.write('  <h4 class="text-center">Maximum values for each norm are highlighted</h2>' + '\n')
        
        # Channel - Relative Norm - Max Norm
        data = [('<a href="#{0}">{0}</a>'.format(plot), relativeNorm[i], maxNorm[i]) for i,plot in enumerate(plotList)]    
        maxRelNorm = max(relativeNorm)
        maxMaxNorm = max(maxNorm)
        table = _tableHead(['Channel', 'Relative Max Norm', 'Infinity Norm'])
        body = '      <tbody>' + '\n'
        for i, d in enumerate(data):
            body += '        <tr>' + '\n'
            body += '          <th scope="row">{}</th>'.format(i+1) + '\n'
            body += '          <td>{0:s}</td>'.format(d[0]) + '\n'
            
            fmt = '{0:0.4e}'
            if d[1] == maxRelNorm:
                body += ('          <td class="cell-highlight">' + fmt + '</td>').format(d[1]) + '\n'
            else:
                body += ('          <td>' + fmt + '</td>').format(d[1]) + '\n'
                    
            if d[2] == maxMaxNorm:
                body += ('          <td class="cell-highlight">' + fmt + '</td>').format(d[2]) + '\n'
            else:
                body += ('          <td>' + fmt + '</td>').format(d[2]) + '\n'
            body += '        </tr>' + '\n'
        body += '      </tbody>' + '\n'
        table += body
        table += '    </table>' + '\n'
        html.write(table)
        
        html.write('    <br>' + '\n')
        html.write('    <div class="row">' + '\n')
        for i,plot in enumerate(plotList):
            html.write('      <div id={} class="col-sm-12 col-md-6 col-lg-6">'.format(plot) + '\n')
            html.write('        <img src="{}" class="center-block img-responsive thumbnail">'.format(plot+".png") + '\n')
            html.write('      </div>' + '\n')
        html.write('    </div>' + '\n')
        html.write('  </div>' + '\n')
        html.write('</body>' + '\n')
        html.write( _htmlTail() )
    html.close()
    
def exportResultsSummary(path, results):
    with open(os.path.join(path, "regression_test_summary.html"), "w") as html:
        
        html.write( _htmlHead("Regression Test Summary") )
        
        html.write('<body>' + '\n')
        html.write('  <h2 class="text-center">{}</h2>'.format("Regression Test Summary") + '\n')
        html.write('  <div class="container">' + '\n')
        
        # Test Case - Pass/Fail - Max Relative Norm            
        data = [('<a href="{0}/{0}.html">{0}</a>'.format(r[0]), r[1]) for i,r in enumerate(results)]
        table = _tableHead(['Test Case', 'Pass/Fail'])
        body = '      <tbody>' + '\n'
        for i, d in enumerate(data):
            body += '        <tr>' + '\n'
            body += '          <th scope="row">{}</th>'.format(i+1) + '\n'
            body += '          <td>{0:s}</td>'.format(d[0]) + '\n'
            
            fmt = '{0:s}'
            if d[1] == "FAIL":
                body += ('          <td class="cell-warning">' + fmt + '</td>').format(d[1]) + '\n'
            else:
                body += ('          <td>' + fmt + '</td>').format(d[1]) + '\n'
                
            body += '        </tr>' + '\n'
        body += '      </tbody>' + '\n'
        table += body
        table += '    </table>' + '\n'
        html.write(table)
            
        html.write('    <br>' + '\n')
        html.write('  </div>' + '\n')
        html.write('</body>' + '\n')
        html.write( _htmlTail() )
    html.close()
    
def exportCaseSummary(path, case, results):
    with open(os.path.join(path, case+".html"), "w") as html:
        
        html.write( _htmlHead(case + " Summary") )
        
        html.write('<body>' + '\n')
        html.write('  <h2 class="text-center">{}</h2>'.format(case + " Summary") + '\n')
        html.write('  <h4 class="text-center"><a href="plots/plots.html">Go To Plots</a></h2>' + '\n')
        html.write('  <h4 class="text-center">Maximum values for each norm are highlighted</h2>' + '\n')
        html.write('  <div class="container">' + '\n')
        
        # Channel - Relative Norm - Max Norm
        data = [(r[0], r[1], r[2]) for i,r in enumerate(results)]
        maxRelNorm = max([r[1] for i,r in enumerate(results)])
        maxMaxNorm = max([r[2] for i,r in enumerate(results)])
        table = _tableHead(['Channel', 'Relative Max Norm', 'Infinity Norm'])
        body = '      <tbody>' + '\n'
        for i, d in enumerate(data):
            body += '        <tr>' + '\n'
            body += '          <th scope="row">{}</th>'.format(i+1) + '\n'
            body += '          <td>{0:s}</td>'.format(d[0]) + '\n'
            
            fmt = '{0:0.4e}'
            if d[1] == maxRelNorm:
                body += ('          <td class="cell-highlight">' + fmt + '</td>').format(d[1]) + '\n'
            else:
                body += ('          <td>' + fmt + '</td>').format(d[1]) + '\n'
                    
            if d[2] == maxMaxNorm:
                body += ('          <td class="cell-highlight">' + fmt + '</td>').format(d[2]) + '\n'
            else:
                body += ('          <td>' + fmt + '</td>').format(d[2]) + '\n'
            body += '        </tr>' + '\n'
        body += '      </tbody>' + '\n'
        table += body
        table += '    </table>' + '\n'
        html.write(table)
        
        html.write('    <br>' + '\n')
        html.write('  </div>' + '\n')
        html.write('</body>' + '\n')
        html.write( _htmlTail() )
    html.close()
