
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

import os
import sys
import shutil
import numpy as np

import rtestlib as rtl
from fast_io import load_output

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

def _plotError(time, test, baseline, xlabel, title1, title2, RTOL_MAGNITUDE, ATOL_MAGNITUDE):
    from bokeh.embed import components
    from bokeh.layouts import gridplot
    from bokeh.plotting import figure
    from bokeh.models.tools import HoverTool

    # Plot the baseline and test channels
    p1 = figure(title=title1)
    p1.title.align = 'center'
    p1.grid.grid_line_alpha=0.3
    p1.xaxis.axis_label = 'Time (s)'
    p1.line(time, baseline, color='green', line_width=3, legend_label='Baseline')
    p1.line(time, test, color='red', line_width=1, legend_label='Local')
    p1.add_tools(HoverTool(tooltips=[('Time','@x'), ('Value', '@y')],mode='vline'))

    # Plot the error and threshold
    p2 = figure(title=title2, x_range=p1.x_range)
    p2.title.align = 'center'
    p2.grid.grid_line_alpha = 0
    p2.xaxis.axis_label = 'Time (s)'
    p2.line(time, abs(baseline - test), color='blue', legend_label="Error")

    # Calculate the threshold
    NUMEPS = 1e-12
    ATOL_MIN = 1e-6
    baseline_offset = baseline - np.min(baseline)
    b_order_of_magnitude = np.floor( np.log10( baseline_offset + NUMEPS ) )
    rtol = 10**(-1 * RTOL_MAGNITUDE)
    atol = 10**(max(b_order_of_magnitude) - ATOL_MAGNITUDE)
    atol = max(atol, ATOL_MIN)
    passfail_line = atol + rtol * abs(baseline)
    p2.line(time, passfail_line, color='red', legend_label="Threshold")
    # p2.cross(xseries, passfail_line)

    p2.add_tools(HoverTool(tooltips=[('Time','@x'), ('Error', '@y')], mode='vline'))

    grid = gridplot([[p1, p2]], width=650, height=375, sizing_mode="scale_both")
    script, div = components(grid)
    
    return script, div

def _replace_id_div(html_string, plot):
    id_start = html_string.find('id=') + 4
    id_end = html_string[id_start:].find('"') + id_start
    html_string = plot.join((html_string[:id_start], html_string[id_end:]))
    return html_string

def _replace_id_script(html_string, plot):
    id_start = html_string.find('const render_items')    
    id_start += html_string[id_start:].find('roots')    
    id_start += html_string[id_start:].find('":"') + 3    
    id_end = html_string[id_start:].find('"') + id_start
    html_string = plot.join((html_string[:id_start], html_string[id_end:]))
    return html_string

def _save_plot(script, div, path, attribute):
    div_class = ' class="col-sm-12 col-md-6 col-lg-6"'

    file_name = "_script".join((attribute, ".txt"))
    with open(os.path.join(path, file_name), 'w') as f:
        script = script.replace('\n', '\n  ')
        script = _replace_id_script(script, attribute)
        f.write(script)
    
    file_name = "_div".join((attribute, ".txt"))
    with open(os.path.join(path, file_name), 'w') as f:
        div = _replace_id_div(div, attribute)
        ix_insert = div.find('></div>')
        div = div_class.join((div[:ix_insert], div[ix_insert:]))
        style = 'style="margin:10 auto"'
        div = div.replace("<div", " ".join(("<div", style)))
        f.write(div)

def plotOpenfastError(testSolution, baselineSolution, attribute, RTOL_MAGNITUDE, ATOL_MAGNITUDE):
    testSolution, baselineSolution, attribute = _validateAndExpandInputs([
        testSolution, baselineSolution, attribute
    ])
    dict1, info1 = _parseSolution(testSolution)
    dict2, info2 = _parseSolution(baselineSolution)

    try:
        channel = info1['attribute_names'].index(attribute)
    except Exception as e:
        rtl.exitWithError("Error: Invalid channel name--{}".format(e))

    title1 = attribute + " (" + info1["attribute_units"][channel] + ")"
    title2 = "abs(Local - Baseline)"
    xlabel = 'Time (s)'

    timevec = dict1[:, 0]
    y1series = np.array(dict1[:, channel], dtype = float)
    y2series = np.array(dict2[:, channel], dtype = float)
    script, div = _plotError(timevec, y1series, y2series, xlabel, title1, title2, RTOL_MAGNITUDE, ATOL_MAGNITUDE)

    basePath = os.path.sep.join(testSolution.split(os.path.sep)[:-1])
    plotPath = os.path.join(basePath, "plots")
    rtl.validateDirOrMkdir(plotPath)
    _save_plot(script, div, plotPath, attribute)
    
def _htmlHead(title):
    from bokeh.resources import CDN

    head  = '<!DOCTYPE html>' + '\n'
    head += '<html>' + '\n'
    head += '<head>' + '\n'
    head += '  <title>{}</title>'.format(title) + '\n'

    # CSS
    head += '  <link href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css" rel="stylesheet" integrity="sha384-BVYiiSIFeK1dGmJRAkycuHAHRg32OmUcww7on3RYdg4Va+PmSTsz/K68vbdEjh4u" crossorigin="anonymous">' + '\n'
    head += '  <link href="https://cdn.pydata.org/bokeh/release/bokeh-widgets-1.2.0.min.css" rel="stylesheet" type="text/css">' + '\n'
    head += '  <link href="https://cdn.pydata.org/bokeh/release/bokeh-1.2.0.min.css" rel="stylesheet" type="text/css">' + '\n'

    # JS
    head += '  <script src="https://ajax.googleapis.com/ajax/libs/jquery/1.12.4/jquery.min.js"></script>' + '\n'
    head += '  <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js" integrity="sha384-Tc5IQib027qvyjSMfHjOMaLkfuWVxZxUPnCJA7l2mCWNIpG9mGCD8wGNIcPD7Txa" crossorigin="anonymous"></script>' + '\n'

    # JS - Bokeh
    head += f'  <script src="{CDN.js_files[0]}"></script>' + '\n'
    head += f'  <script src="{CDN.js_files[2]}"></script>' + '\n'
    head += '  <script type="text/javascript"> Bokeh.set_log_level("info"); </script>' + '\n'
    
    head += '  <style media="screen" type="text/css">'
    head += '    .cell-warning {'
    head += '      background-color: #efc15c;'
    head += '    }'
    head += '    .cell-highlight {'
    head += '      background-color: #f5ed86 ;'
    head += '    }'
    head += '  </style>'
    head += '</head>' + '\n'
    return head

def _htmlTail():
    tail = '</html>' + '\n'
    return tail

def _tableHead(columns):
    head  = '    <table class="table table-bordered table-hover table-sm" style="margin: auto; width: 100%; font-size:80%">' + '\n'
    head += '      <thead>' + '\n'
    head += '        <tr>' + '\n'
    head += '          <th>#</th>' + '\n'
    for column in columns:
        head += f'          <th>{column}</th>' + '\n'
    head += '        </tr>' + '\n'
    head += '      </thead>' + '\n'
    return head

def finalizePlotDirectory(test_solution, plot_list, case):
    base_path = os.path.sep.join(test_solution.split(os.path.sep)[:-1])
    plot_path = os.path.join(base_path, "plots")

    with open(os.path.join(base_path, '.'.join((case, 'html'))), 'r') as html:
        html = html.read()

        # Since this routine adds to an existing html file,
        # remove the existing endings and they will be added later.
        html = html.replace("  </div>\n</body>\n</html>", "")

        script_ix = html.rfind('</script>\n') + len('</script>\n')

        for i, plot in enumerate(plot_list):
            _path = os.path.join(plot_path, plot + '_div.txt')
            try:
                with open(_path, 'r') as f:
                    div = f.read().strip().join(('      ', '\n'))
                html = ''.join((html, div))
            except FileNotFoundError:
                print("The file for ", plot ," does not exist.")

        html = ''.join((html, '    </div>' + '\n'))
        html = ''.join((html, '  </div>' + '\n'))
        html = ''.join((html, '</body>' + '\n'))
        html = ''.join((html, _htmlTail()))

    script = "" # initialize in case plot_list is empty
    for i, plot in enumerate(plot_list):
        _path = os.path.join(plot_path, f'{plot}_script.txt')
        try:
            with open(_path, "r") as f:
                _s = f.read()
            if i == 0:
                script = _s
            else:
                script = ''.join((script, _s))
        except FileNotFoundError:
            msg = "The file does not exist."
        
    shutil.rmtree(plot_path, ignore_errors=True)

    script = ''.join((script, '\n'))
    html = script.join((html[:script_ix], html[script_ix:]))
    with open(os.path.join(base_path, '.'.join((case, 'html'))), 'w') as f:
        f.write(html)
    
def exportResultsSummary(path, results):
    with open(os.path.join(path, "regression_test_summary.html"), "w") as html:
        
        html.write( _htmlHead("Regression Test Summary") )
        
        html.write('<body>' + '\n')
        html.write('  <h2 class="text-center">{}</h2>'.format("Regression Test Summary") + '\n')
        html.write('  <div class="container">' + '\n')
        
        # Test Case - Pass/Fail - Max Relative Norm            
        data = [('<a href="{0}/{0}.html">{0}</a>'.format(r[0]), r[1],r[2],'<a href="{0}/{0}.log">{0}.log</a>'.format(r[0])) for i,r in enumerate(results)]
        table = _tableHead(['Test Case', 'Pass/Fail', 'Completion Code', 'Screen Output'])
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

            if d[2] != 0:
                body += ('          <td class="cell-warning">{}</td>').format(d[2]) + '\n'
            else:
                body += ('          <td>{}</td>').format(d[2]) + '\n'
                
            body += '          <td>{0:s}</td>'.format(d[3]) + '\n'

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
    
def exportCaseSummary(path, case, channel_names, passing_channels, norms):
    """
    norms: first dimension is the channel and second dimension contains the norms
    passing_channels: 1d boolean array for whether a channel passed the comparison
    """
    
    with open(os.path.join(path, case+".html"), "w") as html:
        html.write( _htmlHead(case + " Summary") )
        
        html.write('<body>\n')
        html.write('  <h2 class="text-center">{}</h2>\n'.format(case + " Summary"))
        html.write('  <div class="container">\n')
        
        channel_tags = [ '<a href="#{0}">{0}</a>'.format(channel) for channel in channel_names ]

        cols = [
            'Channel',
            'Relative Max Norm',
            'Relative L2 Norm',
            'Infinity Norm'
        ]
        table = _tableHead(cols)
        
        body = '      <tbody>' + '\n'

        for i, channel_tag in enumerate(channel_tags):
            body += '        <tr>' + '\n'
            body += '          <th scope="row">{}</th>'.format(i+1) + '\n'
            body += '          <td>{0:s}</td>'.format(channel_tag) + '\n'
            
            fmt = '{0:0.4e}'
            for val in norms[i]:
                if not passing_channels[i]:
                    body += ('          <td class="cell-highlight">' + fmt + '</td>\n').format(val)
                else:
                    body += ('          <td>' + fmt + '</td>\n').format(val)

            body += '        </tr>' + '\n'
        body += '      </tbody>' + '\n'
        table += body
        table += '    </table>' + '\n'
        html.write(table)
        
        html.write('    <br>' + '\n')
        html.write('  </div>' + '\n')
        html.write('</body>' + '\n')
        html.write( _htmlTail() )
