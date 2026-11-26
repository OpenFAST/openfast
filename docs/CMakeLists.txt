#*********************************************************************************************************************************
# LICENSING
# Copyright (C) 2017-2025  National Renewable Energy Laboratory
# Copyright (C) 2025-2026  National Laboratory of the Rockies
#
#    This file is part of OpenFAST.
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
#*********************************************************************************************************************************
find_package(Doxygen)
find_package(Sphinx)

if(NOT DOXYGEN_FOUND)
    message(WARNING "Skipping Doxygen documentation due to missing dependencies.")
else()
    configure_file(Doxyfile.in Doxyfile @ONLY)
    file(COPY "DoxygenLayout.xml" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    
    add_custom_target(doxygen
        COMMAND ${DOXYGEN_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile" )
endif()

if(NOT SPHINX_FOUND)
    message(WARNING "Skipping Sphinx documentation due to missing dependencies.")
else()
    set(SPHINX_GENERATOR html)
    file(COPY "conf.py" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    file(COPY "_static" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

    add_custom_target(sphinx
        COMMAND ${SPHINX_EXECUTABLE} -b ${SPHINX_GENERATOR}
        -c ${CMAKE_CURRENT_BINARY_DIR}
        "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/${SPHINX_GENERATOR}"
        )
    add_dependencies(sphinx doxygen)
    add_custom_target(docs DEPENDS sphinx)

    # Add a sphinx-only HTML target to avoid building doxygen while developing documentation
    add_custom_target(sphinx-html
        COMMAND ${SPHINX_EXECUTABLE} -M html
        "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}")

    # PDF documentation generation using sphinx -> latex -> pdflatex
    find_package(LATEX COMPONENTS PDFLATEX BIBTEX)
    if (LATEX_FOUND)
        add_custom_target(sphinx-pdf
            COMMAND ${SPHINX_EXECUTABLE} -M latexpdf
            "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}")
    endif()
endif()
