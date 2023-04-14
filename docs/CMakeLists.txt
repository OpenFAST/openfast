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
