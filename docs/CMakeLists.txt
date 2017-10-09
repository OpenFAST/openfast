find_package(Doxygen)
if(DOXYGEN_FOUND)

  if (DOXYGEN_FOUND_DOT)
    set(DOXYGEN_HAVE_DOT "YES"
      CACHE STRING "Flag indicating whether Graphviz is installed")
  else()
    set(DOXYGEN_HAVE_DOT "NO"
      CACHE STRING "Flag indicating whether Graphviz is installed")
  endif()
  

  if (ENABLE_SPHINX_API_DOCS)
    set(DOXYGEN_GENERATE_XML "YES"
      CACHE STRING "Create XML tags for Sphinx Breathe")
  else()
    set(DOXYGEN_GENERATE_XML "NO"
      CACHE STRING "Create XML tags for Sphinx Breathe")
  endif()
  configure_file(Doxyfile.in Doxyfile @ONLY)
  add_custom_target(doxygen
    COMMAND ${DOXYGEN_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile" )
  
else()
  
  message(WARNING "Skipping Doxygen documentation due to missing dependencies.")
  
endif()

find_package(Sphinx)
if(SPHINX_FOUND)
    set(SPHINX_GENERATOR html)
    file(COPY "conf.py" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    file(COPY "_static" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    #file(COPY "_static/openfastlogo.jpg" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


    if (ENABLE_SPHINX_API_DOCS)
      set(SPHINX_OPTS "-t" "use_breathe")
    else()
      set(SPHINX_OPTS " ")
    endif()

    message("SPHINX_OPTS" ${SPHINX_OPTS})
    
    add_custom_target(sphinx-html
        COMMAND ${SPHINX_EXECUTABLE} -M html
        "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}"
        -c ${CMAKE_CURRENT_BINARY_DIR} ${SPHINX_OPTS})
    add_dependencies(sphinx-html doxygen)
    add_custom_target(docs DEPENDS sphinx-html)
    

else()
    message(WARNING "Skipping Sphinx documentation due to missing dependencies.")
endif()
