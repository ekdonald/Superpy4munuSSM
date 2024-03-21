add_subdirectory(examples)

find_package(Doxygen)
if(DOXYGEN_FOUND)
  if(${PROJECT_NAME}_INTERNAL_DOC)
    set(DOXYGEN_INTERNAL_DOCS YES)
    set(INTERNALS "${PROJECT_SOURCE_DIR}/src/")
  endif()
  set(DOXYGEN_USE_MATHJAX YES)
  set(DOXYGEN_STRIP_FROM_PATH "${PROJECT_SOURCE_DIR}/include" ${INTERNALS})
  set(DOXYGEN_STRIP_FROM_INC_PATH "${PROJECT_SOURCE_DIR}/include" ${INTERNALS})
  set(DOXYGEN_FILE_PATTERNS "*.hpp;*.md;*.dox")
  set(DOXYGEN_USE_MDFILE_AS_MAINPAGE "${PROJECT_SOURCE_DIR}/README.md")
  set(DOXYGEN_JAVADOC_AUTOBRIEF YES)
  set(DOXYGEN_DISTRIBUTE_GROUP_DOC YES)
  set(DOXYGEN_ENABLE_PREPROCESSING YES)
  set(DOXYGEN_MACRO_EXPANSION YES)
  set(DOXYGEN_EXPAND_ONLY_PREDEF YES)
  set(DOXYGEN_GENERATE_XML YES)
  set(DOXYGEN_ALIASES "rststar=@verbatim embed:rst:leading-asterisk"
                      "endrststar=@endverbatim")
  set(DOXYGEN_PREDEFINED "REQUIRES(...)=" "HIGGSBOUNDS_HAS_CONCEPTS" "HIGGSTOOLS_EXPORT")
  set(DOXYGEN_EXTRACT_STATIC YES)

  doxygen_add_docs(
    dox-doc "${PROJECT_SOURCE_DIR}/include/Higgs/"
    "${PROJECT_SOURCE_DIR}/README.md"
    "${CMAKE_CURRENT_SOURCE_DIR}/Namespaces.dox" ${INTERNALS})
endif()

find_program(
  JSON_SCHEMA_FOR_HUMANS_EXECUTABLE
  NAMES generate-schema-doc
  DOC "Path to generate-schema-doc executable")
if(JSON_SCHEMA_FOR_HUMANS_EXECUTABLE)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/json-schema)
  add_custom_target(
    json-doc
    COMMAND
      ${JSON_SCHEMA_FOR_HUMANS_EXECUTABLE} --expand-buttons --config
      collapse_long_descriptions=false
      ${CMAKE_CURRENT_SOURCE_DIR}/../json/Limit.schema.json
      ${CMAKE_CURRENT_SOURCE_DIR}/_static/Limit_schema.html)
endif()

find_program(
  SPHINX_EXECUTABLE
  NAMES sphinx-build
  DOC "Path to sphinx-build executable")
if(SPHINX_EXECUTABLE AND DOXYGEN_FOUND)
  add_custom_target(
    sphinx-doc
    COMMAND
      ${SPHINX_EXECUTABLE} -b html
      -Dbreathe_projects.${PROJECT_NAME}=${CMAKE_CURRENT_BINARY_DIR}/xml
      ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/sphinx
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  add_dependencies(sphinx-doc dox-doc doc-examples)
  if(JSON_SCHEMA_FOR_HUMANS_EXECUTABLE)
    add_dependencies(sphinx-doc json-doc)
  endif()
endif()

