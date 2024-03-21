cmake_minimum_required(VERSION 3.17)
project(
  HiggsTools
  VERSION 1.1.3
  LANGUAGES CXX)

# -------------------------------- setup cmake --------------------------------
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
                                               "MinSizeRel" "RelWithDebInfo")
endif()

if(PROJECT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  set(${PROJECT_NAME}_IS_TOPLEVEL TRUE)
else()
  set(${PROJECT_NAME}_IS_TOPLEVEL FALSE)
endif()
# ---------------------------------- options ----------------------------------
option(${PROJECT_NAME}_BUILD_EXAMPLES "Build the ${PROJECT_NAME} examples."
       ${${PROJECT_NAME}_IS_TOPLEVEL})
option(${PROJECT_NAME}_BUILD_TESTING "Build the ${PROJECT_NAME} tests."
       ${${PROJECT_NAME}_IS_TOPLEVEL})
option(${PROJECT_NAME}_BUILD_PYTHON_MODULE
       "Build the ${PROJECT_NAME} python module (requires python)." ${SKBUILD})
option(${PROJECT_NAME}_BUILD_MATHEMATICA_INTERFACE
       "Build the WSTP-based ${PROJECT_NAME}-Mathematica interface" OFF)
option(${PROJECT_NAME}_INTERNAL_DOC
       "include documentation for ${PROJECT_NAME} internals" OFF)
option(${PROJECT_NAME}_RUN_IWYU "Run include-what-you-use with the compiler."
       OFF)
option(${PROJECT_NAME}_STATIC_CXX_STDLIB "Statically link the C++ std library"
       OFF)
option(${PROJECT_NAME}_BUILD_SHARED_LIBS
       "Build libHiggsTools as a shared library" ${${PROJECT_NAME}_IS_TOPLEVEL})
# ----------------  compiler versions and C++ standard ------------------------
if(NOT "${CMAKE_CXX_STANDARD}" OR "${CMAKE_CXX_STANDARD}" LESS 17)
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED)
endif()
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL
                                             "AppleClang")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "5")
    message(FATAL_ERROR "${PROJECT_NAME} requires clang>=5 or gcc>=9.")
  endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9")
    message(FATAL_ERROR "${PROJECT_NAME} requires clang>=5 or gcc>=9.")
  endif()
  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL "13")
    list(APPEND HiggsTools_EXTRA_COMPILE_ARGS "-Wno-dangling-reference")
  endif()
else()
  message(WARNING "Unknown compiler, proceed at your own risk.")
endif()
# ---------------------   shared library settings -----------------------------
if(${PROJECT_NAME}_BUILD_PYTHON_MODULE OR ${PROJECT_NAME}_BUILD_SHARED_LIBS)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN TRUE)
if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
  include(CheckIPOSupported)
  check_ipo_supported(
    RESULT IPO_SUPPORTED
    OUTPUT IPO_CHECK_OUTPUT
    LANGUAGES CXX)
  if(IPO_SUPPORTED)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
  else()
    message(WARNING "IPO is not supported: ${IPO_CHECK_OUTPUT}")
  endif()
endif()
include(GenerateExportHeader)
# -------------------------------  dependencies -------------------------------
add_subdirectory(external EXCLUDE_FROM_ALL)

if(${PROJECT_NAME}_BUILD_MATHEMATICA_INTERFACE)
  find_package(Mathematica REQUIRED WSTP)
endif()

# ------------------------------- documentation -------------------------------
if(PROJECT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  add_subdirectory(doc)
endif()
# ---------------------------- include-what-you-use ---------------------------
if(${PROJECT_NAME}_RUN_IWYU)
  find_program(IWYU_COMMAND NAMES include-what-you-use iwyu REQUIRED)
  set(CMAKE_CXX_INCLUDE_WHAT_YOU_USE
      "${IWYU_COMMAND};-Xiwyu;--mapping_file=${PROJECT_SOURCE_DIR}/.iwyu.imp")
endif()
# ------------------------------ subdirectories -------------------------------
add_subdirectory(src)
if(${PROJECT_NAME}_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()
if(${PROJECT_NAME}_BUILD_PYTHON_MODULE)
  unset(CMAKE_CXX_INCLUDE_WHAT_YOU_USE)
  add_subdirectory(python/Higgs)
endif()
if(${PROJECT_NAME}_BUILD_MATHEMATICA_INTERFACE)
  add_subdirectory(wstp)
endif()

# ---------------------------------- testing ----------------------------------
if(${PROJECT_NAME}_BUILD_TESTING)
  unset(CMAKE_CXX_INCLUDE_WHAT_YOU_USE)
  FetchContent_Declare(
    catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG devel)
  FetchContent_MakeAvailable(catch2)
  list(APPEND CMAKE_MODULE_PATH "${catch2_SOURCE_DIR}/extras")
  include(CTest)
  add_subdirectory(tests)
endif()

# ----------------------- alias for use as subdirectory -----------------------
add_library(HiggsTools::HiggsTools ALIAS HiggsTools)
# -------------------------- export and installation --------------------------
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  HiggsToolsConfigVersion.cmake
  VERSION ${PACKAGE_VERSION}
  COMPATIBILITY SameMajorVersion)

export(
  TARGETS HiggsTools
  NAMESPACE HiggsTools::
  FILE HiggsToolsConfig.cmake)
if(NOT SKBUILD)
  install(TARGETS HiggsTools EXPORT HiggsToolsConfig)
  install(
    DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
    PATTERN "*.hpp")
  install(
    EXPORT HiggsToolsConfig
    FILE HiggsToolsConfig.cmake
    NAMESPACE HiggsTools::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/HiggsTools)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/HiggsToolsConfigVersion.cmake
          DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/HiggsTools)
endif()
