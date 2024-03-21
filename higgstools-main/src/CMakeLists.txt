file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/empty.cpp" )
if(${PROJECT_NAME}_BUILD_SHARED_LIBS)
  add_library(HiggsTools SHARED ${CMAKE_CURRENT_BINARY_DIR}/empty.cpp)
else()
  add_library(HiggsTools STATIC ${CMAKE_CURRENT_BINARY_DIR}/empty.cpp)
endif()
generate_export_header(HiggsTools EXPORT_FILE_NAME
                       ${PROJECT_BINARY_DIR}/include/Higgs/HiggsTools_export.h)

add_subdirectory(utilities)
add_subdirectory(predictions)
add_subdirectory(bounds)
add_subdirectory(signals)

target_link_libraries(HiggsTools PRIVATE HiggsBounds HiggsSignals
                                         HiggsPredictions HiggsUtilities)
target_include_directories(
  HiggsTools
  PUBLIC $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
         $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
         $<INSTALL_INTERFACE:include>)
target_compile_features(HiggsUtilities PUBLIC cxx_std_17)

if(NOT SKBUILD)
  install(FILES ${PROJECT_BINARY_DIR}/include/Higgs/HiggsTools_export.h
          DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Higgs)
endif()
