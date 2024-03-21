mathematica_wstp_add_executable(MHiggsTools MHiggsTools.tmpp MHiggsTools.cpp)
target_link_libraries(MHiggsTools PRIVATE HiggsPredictions HiggsBounds
                                          HiggsSignals HiggsUtilities)
mathematica_absolutize_library_dependencies(MHiggsTools)
target_include_directories(MHiggsTools PRIVATE ${Mathematica_WSTP_INCLUDE_DIR})
