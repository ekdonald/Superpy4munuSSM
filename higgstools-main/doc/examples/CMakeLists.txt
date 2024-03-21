add_library(doc-examples OBJECT addingParticles.cpp settingBRs.cpp
                                settingDecayWidths.cpp usingHB.cpp usingHS.cpp)
target_link_libraries(doc-examples PRIVATE HiggsTools)
