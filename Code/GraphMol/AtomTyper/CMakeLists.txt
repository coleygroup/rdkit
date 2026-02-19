rdkit_library(AtomTyper
              atom_typer.cpp
                            atom_typer_query_reorder.cpp
              smarts_analyzer.cpp
              expression_builder.cpp
            #   parser.cpp
              LINK_LIBRARIES GraphMol SmilesParse ChemTransforms)
target_compile_definitions(AtomTyper PRIVATE RDKIT_ATOMTYPER_BUILD)

option(RDK_BUILD_TRACY "Enable Tracy profiling for AtomTyper targets" OFF)
if(RDK_BUILD_TRACY)
    set(_TRACY_SOURCE_DIR "${CMAKE_SOURCE_DIR}/../tracy")
    if(EXISTS "${_TRACY_SOURCE_DIR}/CMakeLists.txt")
        add_subdirectory("${_TRACY_SOURCE_DIR}" "${CMAKE_BINARY_DIR}/tracy" EXCLUDE_FROM_ALL)
        target_link_libraries(AtomTyper PUBLIC TracyClient)
        target_compile_definitions(AtomTyper PUBLIC TRACY_ENABLE)
    else()
        message(WARNING "RDK_BUILD_TRACY=ON but Tracy source not found at ${_TRACY_SOURCE_DIR}")
    endif()
endif()

rdkit_headers(atom_typer.hpp
                            atom_typer_query_reorder.hpp
              smarts_analyzer.hpp
              expression_builder.hpp
            #   parser.hpp
              DEST GraphMol/AtomTyper)

rdkit_catch_test(testAtomTyper test_atom_typer.cpp
                 LINK_LIBRARIES AtomTyper SmilesParse)

rdkit_catch_test(testSmartsAnalyzer test_smarts_analyzer.cpp
                 LINK_LIBRARIES AtomTyper SmilesParse)

rdkit_catch_test(testExpressionBuilder test_expression_builder.cpp
                 LINK_LIBRARIES AtomTyper SmilesParse)

option(RDK_BUILD_ATOMTYPER_SMOKE "Build AtomTyper smoke executables" ON)
if(RDK_BUILD_ATOMTYPER_SMOKE)
    add_executable(smartsAnalyzerSmoke smarts_analyzer_smoke.cpp)
    target_link_libraries(smartsAnalyzerSmoke PRIVATE AtomTyper SmilesParse SubstructMatch)
    if(RDK_BUILD_TRACY)
        target_link_libraries(smartsAnalyzerSmoke PRIVATE TracyClient)
    endif()

    add_executable(smartsAnalyzerMatchsetSmoke smarts_analyzer_matchset_smoke.cpp)
    target_link_libraries(smartsAnalyzerMatchsetSmoke PRIVATE AtomTyper SmilesParse SubstructMatch)
    if(RDK_BUILD_TRACY)
        target_link_libraries(smartsAnalyzerMatchsetSmoke PRIVATE TracyClient)
    endif()

    add_executable(smartsAnalyzerMechsSmoke smarts_analyzer_mechs_smoke.cpp)
    target_link_libraries(smartsAnalyzerMechsSmoke PRIVATE AtomTyper SmilesParse SubstructMatch)
    if(RDK_BUILD_TRACY)
        target_link_libraries(smartsAnalyzerMechsSmoke PRIVATE TracyClient)
    endif()

    add_executable(smartsAnalyzerEquivalenceTest unitTests.cpp)
    target_link_libraries(smartsAnalyzerEquivalenceTest PRIVATE AtomTyper SmilesParse GraphMol RDGeneral SubstructMatch)

    # add_executable(smartsAnalyzerPrimitivesSmoke smarts_analyzer_primitives_smoke.cpp)
    # target_link_libraries(smartsAnalyzerPrimitivesSmoke PRIVATE AtomTyper SmilesParse)
    # if(RDK_BUILD_TRACY)
    #     target_link_libraries(smartsAnalyzerPrimitivesSmoke PRIVATE TracyClient)
    # endif()
endif()

if(RDK_BUILD_PYTHON_WRAPPERS)
    add_subdirectory(Wrap)
endif()
