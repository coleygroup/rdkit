rdkit_library(AtomTyper
              atom_typer.cpp
              smarts_analyzer.cpp
              expression_builder.cpp
              parser.cpp
              LINK_LIBRARIES GraphMol SmilesParse)
target_compile_definitions(AtomTyper PRIVATE RDKIT_ATOMTYPER_BUILD)

rdkit_headers(atom_typer.hpp
              smarts_analyzer.hpp
              expression_builder.hpp
              parser.hpp
              DEST GraphMol/AtomTyper)

rdkit_catch_test(testAtomTyper test_atom_typer.cpp
                 LINK_LIBRARIES AtomTyper SmilesParse)

rdkit_catch_test(testSmartsAnalyzer test_smarts_analyzer.cpp
                 LINK_LIBRARIES AtomTyper SmilesParse)

rdkit_catch_test(testExpressionBuilder test_expression_builder.cpp
                 LINK_LIBRARIES AtomTyper SmilesParse)

add_executable(smartsAnalyzerSmoke smarts_analyzer_smoke.cpp)
target_link_libraries(smartsAnalyzerSmoke PRIVATE AtomTyper SmilesParse)

if(RDK_BUILD_PYTHON_WRAPPERS)
    add_subdirectory(Wrap)
endif()
