rdkit_library(AtomTyper
              atom_typer.cpp
              smarts_analyzer.cpp
              expression_builder.cpp
              LINK_LIBRARIES GraphMol)
target_compile_definitions(AtomTyper PRIVATE RDKIT_ATOMTYPER_BUILD)

rdkit_headers(atom_typer.hpp
              smarts_analyzer.hpp
              expression_builder.hpp
              DEST GraphMol/AtomTyper)

rdkit_test(testAtomTyper test_atom_typer.cpp
           LINK_LIBRARIES AtomTyper SmilesParse)

rdkit_test(testSmartsAnalyzer test_smarts_analyzer.cpp
           LINK_LIBRARIES AtomTyper SmilesParse)

rdkit_test(testExpressionBuilder test_expression_builder.cpp
           LINK_LIBRARIES AtomTyper SmilesParse)

if(RDK_BUILD_PYTHON_WRAPPERS)
    add_subdirectory(Wrap)
endif()
