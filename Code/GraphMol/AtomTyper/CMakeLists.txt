rdkit_library(AtomTyper
              atom_typer.cpp
                            atom_typer_query_reorder.cpp
              expression_builder.cpp
                            reaction_extractor.cpp
                            LINK_LIBRARIES GraphMol SmilesParse ChemTransforms ChemReactions)
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
              expression_builder.hpp
                            reaction_extractor.hpp
              DEST GraphMol/AtomTyper)

rdkit_catch_test(testAtomTyper unitTests.cpp
                 LINK_LIBRARIES AtomTyper SmilesParse ChemReactions GraphMol RDGeneral SubstructMatch)

if(RDK_BUILD_PYTHON_WRAPPERS)
    add_subdirectory(Wrap)
endif()
