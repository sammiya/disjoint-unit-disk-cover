cmake_minimum_required(VERSION 3.25)
if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)
endif()
project(disjoint_unit_disk_cpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(DISJOINT_DISK_ENABLE_CLANG_TIDY "Run clang-tidy during builds" OFF)
option(DISJOINT_DISK_ENABLE_IWYU "Run include-what-you-use as a custom target" OFF)
option(DISJOINT_DISK_BUILD_TESTS "Build C++ tests" ON)

set(THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third_party)
find_package(Threads REQUIRED)
set(CGAL_DIR ${THIRD_PARTY_DIR}/cgal CACHE PATH "Path to vendored CGAL" FORCE)
list(APPEND CMAKE_PREFIX_PATH ${CGAL_DIR})
list(APPEND CMAKE_MODULE_PATH ${CGAL_DIR}/cmake ${CGAL_DIR}/scripts/cmake/modules ${CGAL_DIR}/Installation/cmake/modules)
find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)
find_package(CGAL 6.0 CONFIG REQUIRED COMPONENTS Core)
find_package(Threads REQUIRED)

find_program(CLANG_TIDY_EXE NAMES clang-tidy)
if(DISJOINT_DISK_ENABLE_CLANG_TIDY AND CLANG_TIDY_EXE)
    set(DISJOINT_DISK_CLANG_TIDY_ARGS
        "${CLANG_TIDY_EXE};-warnings-as-errors=*;-extra-arg=-D_LIBCPP_REMOVE_TRANSITIVE_INCLUDES=1;-extra-arg=-Wno-unknown-warning-option"
    )
elseif(DISJOINT_DISK_ENABLE_CLANG_TIDY AND NOT CLANG_TIDY_EXE)
    message(WARNING "DISJOINT_DISK_ENABLE_CLANG_TIDY=ON but clang-tidy not found")
endif()

find_program(IWYU_EXE NAMES include-what-you-use iwyu)
find_program(IWYU_TOOL NAMES iwyu_tool.py iwyu_tool)
find_package(Python3 COMPONENTS Interpreter)

add_library(disjoint_unit_disk_cover_core
    src/enumeration.cpp
    src/search.cpp
    src/workflow.cpp
    src/logging.cpp
    src/cli_io.cpp
)

target_include_directories(disjoint_unit_disk_cover_core
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${THIRD_PARTY_DIR}/toml++/include
)

target_compile_options(disjoint_unit_disk_cover_core PRIVATE -Wall -Wextra -Wpedantic)
target_link_libraries(disjoint_unit_disk_cover_core PUBLIC CGAL::CGAL)

add_executable(disjoint_unit_disk_cover_cli
    src/main.cpp
)

target_include_directories(disjoint_unit_disk_cover_cli
    PRIVATE
        ${THIRD_PARTY_DIR}/CLI11/include
)

target_link_libraries(disjoint_unit_disk_cover_cli
    PRIVATE
        disjoint_unit_disk_cover_core
        Threads::Threads
)
if(DISJOINT_DISK_ENABLE_CLANG_TIDY AND DISJOINT_DISK_CLANG_TIDY_ARGS)
    set_target_properties(disjoint_unit_disk_cover_core disjoint_unit_disk_cover_cli
        PROPERTIES CXX_CLANG_TIDY "${DISJOINT_DISK_CLANG_TIDY_ARGS}")
endif()

target_compile_options(disjoint_unit_disk_cover_cli PRIVATE -Wall -Wextra -Wpedantic)

set_target_properties(disjoint_unit_disk_cover_cli PROPERTIES OUTPUT_NAME disjoint_unit_disk_cover)

find_program(CLANG_FORMAT_EXE NAMES clang-format)
if(CLANG_FORMAT_EXE)
    file(GLOB_RECURSE FORMAT_CXX_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp
    )
    add_custom_target(format
        COMMAND ${CLANG_FORMAT_EXE} -i ${FORMAT_CXX_SOURCES}
        COMMENT "Running clang-format"
        VERBATIM
    )
endif()

file(GLOB_RECURSE PUBLIC_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h
)
add_custom_target(check-headers
    COMMENT "Checking public headers are self-contained"
)
foreach(_hdr IN LISTS PUBLIC_HEADERS)
    add_custom_command(TARGET check-headers
        POST_BUILD
        COMMAND ${CMAKE_CXX_COMPILER}
                -std=c++20
                -fsyntax-only
                -D_LIBCPP_REMOVE_TRANSITIVE_INCLUDES=1
                -I${CMAKE_CURRENT_SOURCE_DIR}/include
                -isystem ${THIRD_PARTY_DIR}/toml++/include
                -include ${_hdr}
                -xc++
                /dev/null
        COMMENT "Checking header ${_hdr}"
        VERBATIM
    )
endforeach()

if(DISJOINT_DISK_ENABLE_IWYU)
    if(IWYU_EXE AND IWYU_TOOL AND Python3_Interpreter_FOUND)
        add_custom_target(iwyu
            COMMAND ${Python3_EXECUTABLE} ${IWYU_TOOL} -p ${CMAKE_BINARY_DIR} -e ${THIRD_PARTY_DIR} --
                    -Xiwyu --no_fwd_decls
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            COMMENT "Running include-what-you-use (iwyu) using compile_commands.json"
            VERBATIM
        )
    else()
        message(WARNING "DISJOINT_DISK_ENABLE_IWYU=ON but include-what-you-use or Python3 was not found")
    endif()
endif()

if(DISJOINT_DISK_BUILD_TESTS)
    include(CTest)
    enable_testing()

    set(GTEST_ROOT ${THIRD_PARTY_DIR}/googletest CACHE PATH "")

    if(NOT EXISTS ${GTEST_ROOT})
        message(FATAL_ERROR "googletest submodule not found at ${GTEST_ROOT}. Please add it as a submodule.")
    endif()

    add_subdirectory(${GTEST_ROOT} ${CMAKE_BINARY_DIR}/googletest)

    add_executable(disjoint_unit_disk_cover_core_test
        tests/test_primitives.cpp
        tests/test_enumeration.cpp
        tests/test_search.cpp
        tests/test_workflow.cpp
        tests/test_cli_io.cpp
    )

    target_link_libraries(disjoint_unit_disk_cover_core_test
        PRIVATE
            disjoint_unit_disk_cover_core
            gtest_main
    )

    include(GoogleTest)
    gtest_discover_tests(disjoint_unit_disk_cover_core_test)
endif()

add_custom_target(clang-tidy
    COMMENT "Run clang-tidy across build targets (enable with -DDISJOINT_DISK_ENABLE_CLANG_TIDY=ON)"
)
add_dependencies(clang-tidy disjoint_unit_disk_cover_core disjoint_unit_disk_cover_cli)
