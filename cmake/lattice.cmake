option(LATTICE_INSTALL OFF)
option(LATTICE_BUILD_TESTS OFF)
option(LATTICE_BUILD_EXAMPLES OFF)

include(FetchContent)
FetchContent_Declare(
  lattice
  GIT_REPOSITORY https://github.com/todo-group/lattice.git
)
list(APPEND FetchContent_includes "${PROJECT_BINARY_DIR}/_deps/lattice-src")
list(APPEND FetchContents lattice)
