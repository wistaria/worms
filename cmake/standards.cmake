option(STANDARDS_INSTALL OFF)
option(STANDARDS_BUILD_TESTS OFF)
option(STANDARDS_BUILD_EXAMPLES OFF)

include(FetchContent)
FetchContent_Declare(
  standards
  GIT_REPOSITORY https://github.com/todo-group/standards.git
)
list(APPEND FetchContent_includes "${PROJECT_BINARY_DIR}/_deps/standards-src")
list(APPEND FetchContents standards)
