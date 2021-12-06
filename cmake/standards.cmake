option(STANDARDS_INSTALL OFF)
option(STANDARDS_BUILD_TESTS OFF)
option(STANDARDS_BUILD_EXAMPLES OFF)

include(FetchContent)
FetchContent_Declare(
  standards
  GIT_REPOSITORY git@github.com:todo-group/standards.git
)
