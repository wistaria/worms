option(BCL_INSTALL OFF)
option(BCL_BUILD_TESTS OFF)
option(BCL_BUILD_EXAMPLES OFF)

include(FetchContent)
FetchContent_Declare(
  bcl
  GIT_REPOSITORY https://github.com/cmsi/bcl.git
)
