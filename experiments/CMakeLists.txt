include(FetchContent)
include(CMakePrintHelpers)

add_library(experiments INTERFACE)
target_include_directories(experiments INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(experiments INTERFACE fmt json)
if (ENABLE_MATPLOTLIB)
  target_link_libraries(experiments INTERFACE matplot)
endif()

# check for git revision
if(EXISTS ${PROJECT_SOURCE_DIR}/.git)
  find_package(Git)
  if(GIT_FOUND)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
      WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
      OUTPUT_VARIABLE "GIT_SHORT_REVISION"
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    target_compile_definitions(experiments INTERFACE "GIT_SHORT_REVISION=\"${GIT_SHORT_REVISION}\"")
  endif()
endif()

# path to the experiments cpp files and EPFL benchmarks
target_compile_definitions(experiments INTERFACE "EXPERIMENTS_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}/\"")

# get the mtKaHyPar lib as a library
FetchContent_Declare(
  MtKaHyPar EXCLUDE_FROM_ALL
  GIT_REPOSITORY https://github.com/kahypar/mt-kahypar
  GIT_TAG        v1.5.2
)
FetchContent_MakeAvailable(MtKaHyPar)

file(GLOB FILENAMES *.cpp)
set(EXPFILE ${FILENAMES})

# TODO: @Jingren Filter out the target file
foreach(filename ${EXPFILE})
  get_filename_component(basename ${filename} NAME_WE)
  add_executable(${basename} ${filename})
  target_link_libraries(${basename} PUBLIC experiments mockturtle MtKaHyPar::mtkahypar)
endforeach()
