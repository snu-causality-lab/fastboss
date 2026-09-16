cmake_minimum_required(VERSION 3.16)
project(fastboss_core VERSION 0.1.0 LANGUAGES C CXX)

option(FASTBOSS_BUILD_SHARED "Build shared library" ON)
option(FASTBOSS_BUILD_EXAMPLES "Build examples" ON)
option(FASTBOSS_BUILD_TESTS "Build and register tests" ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN YES)

add_library(fastboss_core_obj OBJECT src/fastboss_core.cpp)
target_include_directories(fastboss_core_obj PUBLIC include)
target_compile_features(fastboss_core_obj PUBLIC cxx_std_17)
set_target_properties(fastboss_core_obj PROPERTIES POSITION_INDEPENDENT_CODE ON)

if(FASTBOSS_BUILD_SHARED)
  target_compile_definitions(fastboss_core_obj PRIVATE FASTBOSS_BUILD_SHARED)
  add_library(fastboss SHARED $<TARGET_OBJECTS:fastboss_core_obj>)
  target_include_directories(fastboss PUBLIC include)
  set_target_properties(fastboss PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 0)
else()
  target_compile_definitions(fastboss_core_obj PRIVATE FASTBOSS_STATIC)
  add_library(fastboss STATIC $<TARGET_OBJECTS:fastboss_core_obj>)
  target_include_directories(fastboss PUBLIC include)
  target_compile_definitions(fastboss PUBLIC FASTBOSS_STATIC)
endif()

install(TARGETS fastboss
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin)
install(FILES include/fastboss.h DESTINATION include)

if(FASTBOSS_BUILD_EXAMPLES)
  add_executable(fastboss_smoke examples/smoke.c)
  target_link_libraries(fastboss_smoke PRIVATE fastboss)
  set_target_properties(fastboss_smoke PROPERTIES LINKER_LANGUAGE CXX)

  add_executable(fastboss_profile examples/profile_boss.cpp)
  target_link_libraries(fastboss_profile PRIVATE fastboss)
endif()

if(FASTBOSS_BUILD_TESTS)
  enable_testing()
  add_executable(fastboss_core_tests tests/core_tests.c)
  target_link_libraries(fastboss_core_tests PRIVATE fastboss)
  set_target_properties(fastboss_core_tests PROPERTIES LINKER_LANGUAGE CXX)
  add_test(NAME fastboss_core_tests COMMAND fastboss_core_tests)

  add_executable(fastboss_determinism_tests tests/determinism_tests.cpp)
  target_include_directories(fastboss_determinism_tests PRIVATE src)
  add_test(NAME fastboss_determinism_tests COMMAND fastboss_determinism_tests)

  if(FASTBOSS_BUILD_EXAMPLES)
    add_test(NAME fastboss_smoke COMMAND fastboss_smoke)
  endif()
endif()
