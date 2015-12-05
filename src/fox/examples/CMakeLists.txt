add_executable(
  wxml_example
  wxml_example.f90)
target_link_libraries(
  wxml_example
  ${PROJECT_NAME}_wxml)
ADD_TEST(wxml_example_TEST ${EXECUTABLE_OUTPUT_PATH}/wxml_example)


add_executable(
  wcml_example
  wcml_example.f90)
target_link_libraries(
  wcml_example
  ${PROJECT_NAME}_wcml)
ADD_TEST(wcml_example_TEST ${EXECUTABLE_OUTPUT_PATH}/wcml_example)

add_executable(
  wkml_example
  wkml_example.f90)
target_link_libraries(
  wkml_example
  ${PROJECT_NAME}_wkml)
ADD_TEST(wkml_example_TEST ${EXECUTABLE_OUTPUT_PATH}/wkml_example)

add_executable(
  sax_example
  sax_example.f90)
target_link_libraries(
  sax_example
  ${PROJECT_NAME}_sax)
ADD_CUSTOM_COMMAND(TARGET sax_example PRE_BUILD COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/examples/staffNS.xml ${PROJECT_BINARY_DIR}/examples)
ADD_TEST(sax_example_TEST ${EXECUTABLE_OUTPUT_PATH}/sax_example)

add_executable(
  sax_example_2
  sax_example_2.f90)
target_link_libraries(
  sax_example_2
  ${PROJECT_NAME}_sax)
ADD_TEST(sax_example_2_TEST ${EXECUTABLE_OUTPUT_PATH}/sax_example_2)

add_executable(
  dom_example_2
  dom_example_2.f90)
target_link_libraries(
  dom_example_2
  ${PROJECT_NAME}_dom)
ADD_CUSTOM_COMMAND(TARGET dom_example_2 PRE_BUILD COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/examples/h2o.xml ${PROJECT_BINARY_DIR}/examples)
ADD_TEST(dom_example_2_TEST ${EXECUTABLE_OUTPUT_PATH}/dom_example_2)

add_executable(
  dom_example_3
  dom_example_3.f90)
target_link_libraries(
  dom_example_3
  ${PROJECT_NAME}_dom)
ADD_CUSTOM_COMMAND(TARGET dom_example_3 PRE_BUILD COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/examples/h2o.xml ${PROJECT_BINARY_DIR}/examples)
ADD_TEST(dom_example_2_TEST ${EXECUTABLE_OUTPUT_PATH}/dom_example_2)


include_directories(
  ${${PROJECT_NAME}_BINARY_DIR}/dom
  ${${PROJECT_NAME}_BINARY_DIR}/wxml
  ${${PROJECT_NAME}_BINARY_DIR}/wcml
  ${${PROJECT_NAME}_BINARY_DIR}/wkml
  ${${PROJECT_NAME}_BINARY_DIR}/common
  ${${PROJECT_NAME}_BINARY_DIR}/utils
  ${${PROJECT_NAME}_BINARY_DIR}/fsys
  ${${PROJECT_NAME}_BINARY_DIR}/sax
)



