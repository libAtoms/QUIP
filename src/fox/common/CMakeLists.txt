set(fox_common_srcs
  FoX_common.F90
)

set(fox_common_srcs_pp
  m_common_attrs.F90
  m_common_buffer.F90
  m_common_charset.F90
  m_common_content_model.F90
  m_common_element.F90
  m_common_elstack.F90
  m_common_entities.F90
  m_common_entity_expand.F90
  m_common_error.F90
  m_common_io.F90
  m_common_namecheck.F90
  m_common_namespaces.F90
  m_common_notations.F90
  m_common_struct.F90
)

add_library(
  ${PROJECT_NAME}_common STATIC
  ${fox_common_srcs}
  ${fox_common_srcs_pp}
)

target_link_libraries(
  ${PROJECT_NAME}_common
  ${PROJECT_NAME}_utils
  ${PROJECT_NAME}_fsys
)

