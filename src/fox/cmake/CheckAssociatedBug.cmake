#=====================================================#
# check for 'associated in restricted expression' bug #
#=====================================================#

MESSAGE(STATUS "Checking for 'associated in restricted expression' bug")
TRY_COMPILE(TEST_ASSOCIATED ${${PROJECT_NAME}_BINARY_DIR}/associated ${${PROJECT_NAME}_SOURCE_DIR}/cmake
  fox_config expr_bug
  OUTPUT_VARIABLE BUILD_OUTPUT 
)

IF(${TEST_ASSOCIATED} MATCHES FALSE)
  MESSAGE("   -> yes")
  ADD_DEFINITIONS(-DRESTRICTED_ASSOCIATED_BUG)
ENDIF(${TEST_ASSOCIATED} MATCHES FALSE)


