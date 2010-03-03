find_path(ANN_INCLUDE_DIRS ANN.h  "$ENV{ANN_ROOT}/include/ANN/")
find_library(ANN_LIBRARIES ANN  "$ENV{ANN_ROOT}/lib")


MESSAGE(STATUS "${ANN_INCLUDE_DIRS}")
MESSAGE(STATUS "${ANN_LIBRARIES}")
set(ANN_FOUND TRUE)

if (NOT ANN_INCLUDE_DIRS)
set(ANN_FOUND FALSE)
else(NOT ANN_INCLUDE_DIRS)
set(ANN_INCLUDE_DIRS "$ENV{ANN_ROOT}/include")
endif (NOT ANN_INCLUDE_DIRS)
if (NOT ANN_LIBRARIES)
set(ANN_FOUND FALSE)
endif (NOT ANN_LIBRARIES)
# ==========================================
IF(NOT ANN_FOUND)
  # make FIND_PACKAGE friendly
  IF(NOT ANN_FIND_QUIETLY)
    IF(ANN_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "ANN required, please specify it's location by setting ENV variable ANN_ROOT.")
    ELSE(ANN_FIND_REQUIRED)
      MESSAGE(STATUS       "ERROR: ANN was not found  please specify it's location by setting ENV variable ANN_ROOT.")
    ENDIF(ANN_FIND_REQUIRED)
  ENDIF(NOT ANN_FIND_QUIETLY)
ENDIF(NOT ANN_FOUND)

