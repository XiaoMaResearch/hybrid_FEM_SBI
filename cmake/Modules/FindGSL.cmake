find_path(GSL_INCLUDE_PATH gsl_math.h
  PATHS ${GSL_DIR} ENV C_INCLUDE_PATH
  PATH_SUFFIXES gsl
  )

find_library(GSL_MAIN_LIBRARY NAME gsl
  PATHS ${GSL_DIR} ENV LIBRARY_PATH
  PATH_SUFFIXES lib
  )

find_library(GSL_BLAS_LIBRARY NAME gslcblas
  PATHS ${GSL_DIR} ENV LIBRARY_PATH
  PATH_SUFFIXES lib
)

mark_as_advanced(GSL_INCLUDE_PATH)
mark_as_advanced(GSL_MAIN_LIBRARY NAME)
mark_as_advanced(GSL_BLAS_LIBRARY NAME)

set(GSL_LIBRARIES ${GSL_MAIN_LIBRARY} ${GSL_BLAS_LIBRARY})


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GSL DEFAULT_MSG
  GSL_LIBRARIES GSL_INCLUDE_PATH)
