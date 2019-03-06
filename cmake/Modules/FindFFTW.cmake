set(FFTW_VERSION "3" CACHE INTEGER "Version of FFTW required")

if (FFTW_FIND_VERSION)
  set(FFTW_VERSION ${FFTW_FIND_VERSION} CACHE INTEGER "Version of FFTW required")
endif()

if (FFTW_VERSION EQUAL "2")
  find_library(FFTW_LIBRARIES fftw
    PATHS ${FFTW_DIR}
    PATH_SUFFIXES fftw/.libs/ lib
    )
  find_path(FFTW_INCLUDE_PATH fftw.h
    PATHS ${FFTW_DIR} $ENV{INCLUDE_PATH}
    PATH_SUFFIXES include fftw
    )
else()

  find_library(FFTW_LIBRARIES fftw3
    PATHS ENV LD_LIBRARY_PATH
  )
  find_library(FFTW_THREAD_LIBRARY fftw3_threads
    PATHS ENV LD_LIBRARY_PATH
    )   
  find_library(FFTW_OPENMP_LIBRARY fftw3_omp
    PATHS ENV LD_LIBRARY_PATH
    )
  find_path(FFTW_INCLUDE_PATH fftw3.h
    PATHS ENV INCLUDE_PATH
    PATH_SUFFIXES include fftw
    )   
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW DEFAULT_MSG
  FFTW_LIBRARIES FFTW_INCLUDE_PATH)


if(NOT FFTW_FOUND)
  set(FFTW_DIR "" CACHE PATH "Location of FFTW library.")
endif(NOT FFTW_FOUND)

mark_as_advanced(FFTW_LIBRARIES FFTW_OPENMP_LIBRARY FFTW_THREAD_LIBRARY FFTW_INCLUDE_PATH FFTW_VERSION)
