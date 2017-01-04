#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([Allink], [1.0], [sabeiro@virgilio.it])
AC_CONFIG_SRCDIR([include/MatematicaStruct.h])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_AWK
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lCGAL':
AC_CHECK_LIB([CGAL], [main])
# FIXME: Replace `main' with a function in `-lDraw':
AC_CHECK_LIB([Draw], [main])
# FIXME: Replace `main' with a function in `-lGL':
AC_CHECK_LIB([GL], [main])
# FIXME: Replace `main' with a function in `-lGLU':
AC_CHECK_LIB([GLU], [main])
# FIXME: Replace `main' with a function in `-lMatematica':
AC_CHECK_LIB([Matematica], [main])
# FIXME: Replace `main' with a function in `-lQt3Support':
AC_CHECK_LIB([Qt3Support], [main])
# FIXME: Replace `main' with a function in `-lQtCore':
AC_CHECK_LIB([QtCore], [main])
# FIXME: Replace `main' with a function in `-lQtGui':
AC_CHECK_LIB([QtGui], [main])
# FIXME: Replace `main' with a function in `-lQtNetwork':
AC_CHECK_LIB([QtNetwork], [main])
# FIXME: Replace `main' with a function in `-lQtOpenGL':
AC_CHECK_LIB([QtOpenGL], [main])
# FIXME: Replace `main' with a function in `-lQtSql':
AC_CHECK_LIB([QtSql], [main])
# FIXME: Replace `main' with a function in `-lQtXml':
AC_CHECK_LIB([QtXml], [main])
# FIXME: Replace `main' with a function in `-lVarData':
AC_CHECK_LIB([VarData], [main])
# FIXME: Replace `main' with a function in `-lX11':
AC_CHECK_LIB([X11], [main])
# FIXME: Replace `main' with a function in `-lXi':
AC_CHECK_LIB([Xi], [main])
# FIXME: Replace `main' with a function in `-lXmu':
AC_CHECK_LIB([Xmu], [main])
# FIXME: Replace `main' with a function in `-lXtst':
AC_CHECK_LIB([Xtst], [main])
# FIXME: Replace `main' with a function in `-lconfuse':
AC_CHECK_LIB([confuse], [main])
# FIXME: Replace `main' with a function in `-lfftw3':
AC_CHECK_LIB([fftw3], [main])
# FIXME: Replace `main' with a function in `-lfreetype':
AC_CHECK_LIB([freetype], [main])
# FIXME: Replace `main' with a function in `-lglut':
AC_CHECK_LIB([glut], [main])
# FIXME: Replace `main' with a function in `-lgmp':
AC_CHECK_LIB([gmp], [main])
# FIXME: Replace `main' with a function in `-lgsl':
AC_CHECK_LIB([gsl], [main])
# FIXME: Replace `main' with a function in `-lgslcblas':
AC_CHECK_LIB([gslcblas], [main])
# FIXME: Replace `main' with a function in `-lguide':
AC_CHECK_LIB([guide], [main])
# FIXME: Replace `main' with a function in `-liomp5':
AC_CHECK_LIB([iomp5], [main])
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [main])
# FIXME: Replace `main' with a function in `-lmkl':
AC_CHECK_LIB([mkl], [main])
# FIXME: Replace `main' with a function in `-lmkl_core':
AC_CHECK_LIB([mkl_core], [main])
# FIXME: Replace `main' with a function in `-lmkl_intel_lp64':
AC_CHECK_LIB([mkl_intel_lp64], [main])
# FIXME: Replace `main' with a function in `-lmkl_sequential':
AC_CHECK_LIB([mkl_sequential], [main])
# FIXME: Replace `main' with a function in `-lmkl_solver':
AC_CHECK_LIB([mkl_solver], [main])
# FIXME: Replace `main' with a function in `-lmpfr':
AC_CHECK_LIB([mpfr], [main])
# FIXME: Replace `main' with a function in `-lpng':
AC_CHECK_LIB([png], [main])
# FIXME: Replace `main' with a function in `-lpngwriter':
AC_CHECK_LIB([pngwriter], [main])
# FIXME: Replace `main' with a function in `-lpthread':
AC_CHECK_LIB([pthread], [main])
# FIXME: Replace `main' with a function in `-ltiff':
AC_CHECK_LIB([tiff], [main])
# FIXME: Replace `main' with a function in `-lxosd':
AC_CHECK_LIB([xosd], [main])
# FIXME: Replace `main' with a function in `-lz':
AC_CHECK_LIB([z], [main])

# Checks for header files.
AC_CHECK_HEADERS([fenv.h limits.h stdlib.h string.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_C_RESTRICT
AC_TYPE_SIZE_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_ALLOCA
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([atexit clock_gettime floor memmove memset pow sqrt strcasecmp strchr strcspn strerror strpbrk strstr])

AC_CONFIG_FILES([Debug/makefile
                 Macro/openFrameworks/billboardAudio/Makefile
                 Macro/openFrameworks/continuumPiano/Makefile
                 Macro/openFrameworks/detection/Makefile
                 Macro/openFrameworks/ofxMSAPhysics/example/Makefile
                 Makefile
                 src/Addraw/Makefile
                 src/Addyn/Makefile
                 src/DrEffect/Makefile
                 src/ElPoly/Makefile
                 src/KinMC/Makefile
                 src/Makefile
                 src/Matematica/Makefile
                 src/PathOpt/Makefile
                 src/Segnali/Makefile
                 src/VarData/Makefile
                 src/fdmdpd2/src/autocorr/Makefile
                 src/fdmdpd2/src/corr/Makefile
                 src/fdmdpd2/src/dwham/Makefile
                 src/fdmdpd2/src/fdmdpd2/Makefile
                 src/fdmdpd2/src/forcematch/Makefile
                 src/fdmdpd2/src/ideal-chain/Makefile
                 src/fdmdpd2/src/imf-ev/Makefile
                 src/fdmdpd2/src/line-tension/Makefile
                 src/fdmdpd2/src/normpot/Makefile
                 src/fdmdpd2/src/undu/Makefile
                 src/fdmdpd2/src/viscosity/Makefile])
AC_CONFIG_SUBDIRS([src/Avvis])
AC_OUTPUT
