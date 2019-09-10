#ifndef __CONFIG_H__
#define __CONFIG_H__

// This is a configuration template header file. The actual configuration is in the CMakeLists.txt files.
// Normally most of the configuration is done automatically by CMake or should be done using CMake.

#define VERSION_MAJOR 4
#define VERSION_MINOR 0
#define VERSION_PATCH 0
#define VERSION "4.0.0"

#define DATA_PATH "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/../Data"

#define SCANMAN_TEST_MECH_INPUT_PATH_CH4 "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Repository/examples/ScanMan/CH4/CH4.72.mech"
#define SCANMAN_TEST_MECH_INPUT_PATH "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/ScanManNew/sandiego20141004_mechCK.txt"
#define SCANMAN_TEST_THERM_INPUT_PATH "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/ScanManNew/sandiego20120907_therm.txt"
#define SCANMAN_TEST_TRANS_INPUT_PATH "/Users/ruixu/Documents/Flamelet/flamemaster/FlameMaster/Build/src/ScanManNew/sandiego20120907_trans.txt"

/* #undef COMP_MSVC */
#define COMP_CLANG
/* #undef COMP_INTEL */
/* #undef COMP_GNU */
/* #undef COMP_MINGW */
/* #undef COMP_UNSPEC */

#ifndef SYSTEM_WINDOWS
/* #undef SYSTEM_WINDOWS */
#endif

#ifndef SYSTEM_DARWIN
#define SYSTEM_DARWIN
#endif

#ifndef SYSTEM_LINUX
/* #undef SYSTEM_LINUX */
#endif

#ifndef SYSTEM_AIX
/* #undef SYSTEM_AIX */
#endif

#ifndef SYSTEM_CYGWIN
/* #undef SYSTEM_CYGWIN */
#endif

#define HAVE_UNISTD_H
/* #undef HAVE_DIRECT_H */
#define HAVE_ALLOCA_H
#define HAVE_STRING_H
#define HAVE_GETOPT_H
#define HAVE_STDLIB_H

/* #undef COMPILE_FORTRAN_SRC */
/* #undef OPTIMAPP */
/* #undef BILIN_OMEGA */
/* #undef NEWTON_PERFORMANCE */
/* #undef WRITE_LOCAL_LEWIS */
#define PREM_UPWIND_CONVEC

/* we always use the custom data type Double */
typedef double Double;

#undef NOXPRODRATES

#ifndef __MYIOCONST__
#define __MYIOCONST__
#define TAB "\t"
#define NEWL "\n"
#endif // __MYIOCONST__

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

typedef const char *const ConstString;
typedef char *String;
typedef char Flag;

#endif // __CONFIG_H_

