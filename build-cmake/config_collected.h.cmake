/* config.h.  Generated from config_collected.h.cmake by CMake.
   It was generated from config_collected.h.cmake which in turn is generated automatically
   from the config.h.cmake files of modules this module depends on. */

/* Define to 1 if you have module dune-PFASST available */
#cmakedefine01 HAVE_DUNE_PFASST


/* Define to 1 if you have module dune-common available */
#cmakedefine01 HAVE_DUNE_COMMON


/* Define to 1 if you have module dune-uggrid available */
#cmakedefine01 HAVE_DUNE_UGGRID


/* Define to 1 if you have module dune-geometry available */
#cmakedefine01 HAVE_DUNE_GEOMETRY


/* Define to 1 if you have module dune-istl available */
#cmakedefine01 HAVE_DUNE_ISTL


/* Define to 1 if you have module dune-typetree available */
#cmakedefine01 HAVE_DUNE_TYPETREE


/* Define to 1 if you have module dune-grid available */
#cmakedefine01 HAVE_DUNE_GRID


/* Define to 1 if you have module dune-localfunctions available */
#cmakedefine01 HAVE_DUNE_LOCALFUNCTIONS


/* Define to 1 if you have module dune-alugrid available */
#cmakedefine01 HAVE_DUNE_ALUGRID


/* Define to 1 if you have module dune-matrix-vector available */
#cmakedefine01 HAVE_DUNE_MATRIX_VECTOR


/* Define to 1 if you have module dune-functions available */
#cmakedefine01 HAVE_DUNE_FUNCTIONS


/* Define to 1 if you have module dune-fufem available */
#cmakedefine01 HAVE_DUNE_FUFEM


/* Define to the version of dune-common */
#define DUNE_COMMON_VERSION "${DUNE_COMMON_VERSION}"

/* Define to the major version of dune-common */
#define DUNE_COMMON_VERSION_MAJOR ${DUNE_COMMON_VERSION_MAJOR}

/* Define to the minor version of dune-common */
#define DUNE_COMMON_VERSION_MINOR ${DUNE_COMMON_VERSION_MINOR}

/* Define to the revision of dune-common */
#define DUNE_COMMON_VERSION_REVISION ${DUNE_COMMON_VERSION_REVISION}

/* Standard debug streams with a level below will collapse to doing nothing */
#define DUNE_MINIMAL_DEBUG_LEVEL ${DUNE_MINIMAL_DEBUG_LEVEL}

/* does the compiler support __attribute__((deprecated))? */
#cmakedefine HAS_ATTRIBUTE_DEPRECATED 1

/* does the compiler support __attribute__((deprecated("message"))? */
#cmakedefine HAS_ATTRIBUTE_DEPRECATED_MSG 1

/* does the compiler support __attribute__((unused))? */
#cmakedefine HAS_ATTRIBUTE_UNUSED 1

/* Define if you have a BLAS library. */
#cmakedefine HAVE_BLAS 1

/* does the compiler support abi::__cxa_demangle */
#cmakedefine HAVE_CXA_DEMANGLE 1

/* Define if you have LAPACK library. */
#cmakedefine HAVE_LAPACK 1

/* Define to 1 if you have the <malloc.h> header file. */
// Not used! #cmakedefine01 HAVE_MALLOC_H

/* Define if you have the MPI library.  */
#cmakedefine HAVE_MPI ENABLE_MPI

/* Define if you have the GNU GMP library. The value should be ENABLE_GMP
   to facilitate activating and deactivating GMP using compile flags. */
#cmakedefine HAVE_GMP ENABLE_GMP

/* Define if you have the Vc library. The value should be ENABLE_VC
   to facilitate activating and deactivating Vc using compile flags. */
#cmakedefine HAVE_VC ENABLE_VC

/* Define to 1 if you have the symbol mprotect. */
#cmakedefine HAVE_MPROTECT 1

/* Define to 1 if you have the <stdint.h> header file. */
#cmakedefine HAVE_STDINT_H 1

/* Define to 1 if you have <sys/mman.h>. */
#cmakedefine HAVE_SYS_MMAN_H 1



/* old feature support macros which were tested until 2.4, kept around for one more release */
/* As these are now always supported due to the new compiler requirements, they are directly */
/* defined without an explicit test. */
#define HAVE_NULLPTR 1
#define HAVE_CONSTEXPR 1
#define HAVE_RANGE_BASED_FOR 1
#define HAVE_NOEXCEPT_SPECIFIER 1
#define HAVE_STD_DECLVAL 1
#define HAVE_KEYWORD_FINAL 1
#define MPI_2 1

/* Define to 1 if the compiler properly supports testing for operator[] */
#cmakedefine HAVE_IS_INDEXABLE_SUPPORT 1

/* Define to ENABLE_UMFPACK if the UMFPack library is available */
#cmakedefine HAVE_UMFPACK ENABLE_SUITESPARSE

/* Define to ENABLE_SUITESPARSE if the SuiteSparse library is available */
#cmakedefine HAVE_SUITESPARSE ENABLE_SUITESPARSE

/* Define to ENABLE_SUITESPARSE if the SuiteSparse's AMD library is available */
#cmakedefine HAVE_SUITESPARSE_AMD ENABLE_SUITESPARSE

/* Define to ENABLE_SUITESPARSE if the SuiteSparse's BTF library is available */
#cmakedefine HAVE_SUITESPARSE_BTF ENABLE_SUITESPARSE

/* Define to ENABLE_SUITESPARSE if the SuiteSparse's CAMD library is available */
#cmakedefine HAVE_SUITESPARSE_CAMD ENABLE_SUITESPARSE

/* Define to ENABLE_SUITESPARSE if the SuiteSparse's CCOLAMD library is available */
#cmakedefine HAVE_SUITESPARSE_CCOLAMD ENABLE_SUITESPARSE

/* Define to ENABLE_SUITESPARSE if the SuiteSparse's CHOLMOD library is available */
#cmakedefine HAVE_SUITESPARSE_CHOLMOD ENABLE_SUITESPARSE

/* Define to ENABLE_SUITESPARSE if the SuiteSparse's COLAMD library is available */
#cmakedefine HAVE_SUITESPARSE_COLAMD ENABLE_SUITESPARSE

/* Define to ENABLE_SUITESPARSE if the SuiteSparse's CXSPARSE library is available */
#cmakedefine HAVE_SUITESPARSE_CXSPARSE ENABLE_SUITESPARSE

/* Define to ENABLE_SUITESPARSE if the SuiteSparse's KLU library is available */
#cmakedefine HAVE_SUITESPARSE_KLU ENABLE_SUITESPARSE

/* Define to ENABLE_SUITESPARSE if the SuiteSparse's LDL library is available */
#cmakedefine HAVE_SUITESPARSE_LDL ENABLE_SUITESPARSE

/* Define to ENABLE_SUITESPARSE if the SuiteSparse's RBIO library is available */
#cmakedefine HAVE_SUITESPARSE_RBIO ENABLE_SUITESPARSE

/* Define to ENABLE_SUITESPARSE if the SuiteSparse's SPQR library is available
   and if it's version is at least 4.3 */
#cmakedefine HAVE_SUITESPARSE_SPQR ENABLE_SUITESPARSE

/* Define to ENABLE_SUITESPARSE if the SuiteSparse's UMFPACK library is available */
#cmakedefine HAVE_SUITESPARSE_UMFPACK ENABLE_SUITESPARSE

/* Define to ENABLE_PARMETIS if you have the Parmetis library.
   This is only true if MPI was found
   by configure _and_ if the application uses the PARMETIS_CPPFLAGS */
#cmakedefine HAVE_PARMETIS ENABLE_PARMETIS

/* Define to 1 if PT-Scotch is available */
#cmakedefine HAVE_PTSCOTCH 1

/* Include always useful headers */
#include "FC.h"
#define FC_FUNC FC_GLOBAL_





/* Define to the version of dune-geometry */
#define DUNE_GEOMETRY_VERSION "${DUNE_GEOMETRY_VERSION}"

/* Define to the major version of dune-geometry */
#define DUNE_GEOMETRY_VERSION_MAJOR ${DUNE_GEOMETRY_VERSION_MAJOR}

/* Define to the minor version of dune-geometry */
#define DUNE_GEOMETRY_VERSION_MINOR ${DUNE_GEOMETRY_VERSION_MINOR}

/* Define to the revision of dune-geometry */
#define DUNE_GEOMETRY_VERSION_REVISION ${DUNE_GEOMETRY_VERSION_REVISION}






/* Define to ENABLE_SUPERLU if the SuperLU library is available */
#cmakedefine HAVE_SUPERLU ENABLE_SUPERLU

/* Define to the integer type that SuperLU was compiled for
   See e.g. what int_t is defined to in slu_sdefs.h */
#cmakedefine SUPERLU_INT_TYPE @SUPERLU_INT_TYPE@

/* Define to 1 if header slu_sdefs.h is there. */
#cmakedefine01 HAVE_SLU_SDEFS_H

/* Define to 1 if header slu_ddefs.h is there. */
#cmakedefine01 HAVE_SLU_DDEFS_H

/* Define to 1 if header slu_cdefs.h is there. */
#cmakedefine01 HAVE_SLU_CDEFS_H

/* Define to 1 if header slu_zdefs.h is there. */
#cmakedefine01 HAVE_SLU_ZDEFS_H

/* Define to ENABLE_ARPACKPP if the ARPACK++ library is available */
#cmakedefine HAVE_ARPACKPP ENABLE_ARPACKPP

/* Define to 0 as all versions since SuperLu 4.0 do no longer provide it that way. */
#define HAVE_MEM_USAGE_T_EXPANSIONS 1

/* define to 1 if SuperLU header slu_ddefs.h contains SLU_DOUBLE */
#cmakedefine SUPERLU_MIN_VERSION_4_3 @SUPERLU_MIN_VERSION_4_3@

/* define to 1 if SuperLU dgssvx takes a GlobalLU_t parameter */
#cmakedefine SUPERLU_MIN_VERSION_5 @SUPERLU_MIN_VERSION_5@

/* Define to the version of dune-istl */
#define DUNE_ISTL_VERSION "${DUNE_ISTL_VERSION}"

/* Define to the major version of dune-istl */
#define DUNE_ISTL_VERSION_MAJOR ${DUNE_ISTL_VERSION_MAJOR}

/* Define to the minor version of dune-istl */
#define DUNE_ISTL_VERSION_MINOR ${DUNE_ISTL_VERSION_MINOR}

/* Define to the revision of dune-istl */
#define DUNE_ISTL_VERSION_REVISION ${DUNE_ISTL_VERSION_REVISION}





/* Define to the version of dune-typetree */
#define DUNE_TYPETREE_VERSION "${DUNE_TYPETREE_VERSION}"

/* Define to the major version of dune-typetree */
#define DUNE_TYPETREE_VERSION_MAJOR ${DUNE_TYPETREE_VERSION_MAJOR}

/* Define to the minor version of dune-typetree */
#define DUNE_TYPETREE_VERSION_MINOR ${DUNE_TYPETREE_VERSION_MINOR}

/* Define to the revision of dune-typetree */
#define DUNE_TYPETREE_VERSION_REVISION ${DUNE_TYPETREE_VERSION_REVISION}

/* Define to 1 if std::initializer_list is supported. */
#cmakedefine HAVE_INITIALIZER_LIST 1

/* Define to 1 if template aliases are supported. */
#cmakedefine HAVE_TEMPLATE_ALIASES 1

/* Define to 1 if decltype if supported. */
#cmakedefine HAVE_STD_DECLTYPE 1

/* Define to 1 if GCC's __typeof__ extension is supported. */
#cmakedefine HAVE_GCC___TYPEOF__ 1





/* Define to the version of dune-grid */
#define DUNE_GRID_VERSION "${DUNE_GRID_VERSION}"

/* Define to the major version of dune-grid */
#define DUNE_GRID_VERSION_MAJOR ${DUNE_GRID_VERSION_MAJOR}

/* Define to the minor version of dune-grid */
#define DUNE_GRID_VERSION_MINOR ${DUNE_GRID_VERSION_MINOR}

/* Define to the revision of dune-grid */
#define DUNE_GRID_VERSION_REVISION ${DUNE_GRID_VERSION_REVISION}

/* If this is set, public access to the implementation of facades like Entity,
   Geometry, etc. is granted. */
#cmakedefine DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS 1

/* Define to 1 if psurface library is found */
#cmakedefine HAVE_PSURFACE 1

/* Define to 1 if AmiraMesh library is found */
#cmakedefine HAVE_AMIRAMESH 1

/* The namespace prefix of the psurface library (deprecated) */
#define PSURFACE_NAMESPACE psurface::

/* Define to 1 if you have at least psurface version 2.0 */
#cmakedefine HAVE_PSURFACE_2_0 1

/* Alberta version found by configure, either 0x200 for 2.0 or 0x300 for 3.0 */
#cmakedefine DUNE_ALBERTA_VERSION @DUNE_ALBERTA_VERSION@

/* This is only true if alberta-library was found by configure _and_ if the
   application uses the ALBERTA_CPPFLAGS */
#cmakedefine HAVE_ALBERTA ENABLE_ALBERTA

/* This is only true if UG was found by configure _and_ if the application
   uses the UG_CPPFLAGS */
#cmakedefine HAVE_UG ENABLE_UG

/* Define to 1 if you have mkstemp function */
#cmakedefine01 HAVE_MKSTEMP







/* Define to the version of dune-localfunctions */
#define DUNE_LOCALFUNCTIONS_VERSION "${DUNE_LOCALFUNCTIONS_VERSION}"

/* Define to the major version of dune-localfunctions */
#define DUNE_LOCALFUNCTIONS_VERSION_MAJOR ${DUNE_LOCALFUNCTIONS_VERSION_MAJOR}

/* Define to the minor version of dune-localfunctions */
#define DUNE_LOCALFUNCTIONS_VERSION_MINOR ${DUNE_LOCALFUNCTIONS_VERSION_MINOR}

/* Define to the revision of dune-localfunctions */
#define DUNE_LOCALFUNCTIONS_VERSION_REVISION ${DUNE_LOCALFUNCTIONS_VERSION_REVISION}






/* Define to the version of dune-matrix-vector */
#define DUNE_MATRIX_VECTOR_VERSION "@DUNE_MATRIX_VECTOR_VERSION@"

/* Define to the major version of dune-matrix-vector */
#define DUNE_MATRIX_VECTOR_VERSION_MAJOR @DUNE_MATRIX_VECTOR_VERSION_MAJOR@

/* Define to the minor version of dune-matrix-vector */
#define DUNE_MATRIX_VECTOR_VERSION_MINOR @DUNE_MATRIX_VECTOR_VERSION_MINOR@

/* Define to the revision of dune-matrix-vector */
#define DUNE_MATRIX_VECTOR_VERSION_REVISION @DUNE_MATRIX_VECTOR_VERSION_REVISION@






/* Define to the version of dune-functions */
#define DUNE_FUNCTIONS_VERSION "@DUNE_FUNCTIONS_VERSION@"

/* Define to the major version of dune-functions */
#define DUNE_FUNCTIONS_VERSION_MAJOR @DUNE_FUNCTIONS_VERSION_MAJOR@

/* Define to the minor version of dune-functions */
#define DUNE_FUNCTIONS_VERSION_MINOR @DUNE_FUNCTIONS_VERSION_MINOR@

/* Define to the revision of dune-functions */
#define DUNE_FUNCTIONS_VERSION_REVISION @DUNE_FUNCTIONS_VERSION_REVISION@






/* Define to the version of dune-fufem */
#define DUNE_FUFEM_VERSION "@DUNE_FUFEM_VERSION@"

/* Define to the major version of dune-fufem */
#define DUNE_FUFEM_VERSION_MAJOR @DUNE_FUFEM_VERSION_MAJOR@

/* Define to the minor version of dune-fufem */
#define DUNE_FUFEM_VERSION_MINOR @DUNE_FUFE;_VERSION_MINOR@

/* Define to the revision of dune-fufem */
#define DUNE_FUFEM_VERSION_REVISION @DUNE_FUFEM_VERSION_REVISION@

/* Define to 1 if adolc library is found */
#cmakedefine HAVE_ADOLC 1

/* Define to 1 if Boost serialization library is found */
#cmakedefine HAVE_BOOST_SERIALIZATION 1



/* begin dune-PFASST
   put the definitions for config.h specific to
   your project here. Everything above will be
   overwritten
*/

/* begin private */
/* Name of package */
#define PACKAGE "@DUNE_MOD_NAME@"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "@DUNE_MAINTAINER@"

/* Define to the full name of this package. */
#define PACKAGE_NAME "@DUNE_MOD_NAME@"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "@DUNE_MOD_NAME@ @DUNE_MOD_VERSION@"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "@DUNE_MOD_NAME@"

/* Define to the home page for this package. */
#define PACKAGE_URL "@DUNE_MOD_URL@"

/* Define to the version of this package. */
#define PACKAGE_VERSION "@DUNE_MOD_VERSION@"

/* end private */

/* Define to the version of dune-PFASST */
#define DUNE_PFASST_VERSION "@DUNE_PFASST_VERSION@"

/* Define to the major version of dune-PFASST */
#define DUNE_PFASST_VERSION_MAJOR @DUNE_PFASST_VERSION_MAJOR@

/* Define to the minor version of dune-PFASST */
#define DUNE_PFASST_VERSION_MINOR @DUNE_PFASST_VERSION_MINOR@

/* Define to the revision of dune-PFASST */
#define DUNE_PFASST_VERSION_REVISION @DUNE_PFASST_VERSION_REVISION@

/* end dune-PFASST
   Everything below here will be overwritten
*/ 

/* Grid type magic for DGF parser */
@GRID_CONFIG_H_BOTTOM@

