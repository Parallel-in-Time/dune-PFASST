/* config.h.  Generated from config_collected.h.cmake by CMake.
   It was generated from config_collected.h.cmake which in turn is generated automatically
   from the config.h.cmake files of modules this module depends on. */

/* Define to 1 if you have module dune-PFASST available */
#define HAVE_DUNE_PFASST 1


/* Define to 1 if you have module dune-common available */
#define HAVE_DUNE_COMMON 1


/* Define to 1 if you have module dune-uggrid available */
#define HAVE_DUNE_UGGRID 0


/* Define to 1 if you have module dune-geometry available */
#define HAVE_DUNE_GEOMETRY 1


/* Define to 1 if you have module dune-istl available */
#define HAVE_DUNE_ISTL 1


/* Define to 1 if you have module dune-typetree available */
#define HAVE_DUNE_TYPETREE 1


/* Define to 1 if you have module dune-grid available */
#define HAVE_DUNE_GRID 1


/* Define to 1 if you have module dune-localfunctions available */
#define HAVE_DUNE_LOCALFUNCTIONS 1


/* Define to 1 if you have module dune-alugrid available */
#define HAVE_DUNE_ALUGRID 0


/* Define to 1 if you have module dune-matrix-vector available */
#define HAVE_DUNE_MATRIX_VECTOR 1


/* Define to 1 if you have module dune-functions available */
#define HAVE_DUNE_FUNCTIONS 1


/* Define to 1 if you have module dune-fufem available */
#define HAVE_DUNE_FUFEM 1


/* Define to the version of dune-common */
#define DUNE_COMMON_VERSION "2.6-git"

/* Define to the major version of dune-common */
#define DUNE_COMMON_VERSION_MAJOR 2

/* Define to the minor version of dune-common */
#define DUNE_COMMON_VERSION_MINOR 6

/* Define to the revision of dune-common */
#define DUNE_COMMON_VERSION_REVISION 0

/* Standard debug streams with a level below will collapse to doing nothing */
#define DUNE_MINIMAL_DEBUG_LEVEL 4

/* does the compiler support __attribute__((deprecated))? */
#define HAS_ATTRIBUTE_DEPRECATED 1

/* does the compiler support __attribute__((deprecated("message"))? */
#define HAS_ATTRIBUTE_DEPRECATED_MSG 1

/* does the compiler support __attribute__((unused))? */
#define HAS_ATTRIBUTE_UNUSED 1

/* Define if you have a BLAS library. */
#define HAVE_BLAS 1

/* does the compiler support abi::__cxa_demangle */
#define HAVE_CXA_DEMANGLE 1

/* Define if you have LAPACK library. */
#define HAVE_LAPACK 1

/* Define to 1 if you have the <malloc.h> header file. */
// Not used! #define HAVE_MALLOC_H 0

/* Define if you have the MPI library.  */
#define HAVE_MPI ENABLE_MPI

/* Define if you have the GNU GMP library. The value should be ENABLE_GMP
   to facilitate activating and deactivating GMP using compile flags. */
/* #undef HAVE_GMP */

/* Define if you have the Vc library. The value should be ENABLE_VC
   to facilitate activating and deactivating Vc using compile flags. */
/* #undef HAVE_VC */

/* Define to 1 if you have the symbol mprotect. */
#define HAVE_MPROTECT 1

/* Define to 1 if you have the <stdint.h> header file. */
/* #undef HAVE_STDINT_H */

/* Define to 1 if you have <sys/mman.h>. */
#define HAVE_SYS_MMAN_H 1



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
#define HAVE_IS_INDEXABLE_SUPPORT 1

/* Define to ENABLE_UMFPACK if the UMFPack library is available */
/* #undef HAVE_UMFPACK */

/* Define to ENABLE_SUITESPARSE if the SuiteSparse library is available */
/* #undef HAVE_SUITESPARSE */

/* Define to ENABLE_SUITESPARSE if the SuiteSparse's AMD library is available */
/* #undef HAVE_SUITESPARSE_AMD */

/* Define to ENABLE_SUITESPARSE if the SuiteSparse's BTF library is available */
/* #undef HAVE_SUITESPARSE_BTF */

/* Define to ENABLE_SUITESPARSE if the SuiteSparse's CAMD library is available */
/* #undef HAVE_SUITESPARSE_CAMD */

/* Define to ENABLE_SUITESPARSE if the SuiteSparse's CCOLAMD library is available */
/* #undef HAVE_SUITESPARSE_CCOLAMD */

/* Define to ENABLE_SUITESPARSE if the SuiteSparse's CHOLMOD library is available */
/* #undef HAVE_SUITESPARSE_CHOLMOD */

/* Define to ENABLE_SUITESPARSE if the SuiteSparse's COLAMD library is available */
/* #undef HAVE_SUITESPARSE_COLAMD */

/* Define to ENABLE_SUITESPARSE if the SuiteSparse's CXSPARSE library is available */
/* #undef HAVE_SUITESPARSE_CXSPARSE */

/* Define to ENABLE_SUITESPARSE if the SuiteSparse's KLU library is available */
/* #undef HAVE_SUITESPARSE_KLU */

/* Define to ENABLE_SUITESPARSE if the SuiteSparse's LDL library is available */
/* #undef HAVE_SUITESPARSE_LDL */

/* Define to ENABLE_SUITESPARSE if the SuiteSparse's RBIO library is available */
/* #undef HAVE_SUITESPARSE_RBIO */

/* Define to ENABLE_SUITESPARSE if the SuiteSparse's SPQR library is available
   and if it's version is at least 4.3 */
/* #undef HAVE_SUITESPARSE_SPQR */

/* Define to ENABLE_SUITESPARSE if the SuiteSparse's UMFPACK library is available */
/* #undef HAVE_SUITESPARSE_UMFPACK */

/* Define to ENABLE_PARMETIS if you have the Parmetis library.
   This is only true if MPI was found
   by configure _and_ if the application uses the PARMETIS_CPPFLAGS */
/* #undef HAVE_PARMETIS */

/* Define to 1 if PT-Scotch is available */
/* #undef HAVE_PTSCOTCH */

/* Include always useful headers */
#include "FC.h"
#define FC_FUNC FC_GLOBAL_





/* Define to the version of dune-geometry */
#define DUNE_GEOMETRY_VERSION "2.6-git"

/* Define to the major version of dune-geometry */
#define DUNE_GEOMETRY_VERSION_MAJOR 2

/* Define to the minor version of dune-geometry */
#define DUNE_GEOMETRY_VERSION_MINOR 6

/* Define to the revision of dune-geometry */
#define DUNE_GEOMETRY_VERSION_REVISION 0






/* Define to ENABLE_SUPERLU if the SuperLU library is available */
/* #undef HAVE_SUPERLU */

/* Define to the integer type that SuperLU was compiled for
   See e.g. what int_t is defined to in slu_sdefs.h */
#define SUPERLU_INT_TYPE int

/* Define to 1 if header slu_sdefs.h is there. */
#define HAVE_SLU_SDEFS_H 0

/* Define to 1 if header slu_ddefs.h is there. */
#define HAVE_SLU_DDEFS_H 0

/* Define to 1 if header slu_cdefs.h is there. */
#define HAVE_SLU_CDEFS_H 0

/* Define to 1 if header slu_zdefs.h is there. */
#define HAVE_SLU_ZDEFS_H 0

/* Define to ENABLE_ARPACKPP if the ARPACK++ library is available */
/* #undef HAVE_ARPACKPP */

/* Define to 0 as all versions since SuperLu 4.0 do no longer provide it that way. */
#define HAVE_MEM_USAGE_T_EXPANSIONS 1

/* define to 1 if SuperLU header slu_ddefs.h contains SLU_DOUBLE */
/* #undef SUPERLU_MIN_VERSION_4_3 */

/* define to 1 if SuperLU dgssvx takes a GlobalLU_t parameter */
/* #undef SUPERLU_MIN_VERSION_5 */

/* Define to the version of dune-istl */
#define DUNE_ISTL_VERSION "2.6-git"

/* Define to the major version of dune-istl */
#define DUNE_ISTL_VERSION_MAJOR 2

/* Define to the minor version of dune-istl */
#define DUNE_ISTL_VERSION_MINOR 6

/* Define to the revision of dune-istl */
#define DUNE_ISTL_VERSION_REVISION 0





/* Define to the version of dune-typetree */
#define DUNE_TYPETREE_VERSION "2.5-dev"

/* Define to the major version of dune-typetree */
#define DUNE_TYPETREE_VERSION_MAJOR 2

/* Define to the minor version of dune-typetree */
#define DUNE_TYPETREE_VERSION_MINOR 5

/* Define to the revision of dune-typetree */
#define DUNE_TYPETREE_VERSION_REVISION 0

/* Define to 1 if std::initializer_list is supported. */
/* #undef HAVE_INITIALIZER_LIST */

/* Define to 1 if template aliases are supported. */
#define HAVE_TEMPLATE_ALIASES 1

/* Define to 1 if decltype if supported. */
#define HAVE_STD_DECLTYPE 1

/* Define to 1 if GCC's __typeof__ extension is supported. */
#define HAVE_GCC___TYPEOF__ 1





/* Define to the version of dune-grid */
#define DUNE_GRID_VERSION "2.6-git"

/* Define to the major version of dune-grid */
#define DUNE_GRID_VERSION_MAJOR 2

/* Define to the minor version of dune-grid */
#define DUNE_GRID_VERSION_MINOR 6

/* Define to the revision of dune-grid */
#define DUNE_GRID_VERSION_REVISION 0

/* If this is set, public access to the implementation of facades like Entity,
   Geometry, etc. is granted. */
/* #undef DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS */

/* Define to 1 if psurface library is found */
/* #undef HAVE_PSURFACE */

/* Define to 1 if AmiraMesh library is found */
/* #undef HAVE_AMIRAMESH */

/* The namespace prefix of the psurface library (deprecated) */
#define PSURFACE_NAMESPACE psurface::

/* Define to 1 if you have at least psurface version 2.0 */
#define HAVE_PSURFACE_2_0 1

/* Alberta version found by configure, either 0x200 for 2.0 or 0x300 for 3.0 */
/* #undef DUNE_ALBERTA_VERSION */

/* This is only true if alberta-library was found by configure _and_ if the
   application uses the ALBERTA_CPPFLAGS */
/* #undef HAVE_ALBERTA */

/* This is only true if UG was found by configure _and_ if the application
   uses the UG_CPPFLAGS */
/* #undef HAVE_UG */

/* Define to 1 if you have mkstemp function */
#define HAVE_MKSTEMP 1







/* Define to the version of dune-localfunctions */
#define DUNE_LOCALFUNCTIONS_VERSION "2.6-git"

/* Define to the major version of dune-localfunctions */
#define DUNE_LOCALFUNCTIONS_VERSION_MAJOR 2

/* Define to the minor version of dune-localfunctions */
#define DUNE_LOCALFUNCTIONS_VERSION_MINOR 6

/* Define to the revision of dune-localfunctions */
#define DUNE_LOCALFUNCTIONS_VERSION_REVISION 0






/* Define to the version of dune-matrix-vector */
#define DUNE_MATRIX_VECTOR_VERSION "0.1"

/* Define to the major version of dune-matrix-vector */
#define DUNE_MATRIX_VECTOR_VERSION_MAJOR 0

/* Define to the minor version of dune-matrix-vector */
#define DUNE_MATRIX_VECTOR_VERSION_MINOR 1

/* Define to the revision of dune-matrix-vector */
#define DUNE_MATRIX_VECTOR_VERSION_REVISION 0






/* Define to the version of dune-functions */
#define DUNE_FUNCTIONS_VERSION "2.5-dev"

/* Define to the major version of dune-functions */
#define DUNE_FUNCTIONS_VERSION_MAJOR 2

/* Define to the minor version of dune-functions */
#define DUNE_FUNCTIONS_VERSION_MINOR 5

/* Define to the revision of dune-functions */
#define DUNE_FUNCTIONS_VERSION_REVISION 0






/* Define to the version of dune-fufem */
#define DUNE_FUFEM_VERSION "0.1"

/* Define to the major version of dune-fufem */
#define DUNE_FUFEM_VERSION_MAJOR 0

/* Define to the minor version of dune-fufem */
#define DUNE_FUFEM_VERSION_MINOR @DUNE_FUFE;_VERSION_MINOR@

/* Define to the revision of dune-fufem */
#define DUNE_FUFEM_VERSION_REVISION 0

/* Define to 1 if adolc library is found */
/* #undef HAVE_ADOLC */

/* Define to 1 if Boost serialization library is found */
#define HAVE_BOOST_SERIALIZATION 1



/* begin dune-PFASST
   put the definitions for config.h specific to
   your project here. Everything above will be
   overwritten
*/

/* begin private */
/* Name of package */
#define PACKAGE "dune-PFASST"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "r.schoebel@fz-juelich.de"

/* Define to the full name of this package. */
#define PACKAGE_NAME "dune-PFASST"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "dune-PFASST 1.0"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "dune-PFASST"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "1.0"

/* end private */

/* Define to the version of dune-PFASST */
#define DUNE_PFASST_VERSION "1.0"

/* Define to the major version of dune-PFASST */
#define DUNE_PFASST_VERSION_MAJOR 1

/* Define to the minor version of dune-PFASST */
#define DUNE_PFASST_VERSION_MINOR 0

/* Define to the revision of dune-PFASST */
#define DUNE_PFASST_VERSION_REVISION 0

/* end dune-PFASST
   Everything below here will be overwritten
*/ 

/* Grid type magic for DGF parser */

/* ONEDGRID not available, enable with cmake variable DUNE_GRID_GRIDTYPE_SELECTOR=ON */
/* YASPGRID not available, enable with cmake variable DUNE_GRID_GRIDTYPE_SELECTOR=ON */

