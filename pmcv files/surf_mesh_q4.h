#ifndef __CM2_SURF_MESH_Q4_H__
#define __CM2_SURF_MESH_Q4_H__

/**
   \file       surf_mesh_q4.h
   \brief      File to include the cm2::surf_mesh_q4 classes & routines.
   \copyright  (C)1999-2017, Computing Objects, France. info@computing-objects.com

   $Rev: 1245 $
   $Date: 2013-02-02 14:43:39 +0100 (sam., 02 févr. 2013) $ 
   */

/******************************************************************************************
   This program is not free software. It is the sole property of Computing Objects, France.
   You are allowed to modify it, to recompile it, to port it on several platforms,
   but not to redistribute it, even the modified version.
 ******************************************************************************************/


// CONFIGURATION
#include "cm2_api_config.h"
#include "cm2_api_stl.h"


// Define CM2SURFMESHQ4IMPORT, CM2SURFMESHQ4EXPORT, CM2SURFMESHQ4LOCAL and CM2SURFMESHQ4PUBLIC.
#if (defined(WIN32) | defined(WIN64))
   #if defined(_USRDLL)
      #define CM2SURFMESHQ4EXPORT __declspec(dllexport)
   #else
      #define CM2SURFMESHQ4EXPORT
   #endif
   #if defined(CM2_LIBIMPORT)
      #define CM2SURFMESHQ4IMPORT
   #else
      #define CM2SURFMESHQ4IMPORT __declspec(dllimport)
   #endif
   #define CM2SURFMESHQ4LOCAL
   #define CM2SURFMESHQ4PUBLIC
#elif defined(CM2_USE_GCC_VISIBILITY_ATTRIBS)
   #define CM2SURFMESHQ4IMPORT
   #define CM2SURFMESHQ4EXPORT __attribute__ ((visibility("default")))
   #define CM2SURFMESHQ4LOCAL __attribute__ ((visibility("hidden")))
   #define CM2SURFMESHQ4PUBLIC __attribute__ ((visibility("default")))
#else
   #define CM2SURFMESHQ4IMPORT
   #define CM2SURFMESHQ4EXPORT
   #define CM2SURFMESHQ4LOCAL
   #define CM2SURFMESHQ4PUBLIC
#endif

// Define CM2SURFMESHQ4_API for DLL builds
#ifdef CM2SURFMESHQ4_EXPORTS
   #define CM2SURFMESHQ4_API CM2SURFMESHQ4EXPORT
#else
   #define CM2SURFMESHQ4_API CM2SURFMESHQ4IMPORT
#endif


// SURF_MESH_Q4 API
#include "surf_mesher_q4.h"


/**
   \namespace cm2::surf_mesh_q4
   \brief Namespace for the 3d-surface triangle mesher.
   */
namespace cm2 {
namespace surf_mesh_q4 {

/// Versioning.
CM2SURFMESHQ4_API const char*
version();

/**
   Function to unlock the DLL.
   The user must provide two strings in order to use this DLL.
   To get or update a license, please contact license@computing-objects.com

   \param[in]   agreement     A short user agreement in human readable form.
   \param[in]   user_key      A secure encoded string.
   */
CM2SURFMESHQ4_API void
registration (const char* agreement, const char* user_key);

}  // namespace surf_mesh_q4
}  // namespace cm2

#endif

