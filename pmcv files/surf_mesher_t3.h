#ifndef __SURF_MESHER_T3_H__
#define __SURF_MESHER_T3_H__

/**
   \file       surf_mesher_t3.h
   \brief      Definition of the triangle 3d-surface mesher.
   \copyright  (C)1999-2017, Computing Objects, France. info@computing-objects.com

   $Rev: 1294 $
   $Date: 2013-04-22 14:32:52 +0200 (lun., 22 avr. 2013) $ 
   */

/******************************************************************************************
   This program is not free software. It is the sole property of Computing Objects, France.
   You are allowed to modify it, to recompile it, to port it on several platforms,
   but not to redistribute it, even the modified version.
 ******************************************************************************************/


#include "fe_defs.h"
#include "math1.h"
#include "misc.h"


namespace cm2 {
namespace surf_mesh_t3 {

/**
   3-D surface triangle mesher.

   Parametric surface mesher (direct meshing) based on the OpenCascade kernel. 

   Runs on STEP and IGES files.

   \image html elephant_t.png "Surface mesh generated by CM2 SurfMesh T3. CAD model made of 57 surfaces (colorized)."

   \note    Meshing failures can be caused by pathological 
            surfaces (for instance with foldings and self-intersections).
            In these cases, check your model and fix.
   \note    At present, the mesh-size gradation control doesn't work between curves nor surfaces
            (only within each curve, and each surface).
   \note    Bad curve definition can cause accumulation points. 
            A post-meshing correction step may be needed (from simple node merging to local meshings).
   \note    Non-connected surfaces (for instance surfaces with gaps between them) give unconnected meshes
            (this tool doesn't fix these kinds of topological error).
   \warning Still experimental.
   */
class CM2SURFMESHT3_API mesher
{
public:

   /**
      This generator is for triangle meshes.
      */
   enum
   {
      generated_FE_type = CM2_FACET3         //!< Generates all-triangle meshes (CM2_FACET3).
   };


   /**
      The settings of surf_mesh_t3::mesher.

      This structure gathers flags and parameters to control the way the %mesher works.
      */
   struct CM2SURFMESHT3_API settings_type
   {

   /**
      Typedef for the display handler.

      \param[in]  pass_thru   The \ref settings_type::pass_thru pointer (user provided).
      \param[in]  level       The importance of the message: 
                                 -  0 = Basic message.
                                 -  1 = Somewhat detailed message.
                                 -  2 = Very detailed message (DEBUG mode only).
      \param[in]  msg         The message.
      */
   typedef void (*display_handler_type) (void* pass_thru, unsigned level, const char* msg);

   /**
      Typedef for the interrupt handler.

      The interrupt handler, if any, is called internally by the %mesher from time to time 
      to query the abortion of process (when the handler returns true) or to continue 
      the meshing process (when the handler returns false).

      \param[in]  pass_thru   The \ref settings_type::pass_thru pointer (user provided).
      \param[in]  progress    Value (between 0 = just starting and 1 = meshing is complete) 
                              to give an estimate of the progress of the meshing.
      */
   typedef bool (*interrupt_handler_type) (void* pass_thru, double progress);


   /**
      Default constructor.

      The default values are:
         - \ref target_h = 0.0
         - \ref max_chordal_error = -0.02
         - \ref min_h = -0.001
         - \ref chordal_control_type = 4
         - \ref max_gradation = 0.5
         - \ref force_even_flag = `false`
         - \ref optim_level = 3
         - \ref display_hdl = `NULL`
         - \ref interrupt_hdl = `NULL`
         - \ref pass_thru = `NULL`

      The most useful parameters are:
         - \ref target_h to control the size of the elements in the mesh.
         - \ref max_chordal_error to control the chordal error between the mesh and the surface.
      */
   settings_type()
      : target_h(0.0), max_chordal_error(-0.02), min_h(-0.001),
        chordal_control_type(4), max_gradation(0.5), 
        force_even_flag(false), optim_level(3), 
        display_hdl(NULL), interrupt_hdl(NULL), pass_thru(NULL)
      { }

   /**
      Returns true if the sanity check is successful.
      */
   bool
   valid() const        { return this->check() == 0; }

   /**
      Sanity check of the settings.

      \return        Returns a negative value in case of error. \n
                     Returned value = -k => the k-th field of the structure is illegal
                     (\p target_h is the first field, \p max_chordal_error is the second...)
      */
   int
   check() const;

   /**
      Resets to the \link cm2::surf_mesh_t3::mesher::settings_type::settings_type() default values \endlink.
      */
   void
   reset();

   /**
      The mesh size. 

      Default = 0.

      \pre        target_h >= 0.
      */
   double                  target_h;

   /**
      Maximum chordal error allowed.

      The mesh size will be reduced locally to limit the chordal error between the mesh and the surface:
         - If negative, this value is relative to the local radii (for instance -0.01 => max chordal error = 1% of local radii).
         - If positive, this value is absolute (for instance 0.1 => max chordal error = 0.1).

      With max_chordal_error = -0.05 (5%), a complete circle is meshed with approximatively 10 elements. \n
      With max_chordal_error = -0.03 (3%), a complete circle is meshed with approximatively 13 elements. \n
      With max_chordal_error = -0.02 (2%), a complete circle is meshed with approximatively 16 elements. \n
      With max_chordal_error = -0.01 (1%), a complete circle is meshed with approximatively 22 elements.

      Default = -0.02 (max error = 2% of the local radii).
      */
   double                  max_chordal_error;

   /**
      Minimum mesh size (security). 
                                          
      The mesh size will not be reduced (by the \ref max_chordal_error criterion) beyond this limit:
         - If negative, this value is relative to \ref target_h (min size = - min_h * \ref target_h).
         - If positive, this value is absolute.

      Default = -0.001 (1E-3 * \ref target_h).
      */
   double                  min_h;

   /**
      Chordal control type:
         -  0: No curvature is computed and the chordal error control is disabled.
         -  1: Approximate curvature (least-square method using the local bases). \n
               Isotropic chordal control: element sizes are limited by the same (isotropic) chordal value.
         -  2: Approximate curvature (least-square method using the local bases). \n
               Anisotropic chordal control: element sizes are limited by two different (anisotropic) chordal values
               along the two principal curvature directions.
         -  3: Exact curvature (provided by `Surface::get_local_curvatures`). \n
               Isotropic chordal control: element sizes are limited by the same (isotropic) chordal value.
         -  4: Exact curvature (provided by `Surface::get_local_curvatures`). \n
               Anisotropic chordal control: element sizes are limited by two different (anisotropic) chordal values
               along the two principal curvature directions.

      Default = 4
      */
   unsigned                chordal_control_type;

   /**
      Gradation of the elements size.

      This parameter controls the gradation of the elements size along the skeleton lines
      and inside the patches towards the \ref target_h mesh size.
      
      A value close to 0 leads to a more progressive variation of mesh size (smoother).

      \note    max_gradation >= 0

      Default = 0.5
      */
   double                  max_gradation;

   /**
      Flag to force the number of edges along each curve to be even.
      
      Default = false.
      */
   bool                    force_even_flag;

   /**
      Optimization level between 0 and 10.

      Controls the trade-off between speed (low values) and high element quality (high values). \n
      The higher the level, the more aggressive the algorithms.

      \note    Level 0 disables any optimization.

      Default = 3
      */
   unsigned                optim_level;

   /**
      The user-provided display function for verbose.

      Level of message: 
         -  0 = Basic message.
         -  1 = Somewhat detailed message.
         -  2 = Very detailed message (debug mode only).

      Default = `NULL`.
      */
   display_handler_type    display_hdl;

   /**
      The user-provided interrupt function.

      When this function returns true, the process is stopped and the meshing aborts.

      The parameter `progress` (between 0 and 1) gives a hint about the progress
      of the meshing.

      Default = `NULL`.
      */
   interrupt_handler_type  interrupt_hdl;

   /**
      The pointer that will be passed to the \ref display_hdl and \ref interrupt_hdl
      (if any) when called by the mesher.

      Example:
      \code
         void my_display_hdl (void* pass_thru, unsigned level, const char* msg)
         {
            window_type*   my_window = static_cast<window_type*>(pass_thru);
            my_window->show(msg);
         }

         mesher         my_mesher;
         data_type      my_mesh_data(pos, connectB);
         window_type    my_window;                    // A window instance.
         my_mesher.settings.display_hdl = &my_display_hdl;
         my_mesher.settings.pass_thru = static_cast<void*>(&my_window);   
         my_mesher.run(my_mesh_data);                 // Will call my_display_hdl.
      \endcode
      */
   void*                   pass_thru;

   };    // settings_type




   /**
      The data type for surf_mesh_t3::mesher.

      This structure is used to gather all the input and output data of the function 
      surf_mesh_t3::mesher::run.
      */
   struct CM2SURFMESHT3_API data_type
   {
   public:

   ///
   typedef settings_type::display_handler_type     display_handler_type;

   /**
      Warning enums for the mesher.

      Non fatal warnings.

      If such a warning is raised (returned in the field \ref warning_code),
      the mesher is not stopped but the user shall check its input and output data for a potential mistake.
      */
   enum warning_type
   { 
      CM2_NO_WARNING                =   0,   //!< OK, no warning.
      CM2_INTERRUPTION              = -10,   //!< The process has been interrupted (by the user through the interrupt handler).
      CM2_FAILED_SURFACES_WARNING   = -11,   //!< One or several surfaces could not be meshed.
      CM2_SHAPE_QUALITY_WARNING     = -12    //!< One or several elements have a bad quality (Qs < 0.01).
   };

   /**
      Error enums for the mesher.
       
      Fatal errors.

      If such an error is raised (returned in the field \ref error_code),
      the mesher is stopped and no final mesh is returned.
      */
   enum error_type
   { 
      CM2_NO_ERROR                  =    0,  //!< OK, no error.
      CM2_LICENSE_ERROR             = -100,  //!< Invalid or expired license.
      CM2_MODE_ERROR                = -101,  //!< Error in settings.
      CM2_FILE_ERROR                = -102,  //!< Error in input file. Missing file or OpenCascade could not read it properly.
      CM2_DEGENERATED_ELEMENT       = -103,  //!< The process leads to some degenerated elements and therefore the mesher cannot build a valid mesh.
      CM2_SYSTEM_MEMORY_ERROR       = -199,  //!< Insufficient memory available.
      CM2_INTERNAL_ERROR            = -200   //!< Something went wrong but don't know what (please contact support@computing-objects.com).
   };


   /// Default constructor.
   data_type()                            { this->clear(); }

   /**
      Copy constructor (hard copy).
      All matrices are hard copied.
      */
   data_type (const data_type& cpy)       { (*this) = cpy; }


   /**@name Copy Members */
   //@{
   /// Copy operator (hard copy).
   const data_type&
   operator= (const data_type& data);
   /// Shallow copy.
   void
   shallow_copy (const data_type& data);
   //@}


   /// Shallow copy.
   void
   extract (DoubleMat& P, UIntMat& connectT3) const
      { P = this->pos; connectT3 = this->connectM; }

   /// Shallow copy. Provided for compatibility with CM2 SurfMesh Q4.
   void
   extract (DoubleMat& P, UIntMat& connectQ4, UIntMat& connectT3) const
      { P = this->pos; connectQ4.clear(); connectT3 = this->connectM; }

   /// Clears all fields and reinit to default values.
   void
   clear();

   /**
      Prints information about the mesh using a display handler.

      \param   hdl   The display handler.

      \note    The first parameter of the handler (param `pass_thru`) is irrelevant (`NULL`).
      */
   void
   print_info (display_handler_type hdl) const;

   /**
      Saves the data into a file.

      \param   filename       The file name for the output (overwrite mode).
      \param   precision      The precision of the output for the floating-point values.
      */
   void
   save (const char* filename, unsigned precision = 16) const;

   /// Checks if an error has been raised.
   inline bool
   error_raised() const       { return error_code != CM2_NO_ERROR; }

   /// Checks if a warning has been raised.
   inline bool
   warning_raised() const     { return warning_code != CM2_NO_WARNING; }


   /**
      The coordinates matrix (column oriented).

      Node `i` has coordinates in column `pos.col(i)`.

      The dimensions of this matrix are 3 x \ref nods.

      Mode = OUT.
      */
   DoubleMat               pos;

   /**
      The connectivity of the triangle mesh (column oriented). 

      The indices refer to columns in the coordinates matrix \ref pos.

      The dimensions of this matrix are 3 x \ref nefs.

      Mode = OUT.
      */
   UIntMat                 connectM;

   /**
      The coloring of the mesh.

      Each element is assigned a color (from 0 to \ref nbr_surfaces - 1) corresponding
      to the ID of its surface.

      The size of this array equals to the number of columns in \ref connectM.

      Mode = OUT.
      */
   UIntVec                 colors;

   /**
      The IDs of the failed surfaces. 

      Mode = OUT.
      */
   UIntMat                 failed;

   /**
      The histogram of the element shape qualities.

      Mode = OUT.
      */ 
   cm2::misc::histogram    histo_Qs; 

   /**
      Number of curves in the model. 

      Mode = OUT.
      */ 
   unsigned                nbr_curves;

   /**
      Number of surfaces in the model. 

      Mode = OUT.
      */ 
   unsigned                nbr_surfaces;

   /**
      Number of elements in the final mesh. 

      Mode = OUT.
      */ 
   unsigned                nefs;

   /**
      Number of quads in the final mesh (always null). 

      Provided for compatibility with surf_mesh_q4::mesher.

      Mode = OUT.
      */ 
   unsigned                nefs_Q4;

   /**
      Number of triangles in the final mesh. Same as \ref nefs.

      Mode = OUT.
      */ 
   unsigned                nefs_T3;

   /**
      Number of nodes in the final mesh. 

      Mode = OUT.
      */ 
   unsigned                nods;

   /**
      Area of the meshed domain (sum of the areas of all elements).

      Mode = OUT.
      */ 
   double                  area;

   /**
      Area of the quadrangles. Always null.

      Provided for compatibility with surf_mesh_q4::mesher.

      Mode = OUT.
      */ 
   double                  area_Q4;

   /**
      Area of the triangles. Same as \ref area.

      Provided for compatibility with surf_mesh_q4::mesher.

      Mode = OUT.
      */ 
   double                  area_T3;

   /**
      Minimum shape quality of the elements upon exit. 

      Mode = OUT.
      */ 
   double                  Qmin;

   /**
      Time spent in the OCC kernel for reading the model. 

      Mode = OUT.
      */ 
   double                  read_time;

   /**
      Time spent for meshing all the curves. 
    
      Mode = OUT.
      */ 
   double                  meshC_time;

   /**
      Time spent for meshing all the surfaces. 
    
      Mode = OUT.
      */ 
   double                  meshS_time;

   /**
      Time spent for meshing the CAD model (curves + surfaces). 
    
      Mode = OUT.
      */ 
   double                  mesh_time;

   /**
      Total time spent for the reading the model and meshing the model. 

      Mode = OUT.
      */ 
   double                  total_time; 

   /**
      Meshing speed, excluding the OCC-reading time (number of elements per second). 

      Mode = OUT.
      */ 
   double                  speed;

   /**
      Error code. 

      Mode = OUT.
      */ 
   error_type              error_code;

   /**
      Warning code. 

      Mode = OUT.
      */ 
   warning_type            warning_code;

   /**
      Error/Warning message #1. 

      Mode = OUT.
      */ 
   char                    msg1[256];

   /**
      Error/Warning message #2.
      Debug purpose only.

      Mode = OUT.
      */ 
   char                    msg2[256];

   };    // data_type


/**
   Constructor with default settings.

   See \link cm2::surf_mesh_t3::mesher::settings_type::settings_type() default values \endlink.
   */
mesher();

/**
   Constructor with specific settings.
   
   \param[in]   settings_        The settings. 
   */
mesher (const settings_type& settings_);


/**
   Main function of the mesh generator.

   \param[in]    filename  The IGES/STEP file containing the CAD model.
   \param[out]   data      The data.
   */
void 
run (const char* filename, data_type& data) const;


/**
   The settings of the %mesher.

   See \link cm2::surf_mesh_t3::mesher::settings_type::settings_type() default values \endlink.
   */
settings_type         settings;

};

}  // namespace surf_mesh_t3
}  // namespace cm2


/** 
   \example tmsh3ds01.cpp
   This is an example of how to use cm2::surf_mesh_t3::mesher and cm2::surf_mesh_q4::mesher to mesh directly onto
   a CAD model (IGES or STEP).

   \image html elephant_t.png "Surface mesh generated with CM2 SurfMesh T3. CAD model made of 57 surfaces (colorized)."
   \image html elephant_qt.png "Surface mesh generated with CM2 SurfMesh Q4. CAD model made of 57 surfaces (colorized)."
 */

#endif
