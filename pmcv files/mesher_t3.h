#ifndef __MESHER_T3_H__
#define __MESHER_T3_H__

/**
   \file       mesher_t3.h
   \brief      Definition of the unstructured isotropic triangle %mesher.
   \copyright  (C)1999-2017, Computing Objects, France. info@computing-objects.com

   $Rev: 2974 $
   $Date: 2017-12-22 18:10:37 +0100 (ven., 22 déc. 2017) $ 
   */

/******************************************************************************************
   This program is not free software. It is the sole property of Computing Objects, France.
   You are allowed to modify it, to recompile it, to port it on several platforms,
   but not to redistribute it, even the modified version.
 ******************************************************************************************/


#include "fe_defs.h"
#include "math1.h"
#include "meshtools.h"
#include "misc.h"


namespace cm2 {
namespace triamesh { 

/**
   Unstructured triangle %mesher.

   \image html tmsh2d09_t.png

   Instances of this class generate triangle meshes on any 2-D domains defined 
   by some boundary edges and possibly some embedded nodes and edges.
 
   They can also be used as shape optimizers of an already existing triangle mesh.

   All the data to be passed to the %mesher are gathered in a single 
   structure of type triamesh::mesher::data_type.

   Example:
   \verbatim

   Mesh of a square

   INPUT:
                                                 3               2
      data.pos = [0 4 4 0                         *-------------*
                  0 0 4 4]                        |             |
                                                  |             |
      data.connectB = [0 1 2 3                    |             |
                       1 2 3 0]                   |             |
                                                  |             |
                                                  |             |
                                                  |             |
                                                  *-------------*
                                                 0               1


   triamesh::mesher   my_mesher;
   my_mesher.run(data);


   OUTPUT:
                                                 3               2
      data.pos = [0  4  4  0  2                   *-------------*
                  0  0  4  4  2]                  |\           /|
                                                  |  \       /  |
      data.connectM = [0 4 3 0                    |    \ 4 /    |
                       1 1 4 4                    |      *      |
                       4 2 2 3]                   |    /   \    |
                                                  |  /       \  |
          (indices may vary)                      |/           \|
                                                  *-------------*
                                                 0               1
   \endverbatim

   \note       Upon entry, the boundary edges in \p data.connectB need to be ordered the same way 
               in the case of several subdomains, or when internal closed lines (such as holes) are prescribed.
               However, edges don't need to be consecutive in the matrix.
   \note       Upon exit, the elements are always oriented counter-clockwise, i.e. with the normal&nbsp;=&nbsp;Oz.
   \note       Thread safe.
   \warning    The third coordinates in matrix \p data.pos is discarded for all points.
               That means that all points are projected onto the Z=0 plane (even for points
               not referenced in \p data.connectB or \p data.isolated_nodes).
   \see        meshtools2d::mesh_struct_T3, quadmesh::mesher, triamesh_aniso::mesher.
   */
class CM2TRIAMESH_API mesher
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
      The settings of triamesh::mesher.

      This structure gathers flags and parameters to control the way the triamesh::mesher works.
      */
   struct CM2TRIAMESH_API settings_type
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
      Basic mode for the %mesher.
      */
   enum basic_mode_type
   { 
      MESH_MODE,                //!< Normal mode of the %mesher (generates a new mesh).
      REGULARIZE_MODE,          //!< Takes an input mesh and optimizes it (nodes are added, removed, smoothed, connectivity is changed).
      CONVEX_HULL_MODE          //!< Produces an approximate convex hull of a set of points.
   };


   /**
      Default constructor.

      The default values are:
         - \ref basic_mode = \ref MESH_MODE 
         - \ref strict_constraints_flag = `true`
         - \ref subdomains_forcing = 0
         - \ref refine_flag = `true`
         - \ref use_default_background_mesh_flag = `false`
         - \ref node_smoothing_flag = `true`
         - \ref node_inserting_flag = `true`
         - \ref node_removing_flag = `true`
         - \ref shell_remeshing_flag = `true`
         - \ref no_clamped_edges_flag = `false`
         - \ref compute_Qh_flag = `false`
         - \ref structured_pattern = -1
         - \ref multi_structured_flag = `false`
         - \ref nodes_limit = `UINT_MAX`
         - \ref optim_level = 3
         - \ref target_metric = 0.0
         - \ref max_gradation = 0.5
         - \ref shape_quality_weight = 0.60
         - \ref length_upper_bound = 1.414
         - \ref display_hdl = `NULL`
         - \ref interrupt_hdl = `NULL`
         - \ref pass_thru = `NULL`
         - \ref progress_start = 0.0
         - \ref progress_range = 1.0
      */
   settings_type()
      : basic_mode(MESH_MODE), 
        strict_constraints_flag(true), subdomains_forcing(0), refine_flag(true), 
        use_default_background_mesh_flag(false), 
        node_smoothing_flag(true), node_inserting_flag(true), node_removing_flag(true), 
        shell_remeshing_flag(true), no_clamped_edges_flag(false), 
        compute_Qh_flag(false), structured_pattern(-1), multi_structured_flag(false),
        nodes_limit(UINT_MAX), optim_level(3),
        target_metric(0.0), max_gradation(0.5),
        shape_quality_weight(0.60), length_upper_bound(1.414),
        display_hdl(NULL), interrupt_hdl(NULL), pass_thru(NULL),
        progress_start(0.0), progress_range(1.0)
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
                     (\ref basic_mode is the first field, \ref strict_constraints_flag is the second...)
      */
   int
   check() const;

   /**
      Resets to the \link cm2::triamesh::mesher::settings_type::settings_type() default values \endlink.
      */
   void
   reset();


   /**
      Basic parameters of the %mesher: \ref MESH_MODE, \ref REGULARIZE_MODE or \ref CONVEX_HULL_MODE.

      Default = \ref MESH_MODE.
      */
   basic_mode_type         basic_mode;

   /**
      Flag for strict enforcement of hard edges and hard nodes.

      If true, an error is raised in the following cases:
         -  Two hard nodes are coincident (same coordinates).
         -  Two hard edges cut each other.
         -  A hard node is inside a hard edge.

      If false, the above cases raise only a warning and one of the erroneous item is discarded.

      \note    When strict_constraints_flag = false, the lowest edges in 
               data_type::connectB are privileged for enforcement.
      \note    When a hard node (resp. edge) is outside the meshed domain
               only a data_type::CM2_NODE_DISCARDED (resp. data_type::CM2_EDGE_DISCARDED) warning is issued,
               even in strict-constraints mode. 
      \note    Used only in \ref MESH_MODE.

      Default = true.
      */
   bool                    strict_constraints_flag;

   /**
      Forcing mode for intern subdomains.
      
      Possible values are:
         -  -1 : All intern subdomains are considered as holes regardless of the orientation of
                 their enclosing boundary.
         -   0 : Intern subdomains are considered as holes when the orientation of
                 their enclosing boundary is opposite to the orientation of the most enclosing domain.\n
                 Otherwise, they are meshed.
         -  +1 : All intern subdomains are meshed regardless of the orientation of
                 their enclosing boundary.

      \note    Used only in \ref MESH_MODE.

      Default = 0.
      */
   int                     subdomains_forcing;

   /**
      Flag to enable refinement. 

      If false, triangulates with the hard nodes only.

      If true, adds new nodes inside the domain to achieve good shape and size qualities.

      \note    Not used in \ref CONVEX_HULL_MODE.
      \note    In \ref REGULARIZE_MODE, unsetting this flag is equivalent to unset flag 
               \ref node_inserting_flag.

      Default = true.
      */
   bool                    refine_flag;

   /**
      Flag to tell the %mesher to use the default background mesh for the interpolation
      of the metric values. 

      For expert users only.

      The default background mesh is the "front" mesh, i.e. the mesh of the triangulation of the hard nodes
      without any new node.

      If set to false, the background mesh used to interpolate the metrics is the mesh currently processed.

      \warning A user-defined background mesh supersedes the default background mesh, even
               with this flag set.
      \note    Using the current mesh is the fastest mode which indeed gives a similar accuracy,
               so it is recommended to leave this flag to false.
      \note    For expert users only.

      Default = false.
      */
   bool                    use_default_background_mesh_flag;

   /**
      Flag to allow the node-smoothing scheme in the optimization step. 

      It can be useful in \ref REGULARIZE_MODE to set this flag to false
      to freeze the positions of the nodes.

      \note    Node smoothing doesn't change the mesh connectivity, only the coordinates of nodes.
      \note    This option, as part of the optimization step, is taken into account only when \ref optim_level > 0.
      \note    Used in \ref MESH_MODE and \ref REGULARIZE_MODE.

      Default = true.
      */
   bool                    node_smoothing_flag;

   /**
      Flag to allow the node-inserting scheme in the optimization step. 
      
      \note    Node inserting increases the number of nodes, changes the mesh 
               connectivity, but doesn't change the other nodes' coordinates.
      \note    This option, as part of the optimization step, is taken into account only when \ref optim_level > 0.
      \note    Used in \ref MESH_MODE and \ref REGULARIZE_MODE.

      Default = true.
      */
   bool                    node_inserting_flag;

   /**
      Flag to allow the node-removing scheme in the optimization step. 

      \note    Node removing decreases the number of nodes, changes the mesh 
               connectivity, but doesn't change the other nodes' coordinates.
      \note    This option, as part of the optimization step, is taken into account only when \ref optim_level > 0.
      \note    Used in \ref MESH_MODE and \ref REGULARIZE_MODE.

      Default = true.
      */
   bool                    node_removing_flag;

   /**
      Flag to allow the edge-swapping scheme in the optimization step. 

      \note    Edge swapping changes the mesh connectivity,
               but doesn't change the number of nodes nor their coordinates.
      \note    This option, as part of the optimization step, is taken into account only when \ref optim_level > 0.
      \note    Used in \ref MESH_MODE and \ref REGULARIZE_MODE.

      Default = true.
      */
   bool                    shell_remeshing_flag;

   /**
      Flag to force the removing of all free edges with two hard nodes.
      
      This triggers a special node-inserting scheme dedicated to breaking the non-hard edges with two hard nodes.

      \note    This option supersedes \ref node_inserting_flag. \n
               It increases the number of nodes, changes the mesh 
               connectivity, but doesn't change the other nodes' coordinates.
      \note    Setting this flag to true increases the number of generated nodes and elements, and can
               reduce the global quality of the mesh (elements' shape quality and elements' size).
      \note    This option, as part of the optimization step, is taken into account only when \ref optim_level > 0.

      Default = false
      */
   bool                    no_clamped_edges_flag;

   /**
      Flag to compute the edge quality histogram (quality in term of edge length).

      \note    Not used in \ref CONVEX_HULL_MODE.

      Default = false.
      */
   bool                    compute_Qh_flag;

   /**
      Meshing mode on rectangular domains (unstructured or structured meshes).

      Enables structured meshing (grid-like) whenever possible and beneficial.

      If true, generates a structured mesh (grid-like) whenever possible, i.e. when:
         - The domain is topologically rectangular (four sharp corners).
         - The discretizations of the four lines match 
           (same number of edges along opposite lines).
         - The domain is roughly convex.
         - There is no isolated hard edge or node (even located outside the domain).
         - There is no repulsive point (even located far outside the domain).
         - There is no valid \ref target_metric specified.

      If false, uses the unstructured Delaunay-based algorithm.

      Three types of structured meshes (grid-like meshes):
      \verbatim
         - Left pattern:

               *------*------* 
               | \    | \    |
               |   \  |   \  |
               |     \|     \|
               *------*------*
               | \    |\     |
               |   \  |  \   |
               |     \|    \ |
               *------*------*

         - Right pattern:

               *------*------* 
               |     /|     /|
               |   /  |   /  |
               | /    | /    |
               *------*------*
               |     /|     /|
               |   /  |   /  |
               | /    | /    |
               *------*------*


         - "Union Jack" pattern:        

               *------*------* 
               | \    |    / |
               |   \  |  /   |
               |     \|/     |
               *------*------*
               |     /|\     |
               |   /  |  \   |
               | /    |    \ |
               *------*------*
      \endverbatim

      Possible values are:
         -  -1    : Always use the Delaunay-frontal algorithm (generates unstructured meshes).
         -   0    : When possible, generates structured left-oriented meshes (simply oriented pattern).
         -   1    : When possible, generates structured right-oriented meshes (simply oriented pattern).
         -  >= 2  : When possible, generates structured UJ meshes ("Union Jack" pattern).

      \note    Used only when \ref refine_flag = true.
      \note    Used in \ref MESH_MODE only.

      Default = -1 (always unstructured mesh).
      */
   int                     structured_pattern;

   /**
      Flag to trigger the multi-structured-part analysis.

      This analysis performs a more complete search for structurable subdomains than normally
      (at the price of a lower speed).

      With this flag on, several rectangular subdomains can be meshed structurably, whereas it is not
      always the case with this flag down.

      \note    Used only when \ref refine_flag = true and \ref structured_pattern != -1.
      \note    Used in \ref MESH_MODE only.
      \note    Disabled with isolated nodes or background mesh.
      \warning Discard embedded hard edges in some cases 
               (subdomains > 1, dangling edges with reversed edges).

      Default = false.
      */
   bool                    multi_structured_flag;

   /**
      Maximum authorized number of nodes.

      Default = `UINT_MAX`.
      */
   unsigned                nodes_limit;

   /**
      Optimization level between 0 and 10.

      Controls the trade-off between speed (low values) and high element quality (high values). \n
      The higher the level, the more aggressive the algorithms.

      \note    The algorithms involved here aim at optimizing both the shape quality of the elements _and_ their size 
               quality (length of elements' edges). \n
               Hence, increasing the value of optim_level may lead sometimes to a
               reduction in the shape qualities (min, average or max value) whereas the _global_ quality (which takes into 
               account both types of quality) gets better. 
      \note    The \ref shape_quality_weight parameter controls the weight of the shape quality in the global 
               quality measure.
      \note    Level 0 disables any optimization.

      Default = 3
      */
   unsigned                optim_level;

   /**
      Element size inside the domain.

      The elements size tends toward this value as they get away from the hard (boundary) edges. 
      
      \warning Discarded when <= 0.
      \warning Discarded when a structured mesh can be done.
      \warning Discarded when a background mesh is used.
      \warning Discarded when \ref use_default_background_mesh_flag = true.
      \note    This parameter is often used to reduce the total number of elements
               (when target_metric > default sizes based on boundary edge lengths) 
               and to save FEM computation time without loosing much on element shape quality.
      \note    The target_metric may be smaller or bigger than the default sizes based on boundary edge lengths.
      \note    target_metric >= 0.
      \note    Used in \ref MESH_MODE only.
      \see     max_gradation, structured_pattern.

      Default = 0 (disabled)
      */
   double                  target_metric;

   /**
      This parameter controls the gradation of the elements size from the boundary sizes 
      to the size defined by \ref target_metric inside the domain. 

      A value close to 0 leads to a more progressive variation of mesh size (smoother).

      \warning Discarded when \ref target_metric <= 0.
      \warning Discarded when a structured mesh is done (see \ref structured_pattern).
      \warning Discarded when a background mesh is used.
      \warning Discarded when \ref use_default_background_mesh_flag = true.
      \note    max_gradation >= 0.
      \note    Used in \ref MESH_MODE only.
      \see     target_metric, structured_pattern.

      Default = 0.5
      */
   double                  max_gradation;

   /**
      Weight on shape qualities versus size qualities.

      Scalar between 0 and 1 indicating the importance of the shape qualities versus
      the size quality for the regularization step:
         - For values greater than 0.5, the shape quality is privileged.
         - For values lesser than 0.5, the size quality is privileged.

      Default = 0.6 (slight preference for shape quality). Should not normally be changed.

      \pre        shape_quality_weight >= 0.
      \pre        shape_quality_weight <= 1.
      */
   double                  shape_quality_weight;

   /**
      Upper bound value for the edge lengths in the final mesh.

      The mesher tends to generate meshes with no edge longer than this value (normalized length).
      But this is _not_ a strict enforcement. Most of the edges will be shorter than this limit, 
      but some may remain somewhat longer.

      The default value (1.414) gives the optimal meshes with respect to the size qualities.
      With this default value (no edge more than 41% longer than the perfect length), 
      the average edge length tends to be 1.

      Sometimes it can be useful to limit the length of the edges to a shorter value (usually between 1 and 1.414),
      and to accept a mean value smaller than 1 (sub-optimal edge qualities). 

      \note    The value is _relative_ to the perfect edge length.

      Default = 1.414
      */
   double                  length_upper_bound;

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

      When this function returns true, the process is stopped. \n
      This can abort the mesh (if the interruption occurred at an early stage),
      or result in a mesh of poor quality (if the %mesher was in the refinement or optimization 
      stage when the interruption occurs).

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

   /**
      The start value for the progress parameter sent to the interrupt handler.

      Default = 0.
      */
   double                  progress_start;

   /**
      The range for the progress parameter sent to the interrupt handler.

      The end value equals to progress_start + progress_range.

      Default = 1.0
      */
   double                  progress_range;

   };    // settings_type



   /**
      The data type for triamesh::mesher.

      This structure is used to gather all the input and output data of the function 
      triamesh::mesher::run.
      */
   struct CM2TRIAMESH_API data_type
   {

   /// The display handler type.
   typedef settings_type::display_handler_type        display_handler_type;


   /**
      Warning enums for the %mesher.

      Non fatal warnings.

      If such a warning is raised (returned in the field \ref warning_code),
      the %mesher is not stopped but the user shall check its input and output data for a potential
      mistake.
      */
   enum warning_type
   { 
      CM2_NO_WARNING                =   0,   //!< OK, no warning.
      CM2_INTERRUPTION              = -10,   //!< The process has been interrupted (by the user through the interrupt handler).
      CM2_NODES_LIMIT_WARNING       = -11,   //!< The limitation on the number of nodes has been reached.
      CM2_NODE_DISCARDED            = -12,   //!< Two nodes are too close to each other (same position or almost). 
      CM2_EDGE_DISCARDED            = -13,   //!< A hard edge is cut by another hard edge, or has a hard node inside it.  
      CM2_BOUNDARY_WARNING          = -14,   //!< The outermost boundary is not watertight or is not properly oriented.
      CM2_SHAPE_QUALITY_WARNING     = -15    //!< One or several elements have a bad quality (Qs < 0.01).
   };

   /**
      Error enums for the %mesher.

      Fatal errors.

      If such an error is raised (returned in the field \ref error_code),
      the %mesher is stopped and no final mesh is returned.
      */
   enum error_type
   { 
      CM2_NO_ERROR                  =    0,  //!< OK, no error.
      CM2_LICENSE_ERROR             = -100,  //!< Invalid or expired license.
      CM2_MODE_ERROR                = -101,  //!< Error in settings.
      CM2_DATA_ERROR                = -102,  //!< Error in input data (`NaN`, `INF`, error in dimensions of arrays...)
      CM2_NODES_LIMIT_ERROR         = -103,  //!< The limitation on the number of nodes is too tight.
      CM2_NODE_ERROR                = -104,  //!< Two nodes are too close to each other (same position or almost) or some are located outside. 
      CM2_EDGE_ERROR                = -105,  //!< A hard edge is cut by another hard edge, or has a hard node inside it, or is located outside.  
      CM2_BOUNDARY_ERROR            = -106,  //!< Illegal boundary(ies) provided by the user (boundary not closed, boundary Xings...) 
      CM2_DEGENERATED_ELEMENT       = -107,  //!< The process leads to at least one degenerated element and therefore the %mesher cannot build a valid mesh.
      CM2_BACKGROUND_MESH_ERROR     = -108,  //!< The background mesh is not valid (bad dimensions, bad indices, bad associated metrics, degenerated elements...)
      CM2_SYSTEM_MEMORY_ERROR       = -199,  //!< Insufficient memory available.
      CM2_INTERNAL_ERROR            = -200   //!< Something went wrong but don't know what (please contact support@computing-objects.com).
   };


   /**
      The type of metric.

      This mesh generator is _isotropic_ and the metrics are scalars.
     */
   typedef double             metric_type;

   /**
      The array type for the metric values.

      This mesh generator is _isotropic_ and the metrics are represented as single scalars (`double` values).
      */
   typedef DoubleVec          metrics_type;



   /// Default constructor.
   data_type()                            { this->clear(); }

   /**
      Copy constructor (hard copy).
      All matrices are hard copied.
      */
   data_type (const data_type& cpy)       { (*this) = cpy; }

   /**
      Constructor.

      Constructs the data with a coordinates matrix and a boundary matrix of edges (shallow copy).

      This is the simplest way to create the data to exchange with the %mesher.

      \param   P           The coordinates matrix. New nodes are appended by the %mesher.
      \param   connectE2   The connectivity matrix of the boundary mesh.

      \note    After the meshing process, you must get back the new coordinate matrix from the field \ref pos,
               and the connectivity matrix \ref connectM as well, or use the \ref extract function.
      \warning Upon exit, the matrix \ref pos is resized to two rows. \n
               The third coordinate (i.e. Z) is lost for all points.
      */
   data_type (const DoubleMat& P, const UIntMat& connectE2)
   {
      this->clear();
      this->pos = P;
      this->connectB = connectE2;
   }


   /**@name Copy members */
   //@{
   /// Copy operator (hard copy).
   const data_type&
   operator= (const data_type& data);

   /// Shallow copy.
   void
   shallow_copy (const data_type& data);
   //@}


   /**
      Returns true if the sanity check is successful.
      */
   bool
   valid() const        { return this->check() == 0; }

   /**
      Sanity check of the input fields.

      \return        Returns a negative value in case of error. \n
                     Returned value = -k => the k-th field of the structure is illegal
                     (\ref pos is the first field, \ref connectB is the second...)
      */
   int
   check() const;

   /// Shallow copy.
   void
   extract (DoubleMat& P, UIntMat& connectT3) const
      { P = this->pos; connectT3 = this->connectM; }

   /// Shallow copy. Provided for compatibility with quadmesh::mesher.
   void
   extract (DoubleMat& P, UIntMat& connectQ4, UIntMat& connectT3) const
      { P = this->pos; connectQ4.clear(); connectT3 = this->connectM; }

   /// Clears all fields and reinit to default values.
   void
   clear();

   /// Clears only the output fields.
   void
   clear_data_out();

   /**
      Gets the order of elements according to their color.

      Simply calls meshtools::get_sorted_order on \ref colors. \n
      Provided for convenience only.

      Example:
      \code
         triamesh::mesher              my_mesher;
         triamesh::mesher::data_type   my_data(pos, connectB);  
         UIntVec                       color_order, xColor;
         UIntVec                       IDs;
         UIntMat                       connect_7;

         // Do the meshing.
         my_mesher.run(my_data);               

         // Get element IDs sorted by color.
         my_data.IDs_sorted_by_color(color_order, xColor);             

         // Get the elements with color 7.
         IDs = color_order.sub_vector(xColor[7], xColor[8]-xColor[7]);     // All the element IDs with color 7.
         matmat::push_back(my_data.connectM, IDs, connect_7);              // connect_7 = connectivity matrix of the elements with color 7.
      \endcode
      */
   void
   IDs_sorted_by_color (UIntVec& color_order, UIntVec& xColor) const
   {
      meshtools::get_sorted_order(this->colors, color_order, xColor);
   }

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

   /**
      Load the data from a file.
      This member is the reciprocal of the save member above.

      \param   filename       The file name for the input.

      \return  Error code (negative value = error).
      */
   int
   load (const char* filename);

   /// Checks if an error has been raised.
   inline bool
   error_raised() const       { return this->error_code != CM2_NO_ERROR; }

   /// Checks if a warning has been raised.
   inline bool
   warning_raised() const     { return this->warning_code != CM2_NO_WARNING; }


   /**
      The coordinates matrix (column oriented).

      Node `i` has coordinates in column `pos.col(i)`.

      Upon entry, this matrix must contains the coordinates of all the points (i.e. nodes referenced in 
      \ref connectB, \ref isolated_nodes, \ref background_mesh and \ref repulsive_points). 

      It may also contain other unreferenced points.

      Upon exit, the new nodes are appended at the back of the matrix and the indices in the
      final mesh (\ref connectM) refer to column numbers in this matrix.

      Mode = IN-OUT in settings_type::MESH_MODE and settings_type::REGULARIZE_MODE, \n
             IN in settings_type::CONVEX_HULL_MODE.

      \warning    Upon exit, the matrix is resized to two rows. \n
                  The third coordinate (Z) is lost for all points.
      */
   DoubleMat               pos;

   /**
      The hard edges (boundary edges + internal hard edges).

      This connectivity matrix must define a valid and closed edge mesh for the 
      external boundary of the domain.

      The orientation of the boundary edges must be similar (either all CW or all CCW). \n 
      However, mixed orientations are accepted when:
         - There is only a single external contour.
         - settings_type::subdomains_forcing = +1. \n
           In this case, all subdomains are filled with elements (no hole).

      Internal (embedded) hard edges can be added to the boundary mesh. 

      In order to define holes inside the domain, the boundary mesh of the holes must be closed
      and oriented oppositely to the external boundary mesh 
      (and \ref settings_type::subdomains_forcing = 0 or -1).

      The indices refer to columns in the coordinates matrix \ref pos.

      \note    In settings_type::REGULARIZE_MODE, 
               external edges and edges between elements with different colors are added.\n
               These new edges, if any, are also considered as hard entities.
      \note    In settings_type::CONVEX_HULL_MODE, only the nodes in \ref connectB 
               and \ref isolated_nodes are taken into account (the edges are not forced to
               be present in the final mesh).
      \see     \ref settings_type::subdomains_forcing, \ref boundary_sgn.
      \see     meshtools::get_mesh_boundaries, meshtools::get_colors_boundaries.

      \pre     \p connectB is empty or \p connectB.rows() == 2 (edge mesh).

      Mode = IN in settings_type::MESH_MODE and settings_type::CONVEX_HULL_MODE, \n
             IN-OUT in settings_type::REGULARIZE_MODE.
      */
   UIntMat                 connectB;

   /**
      The internal (embedded) hard nodes.

      \note    The indices refer to columns in the coordinates matrix \ref pos.
      \note    Isolated nodes (even located outside of the domain) disable 
               the structured meshing algorithm.

      Mode = IN.
      */
   UIntVec                 isolated_nodes;

   /**
      The nodes of the mesh (that are not hard nodes) will be placed away from
      these repulsive points.

      \note    The indices refer to columns in the coordinates matrix \ref pos.
      \note    Repulsive points disable structured meshing (even when located outside of the domain).
      \note    Only repulsive points associated to a valid metric (i.e. positive value in
               \ref metrics) are taken into account.
      \note    Not used in settings_type::CONVEX_HULL_MODE.

      Mode = IN.
      */
   UIntVec                 repulsive_points;

   /**
      The connectivity matrix of the background mesh, made of triangles (column oriented). 

      The indices of the nodes refer to columns in the coordinates matrix \ref pos. \n
      These indices must also refer to valid values in the \ref metrics array.

      The background mesh doesn't need to match exactly the domain to be meshed. \n
      It is allowed to extend the domain or to cover only a fraction of it. \n
      For subdomains outside of the background mesh, the standard interpolation of the metrics
      defined (or computed) at the hard nodes is used.

      \warning    The use of a background mesh may slow down the generator significantly. \n
                  The finer the background mesh, the slower the meshing.

      Mode = IN in settings_type::MESH_MODE and settings_type::REGULARIZE_MODE.
      */
   UIntMat                 background_mesh;

   /**
      Element sizes specified at nodes.

      For nodes with invalid metric (missing, negative or null), an automatic size is 
      computed using the distances to their neighbours and stored in this vector upon exit.

      \note    Upon exit, this array takes into account the values of settings_type::target_metric 
               and settings_type::max_gradation.
      \note    Not used in settings_type::CONVEX_HULL_MODE.

      Mode = IN-OUT.
      */
   DoubleVec               metrics;

   /**
      The connectivity of the elements (column oriented). 

      The indices refer to columns in the coordinates matrix \ref pos.

      The dimensions of this matrix are 3 x \ref nefs.

      Mode = OUT in settings_type::MESH_MODE and settings_type::CONVEX_HULL_MODE, \n
             IN-OUT in settings_type::REGULARIZE_MODE.
      */
   UIntMat                 connectM;

   /**
      The boundaries that could not be enforced in the mesh.

      This vector gives the column IDs in the \ref connectB matrix that could not be enforced in
      the mesh (i.e. edges intersected by other hard nodes or hard edges, edges outside the domain).

      \note    In strict-constraints mode (settings_type::strict_constraints_flag = true), 
               the generator stops with an error only for intersecting entities, not for
               entities outside the meshed domain.
      \note    Used in settings_type::MESH_MODE only.

      Mode = OUT.
      */
   UIntVec                 unenforced_boundary_IDs;        

   /**
      The nodes that could not be enforced in the mesh.

      This vector gives the indices of the nodes that could not be enforced in
      the final mesh (i.e. nodes coincident with another one, nodes inside hard edges, nodes outside the domain).

      \note    In strict-constraints mode (settings_type::strict_constraints_flag = true), 
               the generator stops with an error only for intersecting entities, not for
               entities outside the meshed domain.
      \note    Used in settings_type::MESH_MODE only.

      Mode = OUT.
      */
   UIntVec                 unenforced_node_IDs;        

   /**
      The boundaries that intersect other boundary(ies).

      This vector gives the column IDs in the \ref connectB matrix of the edges that intersect other hard entities.

      \note    In strict-constraints mode (settings_type::strict_constraints_flag = true), 
               the generator stops with an error when this vector is not empty.
      \note    Used in settings_type::MESH_MODE only.

      Mode = OUT.
      */
   UIntVec                 pathological_boundary_IDs;        

   /**
      The "color" of the elements. 

      Each subdomain is assigned with a color ID (consecutive values starting from 0). \n
      Each element is assigned with the color ID of the subdomain it belongs to.

      A subdomain is a subset of connected elements (not separated by any boundary). \n
      Two elements are in the same subdomain when there is a element path between them 
      which is not crossed by any boundary.
      
      Upon exit, the size of this vector equals to the number of elements in \ref connectM.

      The indices refer to columns in the connectivity matrix \ref connectM.

      In settings_type::REGULARIZE_MODE, this vector is used upon entry to search for additional
      hard boundaries (boundaries between elements of different colors). \n
      The user can also provide an empty color vector. In this case, all elements 
      are assumed to be in the same subdomain (color = 0).

      In settings_type::CONVEX_HULL_MODE, all elements are in the same subdomain (color = 0).

      Mode = OUT in settings_type::MESH_MODE and settings_type::CONVEX_HULL_MODE, \n
             IN-OUT in settings_type::REGULARIZE_MODE.
      */
   UIntVec                 colors;

   /**
      The neighbors matrix. 

      `neighbors(i,j)` is the neighbor of element `j` with edge `i`.

      \post    The dimensions of this matrix are 3 x \ref nefs.

      \see  The User's Manual.

      Mode = OUT.
      */ 
   UIntMat                 neighbors;

   /**
      The ancestors vector. 

      `ancestors[i]` is the index of one of the elements connected to node `i`.

      These indices refer to columns in matrices \ref connectM and \ref neighbors, and to position
      in vector \ref colors.

      If a point is not connected to any element, its ancestors' index equals to `CM2_NONE` (i.e. `unsigned(-1)`).

      \post    The size of this vector equals to `pos.cols()` (i.e. number of points).

      Mode = OUT.
      */ 
   UIntVec                 ancestors;

   /**
      The shape qualities vector. 

      The size of this vector is \ref nefs.

      \see  The User's Manual.

      Mode = OUT.
      */ 
   DoubleVec               shape_qualities;

   /**
      The histogram of the element shape qualities.

      Mode = OUT.
      */ 
   misc::histogram         histo_Qs; 

   /**
      The histogram of the edges qualities (edge length).

      Computed only when settings_type::compute_Qh_flag = true.

      Mode = OUT.
      */ 
   misc::histogram         histo_Qh;

   /**
      Number of elements in the final mesh. 

      Mode = OUT.
      */ 
   unsigned                nefs;

   /**
      Number of quads in the final mesh. Always null. 

      Provided for compatibility with quadmesh::mesher.

      Mode = OUT.
      */ 
   unsigned                nefs_Q4;

   /**
      Number of triangles in the final mesh. Same as \ref nefs.

      Provided for compatibility with quadmesh::mesher.

      Mode = OUT.
      */ 
   unsigned                nefs_T3;

   /**
      Number of nodes in the final mesh. 

      Mode = OUT.
      */ 
   unsigned                nods;

   /**
      Total number of points in the coordinates matrix.

      \post Greater or equal to \ref nods.
      \post Equal to `ancestors.cols()`.
      \post Equal to `metrics.size()`.

      Mode = OUT.
      */ 
   unsigned                total_nods;

   /**
      The number of initial hard nodes.

      Mode = OUT.
      */ 
   unsigned                hard_nodes_in;

   /**
      The number of initial hard edges. 

      Mode = OUT.
      */ 
   unsigned                hard_edges_in;

   /**
      The number of final hard nodes (lesser than or equal to \ref hard_nodes_in). 

      Mode = OUT.
      */ 
   unsigned                hard_nodes_out;

   /**
      The number of final hard edges (lesser than or equal to \ref hard_edges_in). 

      Mode = OUT.
      */ 
   unsigned                hard_edges_out;

   /**
      Area of the meshed domain (sum of the areas of all elements).

      Mode = OUT.
      */ 
   double                  area;

   /**
      Area of the quadrangles. Always null.

      Provided for compatibility with quadmesh::mesher.

      Mode = OUT.
      */ 
   double                  area_Q4;

   /**
      Area of the triangles. Same as \ref area.

      Provided for compatibility with quadmesh::mesher.

      Mode = OUT.
      */ 
   double                  area_T3;

   /**
      Minimum shape quality of the elements upon entry. 

      Mode = OUT in settings_type::REGULARIZE_MODE only.
      */ 
   double                  Qmin_in;

   /**
      Minimum shape quality of the elements upon exit. 

      Mode = OUT.
      */ 
   double                  Qmin;

   /**
      The bounding box's min corner (lower-left). 

      Mode = OUT.
      */ 
   DoubleVec2              min_box;

   /**
      The bounding box's max corner (upper-right). 

      Mode = OUT.
      */ 
   DoubleVec2              max_box;

   /**
      Number of subdomains in the domain.

      Mode = OUT.
      */ 
   unsigned                subdomains;

   /**
      The orientation of the external contour.

      The orientation sign is defined by:
         -  +1 when all edges of the external contour are oriented counter-clockwise.
         -  -1 when all edges of the external contour are oriented clockwise.
         -  0  when some edges are oriented counter-clockwise and some are oriented clockwise 
               (i.e. the external contour orientation is not well defined).

      \note    If the orientation sign is null, the mesher has undefined behaviour with
               respect to any internal/external closed contour (may be considered as hole or not),
               except when forced (by using \ref settings_type::subdomains_forcing).
      \note    This parameter is computed in settings_type::MESH_MODE only.
      \see     \ref settings_type::subdomains_forcing.

      Mode = OUT.
      */ 
   int                     boundary_sgn;

   /**
      Time spent to do the initial front mesh. 

      Mode = OUT.
      */ 
   double                  front_time;

   /**
      Time spent for refining the mesh. 

      Mode = OUT.
      */ 
   double                  refine_time;

   /**
      Time spent for optimizing the mesh.

      Mode = OUT.
      */ 
   double                  optim_time;

   /**
      Total time spent for the mesh.

      Mode = OUT.
      */ 
   double                  total_time; 

   /**
      Number of elements per second. 

      Global speed including all steps of the meshing/optimizing process.

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

   See \link cm2::triamesh::mesher::settings_type::settings_type() default settings \endlink.
   */
mesher();

/**
   Constructor with specific settings.

   \param[in]   settings_        The settings. 
   */
mesher (const settings_type& settings_);

/**
   Is the %mesher authorized to run?

   \see  triamesh::registration.
   */
static bool
authorized();

/**
   Main function of the %mesher.

   Depending on the settings, the %mesher can work
   as a regular %mesher or as a mesh optimizer.

   \param[in,out]   data   The data.
   */
void 
run (data_type& data) const;


/**
   The settings of the %mesher.

   See \link cm2::triamesh::mesher::settings_type::settings_type() default settings \endlink.
   */
settings_type         settings;

};

}  // end namespace triamesh
}  // end namespace cm2



/** 
   \example tmsh2d01.cpp
    This is an example of how to use the 2-D isotropic unstructured meshers cm2::triamesh::mesher and cm2::quadmesh::mesher.

    \image html tmsh2d01_t.png "Simple square with 40 boundary edges, CM2 TriaMesh."
    \image html tmsh2d01_q.png "Simple square with 40 boundary edges, CM2 QuadMesh (same result with all_quad_flag = true or false)." 
 */

/** 
   \example tmsh2d03.cpp
    This is an example of how to prescribe a point with the 2-D meshers cm2::triamesh::mesher and cm2::quadmesh::mesher. \n
    The field `data.metrics` is also used to show how to get a finer mesh near this point.

    \image html tmsh2d03_t.png "Simple square with input 40 boundary edges and one embedded node at centre with small size, CM2 TriaMesh." 
    \image html tmsh2d03_qt.png "Simple square with input 40 boundary edges and one embedded node at centre with small size, CM2 QuadMesh (all_quad_flag = false)."
    \image html tmsh2d03_q.png "Simple square with input 40 boundary edges and one embedded node at centre with small size, CM2 QuadMesh (all_quad_flag = true)."
 */

/** 
   \example tmsh2d34.cpp
    This is an example of how to use the background mesh option with the 2-D meshers cm2::triamesh::mesher and cm2::quadmesh::mesher.

    \image html tmsh2d34_t.png "Square with metrics defined on a background mesh, CM2 TriaMesh." 
    \image html tmsh2d34_qt.png "Square with metrics defined on a background mesh, CM2 QuadMesh (all_quad_flag = false)." 
 */

#endif
