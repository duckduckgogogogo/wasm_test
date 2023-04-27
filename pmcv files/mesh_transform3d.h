/**
   \file       mesh_transform3d.h
   \brief      Declarations of some 3-D mesh transformation routines.
   \copyright  (C)1999-2017, Computing Objects, France. info@computing-objects.com

   $Rev: 2963 $
   $Date: 2017-12-15 12:24:31 +0100 (ven., 15 déc. 2017) $ 
   */

/******************************************************************************************
   This program is not free software. It is the sole property of Computing Objects, France.
   You are allowed to modify it, to recompile it, to port it on several platforms,
   but not to redistribute it, even the modified version.
 ******************************************************************************************/

#ifndef __CM2_MESH_TRANSFORM3D_H__
#define __CM2_MESH_TRANSFORM3D_H__



namespace cm2 {
namespace meshtools3d {


/****************************************************************************************/
/** @name Reorientation of 3-D elements */
//@{

/**
   Changes the orientation of the solid elements in a 3-D mesh 
   (cm2::CM2_TETRA4, cm2::CM2_PYRAMID5, cm2::CM2_WEDGE6, cm2::CM2_HEXA8, cm2::CM2_SOLID_MIX,
    cm2::CM2_TETRA10, cm2::CM2_PYRAMID14, cm2::CM2_WEDGE18, cm2::CM2_HEXA27)

   This flip the orientation (i.e. the volume) of the elements.
   Applied twice does not change anything.

   \param[in]     connect     The connectivity matrix.
   \param[in]     FE_type     The type of the elements (see \ref fe_geo_page).

   \return        Error code. Returns 0 when successful.
                  Returned value = -k => the k-th argument had an illegal value.

   \note          Performed swaps:
                     - cm2::CM2_TETRA4 and cm2::CM2_PYRAMID5: swaps the local nodes #0 and #2.
                     - cm2::CM2_WEDGE6: swaps the local nodes #0 and #2, and #3 and #5.
                     - cm2::CM2_HEXA8: swaps the local nodes #0 and #2, and #4 and #6.
                     - High-order elements (cm2::CM2_TETRA10, cm2::CM2_PYRAMID15, cm2::CM2_WEDGE18, cm2::CM2_HEXA27): swaps also the mid-edge nodes.
   \note          Works also with mixed connectivity matrices 
                  (i.e. \p connect.rows() == 8 and \p FE_type == cm2::CM2_SOLID_MIX) and
                  mixed high-order matrices
                  (i.e. \p connect.rows() == 27 and \p FE_type == cm2::CM2_SOLID_MIX).

   \pre           \p Number of rows in \p connect matches \p FE_type or \p connect is empty.

   \see           \ref fe_geo_page for the definition of the local nodes.
   */
CM2MESHTOOLS3D_API int
flip (UIntMat& connect, cm2::element_type FE_type);

//@}



/****************************************************************************************/
/** @name Special transformations of 3-D meshes */
//@{

/**
   Duplicates the nodes of a triangle mesh embedded inside a tetrahedral mesh,
   and change the connectivity of the tets over the surface so that they are connected
   to the new nodes.
   
   The new nodes are appended at the back of matrix \p pos.
   
   \param[in,out] pos               The coordinates matrix. New points are appended.
   \param[in,out] connectTH4        The connectivity matrix of the tetrahedral mesh (4xNEFS_TH4).
                                    Upon exit, this matrix holds the new mesh with reference to the new 
                                    duplicated nodes. The initial content of the matrix is lost.
   \param[in]     neighborsTH4      The neighbors matrix of the tetrahedral mesh (4xNEFS_TH4).
                                    This matrix can be computed using meshtools::get_neighbors.
   \param[in]     ancestors         The ancestors array (for each node, gives a tetra connected to it) (NODS).
                                    This array can be computed using meshtools::get_ancestors.
   \param[in]     connectT3         The connectivity matrix of the surface mesh (3xNEFS_T3).
   \param[in]     except_nodes      Nodes that must not be duplicated.
   \param[in]     upper_nodes_flag  Flag to tell if the nodes over the \p connectT3 surface should be duplicated.
   \param[in]     lower_nodes_flag  Flag to tell if the nodes under the \p connectT3 surface should be duplicated.
                                    If both flags are `true`, each node of the surface (but except_nodes)
                                    are "split" into three coincident nodes (two new nodes are created). 

   \return        Error code. Returns 0 when successful.
                  Returned value = -k => the k-th argument had an illegal value.

   \pre        pos.rows() == 3
   \pre        connectTH4.rows() == 4
   \pre        neighborsTH4.rows() == 4
   \pre        connectTH4.cols() == neighborsTH4.cols()
   \pre        matscal::max_value (connectT3) < ancestors.size()

   \note       A call to meshtools::merge resets the mesh to its initial state.
   */
CM2MESHTOOLS3D_API int 
duplicate_surface_nodes_TH4 
   (DoubleMat& pos, UIntMat& connectTH4, const UIntMat& neighborsTH4, 
    const UIntVec& ancestors, const UIntMat& connectT3, 
    const UIntVec& except_nodes, bool upper_nodes_flag, bool lower_nodes_flag);    


/**
   Splits the input solid elements having a Q4 face in the input list.

   \warning    Experimental only. Do not use.
   */
CM2MESHTOOLS3D_API int
split_solid_Q4
   (DoubleMat& pos, const UIntMat& connectM, const UIntMat& connectQ4,
    const UIntVec& diagonals, double min_qx, UIntVec& failed_indices, UIntMat& connectOut);

//@}



}  // end namespace meshtools3d
}  // end namespace cm2

#endif
