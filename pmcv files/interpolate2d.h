/**
   \file       interpolate2d.h
   \brief      Declarations of some interpolation routines on 2-D meshes.
   \copyright  (C)1999-2017, Computing Objects, France. info@computing-objects.com

   $Rev: 2963 $
   $Date: 2017-12-15 12:24:31 +0100 (ven., 15 déc. 2017) $ 
   */

/******************************************************************************************
   This program is not free software. It is the sole property of Computing Objects, France.
   You are allowed to modify it, to recompile it, to port it on several platforms,
   but not to redistribute it, even the modified version.
 ******************************************************************************************/

#ifndef __CM2_INTERPOLATE_2D_H__
#define __CM2_INTERPOLATE_2D_H__


namespace cm2 {
namespace meshtools2d {

/****************************************************************************************/
/** @name Fields interpolation on 2-D triangle meshes */
//@{

/**
   Interpolates a scalar field (doubles) defined on the nodes of a 2-D triangle mesh.

   The interpolated values are computed on a specific set of nodes that are not in the mesh.

   \param[in]     pos            The 2-D coordinates matrix (2xNODS). 
   \param[in]     connectT3      The connectivity matrix of the triangle mesh (3xNEFS).
                                 Each node of the mesh should have a valid associated field value.
   \param[in]     neighbors      The neighbors matrix (3xNEFS).
                                 Can be computed by the meshtools::get_neighbors function.
   \param[in]     ancestors      The ancestors vector. For each node referenced in \p connectT3,
                                 the vector gives one triangle connected to it.
                                 For nodes not in \p connectT3, the ancestor should be `CM2_NONE` (i.e. `unsigned(-1)`).
                                 Can be computed by the meshtools::get_ancestors function.
   \param[in,out] field          The scalar field (vector of size NODS).
   \param[in]     nodes          The node IDs for which the field values must be computed.
                                 Values for nodes of the mesh are not recomputed.

   \return        Error code. Returns 0 when successful.
                  Returned value = -k => the k-th argument had an illegal value.
   
   \pre        pos.rows() == 2.
   \pre        connectT3.rows() == 3 (i.e. triangle mesh).
   \pre        field.size() == pos.cols().
   */
CM2MESHTOOLS2D_API int 
interpolate 
   (const DoubleMat& pos, const UIntMat& connectT3, const UIntMat& neighbors, const UIntVec& ancestors,
    DoubleVec& field, const UIntVec& nodes);                          


/**
   Interpolates a vectorial field (doubles) defined on the nodes of a 2-D triangle mesh.

   The interpolated values are computed on a specific set of nodes that are not in the mesh.

   \param[in]     pos            The 2-D coordinates matrix (2xNODS). 
   \param[in]     connectT3      The connectivity matrix of the triangle mesh (3xNEFS).
                                 Each node of the mesh should have a valid associated field vector.
   \param[in]     neighbors      The neighbors matrix (3xNEFS).
                                 Can be computed by the meshtools::get_neighbors function.
   \param[in]     ancestors      The ancestors vector. For each node referenced in \p connectT3,
                                 the vector gives one triangle connected to it.
                                 For nodes not in \p connectT3, the ancestor should be `CM2_NONE` (i.e. `unsigned(-1)`).
                                 Can be computed by the meshtools::get_ancestors function.
   \param[in,out] field          The vectorial field (matrix of dimensions FxNODS).
                                 Each component is interpolated independently from the others.
   \param[in]     nodes          The node IDs for which the field vectors must be computed.
                                 Vectors for nodes already in the mesh are not recomputed.

   \return        Error code. Returns 0 when successful.
                  Returned value = -k => the k-th argument had an illegal value.
   
   \note       Metrics cannot be interpolated with this function as the metric components are not 
               independent from each other. Use meshtools2d::interpolate_metrics instead.
   \pre        pos.rows() == 2.
   \pre        connectT3.rows() == 3 (i.e. triangle mesh).
   \pre        field.cols() == pos.cols().
   */
CM2MESHTOOLS2D_API int 
interpolate 
   (const DoubleMat& pos, const UIntMat& connectT3, const UIntMat& neighbors, const UIntVec& ancestors,
    DoubleMat& field, const UIntVec& nodes);                          
//@}


}  // end namespace meshtools2d
}  // end namespace cm2

#endif
