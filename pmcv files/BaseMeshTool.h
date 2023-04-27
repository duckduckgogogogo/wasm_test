#pragma once

#define VC_EXTRALEAN		// Exclude rarely-used stuff from Windows headers

#ifdef strcasecmp
#undef strcasecmp
#endif

#include <vector>

#ifdef INLINE
#undef INLINE
#endif

#if !defined(_NONE_MESH)
#include "MeshUtil.h"
#include "math1.h"
#include "mesh_packet_2d.h"
#include "meshtools2d.h"

#define D_MESHER_LOOP       0
#define D_MESHER_GRID       1
#define D_MESHER_DELAUNAY   2
#define D_MESHER_QMORPH			3

#define D_MESHTYPE_QUAD     0
#define D_MESHTYPE_TRIAQUAD 1
#define D_MESHTYPE_TRIA     2
#endif

#include "HeaderPre.h"

class CBaseMeshTool  
{
public:
    CBaseMeshTool(); 
    virtual ~CBaseMeshTool();

#if !defined(_NONE_MESH)
    //========<2D>========
    static int MeshGenerator2D(mesh_packet_2d::Mesh_Packet_2D& mp, CArray<T_MESH_ELEM_D,T_MESH_ELEM_D&>& Elems,
        const int Mesher, const int MeshType, bool (*interrupt_hdl)(void* pass_thru, double progress),
        CArray<UINT, UINT>& raFailedInd);

    /*
    static bool MeshGenerator2D_Adaptive(mesh_packet_2d::Mesh_Packet_2D& mp,
    mesh_packet_2d::Mesh_Packet_2D& BGMmp,
    CArray<double, double>& aBGMSize,
    CArray<T_ELEM_D,T_ELEM_D&>& Elems, const int Mesher,
    const int MeshType, bool (*interrupt_hdl)());

    static bool gen_TriaMesh_Adaptive(vector<mesh_packet_2d::Mesh_Packet_2D::I_Node>& vMeshNodes, 
    vector<int>& vIsolatedNodes, vector<double>& vIsolatedSizes,
    CArray<T_ELEM_D,T_ELEM_D&>& Elems, cm2::DoubleMat& pos,
    cm2::UIntMat& connectB, cm2::UIntMat& BMG, FloatVec& BMGsizes,
    const double dRefineFactor,	bool (*interrupt_hdl)());

    static bool CoonsMeshGenerator(mesh_packet_3d_inner::Mesh_Packet_3D_Inner& mp, 
    const int CoonsType, 
    const int MeshType, 
    bool (*interrupt_hdl)());

    static bool GenerateTriaByOnlyPnt(const CArray<T_NODE_D, T_NODE_D&>& aNodeD, CArray<T_ELEM_D, T_ELEM_D&>& raElemD);

    */

    static int gen_TriaMesh(std::vector<mesh_packet_2d::Mesh_Packet_2D::I_Node>& MeshNodes,
        std::vector<int>& vIsolatedNodes, std::vector<double>& vIsolatedSizes,
        CArray<T_MESH_ELEM_D,T_MESH_ELEM_D&>& Elems, 
        cm2::DoubleMat& pos, cm2::UIntMat& connectB, const double dRefineFactor,
        bool (*interrupt_hdl)(void* pass_thru, double progress),
        CArray<UINT, UINT>& raFailedInd);




    static int gen_QuadMesh(std::vector<mesh_packet_2d::Mesh_Packet_2D::I_Node>& MeshNodes,
		std::vector<int>& vIsolatedNodes, std::vector<double>& vIsolatedSizes,
        CArray<T_MESH_ELEM_D,T_MESH_ELEM_D&>& Elems, 
        cm2::DoubleMat& pos, cm2::UIntMat& connectB, const double dRefineFactor,
        const bool dPureQuad,
        bool (*interrupt_hdl)(void* pass_thru, double progress));

    //========<3D>========
    /*
    static int MeshGenerator3D(	const CArray<T_NODE_D,T_NODE_D&>& BdyNodes,
    const CArray<T_ELEM_D,T_ELEM_D&>& BdyElems, 
    CArray<T_NODE_D,T_NODE_D&>& tetraNodes,
    CArray<T_ELEM_D,T_ELEM_D&>& tetraElems, 
    const mesh_packet_3d_inner::Mesh_Packet_3D_Inner& mpinner,
    bool (*interrupt_hdl)(),
    CMap<UINT, UINT, UINT, UINT>&  mMergedInd,
    CMap<UINT, UINT, UINT, UINT>&  mMovedInd,
    const double dRefineCoeff,
    const bool bNoClamped,
    CArray<UINT, UINT>& raFailedInd);
    static int gen_TetraMesh(CArray<T_NODE_D,T_NODE_D&>& tetraNodes,
    CArray<T_ELEM_D,T_ELEM_D&>& Elems, cm2::DoubleMat& pos, 
    cm2::UIntMat& connectB, cm2::UIntMat& connectE,
    const mesh_packet_3d_inner::Mesh_Packet_3D_Inner& mpinner,
    bool (*interrupt_hdl)(),
    const double Coeff,
    const bool bNoClamped,
    CArray<UINT, UINT>& raFailedInd);

    static int gen_TetraMesh(CArray<T_NODE_D,T_NODE_D&>& tetraNodes,
    CArray<T_ELEM_D,T_ELEM_D&>& Elems, cm2::DoubleMat& pos, 
    cm2::UIntMat& connectB, cm2::UIntMat& connectE,
    cm2::UIntMat& BGM,
    FloatVec& sizes,
    const mesh_packet_3d_inner::Mesh_Packet_3D_Inner& mpinner,
    bool (*interrupt_hdl)(),
    const double Coeff,
    const bool bNoClamped,
    CArray<UINT, UINT>& raFailedInd);
    */

    static bool Seek_UnusedIndex(mesh_packet_2d::Mesh_Packet_2D& mp,
		std::vector<int>& v_unused_node_chain_index,
		std::vector<int>& v_unused_node_index,
        const bool bInnerDomain,
        const bool bIncludeInteriorEdge,
        bool (*interrupt_hdl)(void* pass_thru, double progress));

    //========<Check>========
    /*
    static int CheckBdyMesh(const CArray<T_NODE_D,T_NODE_D&>& aBdyNode,
    const CArray<T_ELEM_D,T_ELEM_D&>& aBdyElem, 
    vector<UINT>& v_fault_index);

    //========<Merge>========
    static bool MergeTriaMesh(CArray<T_NODE_D,T_NODE_D&>&   Nodes,
    CArray<T_ELEM_D,T_ELEM_D&>&   Elems,
    CMap<UINT, UINT, UINT, UINT>& mMergedInd,
    CMap<UINT, UINT, UINT, UINT>& mMovedInd,
    const double MergeTol);

    static bool Merge1DMesh(vector<mesh_packet_2d::Mesh_Packet_2D::I_Node>&       preNodes, 
    vector<mesh_packet_2d::Mesh_Packet_2D::I_Node_Chain>& preConnect,
    TColStd_SequenceOfInteger& newIndex,
    const double MergeTol=1.e-6);

    static bool Merge1DMesh(CArray<T_NODE_D,T_NODE_D&>& aNodeD, 
    CArray<T_ELEM_D,T_ELEM_D&>& aElemD, 
    const double dMergeTol);
    static bool MergeTetraMesh(CArray<T_NODE_D,T_NODE_D&>& Nodes,
    CArray<T_ELEM_D,T_ELEM_D&>& Elems,
    CMap<UINT, UINT, UINT, UINT>& mMergedInd,
    CMap<UINT, UINT, UINT, UINT>& mMovedInd,
    const double MergeTol);
    */

    static bool Merge2DMesh(CArray<T_MESH_NODE_D,T_MESH_NODE_D&>& Nodes,
        CArray<T_MESH_ELEM_D,T_MESH_ELEM_D&>& Elems,
        CMap<UINT, UINT, UINT, UINT>& mMergedInd,
        CMap<UINT, UINT, UINT, UINT>& mMovedInd,
        const double MergeTol,
        const bool bHighOrder = false);


    //========<Remesh>========

    // HJK : added hard edge array
    /*
    static bool RemeshTria(const CArray<T_NODE_D,T_NODE_D&>& Nodes,
    const CArray<T_ELEM_D,T_ELEM_D&>& Elems,
    const CArray<UINT, UINT&>& aHardFace,
    const CArray<T_HARDEDGE, T_HARDEDGE&>& aHardEdge,
    const CArray<T_POS_DBL_PAIR, T_POS_DBL_PAIR&>& aInnerVtxSizeInd,
    const double dMinH, const double dMaxH,
    const double dPatchAng,
    const int    nHardFlag, // boundary index를 지켜줘야하기 때문에 반드시 필요함.
    CArray<T_NODE_D,T_NODE_D&>& raNewNodeD,
    CArray<T_ELEM_D,T_ELEM_D&>& raNewElemD,
    bool (*interrupt_hdl)(),
    const bool bMerge = false);

    static bool RemeshQuad(const CArray<T_NODE_D,T_NODE_D&>& Nodes,
    const CArray<T_ELEM_D,T_ELEM_D&>& Elems,
    const CArray<T_HARDEDGE, T_HARDEDGE&>& aHardEdge,
    const CArray<T_POS_DBL_PAIR, T_POS_DBL_PAIR&>& aInnerVtxSizeInd,
    const double dMinH, const double dMaxH,
    const double dPatchAng,
    const int    nMode, // 0:tria+quad, 1:pure quad
    const int    nHardFlag, // boundary index를 지켜줘야하기 때문에 반드시 필요함.
    CArray<T_NODE_D,T_NODE_D&>& raNewNodeD,
    CArray<T_ELEM_D,T_ELEM_D&>& raNewElemD,
    bool (*interrupt_hdl)(),
    const bool bMerge = false);

    //========<Adaptive Mesh>==========
    // HJK
    static int AdaptiveMeshTria(const CArray<T_NODE_D,T_NODE_D&>& aNodes,
    const CArray<T_ELEM_D,T_ELEM_D&>& aElems,
    //const CArray<T_HARDEDGE, T_HARDEDGE&>& aHardEdge,
    const double dSize,
    const CMap<UINT, UINT, double, double>& mSize,
    CArray<T_NODE_D,T_NODE_D&>& raNewNodeD,
    CArray<T_ELEM_D,T_ELEM_D&>& raNewElemD,
    const double dRefineFactor,
    bool (*interrupt_hdl)(),
    CArray<UINT, UINT>& raFailedInd,
    const bool bMerge = false);

    static int AdaptiveMeshTetra(const CArray<T_NODE_D,T_NODE_D&>& aNodes,
    const CArray<T_ELEM_D,T_ELEM_D&>& aBGMElems,
    const CArray<T_ELEM_D,T_ELEM_D&>& aBndrElems,
    const double dSize,
    const CMap<UINT, UINT, double, double>& mBGMSize,
    CArray<T_NODE_D,T_NODE_D&>& raNewNodeD,
    CArray<T_ELEM_D,T_ELEM_D&>& raNewElemD,
    const double dRefineFactor,
    bool (*interrupt_hdl)(),
    CArray<UINT, UINT>& raFailedInd,
    const bool bMerge = false,
    const bool bNoClamped = false);

    static void cube_boundary (double L, unsigned N, cm2::DoubleMat& pos, cm2::UIntMat& connectS);
    static void display_hdl(unsigned level, const char* msg);
    */
#endif // !defined(_NONE_MESH)
};
#include "HeaderPost.h"
