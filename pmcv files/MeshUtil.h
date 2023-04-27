#pragma once

#include <vector>
#include "HeaderPre.h"

#define VTK_EMPTY_CELL                0
#define VTK_VERTEX                    1
#define VTK_POLY_VERTEX               2
#define VTK_LINE                      3
#define VTK_POLY_LINE                 4
#define VTK_TRIANGLE                  5
#define VTK_TRIANGLE_STRIP            6
#define VTK_POLYGON                   7
#define VTK_PIXEL                     8
#define VTK_QUAD                      9
#define VTK_TETRA                    10
#define VTK_VOXEL                    11
#define VTK_HEXAHEDRON               12
#define VTK_WEDGE                    13
#define VTK_PYRAMID                  14
#define VTK_PENTAGONAL_PRISM         15
#define VTK_HEXAGONAL_PRISM          16
// Quadratic, iso-parametric cells
#define VTK_QUADRATIC_EDGE           21
#define VTK_QUADRATIC_TRIANGLE       22
#define VTK_QUADRATIC_QUAD           23
#define VTK_QUADRATIC_TETRA          24
#define VTK_QUADRATIC_HEXAHEDRON     25
#define VTK_QUADRATIC_WEDGE          27
#define VTK_QUADRATIC_PYRAMID        26
// Special class of cells formed by convex group of points
#define VTK_CONVEX_POINT_SET 41
// Higher order cells in parametric form
#define VTK_PARAMETRIC_CURVE         51
#define VTK_PARAMETRIC_SURFACE       52
#define VTK_PARAMETRIC_TRI_SURFACE   53
#define VTK_PARAMETRIC_QUAD_SURFACE  54
#define VTK_PARAMETRIC_TETRA_REGION  55
#define VTK_PARAMETRIC_HEX_REGION    56
// Higher order cells
#define VTK_HIGHER_ORDER_EDGE        60
#define VTK_HIGHER_ORDER_TRIANGLE    61
#define VTK_HIGHER_ORDER_QUAD        62
#define VTK_HIGHER_ORDER_POLYGON     63
#define VTK_HIGHER_ORDER_TETRAHEDRON 64
#define VTK_HIGHER_ORDER_WEDGE       65
#define VTK_HIGHER_ORDER_PYRAMID     66
#define VTK_HIGHER_ORDER_HEXAHEDRON  67

#define MESHSIZETOL 1.0e-5

const double dMeshTolFactor = 0.1; // of the minMeshSize

#ifndef T_MESH_NODE_K
#define T_MESH_NODE_K UINT
typedef CArray<T_MESH_NODE_K, T_MESH_NODE_K> T_MESH_NODE_K_LIST;
#endif

struct T_MESH_NODE_D
{
    double x, y, z;
    UINT  CordK;
    UINT  nNo;
    bool  bPost;
    bool  bHard;

    void Initialize();
    T_MESH_NODE_D() { Initialize(); }
    T_MESH_NODE_D(const T_MESH_NODE_D& src) { *this = src; }
    T_MESH_NODE_D& operator=(const T_MESH_NODE_D& src)
    {
        x     = src.x;
        y     = src.y;
        z     = src.z;
        CordK = src.CordK;
        nNo   = src.nNo;
        bPost = src.bPost;
        bHard = src.bHard;
        return *this;
    }
};

typedef CArray<T_MESH_NODE_D, T_MESH_NODE_D&> T_MESH_NODE_D_LIST;
typedef CMap<INT_PTR, INT_PTR, T_MESH_NODE_D, T_MESH_NODE_D&> T_MESH_NODE_MAPEX;
typedef CMap<INT_PTR, INT_PTR, T_MESH_NODE_D, T_MESH_NODE_D&> T_MESH_NODE_MAP;

struct T_MESH_ELEM_D
{
    int vtktyp;	// 절점 개수와 요소 형상에 따른 타입. vtkCelltype.h참조.
    int eltyp;  // ????_EL (defined in DBStructure.h) // 현재 Post에서만 사용되고 있음.
    int nKind;  // v1006추가, 0:Normal, 1:Interface Elem, 2:Reinforcement Section
    int elmat;  // Material Number // 현재 사용안함.
    int elpro;  // Section or Thickness Number

    UINT nNo;   // v1001 추가, 요소 번호
    T_MESH_NODE_K_LIST aNodeK;
    T_MESH_NODE_K_LIST aInfoK;
    UINT   CordK;
    bool   bPost;

    void Initialize();
    T_MESH_ELEM_D() {}
    T_MESH_ELEM_D(int nType) { vtktyp = nType; }
    T_MESH_ELEM_D(const T_MESH_ELEM_D& src) { *this = src; }
    T_MESH_ELEM_D& operator=(const T_MESH_ELEM_D& src)
    {
        vtktyp = src.vtktyp;
        eltyp = src.eltyp;
        nKind = src.nKind;
        elmat = src.elmat;
        elpro = src.elpro;
        nNo   = src.nNo;
        aNodeK.Copy(src.aNodeK);
        aInfoK.Copy(src.aInfoK);
        CordK = src.CordK;
        bPost = src.bPost;
        return *this;
    }

    bool IsEqualElem(const T_MESH_ELEM_D& Other)
    {
        if ( vtktyp != Other.vtktyp ) return false;
        if ( eltyp  != Other.eltyp )  return false;
        if ( nKind  != Other.nKind )  return false;
        if ( elmat  != Other.elmat )  return false;
        if ( elpro  != Other.elpro )  return false;
        if ( nNo    != Other.nNo )    return false;
        if ( aNodeK.GetSize() != Other.aNodeK.GetSize() ) return false;
        if ( aInfoK.GetSize() != Other.aInfoK.GetSize() ) return false;
        for ( int i=0; i<aNodeK.GetSize(); i++ ) { if ( aNodeK[i] != Other.aNodeK[i] ) return false; }
        for ( int i=0; i<aInfoK.GetSize(); i++ ) { if ( aInfoK[i] != Other.aInfoK[i] ) return false; }
        if ( CordK  != Other.CordK )  return false;
        if ( bPost  != Other.bPost )  return false;

        return true;
    }
};

typedef CArray<T_MESH_ELEM_D, T_MESH_ELEM_D&> T_MESH_ELEM_D_LIST;
typedef CMap<INT_PTR, INT_PTR, T_MESH_ELEM_D, T_MESH_ELEM_D&> T_MESH_ELEM_MAPEX;
typedef CMap<INT_PTR, INT_PTR, T_MESH_ELEM_D, T_MESH_ELEM_D&> T_MESH_ELEM_MAP;

struct T_AUTO_MESH_D
{
    T_MESH_NODE_D_LIST aNode;
    T_MESH_ELEM_D_LIST aElem;

    void Initialize()
    {
        aNode.RemoveAll();
        aElem.RemoveAll();
    }
    T_AUTO_MESH_D() { Initialize(); }
    T_AUTO_MESH_D(const T_AUTO_MESH_D& src) { *this = src; }
    T_AUTO_MESH_D& operator=(const T_AUTO_MESH_D& src)
    {
        aNode.Copy(src.aNode);
        aElem.Copy(src.aElem);
        return *this;
    }
};

enum SEEDING_METHOD
{
    SEED_SIZE_METHOD,
    SEED_DIVISION_METHOD,
};

struct MESHSIZEINFO
{
    int    bEmpty;
    int    nConstraintType; // as SEEDING_METHOD
    int    nDivision;
    double dSize;
    int    bAdapSeed, bSymSeed;
    int    bManualSeed; // DB에 까지는 넣을 필요없다. 필요할 때 만들기 때문.

    double GradingSLen;
    double GradingELen;
    int    GradingDiv;
    double GradingRatio;
    double GradingCParam;

    CArray<double, double> aPosRatio;

    MESHSIZEINFO()
    {
        bEmpty          = 1;
        nConstraintType = 0;
        dSize           = 0;
        nDivision       = 0;
        bAdapSeed       = 0;
        bSymSeed        = 0;
        GradingSLen     = 0.;
        GradingELen     = 0.;
        GradingDiv      = 0;
        GradingRatio    = 0.;
        GradingCParam   = 0.;
        bManualSeed     = 0;
    }
    MESHSIZEINFO(const MESHSIZEINFO& src) { *this = src; }
    MESHSIZEINFO& operator=(const MESHSIZEINFO& src)
    {
        bEmpty          = src.bEmpty;
        nConstraintType = src.nConstraintType;
        dSize           = src.dSize;
        nDivision       = src.nDivision;
        bAdapSeed       = src.bAdapSeed;
        bSymSeed        = src.bSymSeed;
        GradingSLen     = src.GradingSLen;
        GradingELen     = src.GradingELen;
        GradingDiv      = src.GradingDiv;
        GradingRatio    = src.GradingRatio;
        GradingCParam   = src.GradingCParam;
        bManualSeed     = src.bManualSeed;
        aPosRatio.Copy(src.aPosRatio);
        return *this;
    }
};

class CMeshUtil
{
public:
    CMeshUtil();
    virtual ~CMeshUtil();

    static bool GetMidNode(const T_MESH_NODE_D& Node1, const T_MESH_NODE_D& Node2, T_MESH_NODE_D& MidNode);
    static bool checkESC_interruption(void* pass_thru, double progress);

public:
    static double CalcLength(const T_MESH_NODE_D& N1, const T_MESH_NODE_D& N2);
};

#include "HeaderPost.h"
