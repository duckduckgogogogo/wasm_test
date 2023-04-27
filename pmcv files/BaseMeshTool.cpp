// CM2MeshTool.cpp: implementation of the CBaseMeshTool class.
//
//////////////////////////////////////////////////////////////////////

#include "stdAfx.h"

#pragma warning ( disable : 4786 )
#pragma warning ( disable : 4018 )
#pragma warning ( disable : 4996 )

/*
cm2mesh 관련하여 일부 제품과 버전이 맞지 않아,
'_NONE_MESH'를 통해 일부만 수정하여 사용상에 문제가 없도록 하였습니다.
BaseMeshTool.h/.cpp과 AutoMesh_PlanarArea.h/.cpp에만 적용되어 있습니다. - kspark -
(현재 cm2mesh와 관련된 DgnEngine dll은 GEN 계열 dll로 
GEN 계열을 제외한 나머지는 'Debug_NoneMesh' configuration이 추가되어 있습니다.)
*/

#include "BaseMeshTool.h"

//2013/03/20 Soilworks에 출시되는 DgnBase는 None Mesh 로 release 되어 나간다. 
//안과장님이 cm2mesh 쓰는 함수들은 soilworks에서 불일 일이 없다고함. 문제되면 안과장님 책임~!!ㅋㅋㅋㅋ
#if !defined(_NONE_MESH)
#if defined(_X64)
#define CM2MATH_LIBNAME "cm2math1_x64_48.lib"
#define CM2QUADMESH_LIBNAME "cm2quadmesh_x64_48.lib"
#define CM2TRIAMESH_LIBNAME "cm2triamesh_x64_48.lib"
#define CM2MESHTOOL_LIBNAME "cm2meshtools_x64_48.lib"
#else
#define CM2MATH_LIBNAME "cm2math1_win32_48.lib"
#define CM2QUADMESH_LIBNAME "cm2quadmesh_win32_48.lib"
#define CM2TRIAMESH_LIBNAME "cm2triamesh_win32_48.lib"
#define CM2MESHTOOL_LIBNAME "cm2meshtools_win32_48.lib"
#endif

#pragma comment (lib, CM2MATH_LIBNAME)
#pragma comment (lib, CM2QUADMESH_LIBNAME)
#pragma comment (lib, CM2TRIAMESH_LIBNAME)
#pragma comment (lib, CM2MESHTOOL_LIBNAME)

// SYSTEM
#include <iostream>
#include <fstream>


// MISC
#include "fe_defs.h"
#include "misc.h"


#include "domain_maker.h"

// TRIA & TETRA MESHERS
#include "triamesh.h"
#include "quadmesh.h"
//#include "tetramesh.h"
//#include "surf_remesh_t3.h"
//#include "surf_remesh_q4.h"

// Mesher_2D
#include "mesher_2d.h"
//#include "loop_mesher.h"
//#include "delaunay_mesher.h"
//#include "grid_mesher.h"
//#include "coons_mesher.h"
#include <typeinfo>

//#define _PROFILING
//#include "..\MIDAS_lib\ProfileWrite.h"

#define MSG_NODE_ERROR            _LS(IDS_MESHGEN___MSG1)
#define MSG_EDGE_ERROR            _LS(IDS_MESHGEN___MSG2)
#define MSG_BOUNDARY2D_ERROR      _LS(IDS_MESHGEN___MSG4)
#define MSG_BOUNDARY3D_ERROR      _LS(IDS_MESHGEN___MSG5)
#define APPENDMSG_RECOMMEND_SIZE  _LS(IDS_MESHGEN___MSG6)
#endif // !defined(_NONE_MESH)

#if !defined(_NONE_MESH)
/*
//============================================================================
// ExportMctInit
//============================================================================
static void ExportMctInit(ofstream& Outmct)
{
Outmct<<";---------------------------------------------------------------------------"<<endl;
Outmct<<";  MIDAS/Civil Text(MCT) File."<<endl;
Outmct<<";  Date : "<<__DATE__<<endl;
Outmct<<";---------------------------------------------------------------------------"<<endl;
Outmct<<endl<<"*VERSION\n   5.8.0"<<endl;	
Outmct<<endl<<"*UNIT    ; Unit System"<<endl;
Outmct<<"; FORCE, LENGTH"<<endl;
Outmct<<"   TONF , M"<<endl;	

Outmct<<endl<<"*PROJINFO    ; Project Information"<<endl;
Outmct<<"   USER    = MIDAS"<<endl;		
Outmct<<"   ADDRESS = MIDAS IT"<<endl;	
}
*/

//============================================================================
// ExportMct
//============================================================================
/*
static void ExportMct(const CString& FileName, 
const cm2::DoubleMat& dmPos,
const cm2::UIntMat& TriaElem, 
const cm2::UIntMat& LineElem)
{
ofstream Outmct(FileName);
int i;
ExportMctInit(Outmct);

Outmct<<endl<<"*NODE    ; Nodes"<<endl;
for(i=0;i<(int)dmPos.cols();i++){
Outmct.width(6);Outmct<<i+1<<", ";
Outmct.width(15);	Outmct<<dmPos(0,i)<<", ";
Outmct.width(15); Outmct<<dmPos(1,i)<<", ";
Outmct.width(15);	Outmct<<dmPos(2,i)<<endl;
}
Outmct<<endl<<"*ELEMENT    ; Elements"<<endl;
Outmct<<"; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iN5, iN6, iN7, iN8 \t; Solid  Element"<<endl;
Outmct<<"; iEL, TYPE, iMAT, iPRO, iN1, iN2, ANGLE, iSUB, EXVAL              ; Frame  Element"<<endl;
Outmct<<"; iEL, TYPE, iMAT, iPRO, iN1, iN2, iN3, iN4, iSUB, iWID            ; Planar Element"<<endl;
Outmct<<"; iEL, TYPE, iMAT, iPRO, iN1, iN2, REF, RPX, RPY, RPZ, iSUB, EXVAL ; Frame(Ref. Point)"<<endl;

for(i=0;i<(int)TriaElem.cols();i++)
{
Outmct.width(6);Outmct<<i+1<<", PLATE ,    1,     1, ";
Outmct.width(5);Outmct<<TriaElem(0,i)+1<<", ";
Outmct.width(5);Outmct<<TriaElem(1,i)+1<<", ";
Outmct.width(5);Outmct<<TriaElem(2,i)+1;
Outmct<<",     0,     1"<<endl;
}

int nTria = TriaElem.cols();
for(i=0;i<(int)LineElem.cols();i++)
{
Outmct.width(6);Outmct<<nTria + i+1<<", BEAM  ,    1,     1, ";
Outmct.width(5);Outmct<<LineElem(0,i)+1<<", ";
Outmct.width(5);Outmct<<LineElem(1,i)+1;
Outmct<<",     0"<<endl;
}

Outmct<<endl<<"*ENDDATA"<<endl;
}
*/

using namespace cm2;
using namespace std;
//============================================================================
// Simplify
// CM2에서 제공되는 simplify는 이동된 index가 제대로 나오지 못한다.
// 두개의 mesh가 merge될때 internal edge가 앞에 놓이도록 connectE를 먼저한다.
//
// MergedInd는 merge되면서 없어지는 정보를 갖고 있으며,
// MovedInd는 없어지는 node로 인해 index가 댕겨지면서 갱신된 정보를 갖는다.
//============================================================================
static void Simplify(cm2::DoubleMat& pos, 
                     cm2::UIntVec& mergedInd,
                     cm2::UIntMat& connectE, 
                     cm2::UIntMat& connectB,
                     CMap<UINT, UINT, UINT, UINT>&  mMergedInd,
                     CMap<UINT, UINT, UINT, UINT>&  mMovedInd,
                     CMap<UINT, UINT, UINT, UINT>*  pmIsolated = NULL)
{

    int i;
    CMap<UINT, UINT, UINT, UINT> mUsedInd;
    UINT dummyInd;
    if(pmIsolated)
    {
        POSITION position = pmIsolated->GetStartPosition();
        UINT     nIsolated;
        while(position != NULL)
        {
            pmIsolated->GetNextAssoc(position, nIsolated, dummyInd);
            mUsedInd.SetAt(nIsolated, 0);
        }
    }
    for(i=0; i<(int)connectE.cols(); i++)
    {
        for(int j=0; j<(int)connectE.rows(); j++)
        {
            if(!mUsedInd.Lookup(connectE(j,i), dummyInd))
                mUsedInd.SetAt(connectE(j,i), 0);
        }
    }
    for(i=0; i<(int)connectB.cols(); i++)
    {
        for(int j=0; j<(int)connectB.rows(); j++)
        {
            if(!mUsedInd.Lookup(connectB(j,i), dummyInd))
                mUsedInd.SetAt(connectB(j,i), 0);
        }
    }

    // merged index
    int nPreNodeSize = pos.cols();
    for(i=0;i<nPreNodeSize;i++)
    {
        if((int)mergedInd[i]!=i)
        {
            mMergedInd.SetAt(i, mergedInd[i]);
        }
    }

    if(mMergedInd.GetCount() == 0 && (mUsedInd.GetCount() == nPreNodeSize)) return; 

    // simplify & move
    cm2::DoubleMat tmpPos = pos;
    pos.clear_hard();
    pos.resize(3, mUsedInd.GetCount());
    int count = 0;

    // merge되어 대체되는 수가 뒤에 놓여있을 경우에 대비하여 map을 만들어 놓는다.
    CMap<UINT, UINT, UINT, UINT> mMergedOrigin;
    mMergedOrigin.InitHashTable((mMergedInd.GetCount()));

    UINT nMovedInd;
    for(i=0; i<nPreNodeSize; i++)
    {
        if(!mUsedInd.Lookup(i, dummyInd)) 
        {
            UINT nMergedInd;
            if(mMergedInd.Lookup(i, nMergedInd)) 
            {
                if(mMovedInd.Lookup(nMergedInd, nMovedInd))
                    mMergedInd.SetAt(i, nMovedInd);
                mMergedOrigin.SetAt(nMergedInd, i);
            }
            else
            {
                mMergedInd.SetAt(i, count);
            }
            continue;
        }
        pos(0,count) = tmpPos(0,i);
        pos(1,count) = tmpPos(1,i);
        pos(2,count) = tmpPos(2,i);
        if(i != count)
        {
            mMovedInd.SetAt(i, count);
            UINT nRepairInd; 
            if(mMergedOrigin.Lookup(i, nRepairInd))
            {
                if(mMergedInd.Lookup(nRepairInd, dummyInd))
                {
                    mMergedInd.SetAt(nRepairInd, count);
                }
                else ASSERT(0); // call Y.Jee ~!
            }
        }
        count++;
    }

    for(i=0; i<(int)connectE.cols(); i++)
    {
        for(int j=0; j<(int)connectE.rows(); j++)
        {
            if(mMovedInd.Lookup(connectE(j,i), nMovedInd))
                connectE(j,i) = nMovedInd;
        }
    }	
    for(i=0; i<(int)connectB.cols(); i++)
    {
        for(int j=0; j<(int)connectB.rows(); j++)
        {
            if(mMovedInd.Lookup(connectB(j,i), nMovedInd))
                connectB(j,i) = nMovedInd;
        }
    }
}

#endif // !defined(_NONE_MESH)

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CBaseMeshTool::CBaseMeshTool()
{

}

CBaseMeshTool::~CBaseMeshTool()
{

}


#if !defined(_NONE_MESH)
//*****************************************************************
//
//	       2 D   M E S H    G E N E R A T I O N 
//
//*****************************************************************

//==============================================================
// MeshGenerator2D
//==============================================================
int CBaseMeshTool::MeshGenerator2D(mesh_packet_2d::Mesh_Packet_2D& mp, CArray<T_MESH_ELEM_D,T_MESH_ELEM_D&>& Elems,
                                   const int Mesher, const int MeshType, bool (*interrupt_hdl)(void* pass_thru, double progress),
                                   CArray<UINT, UINT>& raFailedInd)
{

    using namespace mesh_packet_2d;

    if (Mesher == D_MESHER_DELAUNAY || (Mesher == D_MESHER_LOOP && MeshType == D_MESHTYPE_TRIA))
    {
        bool bChangeMethod = false;
        int seedSize    = mp.v_node.size();
        int connectSize =0;
        for(int i=0;i<(int)mp.v_node_chain.size();i++)
        {
            if(mp.v_node_chain[i].v_node_index.size() == 0) continue;
            if(mp.v_node_chain[i].t_type != Mesh_Packet_2D::I_Node_Chain::CT_OPEN)
            {
                connectSize += mp.v_node_chain[i].v_node_index.size()-1;
            }else{
                connectSize += 2*(mp.v_node_chain[i].v_node_index.size()-1);
            }
        }

        cm2::DoubleMat pos      = cm2::DoubleMat(2, seedSize);
        cm2::UIntMat   connectB = cm2::UIntMat  (2, connectSize);

        // isolated node를 뒤에 놓을 경우 전체의 mesh size가 node사이즈에 지배된다.
        // 맨앞에 놓아야 grading 되는 듯한 형상의 mesh들을 얻을 수 있다.
        vector<int>    vIsolatedNodes;
        vector<double> vIsolatedSizes;
        for(int i=0; i<seedSize; ++i)
        {
            if(mp.v_node[i].b_island_flag) 
            {
                vIsolatedNodes.push_back(i);
                vIsolatedSizes.push_back(mp.v_node[i].d_mesh_size);
            }
            pos(0,i) = mp.v_node[i].a_uv[0];
            pos(1,i) = mp.v_node[i].a_uv[1];
        }
        int count=0;
        for(int i=0;i<(int)mp.v_node_chain.size();i++)
        {
            if(mp.v_node_chain[i].t_type != Mesh_Packet_2D::I_Node_Chain::CT_OPEN)
            {
                if(MeshType==D_MESHTYPE_QUAD) // 짝수만 가능
                {
                    if((mp.v_node_chain[i].v_node_index.size()-1)%2 ==1 )
                    {
                        // [2009-09-21] Kim, Geun Young (Tel: 2042, gykim@midasit.com) : puad + tri로 변경해준다.
                        //MessageBox(NULL, _LS(IDS_MESHGEN___MSG7), _LS(IDS_MESHGEN___MSG8), MB_ICONINFORMATION);            
                        //return 0;
                        bChangeMethod = true;            
                    }
                }
                for(int j=0;j<(int)mp.v_node_chain[i].v_node_index.size()-1;j++,count++)
                {
                    connectB(0,count) = mp.v_node_chain[i].v_node_index[j];
                    connectB(1,count) = mp.v_node_chain[i].v_node_index[j+1];
                }
            }
            else
            {
                for(int j=0 ; j<(int)mp.v_node_chain[i].v_node_index.size()-1 ; j++,count++)
                {
                    connectB(0,count) = mp.v_node_chain[i].v_node_index[j];
                    connectB(1,count) = mp.v_node_chain[i].v_node_index[j+1];
                }
                for(int j=mp.v_node_chain[i].v_node_index.size()-1; j>0 ; j--,count++)
                {
                    connectB(0,count) = mp.v_node_chain[i].v_node_index[j];
                    connectB(1,count) = mp.v_node_chain[i].v_node_index[j-1];
                }
            }
        }

        // 항상 Mesher == D_MESHER_DELAUNAY일때이므로
        switch(MeshType)
        {
        case D_MESHTYPE_TRIA:
            return gen_TriaMesh(mp.v_node, vIsolatedNodes, vIsolatedSizes, 
                Elems, pos, connectB, mp.d_refinement_ratio, interrupt_hdl, raFailedInd);
        case D_MESHTYPE_TRIAQUAD:      
            return gen_QuadMesh(mp.v_node, vIsolatedNodes, vIsolatedSizes,
                Elems, pos, connectB, mp.d_refinement_ratio, false, interrupt_hdl);        
            break;
        case D_MESHTYPE_QUAD:
            if(bChangeMethod)
            {
                int nRetVal = gen_QuadMesh(mp.v_node, vIsolatedNodes, vIsolatedSizes,
                    Elems, pos, connectB, mp.d_refinement_ratio, false, interrupt_hdl);        
                if(nRetVal != 0)        
                    return 2;
                else
                    return nRetVal;
            }
            else
                return gen_QuadMesh(mp.v_node, vIsolatedNodes, vIsolatedSizes,
                Elems, pos, connectB, mp.d_refinement_ratio, true, interrupt_hdl);
            break;
        default:
            break;
        }

    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    // 단위계로 인하여 매우 작은 수치로 인한 버그를 방지하기 위한 scale factor를 넣어준다.
    //CUnitCtrl* pUnitCtrl = CDBDoc::GetDocPoint()->m_pUnitCtrl; 
    double dFactor=0.0;
    if (dFactor<1.0e-6) { ASSERT(0); dFactor=1.0; }
    mp.d_scale_factor = 1./dFactor;
    ///////////////////////////////////////////////////////////////////////////////////////////


    //-----------------------------------------------------------------------------
    mesher_2d::Mesher_2D* p_mesher_2d = 0;
    /*
    switch (Mesher)
    {
    case D_MESHER_LOOP: p_mesher_2d = new loop_mesher::Loop_Mesher; break;
    case D_MESHER_GRID: 
    p_mesher_2d = new grid_mesher::Grid_Mesher;
    break;
    }
    */
    switch (MeshType)
    {
    case D_MESHTYPE_QUAD:     mp.t_mesh_type = Mesh_Packet_2D::MT_QUAD; break;
    case D_MESHTYPE_TRIAQUAD: mp.t_mesh_type = Mesh_Packet_2D::MT_COMB; break;
    case D_MESHTYPE_TRIA:     mp.t_mesh_type = Mesh_Packet_2D::MT_TRIA; break;
    }
    /////////////////////////////////////////////
    while (true)
    {
        try
        {
			bool(*temp_interrupt)(void) = nullptr;
			p_mesher_2d->generate_mesh(mp, temp_interrupt);
            break;
        }
        catch (mesher_2d::interrupt_exception)
        {
            delete p_mesher_2d;			
            return 0;
        }
        catch (exception& e)
        {	
            /*
            if (typeid(*p_mesher_2d) == typeid(loop_mesher::Loop_Mesher))
            {
            delete p_mesher_2d;
            p_mesher_2d    = new grid_mesher::Grid_Mesher;
            mp.t_mesh_type = Mesh_Packet_2D::MT_COMB;
            continue;
            }
            */
            delete p_mesher_2d;
            CString csMsg = e.what();
            return 0;
        }
    }	
    delete p_mesher_2d;

    /////////////////////////////////////////////////////////////////////
    // 만들어진 node정보는 mp.v_node에 담아놓고 3D로 다시 변환해야한다.
    //
    // element정보는 3D로 변해도 상관없으므로 Elems로 변환해서 내보낸다.
    //
    int nElem = mp.v_element.size();
    Elems.RemoveAll();
    Elems.SetSize(nElem);
    T_MESH_ELEM_D elem;
    for (int i=0; i<nElem; ++i)
    {
        const Mesh_Packet_2D::I_Element& i_element = mp.v_element[i];

        if (i_element.v_node_index.size() == 3)
        {
            elem.Initialize();
            elem.vtktyp = VTK_TRIANGLE;
            elem.aNodeK.SetSize(3);
            for (int j=0; j<3; ++j) elem.aNodeK.SetAt(j, i_element.v_node_index[j]);
            Elems.SetAt(i, elem);
        }
        else
        {
            elem.Initialize();
            elem.vtktyp = VTK_QUAD;
            elem.aNodeK.SetSize(4);
            for (int j=0; j<4; ++j) elem.aNodeK.SetAt(j, i_element.v_node_index[j]);
            Elems.SetAt(i, elem);		
        }
    }

    return 1;
}

//==============================================================
// MeshGenerator2D
//==============================================================
/*
bool CBaseMeshTool::MeshGenerator2D_Adaptive(mesh_packet_2d::Mesh_Packet_2D& mp,
mesh_packet_2d::Mesh_Packet_2D& BGMmp,
CArray<double, double>& aBGMSize,
CArray<T_ELEM_D,T_ELEM_D&>& Elems, const int Mesher,
const int MeshType, bool (*interrupt_hdl)())
{
M_PROFILE("bool CBaseMeshTool::MeshGenerator2D_Adaptive()");

int i,j;
using namespace mesh_packet_2d;

if (Mesher == D_MESHER_DELAUNAY || (Mesher == D_MESHER_LOOP && MeshType == D_MESHTYPE_TRIA))
{
int seedSize    = mp.v_node.size();
int connectSize =0;
for(i=0;i<(int)mp.v_node_chain.size();i++)
{
if(mp.v_node_chain[i].v_node_index.size() == 0) continue;
if(mp.v_node_chain[i].t_type != Mesh_Packet_2D::I_Node_Chain::CT_OPEN)
{
connectSize += mp.v_node_chain[i].v_node_index.size()-1;
}else{
connectSize += 2*(mp.v_node_chain[i].v_node_index.size()-1);
}
}

int nBGMNode = BGMmp.v_node.size();
int nAllNode = seedSize + nBGMNode;
int nBGMElem = BGMmp.v_element.size();
cm2::DoubleMat pos      = cm2::DoubleMat(2, nAllNode);
cm2::UIntMat   connectB = cm2::UIntMat  (2, connectSize);
cm2::UIntMat   BGM			= cm2::UIntMat  (3, nBGMElem);
cm2::cm2::UIntMat   indices;
cm2::FloatVec  BGMSizes;

// isolated node를 뒤에 놓을 경우 전체의 mesh size가 node사이즈에 지배된다.
// 맨앞에 놓아야 grading 되는 듯한 형상의 mesh들을 얻을 수 있다.
vector<int>    vIsolatedNodes;
vector<double> vIsolatedSizes;
for(i=0; i<seedSize; ++i)
{
if(mp.v_node[i].b_island_flag) 
{
vIsolatedNodes.push_back(i);
vIsolatedSizes.push_back(mp.v_node[i].d_mesh_size);
}
pos(0,i) = mp.v_node[i].a_uv[0];
pos(1,i) = mp.v_node[i].a_uv[1];
}
for (i = 0; i < nBGMNode; i++)
{
int index = i+seedSize;
pos(0,index) = BGMmp.v_node[i].a_uv[0];
pos(1,index) = BGMmp.v_node[i].a_uv[1];
}

int count=0;
for(i=0;i<(int)mp.v_node_chain.size();i++)
{
if(mp.v_node_chain[i].t_type != Mesh_Packet_2D::I_Node_Chain::CT_OPEN)
{
if(MeshType==D_MESHTYPE_QUAD) // 짝수만 가능
{
if((mp.v_node_chain[i].v_node_index.size()-1)%2 ==1 )
{
MessageBox(NULL, _LS(IDS_MESHGEN___MSG7), _LS(IDS_MESHGEN___MSG8), MB_ICONINFORMATION);
return false;
}
}
for(j=0;j<(int)mp.v_node_chain[i].v_node_index.size()-1;j++,count++)
//for(j=(int)mp.v_node_chain[i].v_node_index.size()-2; j>=0; j--,count++)
{
connectB(0,count) = mp.v_node_chain[i].v_node_index[j];
connectB(1,count) = mp.v_node_chain[i].v_node_index[j+1];
//connectB(1,count) = mp.v_node_chain[i].v_node_index[j];
//connectB(0,count) = mp.v_node_chain[i].v_node_index[j+1];
}
}
else
{
for(j=0;j<(int)mp.v_node_chain[i].v_node_index.size()-1;j++,count++)
//for(j=(int)mp.v_node_chain[i].v_node_index.size()-2; j>=0; j--,count++)
{
connectB(0,count) = mp.v_node_chain[i].v_node_index[j];
connectB(1,count) = mp.v_node_chain[i].v_node_index[j+1];
//connectB(1,count) = mp.v_node_chain[i].v_node_index[j];
//connectB(0,count) = mp.v_node_chain[i].v_node_index[j+1];
}
for(j=0;j<(int)mp.v_node_chain[i].v_node_index.size()-1;j++,count++)
//for(j=(int)mp.v_node_chain[i].v_node_index.size()-2; j>=0; j--,count++)
{
connectB(0,count) = mp.v_node_chain[i].v_node_index[j];
connectB(1,count) = mp.v_node_chain[i].v_node_index[j-1];
//connectB(1,count) = mp.v_node_chain[i].v_node_index[j];
//connectB(0,count) = mp.v_node_chain[i].v_node_index[j-1];
}
}
}

int nElem = BGMmp.v_element.size();
for (i = 0; i < nElem; i++)
{
int nNode = BGMmp.v_element[i].v_node_index.size();
for (j = 0; j < nNode; j++)
{
BGM(j, i) = BGMmp.v_element[i].v_node_index[j] + seedSize;
//TRACE("BMG %d %d : %d\n", j, i, BGM(j, i));
}	
}

meshtools::unique_indices(indices, BGM);
BGMSizes.resize(nAllNode, 0.0F);
for (i = 0; i < indices.size(); i++)
{
int n = indices[i];
BGMSizes[n] = aBGMSize[n-seedSize];
}

return gen_TriaMesh_Adaptive(mp.v_node, vIsolatedNodes, vIsolatedSizes, 
Elems, pos, connectB, BGM, BGMSizes, mp.d_refinement_ratio, interrupt_hdl);
}

///////////////////////////////////////////////////////////////////////////////////////////
// 단위계로 인하여 매우 작은 수치로 인한 버그를 방지하기 위한 scale factor를 넣어준다.
CUnitCtrl* pUnitCtrl = CMECDocBase::GetCurDoc()->GetUnitCtrl(); 
mp.d_scale_factor = 1./pUnitCtrl->GetConvertFactorCurrent(D_UNITSYS_BASE_LENGTH);
///////////////////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
mesher_2d::Mesher_2D* p_mesher_2d = 0;
switch (Mesher)
{
case D_MESHER_LOOP: p_mesher_2d = new loop_mesher::Loop_Mesher; break;
case D_MESHER_GRID: 
p_mesher_2d = new grid_mesher::Grid_Mesher;
break;
}
switch (MeshType)
{
case D_MESHTYPE_QUAD:     mp.t_mesh_type = Mesh_Packet_2D::MT_QUAD; break;
case D_MESHTYPE_TRIAQUAD: mp.t_mesh_type = Mesh_Packet_2D::MT_COMB; break;
case D_MESHTYPE_TRIA:     mp.t_mesh_type = Mesh_Packet_2D::MT_TRIA; break;
}
/////////////////////////////////////////////
while (true)
{
try
{
p_mesher_2d->generate_mesh(mp, interrupt_hdl);
break;
}
catch (mesher_2d::interrupt_exception)
{
delete p_mesher_2d;
MessageBox(NULL, _LS(IDS_MESHGEN___MSG9), _LS(IDS_MESHGEN___MSG8), MB_ICONINFORMATION);
return false;
}
catch (exception& e)
{			
if (typeid(*p_mesher_2d) == typeid(loop_mesher::Loop_Mesher))
{
delete p_mesher_2d;
p_mesher_2d    = new grid_mesher::Grid_Mesher;
mp.t_mesh_type = Mesh_Packet_2D::MT_COMB;
continue;
}
delete p_mesher_2d;
CString csMsg = e.what();
MessageBox(NULL, csMsg, _LS(IDS_MESHGEN___MSG8), MB_ICONWARNING);
return false;
}
}	
delete p_mesher_2d;

/////////////////////////////////////////////////////////////////////
// 만들어진 node정보는 mp.v_node에 담아놓고 3D로 다시 변환해야한다.
//
// element정보는 3D로 변해도 상관없으므로 Elems로 변환해서 내보낸다.
//
int nElem = mp.v_element.size();
Elems.RemoveAll();
Elems.SetSize(nElem);
T_ELEM_D elem;
for (i=0; i<nElem; ++i)
{
const Mesh_Packet_2D::I_Element& i_element = mp.v_element[i];

if (i_element.v_node_index.size() == 3)
{
elem.Initialize();
elem.vtktyp = VTK_TRIANGLE;
elem.aNodeK.SetSize(3);
for (int j=0; j<3; ++j) elem.aNodeK.SetAt(j, i_element.v_node_index[j]);
Elems.SetAt(i, elem);
}
else
{
elem.Initialize();
elem.vtktyp = VTK_QUAD;
elem.aNodeK.SetSize(4);
for (int j=0; j<4; ++j) elem.aNodeK.SetAt(j, i_element.v_node_index[j]);
Elems.SetAt(i, elem);		
}
}

return true;
}
*/

int CBaseMeshTool::gen_TriaMesh(vector<mesh_packet_2d::Mesh_Packet_2D::I_Node>& vMeshNodes, 
                                vector<int>& vIsolatedNodes, vector<double>& vIsolatedSizes,
                                CArray<T_MESH_ELEM_D,T_MESH_ELEM_D&>& Elems, 
                                cm2::DoubleMat& pos, cm2::UIntMat& connectB, const double dRefineFactor,
                                bool (*interrupt_hdl)(void* pass_thru, double progress),
                                CArray<UINT, UINT>& raFailedInd)
{

    int i,j;
    triamesh::mesher triamesher;
    triamesh::mesher::data_type data(pos, connectB);

#if defined(_X64)
    triamesh::registration("Licensed to MIDAS_IT.", "0D47869FF95E"); // 20030614 license update  
#else
    triamesh::registration("Licensed to MIDAS_IT.", "7AA782C37EB9"); // 20030614 license update  
#endif


    // [2009-10-09] Kim, Geun Young (Tel: 2042, gykim@midasit.com) : Gen에서는 Mesh Pref.가 없으므로 default 값 사용.    
    //T_PREFERENCE pref;
    //CDBDoc::GetDocPoint()->m_pInitCtrl->GetPreference(pref);
    //triamesher.settings.optim_level          = pref.Mesh.nOptLev2D;
    //triamesher.settings.shape_quality_weight = float(pref.Mesh.nWeightQlt2D)*.1;
    triamesher.settings.optim_level          = 7;
    triamesher.settings.shape_quality_weight = 0.7;	
    //triamesher.settings.sizes_mult_coeff     = dRefineFactor;

    for(i=0;i<(int)vIsolatedNodes.size();i++)
    {
        data.isolated_nodes.push_back(vIsolatedNodes[i]);
        data.metrics.resize(vIsolatedNodes[i]+1, vIsolatedSizes[i]);
    }
    triamesher.settings.interrupt_hdl = interrupt_hdl;

    clock_t  StartTime = clock();
    triamesher.run (data);
    /*
    if (CProduct::GetTestEnvValue(_T("ProfilingMesh")) == _T("yes")) 
    {
    CString sTimeMsg;
    double dTime = (double(clock() - StartTime) / CLOCKS_PER_SEC);
    sTimeMsg.Format(_T("[PROFILING] cm2::Tria-Mesh Core Time : %.5g[sec]"), dTime);
    CDBDoc::GetDocPoint()->OutputMessage(sTimeMsg, 0);	
    }
    */
    _tremove(_T("WIN32 TRIAMESH.log"));

    // interrupted state
    if (data.warning_code == triamesh::mesher::data_type::CM2_INTERRUPTION) return 0;

    if (data.error_code != triamesh::mesher::data_type::CM2_NO_ERROR)
    {
        CString str;
        switch(data.error_code)
        {
        //case triamesh::mesher::data_type::CM2_DATA_ERROR:
        //    str = APPENDMSG_RECOMMEND_SIZE;
        //    break;
        //case triamesh::mesher::data_type::CM2_BOUNDARY_ERROR:
        //    str = MSG_BOUNDARY2D_ERROR;
        //    break;
        //case triamesh::mesher::data_type::CM2_NODE_ERROR:
        //    str = MSG_NODE_ERROR;
        //    str += APPENDMSG_RECOMMEND_SIZE;
        //    break;
        //case triamesh::mesher::data_type::CM2_EDGE_ERROR:
        //    str = MSG_EDGE_ERROR;
        //    str += APPENDMSG_RECOMMEND_SIZE;
        //    break;
            /*
            case triamesh::mesher::data_type::CM2_SYSTEM_MEMORY_ERROR:
            str = _LS(IDS_MESHGEN___MSG10);
            break;
            */  
        default:
            str.Format(_T("1. %s\n2. %s"), data.msg1, data.msg2);
            break;
        }
        /*#ifdef _DEBUG
        CString sDebugFile;
        sDebugFile.Format(_T("c:/temp/Tria%s.txt"),__TIME__);
        data.save(LPCTSTR(sDebugFile));
        #endif*/

        cm2::UIntVec vFailedBdy = data.unenforced_boundary_IDs;
        raFailedInd.RemoveAll();
        raFailedInd.SetSize(vFailedBdy.size());
        for(i=0; i<vFailedBdy.size(); i++)
        {
            raFailedInd.SetAt(i, vFailedBdy[i]);
        }   

        return -1;
    }

    int nOutNode = data.pos.cols();
    using namespace mesh_packet_2d;
    vMeshNodes = vector<mesh_packet_2d::Mesh_Packet_2D::I_Node>(nOutNode);

    for(i=0; i<nOutNode; i++)
    {
        Mesh_Packet_2D::I_Node i_node;
        i_node.b_hard_flag   = 1;
        i_node.a_uv[0]       = data.pos(0,i);
        i_node.a_uv[1]       = data.pos(1,i);
        vMeshNodes[i] = i_node;
    }

    int nElem = data.connectM.cols();
    Elems.RemoveAll();
    Elems.SetSize(nElem);

    T_MESH_ELEM_D elem;
    for(i=0; i<nElem; i++)
    {
        elem.Initialize();
        elem.vtktyp = VTK_TRIANGLE;
        elem.aNodeK.SetSize(3);
        for(j=0; j<3; j++) elem.aNodeK.SetAt(j, data.connectM(j,i));
        Elems.SetAt(i, elem);
    }

    data.clear();

    return 1;
}

/*
bool CBaseMeshTool::gen_TriaMesh_Adaptive(vector<mesh_packet_2d::Mesh_Packet_2D::I_Node>& vMeshNodes, 
vector<int>& vIsolatedNodes, vector<double>& vIsolatedSizes,
CArray<T_ELEM_D,T_ELEM_D&>& Elems, cm2::DoubleMat& pos,
cm2::UIntMat& connectB, cm2::UIntMat& BGM, FloatVec& BGMSizes,
const double dRefineFactor,	bool (*interrupt_hdl)())
{
M_PROFILE(_T("bool CBaseMeshTool::gen_TriaMesh_Adaptive()"));

int i,j;
triamesh::mesher triamesher;
triamesh::mesher::data_type data(pos, connectB);
//triamesh::registration(_LSX(Licensed to MIDAS_IT.), _T("23A165BC19BF")); 
triamesh::registration(_LSX(Licensed to MIDAS_IT.), _T("7AA782C37EB9")); // 20030614 license update  
T_PREFERENCE pref;
CDBDoc::GetDocPoint()->GetInitCtrl()->GetPreference(pref);
//triamesher.settings.optim_level          = pref.Mesh.nOptLev2D;
//triamesher.settings.shape_quality_weight = float(pref.Mesh.nWeightQlt2D)*.1;
//triamesher.settings.sizes_mult_coeff     = dRefineFactor;
data.background_mesh = BGM;
data.metrics = BGMSizes;

#ifdef _DEBUG
data.save(_T("c:/temp/Adaptive_data.txt"));
#endif

//meshtools::vtk_output(_T("c:/vtk_win32/cm2/vtk/win32/Adaptive_bndr.vtk"), pos, data.connectB, CM2_EDGE2);
//meshtools::vtk_output(_T("c:/vtk_win32/cm2/vtk/win32/Adaptive_BGM.vtk"), pos, data.background_mesh, CM2_FACET3);

for(i=0;i<(int)vIsolatedNodes.size();i++)
{
data.isolated_nodes.push_back(vIsolatedNodes[i]);
data.metrics.resize(vIsolatedNodes[i]+1, vIsolatedSizes[i]);
}
triamesher.settings.interrupt_hdl = interrupt_hdl;

clock_t  StartTime = clock();
triamesher.run (data);
if (CProduct::GetTestEnvValue(_T("ProfilingMesh")) == _T("yes")) 
{
CString sTimeMsg;
double dTime = (double(clock() - StartTime) / CLOCKS_PER_SEC);
sTimeMsg.Format(_T("[PROFILING] cm2::Tria-Mesh Core Time : %.5g[sec]"), dTime);
CMECDocBase::GetCurDoc()->OutputMessage(sTimeMsg, 0);	
}

remove(_T("WIN32 TRIAMESH.log"));

// interrupted state
if (data.warning_code == triamesh::mesher::data_type::CM2_INTERRUPTION) return 0;

if (data.error_code != triamesh::mesher::data_type::CM2_NO_ERROR)
{
CString str;
switch(data.error_code)
{
case triamesh::mesher::data_type::CM2_DATA_ERROR:
str = APPENDMSG_RECOMMEND_SIZE;
break;
case triamesh::mesher::data_type::CM2_BOUNDARY_ERROR:
str = MSG_BOUNDARY2D_ERROR;
break;
case triamesh::mesher::data_type::CM2_NODE_ERROR:
str = MSG_NODE_ERROR;
str += APPENDMSG_RECOMMEND_SIZE;
break;
case triamesh::mesher::data_type::CM2_EDGE_ERROR:
str = MSG_EDGE_ERROR;
str += APPENDMSG_RECOMMEND_SIZE;
break;
case triamesh::mesher::data_type::CM2_SYSTEM_MEMORY_ERROR:
str = _LS(IDS_MESHGEN___MSG10);
break;
default:
str.Format(_T("1. %s\n2. %s"), data.msg1, data.msg2);
break;
}

MessageBox(NULL, str, _LS(IDS_MESHGEN___MSG11), MB_ICONWARNING);
return false;
}

int nOutNode = data.pos.cols();
using namespace mesh_packet_2d;
vMeshNodes = vector<mesh_packet_2d::Mesh_Packet_2D::I_Node>(nOutNode);

for(i=0; i<nOutNode; i++)
{
Mesh_Packet_2D::I_Node i_node;
i_node.b_hard_flag   = 1;
i_node.a_uv[0]       = data.pos(0,i);
i_node.a_uv[1]       = data.pos(1,i);
vMeshNodes[i] = i_node;
}

int nElem = data.connectM.cols();
Elems.RemoveAll();
Elems.SetSize(nElem);

T_ELEM_D elem;
for(i=0; i<nElem; i++)
{
elem.Initialize();
elem.vtktyp = VTK_TRIANGLE;
elem.aNodeK.SetSize(3);
for(j=0; j<3; j++) elem.aNodeK.SetAt(j, data.connectM(j,i));
Elems.SetAt(i, elem);
}

data.clear();

return true;
}
*/

int CBaseMeshTool::gen_QuadMesh(vector<mesh_packet_2d::Mesh_Packet_2D::I_Node>& vMeshNodes, 
                                vector<int>& vIsolatedNodes, vector<double>& vIsolatedSizes,
                                CArray<T_MESH_ELEM_D,T_MESH_ELEM_D&>& Elems, 
                                cm2::DoubleMat& pos, cm2::UIntMat& connectB, const double dRefineFactor,
                                const bool dPureQuad,
                                bool (*interrupt_hdl)(void* pass_thru, double progress))
{

    quadmesh::mesher quadmesher;
    quadmesh::mesher::data_type data(pos, connectB);	

#if defined(_X64)
    quadmesh::registration("Licensed to MIDAS_IT.", "AC3841C9B2A1"); // 20030614 license update  
#else
    quadmesh::registration("Licensed to MIDAS_IT.", "4590C38C678E"); // 20030614 license update
#endif

    // [2009-10-09] Kim, Geun Young (Tel: 2042, gykim@midasit.com) : Gen에서는 Mesh Pref.가 없으므로 default 값 사용.  
    //T_PREFERENCE pref;
    //CDBDoc::GetDocPoint()->m_pInitCtrl->GetPreference(pref);
    //quadmesher.settings.optim_level          = pref.Mesh.nOptLev2D;
    //quadmesher.settings.shape_quality_weight = float(pref.Mesh.nWeightQlt2D)*.1;
    quadmesher.settings.optim_level          = 7;  
    quadmesher.settings.shape_quality_weight = 0.7;
    //quadmesher.settings.sizes_mult_coeff     = dRefineFactor;
    quadmesher.settings.all_quad_flag = dPureQuad;
    for(int i=0;i<(int)vIsolatedNodes.size();i++)
    {
        data.isolated_nodes.push_back(vIsolatedNodes[i]);
        data.metrics.resize(vIsolatedNodes[i]+1, vIsolatedSizes[i]);
    }
    quadmesher.settings.interrupt_hdl = interrupt_hdl;

    clock_t  StartTime = clock();
    quadmesher.run (data);
    /*
    if (CProduct::GetTestEnvValue(_T("ProfilingMesh")) == _T("yes")) 
    {
    CString sTimeMsg;
    double dTime = (double(clock() - StartTime) / CLOCKS_PER_SEC);
    sTimeMsg.Format(_T("[PROFILING] cm2::Quad-Mesh Core Time : %.5g[sec]"), dTime);
    CMECDocBase::GetCurDoc()->OutputMessage(sTimeMsg, 0);	
    }
    remove(_T("WIN32 QUADMESH.log"));
    */  

    // if interrupted state
    if (data.warning_code == quadmesh::mesher::data_type::CM2_INTERRUPTION) return 0;

    bool bChangeMethod = false;
    // pure quad 변경해서 mesh 재수행 해본다.

    //임시 주석.
    //   if (quadmesher.settings.pure_quad_flag == true &&
    //       data.error_code == quadmesh::mesher::data_type::CM2_BOUNDARY_ERROR)
    //   { 
    //     data.clear();
    //     data.pos.shallow_copy(pos);
    //     data.connectB.shallow_copy(connectB);
    //     for(i=0;i<(int)vIsolatedNodes.size();i++)
    //     {
    //       data.isolated_nodes.push_back(vIsolatedNodes[i]);
    //       data.metrics.resize(vIsolatedNodes[i]+1, vIsolatedSizes[i]);
    //     }	  
    //     quadmesher.settings.pure_quad_flag = false;
    //     quadmesher.run (data);    
    //     bChangeMethod = true;
    //   }

    if (data.error_code != quadmesh::mesher::data_type::CM2_NO_ERROR)
    {
        CString str;
        switch(data.error_code)
        {
        //case quadmesh::mesher::data_type::CM2_DATA_ERROR:
        //    str = APPENDMSG_RECOMMEND_SIZE;
        //    break;
        //case quadmesh::mesher::data_type::CM2_BOUNDARY_ERROR:
        //    str = _LS(IDS_MESHGEN___MSG26); // seed 수가 적어도 발생할 수 있음.
        //    break;
        //case quadmesh::mesher::data_type::CM2_NODE_ERROR:
        //    str = MSG_NODE_ERROR;
        //    str += APPENDMSG_RECOMMEND_SIZE;
        //    break;
        //case quadmesh::mesher::data_type::CM2_EDGE_ERROR:
        //    str = MSG_EDGE_ERROR;
        //    str += APPENDMSG_RECOMMEND_SIZE;
        //    break;
            /*  
            case quadmesh::mesher::data_type::CM2_SYSTEM_MEMORY_ERROR:
            str = _LS(IDS_MESHGEN___MSG10);
            break;
            */  
        default:
            str.Format(_T("1. %s\n2. %s"), CString(data.msg1), CString(data.msg2));
            break;
        }
        /*#ifdef _DEBUG
        CString sDebugFile;
        sDebugFile.Format(_T("c:/temp/Quad%s.txt"),__TIME__);
        data.save(LPCTSTR(sDebugFile));
        #endif*/
        return 0;
    }

    using namespace mesh_packet_2d;
    vMeshNodes.clear();
    Elems.RemoveAll();

    UINT nNode = data.pos.cols();
    vMeshNodes = vector<mesh_packet_2d::Mesh_Packet_2D::I_Node>(nNode);
    for(UINT i=0; i<nNode; i++)
    {
        Mesh_Packet_2D::I_Node i_node;
        i_node.b_hard_flag   = 1;
        i_node.a_uv[0]       = data.pos(0,i);
        i_node.a_uv[1]       = data.pos(1,i);
        vMeshNodes[i] = i_node;
    }

    T_MESH_ELEM_D elem;
    UINT nElem = data.connectM.cols();
    Elems.SetSize(nElem);
    for(UINT i=0; i<nElem; i++)
    {
        elem.Initialize();
        elem.vtktyp = VTK_QUAD;
        elem.aNodeK.SetSize(4);
        for(UINT j=0; j<4; j++)
        {
            UINT NodeK = data.connectM(j,i);
            if(j==3 && NodeK>nNode) // triangle
            {
                elem.vtktyp = VTK_TRIANGLE;
                elem.aNodeK.SetSize(3);
                break;
            }
            elem.aNodeK.SetAt(j, data.connectM(j,i));
        }
        Elems.SetAt(i, elem);
    }

    if(bChangeMethod) return 2;

    return 1;
}

/*
bool CBaseMeshTool::CoonsMeshGenerator(mesh_packet_3d_inner::Mesh_Packet_3D_Inner& mp, 
const int CoonsType, 
const int Mesher, 
bool (*interrupt_hdl)())
{
//M_PROFILE(_T("bool CBaseMeshTool::CoonsMeshGenerator()"));

using namespace coons_mesher;
Coons_Mesher CoonsMesher;

if(CoonsType == 3)
CoonsMesher.generate_coons_3_mesh(mp, Mesher, interrupt_hdl);
else if(CoonsType == 4)
CoonsMesher.generate_coons_4_mesh(mp, Mesher, interrupt_hdl);

return true;
}
*/

/*
bool CBaseMeshTool::GenerateTriaByOnlyPnt(const CArray<T_NODE_D, T_NODE_D&>& aNodeD, CArray<T_ELEM_D, T_ELEM_D&>& raElemD)
{
int i,j;
raElemD.RemoveAll();
int nSeed = aNodeD.GetSize();

cm2::DoubleMat pos      = cm2::DoubleMat(2, nSeed);
cm2::UIntMat   connectB = cm2::UIntMat  (2, 0);

for(i=0; i<nSeed; i++)
{
pos(0,i) = aNodeD[i].x;
pos(1,i) = aNodeD[i].y;
}
triamesh::mesher triamesher;
triamesher.settings.basic_mode = cm2::triamesh::mesher::operating_mode_type::CONVEX_HULL_MODE;
triamesh::mesher::data_type data(pos, connectB);  
for(i=0; i<nSeed; i++)
{
data.isolated_nodes.push_back(i);
}
triamesh::registration(_LSX(Licensed to MIDAS_IT.), _T("7AA782C37EB9")); 

clock_t  StartTime = clock();
triamesher.run (data);
remove(_T("WIN32 TRIAMESH.log"));

// interrupted state
if (data.warning_code == triamesh::mesher::data_type::CM2_INTERRUPTION) return 0;

if (data.error_code != triamesh::mesher::data_type::CM2_NO_ERROR)
{
CString str;
switch(data.error_code)
{
case triamesh::mesher::data_type::CM2_DATA_ERROR:
str = APPENDMSG_RECOMMEND_SIZE;
break;
case triamesh::mesher::data_type::CM2_BOUNDARY_ERROR:
str = MSG_BOUNDARY2D_ERROR;
break;
case triamesh::mesher::data_type::CM2_NODE_ERROR:
str = MSG_NODE_ERROR;
str += APPENDMSG_RECOMMEND_SIZE;
break;
case triamesh::mesher::data_type::CM2_EDGE_ERROR:
str = MSG_EDGE_ERROR;
str += APPENDMSG_RECOMMEND_SIZE;
break;
case triamesh::mesher::data_type::CM2_SYSTEM_MEMORY_ERROR:
str = _LS(IDS_MESHGEN___MSG10);
break;
default:
str.Format(_T("1. %s\n2. %s"), data.msg1, data.msg2);
break;
}
MessageBox(NULL, str, _LS(IDS_MESHGEN___MSG11), MB_ICONWARNING);
return false;
}

int nElem = data.connectM.cols();
raElemD.SetSize(nElem);
T_ELEM_D TriaElem; TriaElem.Initialize();
TriaElem.vtktyp = VTK_TRIANGLE;
TriaElem.aNodeK.SetSize(3);
for(i=0; i<nElem; i++)
{
TriaElem.aNodeK.SetAt(0, data.connectM(0,i));
TriaElem.aNodeK.SetAt(1, data.connectM(1,i));
TriaElem.aNodeK.SetAt(2, data.connectM(2,i));
raElemD.SetAt(i, TriaElem);
}

data.clear();
return true;
}
*/


//*****************************************************************
//
//	        3 D   M E S H    G E N E R A T I O N 
//
//*****************************************************************
//==================================================================================
// MeshGenerator3D()
//==================================================================================
// 실패시 bdy mesh만든다.
/*
int CBaseMeshTool::MeshGenerator3D(const CArray<T_NODE_D,T_NODE_D&>& BdyNodes,
const CArray<T_ELEM_D,T_ELEM_D&>& BdyElems, 
CArray<T_NODE_D,T_NODE_D&>& tetraNodes,
CArray<T_ELEM_D,T_ELEM_D&>& tetraElems,
const mesh_packet_3d_inner::Mesh_Packet_3D_Inner& mpinner,
bool (*interrupt_hdl)(),
CMap<UINT, UINT, UINT, UINT>&  mMergedInd,
CMap<UINT, UINT, UINT, UINT>&  mMovedInd,
const double Coeff,
const bool bNoClamped,
CArray<UINT, UINT>& raFailedInd)
{
M_PROFILE(_T("bool CBaseMeshTool::MeshGenerator3D()"));

int i,j,count;
int NbBdyNodeSize    = BdyNodes.GetSize();
int NbBdyConnectSize = BdyElems.GetSize();
int NbInnerEdges     = mpinner.v_node_chain.size();
raFailedInd.RemoveAll();

cm2::DoubleMat pos      = cm2::DoubleMat(3, NbBdyNodeSize    + mpinner.v_node.size() + mpinner.v_chainNodes.size());
cm2::UIntMat   connectB = cm2::UIntMat(3,   NbBdyConnectSize);
cm2::UIntMat   connectE;

// isolated node를 뒤에 놓을 경우 전체의 mesh size가 node사이즈에 지배된다.
// 맨앞에 놓아야 grading 되는 듯한 형상의 mesh들을 얻을 수 있다.
int NbInnerVtxNodes  = mpinner.v_node.size();
int NbInnerEdgeNodes = mpinner.v_chainNodes.size();
int NbInner = NbInnerVtxNodes + NbInnerEdgeNodes;

CMap<UINT, UINT, UINT, UINT> mInnerVtx;
mInnerVtx.InitHashTable(GF_GetHashSizeByCount(NbInnerVtxNodes));
for(i=0;i<NbInnerVtxNodes;i++)
{
pos(0,i) = mpinner.v_node[i].a_coord[0];
pos(1,i) = mpinner.v_node[i].a_coord[1];
pos(2,i) = mpinner.v_node[i].a_coord[2];
mInnerVtx.SetAt(i,0);
}
for(count=NbInnerVtxNodes, i=0;i<NbInnerEdgeNodes;i++, count++)
{
pos(0,count) = mpinner.v_chainNodes[i].a_coord[0];
pos(1,count) = mpinner.v_chainNodes[i].a_coord[1];
pos(2,count) = mpinner.v_chainNodes[i].a_coord[2];
}
count = NbInner;
for(i=0;i<NbBdyNodeSize;i++,count++)
{
pos(0,count) = BdyNodes[i].x;
pos(1,count) = BdyNodes[i].y;
pos(2,count) = BdyNodes[i].z;
}

int NbInnerEdgeElems=0;
for(i=0;i<NbInnerEdges;i++)
{
if(mpinner.v_node_chain[i].v_node_index.size() == 0) continue;
NbInnerEdgeElems += mpinner.v_node_chain[i].v_node_index.size()-1;
if(mpinner.v_node_chain[i].t_type == mesh_packet_3d_inner::Mesh_Packet_3D_Inner::I_Node_Chain::CT_CLOSED) 
NbInnerEdgeElems++;
}
connectE.resize(2,NbInnerEdgeElems);

double dEdgeMergeTol;
if(NbInnerEdges)
{
int ind1 = NbInnerVtxNodes;
int ind2 = ind1+1;
dEdgeMergeTol = pow((pos(0,ind1)-pos(0,ind2)),2) + pow((pos(1,ind1)-pos(1,ind2)),2) + pow((pos(2,ind1)-pos(2,ind2)),2);
dEdgeMergeTol = sqrt(dEdgeMergeTol)*1.e-3;
}		

count=0;
for(i=0;i<NbInnerEdges;i++)
{
if(mpinner.v_node_chain[i].v_node_index.size() == 0) continue;
for(j=0;j<(int)mpinner.v_node_chain[i].v_node_index.size()-1;j++,count++)
{
connectE(0,count) = mpinner.v_node_chain[i].v_node_index[j]   + NbInnerVtxNodes;
connectE(1,count) = mpinner.v_node_chain[i].v_node_index[j+1] + NbInnerVtxNodes;
}
if(mpinner.v_node_chain[i].t_type == mesh_packet_3d_inner::Mesh_Packet_3D_Inner::I_Node_Chain::CT_CLOSED)
{
connectE(0,count) = mpinner.v_node_chain[i].v_node_index[j]   + NbInnerVtxNodes;
connectE(1,count) = mpinner.v_node_chain[i].v_node_index[0] + NbInnerVtxNodes;
count++;
}
}

for(i=0;i<NbBdyConnectSize;i++)
{
for(j=0;j<3;j++) connectB(j,i)=BdyElems[i].aNodeK[j] + NbInner;
}

if(NbInnerEdges != 0) 
{
cm2::UIntMat mergedInd;
meshtools::merge(pos,connectE, connectB, mergedInd, dEdgeMergeTol, 0);
Simplify(pos, mergedInd, connectE, connectB, mMergedInd, mMovedInd, &mInnerVtx);
}

#ifdef _DEBUG
ExportMct(_T("c:/Temp/FXmesh_PreTetra.mct"), pos, connectB, connectE);
#endif

return gen_TetraMesh(tetraNodes, tetraElems, pos, connectB, connectE, mpinner, interrupt_hdl, Coeff, bNoClamped, raFailedInd);
}
*/

//==================================================================================
// gen_TetraMesh()
//==================================================================================
/*
int CBaseMeshTool::gen_TetraMesh(CArray<T_NODE_D,T_NODE_D&>& tetraNodes,
CArray<T_ELEM_D,T_ELEM_D&>& Elems, 
cm2::DoubleMat& pos, 
cm2::UIntMat& connectB,
cm2::UIntMat& connectE,
const mesh_packet_3d_inner::Mesh_Packet_3D_Inner& mpinner,
bool (*interrupt_hdl)(),
const double Coeff,
const bool bNoClamped,
CArray<UINT, UINT>& raFailedInd)
{
M_PROFILE(_T("bool CBaseMeshTool::gen_TetraMesh()"));

int i,j;
tetramesh::mesher tetramesher;
tetramesh::mesher::data_type data(pos, connectB);
//tetramesh::registration(_LSX(Licensed to MIDAS_IT.), _T("AFF5600F5425"));   //old version
tetramesh::registration (_LSX(Licensed to MIDAS_IT.), _T("46EF1F05BB3F"));
tetramesher.settings.default_numbering_flag = false;
tetramesher.settings.no_clamped_tetra_flag = bNoClamped;

T_PREFERENCE pref;
CDBDoc::GetDocPoint()->GetInitCtrl()->GetPreference(pref);
tetramesher.settings.optim_level          = pref.Mesh.nOptLev3D;
tetramesher.settings.shape_quality_weight = float(pref.Mesh.nWeightQlt3D)*.1;
tetramesher.settings.sizes_mult_coeff     = Coeff;
tetramesher.settings.interrupt_hdl        = interrupt_hdl;

for(i=0;i<(int)mpinner.v_node.size();i++)
{
data.isolated_nodes.push_back(i);
data.metrics.resize(i+1, mpinner.v_node[i].d_mesh_size);
}
data.connectE = connectE;

clock_t  StartTime = clock();
tetramesher.run (data);
if (CProduct::GetTestEnvValue(_T("ProfilingMesh")) == _T("yes")) 
{
CString sTimeMsg;
double dTime = (double(clock() - StartTime) / CLOCKS_PER_SEC);
sTimeMsg.Format(_T("[PROFILING] cm2::Tetra-Mesh Core Time : %.5g[sec]"), dTime);
CMECDocBase::GetCurDoc()->OutputMessage(sTimeMsg, 0);	
}
remove(_T("WIN32 TETRAMESH.log"));

//meshtools::vtk_output(_T("c:/vtk_win32/cm2/vtk/win32/gen_TetraMesh1.vtk"), pos, data.connectB, CM2_FACET3);
//meshtools::vtk_output(_T("c:/vtk_win32/cm2/vtk/win32/gen_TetraMesh2.vtk"), pos, data.connectM, CM2_FACET3);

if (data.warning_code == tetramesh::mesher::data_type::CM2_INTERRUPTION) return 0;

if (data.error_code != tetramesh::mesher::data_type::CM2_NO_ERROR)
{
CString str;
switch(data.error_code)
{
case tetramesh::mesher::data_type::CM2_DATA_ERROR:
case tetramesh::mesher::data_type::CM2_INTERNAL_ERROR:
str = APPENDMSG_RECOMMEND_SIZE;
break;
case tetramesh::mesher::data_type::CM2_BOUNDARY_ERROR:
str = MSG_BOUNDARY3D_ERROR;
break;
case tetramesh::mesher::data_type::CM2_NODE_ERROR:
str = MSG_NODE_ERROR;
str += APPENDMSG_RECOMMEND_SIZE;
break;
case tetramesh::mesher::data_type::CM2_EDGE_ERROR:
str = MSG_EDGE_ERROR;
str += APPENDMSG_RECOMMEND_SIZE;
break;
case tetramesh::mesher::data_type::CM2_FACE_ERROR:
str = MSG_FACE_ERROR;
str += APPENDMSG_RECOMMEND_SIZE;
break;
default:
str.Format(_T("1. %s\n2. %s"), data.msg1, data.msg2);
break;
}
MessageBox(NULL, str, _LS(IDS_MESHGEN___MSG11), MB_ICONWARNING);

cm2::UIntMat vFailedBdy = data.unenforced_boundary_IDs;
raFailedInd.RemoveAll();
raFailedInd.SetSize(vFailedBdy.size());
for(i=0; i<vFailedBdy.size(); i++)
{
raFailedInd.SetAt(i, vFailedBdy[i]);
}
return -1;
}

////////////////////////////////////////////////////////////////////////////////////////
// simplify 함수에서 Tetra일때 connectivity를 update하는데에 문제있는것 같음.(20050211)
// 라이브러리가 수정되기 전까지 임시로 걸러내는 루틴을 만들어서 사용한다.
//int nChkSim = meshtools::simplify(data.pos, data.connectM); 
CMap<UINT, UINT, int, int> mNodeInd;
int nIndCount = 0;
int dum;

mNodeInd.InitHashTable(GF_GetHashSizeByCount(data.pos.cols()));
for(i=0; i<(int)data.connectM.cols(); i++)
{
for(j=0; j<4; j++)
{
UINT nNodeInd = data.connectM(j,i);
if(!mNodeInd.Lookup(nNodeInd, dum))
{
mNodeInd.SetAt(nNodeInd, 0);
nIndCount++;
}
}
}
int nPrevNodeCount = data.pos.cols();
if(nPrevNodeCount != nIndCount)
{
bool bChkBdyNodeShift = false;
bool bShift = true;
int nNewCount = 0;
cm2::DoubleMat TmpPos(3, nIndCount);
for(i=0; i<nPrevNodeCount; i++)
{
if(mNodeInd.Lookup(i, dum))
{
if(!bChkBdyNodeShift && i != nNewCount)
{
if(i<pos.cols()) { bShift = false; break; }
bChkBdyNodeShift = true;
}
TmpPos(0, nNewCount) = data.pos(0,i);
TmpPos(1, nNewCount) = data.pos(1,i);
TmpPos(2, nNewCount) = data.pos(2,i);
mNodeInd.SetAt(i, nNewCount);
nNewCount++;
}
}
if(bShift)
{
data.pos.clear_hard();
data.pos.copy(TmpPos);
for(i=0; i<data.connectM.cols(); i++)
{
data.connectM(0,i) = mNodeInd[data.connectM(0,i)];
data.connectM(1,i) = mNodeInd[data.connectM(1,i)];
data.connectM(2,i) = mNodeInd[data.connectM(2,i)];
data.connectM(3,i) = mNodeInd[data.connectM(3,i)];
}
}
}
////////////////////////////////////////////////////////////////////////////////////////

tetraNodes.RemoveAll();
Elems.RemoveAll();

tetraNodes.SetSize(data.pos.cols());

T_NODE_D  node;
for(i=0; i<(int)data.pos.cols(); i++)
{
node.x = data.pos(0,i);
node.y = data.pos(1,i);
node.z = data.pos(2,i);
tetraNodes.SetAt(i, node);
}

Elems.SetSize(data.connectM.cols());
T_ELEM_D elem;
elem.Initialize();
elem.vtktyp = VTK_TETRA;
elem.aNodeK.SetSize(4);
for(i=0; i<(int)data.connectM.cols(); i++)
{
for(j=0; j<4; j++) elem.aNodeK.SetAt(j, data.connectM(j,i));
Elems.SetAt(i, elem);
}
return 1;
}
*/

//==================================================================================
// Seek_UnusedIndex()
//==================================================================================
bool CBaseMeshTool::Seek_UnusedIndex(mesh_packet_2d::Mesh_Packet_2D& mp,
                                     vector<int>&    v_unused_node_chain_index, 
                                     vector<int>&    v_unused_node_index,
                                     const bool      bInnerDomain,
                                     const bool      bIncludeInteriorEdge,
                                     bool (*interrupt_hdl)(void* pass_thru, double progress))
{

	using namespace domain_maker;
	Domain_Maker* p_domain_maker = new Domain_Maker;
    bool b_inner    = (bInnerDomain)? true:false;
    bool b_int_edge = (bIncludeInteriorEdge)? true:false;

    //CUnitCtrl* pUnitCtrl = CDBDoc::GetDocPoint()->m_pUnitCtrl; 
    double dFactor = 1.0;
    if (dFactor < 1.0e-6) {ASSERT(0); dFactor=1.0;}
    mp.d_scale_factor = 1./dFactor;

    /////////////////////////////////////////////
    try
    {
		bool(*interrupt_temp)(void) = nullptr;
        p_domain_maker->make_domain(mp, v_unused_node_chain_index, v_unused_node_index, b_inner, b_int_edge, interrupt_temp);
    }
    catch (exception& e)
    {
        delete p_domain_maker;
        p_domain_maker = NULL;
        CString csMsg = e.what();
        MessageBox(NULL, csMsg, _T("Auto-Mesh Planar Area"), MB_ICONINFORMATION);
        return false;
    }
    delete p_domain_maker;
    p_domain_maker = NULL;

    return true;
}

/*
/////////////////////////////////////////////////
namespace check_bdry_mesh_20040518
{
struct Bound
{
double MaxX;
double MaxY;
double MaxZ;
double MinX;
double MinY;
double MinZ;
};
///////////////////////////////////////////////
inline void cross_product_3d(const double* a_vector_1, const double* a_vector_2, double* a_result_vector)
{
a_result_vector[0] = a_vector_1[1] * a_vector_2[2] - a_vector_2[1] * a_vector_1[2];
a_result_vector[1] = a_vector_2[0] * a_vector_1[2] - a_vector_1[0] * a_vector_2[2];
a_result_vector[2] = a_vector_1[0] * a_vector_2[1] - a_vector_2[0] * a_vector_1[1];
} // end: cross_product_3d()
/////////////////////////////////////////////////////////////////////
bool is_clockwise_rotation_3d(const double* a_tria_xyz_1, const double* a_tria_xyz_2, const double* a_tria_xyz_3, const double* a_normal)
{  
double a_vector[3];
a_vector[0]          = (((a_tria_xyz_2[1] - a_tria_xyz_1[1]) * (a_tria_xyz_3[2] - a_tria_xyz_1[2]))
- ( (a_tria_xyz_3[1] - a_tria_xyz_1[1]) * (a_tria_xyz_2[2] - a_tria_xyz_1[2])));
a_vector[1]          = (((a_tria_xyz_2[2] - a_tria_xyz_1[2]) * (a_tria_xyz_3[0] - a_tria_xyz_1[0]))
- ( (a_tria_xyz_3[2] - a_tria_xyz_1[2]) * (a_tria_xyz_2[0] - a_tria_xyz_1[0])));
a_vector[2]          = (((a_tria_xyz_2[0] - a_tria_xyz_1[0]) * (a_tria_xyz_3[1] - a_tria_xyz_1[1]))
- ( (a_tria_xyz_3[0] - a_tria_xyz_1[0]) * (a_tria_xyz_2[1] - a_tria_xyz_1[1])));
double d_dot_product = a_vector[0] * a_normal[0] + a_vector[1] * a_normal[1] + a_vector[2] * a_normal[2];
return (d_dot_product >= 0.0);
} // end: is_clockwise_rotation_3d()
/////////////////////////////////////////////////////////////////////
bool triangle_line_intersection_3d(const double* a_vector_12,  const double* a_vector_23,
const double* a_tria_xyz_1, const double* a_tria_xyz_2, const double* a_tria_xyz_3,
const double* a_line_xyz_1, const double* a_line_xyz_2)
{
const double ZERO = 1.0e-6;
/////////////////////////////////////////////
//	double a_vector_12[3], a_vector_23[3];
//	transform(a_tria_xyz_2, a_tria_xyz_2+3, a_tria_xyz_1, a_vector_12,   std::minus<double>());
//	transform(a_tria_xyz_3, a_tria_xyz_3+3, a_tria_xyz_2, a_vector_23,   std::minus<double>());
double a_normal[3], a_line_vector[3];
transform(a_line_xyz_2, a_line_xyz_2+3, a_line_xyz_1, a_line_vector, std::minus<double>());
cross_product_3d(a_vector_12, a_vector_23, a_normal);
/////////////////////////////////////////////
double d_dot_product = a_normal[0] * a_line_vector[0] + a_normal[1] * a_line_vector[1] + a_normal[2] * a_line_vector[2];
if (fabs(d_dot_product) < ZERO) return false;
/////////////////////////////////////////////
double a_ip_xyz[3];
double d_t = - (a_normal[0] * (a_line_xyz_1[0] - a_tria_xyz_1[0])
+  a_normal[1] * (a_line_xyz_1[1] - a_tria_xyz_1[1])
+  a_normal[2] * (a_line_xyz_1[2] - a_tria_xyz_1[2]))
/ (a_normal[0] * a_line_vector[0] + a_normal[1] * a_line_vector[1] + a_normal[2] * a_line_vector[2]);
if (d_t < 0.0 || d_t > 1.0) return false;
for (int i=0; i<3; ++i) a_ip_xyz[i] = a_line_xyz_1[i] + a_line_vector[i] * d_t;
/////////////////////////////////////////////
return (is_clockwise_rotation_3d(a_tria_xyz_1, a_tria_xyz_2, a_ip_xyz, a_normal) && 
is_clockwise_rotation_3d(a_tria_xyz_2, a_tria_xyz_3, a_ip_xyz, a_normal) && 
is_clockwise_rotation_3d(a_tria_xyz_3, a_tria_xyz_1, a_ip_xyz, a_normal));
} // end: triangle_line_intersection_3d()
}
*/

/*
/////////////////////////////////////////////////
int CBaseMeshTool::CheckBdyMesh(const T_NODE_D_LIST& aBdyNode,
const T_ELEM_D_LIST& aBdyElem, 
vector<UINT>& v_fault_index)
{
M_PROFILE(_T("bool CBaseMeshTool::CheckBdyMesh()"));

using namespace check_bdry_mesh_20040518;
///////////////////////////////////////////////
v_fault_index.clear();
///////////////////////////////////////////////
if (aBdyNode.GetSize() < 4 || aBdyElem.GetSize() < 4) return 0;
///////////////////////////////////////////////
vector<Bound> vElemBound(aBdyElem.GetSize());
for (int i=0; i<aBdyElem.GetSize(); ++i)
{
T_ELEM_D ElemD  = aBdyElem[i];
T_NODE_D NodeD1 = aBdyNode[ElemD.aNodeK[0]];
T_NODE_D NodeD2 = aBdyNode[ElemD.aNodeK[1]];
T_NODE_D NodeD3 = aBdyNode[ElemD.aNodeK[2]];
Bound    bound;
bound.MaxX      = __max(__max(NodeD1.x, NodeD2.x), NodeD3.x);
bound.MaxY      = __max(__max(NodeD1.y, NodeD2.y), NodeD3.y);
bound.MaxZ      = __max(__max(NodeD1.z, NodeD2.z), NodeD3.z);
bound.MinX      = __min(__min(NodeD1.x, NodeD2.x), NodeD3.x);
bound.MinY      = __min(__min(NodeD1.y, NodeD2.y), NodeD3.y);
bound.MinZ      = __min(__min(NodeD1.z, NodeD2.z), NodeD3.z);
vElemBound[i]   = bound;
}
///////////////////////////////////////////////
for (i=0; i<aBdyElem.GetSize()-1; ++i)
{
const T_ELEM_D& ElemD1     = aBdyElem[i];
int             aNdKey1[3] = { ElemD1.aNodeK[0], ElemD1.aNodeK[1], ElemD1.aNodeK[2] };
const T_NODE_D& NodeD1     = aBdyNode[aNdKey1[0]];
const T_NODE_D& NodeD2     = aBdyNode[aNdKey1[1]];
const T_NODE_D& NodeD3     = aBdyNode[aNdKey1[2]];
double          aXYZ1[3]   = { NodeD1.x, NodeD1.y, NodeD1.z };
double          aXYZ2[3]   = { NodeD2.x, NodeD2.y, NodeD2.z };
double          aXYZ3[3]   = { NodeD3.x, NodeD3.y, NodeD3.z };
double          aVec12[3]  = { aXYZ2[0]-aXYZ1[0], aXYZ2[1]-aXYZ1[1], aXYZ2[2]-aXYZ1[2] };
double          aVec23[3]  = { aXYZ3[0]-aXYZ2[0], aXYZ3[1]-aXYZ2[1], aXYZ3[2]-aXYZ2[2] };
const Bound&    bound1     = vElemBound[i];
for (int j=i+1; j<aBdyElem.GetSize(); ++j)
{
const T_ELEM_D& ElemD2 = aBdyElem[j];
const Bound&    bound2 = vElemBound[j];
///////////////////////////////////////////
if ((bound1.MaxX < bound2.MinX) || (bound1.MinX > bound2.MaxX)) continue;
if ((bound1.MaxY < bound2.MinY) || (bound1.MinY > bound2.MaxY)) continue;
if ((bound1.MaxZ < bound2.MinZ) || (bound1.MinZ > bound2.MaxZ)) continue;
///////////////////////////////////////////
int aNdKey2[3] = { ElemD2.aNodeK[0], ElemD2.aNodeK[1], ElemD2.aNodeK[2] };
bool bCoNode   = false;
for (int k=0; k<3; ++k)
{
if (find(aNdKey2, aNdKey2+3, aNdKey1[k]) != aNdKey2+3)
{
bCoNode = true;
break;
}
}
if (bCoNode) continue;
///////////////////////////////////////////
const T_NODE_D& NodeD4   = aBdyNode[aNdKey2[0]];
const T_NODE_D& NodeD5   = aBdyNode[aNdKey2[1]];
const T_NODE_D& NodeD6   = aBdyNode[aNdKey2[2]];
double          aXYZ4[3] = { NodeD4.x, NodeD4.y, NodeD4.z };
double          aXYZ5[3] = { NodeD5.x, NodeD5.y, NodeD5.z };
double          aXYZ6[3] = { NodeD6.x, NodeD6.y, NodeD6.z };
for (k=0; k<3; ++k)
{
const double* aXYZ_a;
const double* aXYZ_b;
switch (k)
{
case 0: aXYZ_a = aXYZ4; aXYZ_b = aXYZ5; break;
case 1: aXYZ_a = aXYZ5; aXYZ_b = aXYZ6; break;
case 2: aXYZ_a = aXYZ6; aXYZ_b = aXYZ4; break;
}
if (triangle_line_intersection_3d(aVec12, aVec23, aXYZ1, aXYZ2, aXYZ3, aXYZ_a, aXYZ_b))
{
v_fault_index.push_back(i);
v_fault_index.push_back(j);
break;
}
}
}
}
if (v_fault_index.empty()) return 0;
///////////////////////////////////////////////
sort(v_fault_index.begin(), v_fault_index.end());
v_fault_index.erase(unique(v_fault_index.begin(), v_fault_index.end()), v_fault_index.end());
///////////////////////////////////////////////
return v_fault_index.size();
}
*/

//============================================================================
// Merge1DMesh
// cm2 merge에서 internal edge를 closed식으로 만들지 않으면 무한루프빠진다.
//============================================================================
/*
bool CBaseMeshTool::Merge1DMesh(vector<mesh_packet_2d::Mesh_Packet_2D::I_Node>&       preNodes, 
vector<mesh_packet_2d::Mesh_Packet_2D::I_Node_Chain>& preConnect,
TColStd_SequenceOfInteger& newIndex,
const double MergeTol)
{
M_PROFILE(_T("bool CBaseMeshTool::Merge1DMesh()"));

using namespace mesh_packet_2d;
int i,j,k;
int seedSize     = preNodes.size();
int PreConntSize = preConnect.size();
int connectSize =0;
for(i=0;i<PreConntSize;i++)
{
if(preConnect[i].v_node_index.size() == 0) continue;
if(preConnect[i].t_type != Mesh_Packet_2D::I_Node_Chain::CT_OPEN)
connectSize += preConnect[i].v_node_index.size()-1;
else
connectSize += 2*(preConnect[i].v_node_index.size()-1);
}

cm2::DoubleMat pos         = cm2::DoubleMat(2, seedSize);
cm2::UIntMat   connectB    = cm2::UIntMat(2, connectSize);
cm2::UIntMat   mergedInd;
CMap<UINT, UINT, UINT, UINT> mIsolatedNodeK;
mIsolatedNodeK.InitHashTable(1009);

for(i=0; i<seedSize; ++i)
{
pos(0,i) = preNodes[i].a_uv[0];
pos(1,i) = preNodes[i].a_uv[1];
if(preNodes[i].b_island_flag)
mIsolatedNodeK.SetAt(i, 0);
}
int count=0;
for(i=0;i<PreConntSize;i++)
{
if(preConnect[i].t_type != Mesh_Packet_2D::I_Node_Chain::CT_OPEN)
{
for(j=0;j<preConnect[i].v_node_index.size()-1;j++,count++){
connectB(0,count) = preConnect[i].v_node_index[j];
connectB(1,count) = preConnect[i].v_node_index[j+1];
}
}else{
for(j=0 ; j<preConnect[i].v_node_index.size()-1 ; j++,count++){
connectB(0,count) = preConnect[i].v_node_index[j];
connectB(1,count) = preConnect[i].v_node_index[j+1];
}
for(j=preConnect[i].v_node_index.size()-1; j>0 ; j--,count++){
connectB(0,count) = preConnect[i].v_node_index[j];
connectB(1,count) = preConnect[i].v_node_index[j-1];
}
}
}

// Isolated HardNode에 대해서 처리할 것

bool bMerged = meshtools::merge(pos, connectB, mergedInd, MergeTol, 0);

for(i=0;i<mergedInd.size();i++)
newIndex.Append(mergedInd[i]);

CMap<UINT, UINT, UINT, UINT>  mMergedInd;
CMap<UINT, UINT, UINT, UINT>  mMovedInd;
mMergedInd.InitHashTable(1009);
mMovedInd.InitHashTable(100003);

if(bMerged)
Simplify(pos, mergedInd, cm2::UIntMat(), connectB, mMergedInd, mMovedInd, &mIsolatedNodeK); 

UINT nTarK;  
if(bMerged)
{
preNodes.clear();
for(i=0;i<pos.cols(); i++)
{
Mesh_Packet_2D::I_Node i_node;
i_node.b_hard_flag   = 1;
i_node.b_island_flag = (mIsolatedNodeK.Lookup(i, nTarK))? true:false;
i_node.a_uv[0]       = pos(0,i);
i_node.a_uv[1]       = pos(1,i);
preNodes.push_back(i_node);		
}

for(j=0;j<PreConntSize;j++)
{
for(k=0;k<preConnect[j].v_node_index.size();k++)
{
if(mMergedInd.Lookup(preConnect[j].v_node_index[k], nTarK))
{
preConnect[j].v_node_index[k] = nTarK;
}
else if(mMovedInd.Lookup(preConnect[j].v_node_index[k], nTarK))
{
preConnect[j].v_node_index[k] = nTarK;
}
}
}
}

return true;
}
*/

//============================================================================
// Merge1DMesh
//============================================================================
/*
bool CBaseMeshTool::Merge1DMesh(T_NODE_D_LIST& aNodeD, T_ELEM_D_LIST& aElemD, const double MergeTol)
{
M_PROFILE(_T("bool CBaseMeshTool::Merge1DMesh()"));

const int posSize     = aNodeD.GetSize();
const int connectSize = aElemD.GetSize();
cm2::DoubleMat pos         = cm2::DoubleMat(3, posSize);
cm2::UIntMat   connectB    = cm2::UIntMat(2, connectSize);
//cm2::UIntMat   mergedInd;

for(int i=0;i<posSize; i++)
{
pos(0,i) = aNodeD[i].x;
pos(1,i) = aNodeD[i].y;
pos(2,i) = aNodeD[i].z;
}
for(i=0;i<connectSize; i++)
{
connectB(0,i)=aElemD[i].aNodeK[0];
connectB(1,i)=aElemD[i].aNodeK[1];
}
// merge each coincide points. but no effect on pos matrix.
if(!meshtools::merge(pos, connectB, MergeTol, 0)) return false;

// simplify the pos matrix to eliminate the unreferenced pts.
meshtools::simplify(pos, connectB);

aNodeD.RemoveAll();
aNodeD.SetSize(pos.cols());
for(i=0;i<pos.cols(); i++)
{
T_NODE_D aNode;
aNode.x=pos(0,i);
aNode.y=pos(1,i);
aNode.z=pos(2,i);
aNodeD.SetAt(i, aNode);
}
for(i=0;i<connectSize; i++)
{
aElemD[i].aNodeK[0]=connectB(0,i);
aElemD[i].aNodeK[1]=connectB(1,i);
}
return true;
}
*/

//============================================================================
// MergeTriaMesh
//============================================================================
/*
bool CBaseMeshTool::MergeTriaMesh(CArray<T_NODE_D,T_NODE_D&>&    Nodes,
CArray<T_ELEM_D,T_ELEM_D&>&    Elems,
CMap<UINT, UINT, UINT, UINT>&  mMergedInd,
CMap<UINT, UINT, UINT, UINT>&  mMovedInd,
const double MergeTol)
{
M_PROFILE(_T("bool CBaseMeshTool::MergeTriaMesh()"));

const int posSize     = Nodes.GetSize();
const int connectSize = Elems.GetSize();
cm2::DoubleMat pos         = cm2::DoubleMat(3, posSize);
cm2::UIntMat   connectB    = cm2::UIntMat(3, connectSize);
cm2::UIntMat   mergedInd;

for(int i=0;i<posSize; i++)
{
pos(0,i) = Nodes[i].x;
pos(1,i) = Nodes[i].y;
pos(2,i) = Nodes[i].z;
}
for(i=0;i<connectSize; i++)
{
connectB(0,i)=Elems[i].aNodeK[0];
connectB(1,i)=Elems[i].aNodeK[1];
connectB(2,i)=Elems[i].aNodeK[2];
}
// merge each coincide points. but no effect on pos matrix.
if(!meshtools::merge(pos, connectB, mergedInd, MergeTol, 0)) return false;

// simplify the pos matrix to eliminate the unreferenced pts.
//meshtools::simplify(pos, connectB);
Simplify(pos, mergedInd, connectB, cm2::UIntMat(), mMergedInd, mMovedInd);

Nodes.RemoveAll();
Nodes.SetSize(pos.cols());
for(i=0;i<(int)pos.cols(); i++)
{
T_NODE_D aNode;
aNode.x=pos(0,i);
aNode.y=pos(1,i);
aNode.z=pos(2,i);
Nodes.SetAt(i, aNode);
}
for(i=0;i<connectSize; i++)
{
Elems[i].aNodeK[0]=connectB(0,i);
Elems[i].aNodeK[1]=connectB(1,i);
Elems[i].aNodeK[2]=connectB(2,i);
}
return true;
}
*/

//============================================================================
// Merge2DMesh
//============================================================================
bool CBaseMeshTool::Merge2DMesh(T_MESH_NODE_D_LIST&   Nodes,
                                T_MESH_ELEM_D_LIST&   Elems,
                                CMap<UINT, UINT, UINT, UINT>& mMergedInd,
                                CMap<UINT, UINT, UINT, UINT>& mMovedInd,
                                const double MergeTol,
                                const bool bHighOrder)
{

    int nNodeArraySize;
    if(bHighOrder) nNodeArraySize = 8;
    else nNodeArraySize = 4;
    const int posSize     = Nodes.GetSize();
    const int connectSize = Elems.GetSize();

    cm2::DoubleMat pos         = cm2::DoubleMat(3, posSize);
    cm2::UIntMat   connectB    = cm2::UIntMat(nNodeArraySize, connectSize);
    cm2::UIntVec   mergedInd;

    for(int i=0;i<posSize; i++)
    {
        pos(0,i) = Nodes[i].x;
        pos(1,i) = Nodes[i].y;
        pos(2,i) = Nodes[i].z;
    }
    for(int i=0;i<connectSize; i++)
    {
        connectB(0,i) = Elems[i].aNodeK[0];
        connectB(1,i) = Elems[i].aNodeK[1];
        connectB(2,i) = Elems[i].aNodeK[2];

        if(!bHighOrder)
        {
            if(Elems[i].vtktyp == VTK_QUAD) 
            {
                connectB(3,i) = Elems[i].aNodeK[3];
            }
            else if(Elems[i].vtktyp == VTK_TRIANGLE) 
            {
                connectB(3,i) = connectB(0,i);
            }
        }
        else // High order
        {
            if(Elems[i].vtktyp == VTK_TRIANGLE) 
            {
                connectB(3,i) = connectB(4,i) = connectB(5,i) = connectB(6,i) = connectB(7,i) = connectB(0,i);
            }
            else if(Elems[i].vtktyp == VTK_QUAD) 
            {
                connectB(3,i) = Elems[i].aNodeK[3];
                connectB(4,i) = connectB(5,i) = connectB(6,i) = connectB(7,i) = connectB(0,i);
            }
            else
            {
                connectB(3,i) = Elems[i].aNodeK[3];
                connectB(4,i) = Elems[i].aNodeK[4];
                connectB(5,i) = Elems[i].aNodeK[5];
                if(Elems[i].vtktyp == VTK_QUADRATIC_QUAD) 
                {
                    connectB(6,i) = Elems[i].aNodeK[6];
                    connectB(7,i) = Elems[i].aNodeK[7];
                }
                else if(Elems[i].vtktyp == VTK_QUADRATIC_TRIANGLE) 
                {
                    connectB(6,i) = connectB(7,i) = connectB(0,i);
                }
            }
        }
    }

    // merge each coincide points. but no effect on pos matrix.
    int nMergeCount = meshtools::merge(pos, connectB, mergedInd, MergeTol, 0);
    if(nMergeCount < 0) { ASSERT(0); return false; }
    if(nMergeCount == 0) return true;

    // CM2에서 제공되는 simplify method는 newindex가 제대로 나오지 못하기때문에 
    // 따로 만들었다.
    //Simplify(pos, mergedInd, connectB, cm2::UIntMat(), mMergedInd, mMovedInd);

    Nodes.RemoveAll();
    Nodes.SetSize(pos.cols());
    for(int i=0;i<(int)pos.cols(); i++)
    {
        T_MESH_NODE_D aNode;
        aNode.x=pos(0,i);
        aNode.y=pos(1,i);
        aNode.z=pos(2,i);
        Nodes.SetAt(i, aNode);
    }
    for(int i=0;i<connectSize; i++)
    {
        Elems[i].aNodeK[0]=connectB(0,i);
        Elems[i].aNodeK[1]=connectB(1,i);
        Elems[i].aNodeK[2]=connectB(2,i);

        if(Elems[i].vtktyp == VTK_TRIANGLE) 
        {
        }
        if(Elems[i].vtktyp == VTK_QUAD) 
        {
            Elems[i].aNodeK[3] = connectB(3,i);
        }
        else if(bHighOrder)
        {
            Elems[i].aNodeK[3]=connectB(3,i);
            Elems[i].aNodeK[4]=connectB(4,i);
            Elems[i].aNodeK[5]=connectB(5,i);
            if(Elems[i].vtktyp == VTK_QUADRATIC_QUAD)
            {
                Elems[i].aNodeK[6]=connectB(6,i);
                Elems[i].aNodeK[7]=connectB(7,i);
            }
        }
    }

    return true;
}


//============================================================================
// MergeTetraMesh
//============================================================================
/*
bool CBaseMeshTool::MergeTetraMesh(T_NODE_D_LIST&    Nodes,
T_ELEM_D_LIST&    Elems,
CMap<UINT, UINT, UINT, UINT>&  mMergedInd,
CMap<UINT, UINT, UINT, UINT>&  mMovedInd,
const double                   MergeTol)
{
M_PROFILE(_T("bool CBaseMeshTool::MergeTetraMesh()"));

const int posSize     = Nodes.GetSize();
const int connectSize = Elems.GetSize();
cm2::DoubleMat pos         = cm2::DoubleMat(3, posSize);
cm2::UIntMat   connectB    = cm2::UIntMat(4, connectSize);
cm2::UIntMat   mergedInd;

for(int i=0;i<posSize; i++)
{
pos(0,i) = Nodes[i].x;
pos(1,i) = Nodes[i].y;
pos(2,i) = Nodes[i].z;
}
for(i=0;i<connectSize; i++)
{
if(Elems[i].vtktyp != VTK_TETRA) continue;
connectB(0,i)=Elems[i].aNodeK[0];
connectB(1,i)=Elems[i].aNodeK[1];
connectB(2,i)=Elems[i].aNodeK[2];
connectB(3,i)=Elems[i].aNodeK[3];
}
// merge each coincide points. but no effect on pos matrix.
if(!meshtools::merge(pos, connectB, mergedInd, MergeTol, 0)) return false;

// simplify the pos matrix to eliminate the unreferenced pts.
Simplify(pos, mergedInd, connectB, cm2::UIntMat(), mMergedInd, mMovedInd);

Nodes.RemoveAll();
Nodes.SetSize(pos.cols());
for(i=0;i<(int)pos.cols(); i++)
{
T_NODE_D aNode;
aNode.x=pos(0,i);
aNode.y=pos(1,i);
aNode.z=pos(2,i);
Nodes.SetAt(i, aNode);
}
for(i=0;i<connectSize; i++)
{
Elems[i].aNodeK[0] = connectB(0,i);
Elems[i].aNodeK[1] = connectB(1,i);
Elems[i].aNodeK[2] = connectB(2,i);
Elems[i].aNodeK[3] = connectB(3,i);
}
return true;
}
*/

//============================================================================
// RemeshTria
//============================================================================
/*
bool CBaseMeshTool::RemeshTria(const CArray<T_NODE_D,T_NODE_D&>& Nodes,
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
const bool bMerge)
{
M_PROFILE(_T("bool CBaseMeshTool::RemeshTria()"));

int i, j; 

for(i=0;i<Elems.GetSize();i++)
{
if(Elems[i].vtktyp != VTK_TRIANGLE) return false; // 4각형이 있다면 쪼개서 들어올것!
}
//////////////////////////////////////////////////
int nAllNodeCount = Nodes.GetSize();
int nAllElemCount = Elems.GetSize();
int nHardEdge     = aHardEdge.GetSize();

using namespace cm2;
cm2::DoubleMat pos(3, nAllNodeCount);
cm2::UIntMat   connectT3(3, nAllElemCount);
cm2::UIntMat   connectE;

for(i=0; i<nAllNodeCount; i++)
{
pos(0,i) = Nodes[i].x;
pos(1,i) = Nodes[i].y;
pos(2,i) = Nodes[i].z;
}
for(i=0;i<nAllElemCount;i++)
{
for(j=0;j<3;j++) connectT3(j,i)=Elems[i].aNodeK[j];
}

if(nHardEdge > 0) 
{
connectE.resize(2,nHardEdge);
for(i=0; i<nHardEdge; i++)
{
connectE(0, i) = aHardEdge[i].nKey1;
connectE(1, i) = aHardEdge[i].nKey2;
}
}

// HJK : set hard faces
cm2::UIntMat hardfaces;
int nHardFace = aHardFace.GetSize();
hardfaces.resize(nHardFace);
for (i = 0; i < nHardFace; i++) hardfaces[i] = aHardFace[i];

// REMESHING...
surf_remesh_t3::mesher               remesher_t3;
surf_remesh_t3::mesher::data_type    dataT3;
surf_remesh_t3::registration (_T("Licensed to MIDAS_IT."), _T("BEDAEDAF68CC"));
T_PREFERENCE pref;
CDBDoc::GetDocPoint()->GetInitCtrl()->GetPreference(pref);
remesher_t3.settings.optim_level          = pref.Mesh.nOptLev2D;
remesher_t3.settings.shape_quality_weight = float(pref.Mesh.nWeightQlt2D)*.1;
//remesher_t3.settings.sizes_mult_coeff     = dRefineFactor;

dataT3.pos = pos;                   // The initial coordinates matrix.
dataT3.connectM.copy(connectT3);        // The connectivity matrix of the initial mesh.
if (nHardEdge > 0) dataT3.hard_edges = connectE;
if (nHardFace > 0) dataT3.hard_faces = hardfaces;

if(aInnerVtxSizeInd.GetSize() > 0)
{
dataT3.metrics.resize(pos.cols(), 0.);
for(i=0; i<aInnerVtxSizeInd.GetSize(); i++)
{
dataT3.hard_nodes.push_back(aInnerVtxSizeInd[i].nPos);
dataT3.metrics[aInnerVtxSizeInd[i].nPos] = aInnerVtxSizeInd[i].dVal;
}
}
remesher_t3.settings.min_h = dMinH;
remesher_t3.settings.max_h = dMaxH;
remesher_t3.settings.initial_cleanup_flag = false;    // 입력되는 절점 순서를 지켜주기 위해서 필요하다.
remesher_t3.settings.final_optimization_flag = bMerge?true:false;
remesher_t3.settings.merge_tolerance = .05; // merge_tolerance는 상대값
remesher_t3.settings.patch_angle_tolerance = dPatchAng;
//remesher_t3.settings.solid_flag = true; // true시키면 지켜야할 내부 free-edge를 머지시켜버린다.
remesher_t3.settings.interrupt_hdl = interrupt_hdl;

clock_t  StartTime = clock();
remesher_t3.run (dataT3);
if (CProduct::GetTestEnvValue(_T("ProfilingMesh")) == _T("yes")) 
{
CString sTimeMsg;
double dTime = (double(clock() - StartTime) / CLOCKS_PER_SEC);
sTimeMsg.Format(_T("[PROFILING] cm2::Remesh-Tria Core Time : %.5g[sec]"), dTime);
CMECDocBase::GetCurDoc()->OutputMessage(sTimeMsg, 0);	
}

//////////////////////////////////////////////////////////////////////////  
CString sErrMsg;
if (dataT3.warning_code == surf_remesh_t3::mesher::data_type::CM2_INTERRUPTION) 
{
sErrMsg = _LS(IDS_MESH_ERR_INTERRUPTED);// interrupted state
}
else if (dataT3.error_code != surf_remesh_t3::mesher::data_type::CM2_NO_ERROR)
{
if(dataT3.error_code == surf_remesh_t3::mesher::data_type::CM2_REMESHING_ERROR)
{
sErrMsg = _LS(IDS_BASE_FAIL_REMESH_TRIA);
}
else
{
sErrMsg.Format(_T("1. %s\n2. %s"), dataT3.msg1, dataT3.msg2);
}
}
if(!sErrMsg.IsEmpty())
{
AfxMessageBox(sErrMsg);
return false;
}
//////////////////////////////////////////////////////////////////////////

// 아래의 cm2::meshtools::simplify함수를 이용하면 앞부분에 놓여질 지켜야할 hard-edge절점의 index를 지키지 못하므로
// 별도의 simplify함수를 만들어야 한다.
//cm2::meshtools::simplify(dataT3.pos, dataT3.connectM); 
CMap<UINT, UINT, UINT, UINT> mNewNodeKey;
mNewNodeKey.InitHashTable(GF_GetHashSizeByCount(nAllNodeCount));
T_NODE_D NewNodeD;
T_ELEM_D NewElemD;

raNewNodeD.SetSize(dataT3.pos.cols());
int nCount = 0;
for(i=0; i<nHardFlag; i++) 
{
mNewNodeKey.SetAt(i, i);
NewNodeD.Initialize();
NewNodeD.x = dataT3.pos(0,i);
NewNodeD.y = dataT3.pos(1,i);
NewNodeD.z = dataT3.pos(2,i);
NewNodeD.bHard = true;
raNewNodeD.SetAt(nCount++, NewNodeD);
}

int nNewElem = dataT3.connectM.cols();
raNewElemD.SetSize(nNewElem);
NewElemD.Initialize();
NewElemD.vtktyp = VTK_TRIANGLE;
NewElemD.aNodeK.SetSize(3);
for(i=0; i<nNewElem; i++)
{
for(j=0; j<3; j++) 
{
UINT NodeInd = dataT3.connectM(j,i);
UINT NewNodeInd;
if(!mNewNodeKey.Lookup(NodeInd, NewNodeInd))
{
NewNodeInd = nCount;
NewNodeD.x = dataT3.pos(0,NodeInd);
NewNodeD.y = dataT3.pos(1,NodeInd);
NewNodeD.z = dataT3.pos(2,NodeInd);
raNewNodeD.SetAt(nCount++, NewNodeD);
mNewNodeKey.SetAt(NodeInd, NewNodeInd);
}
NewElemD.aNodeK.SetAt(j, NewNodeInd);
}
raNewElemD.SetAt(i, NewElemD);
}
raNewNodeD.SetSize(nCount);
raNewNodeD.FreeExtra();

return true;
}
*/

//============================================================================
// MergeTetraMesh
//============================================================================
/*
bool CBaseMeshTool::RemeshQuad(const CArray<T_NODE_D,T_NODE_D&>& Nodes,
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
const bool bMerge)
{ 
M_PROFILE(_T("bool CBaseMeshTool::RemeshQuad()"));

int i,j;

//////////////////////////////////////////////////
int nAllNodeCount = Nodes.GetSize();
int nAllElemCount = Elems.GetSize();
int nHardEdge     = aHardEdge.GetSize();

using namespace cm2;
cm2::DoubleMat pos(3, nAllNodeCount);
cm2::UIntMat   connectQ4(3, nAllElemCount);
cm2::UIntMat   connectE;

for(i=0; i<nAllNodeCount; i++)
{
pos(0,i) = Nodes[i].x;
pos(1,i) = Nodes[i].y;
pos(2,i) = Nodes[i].z;
}
for(i=0;i<nAllElemCount;i++)
{
for(j=0;j<3;j++) connectQ4(j,i)=Elems[i].aNodeK[j];
}

if(nHardEdge > 0) 
{
connectE.resize(2,nHardEdge);
for(i=0; i<nHardEdge; i++)
{
connectE(0, i) = aHardEdge[i].nKey1;
connectE(1, i) = aHardEdge[i].nKey2;
}
}

// REMESHING...
surf_remesh_q4::mesher               remesher_q4;
surf_remesh_q4::mesher::data_type    dataQ4;
surf_remesh_q4::registration (_T("Licensed to MIDAS_IT."), _T("ACDA93AF66CC")); // 20071214 

T_PREFERENCE pref;
CDBDoc::GetDocPoint()->GetInitCtrl()->GetPreference(pref);
remesher_q4.settings.optim_level          = pref.Mesh.nOptLev2D;
remesher_q4.settings.shape_quality_weight = float(pref.Mesh.nWeightQlt2D)*.1;
//remesher_q4.settings.sizes_mult_coeff     = dRefineFactor;

dataQ4.pos = pos;                   // The initial coordinates matrix.
dataQ4.connectM.copy(connectQ4);        // The connectivity matrix of the initial mesh.
if(nHardEdge>0) dataQ4.hard_edges = connectE;
if(aInnerVtxSizeInd.GetSize() > 0)
{
dataQ4.metrics.resize(pos.cols(), 0.);
for(i=0; i<aInnerVtxSizeInd.GetSize(); i++)
{
dataQ4.hard_nodes.push_back(aInnerVtxSizeInd[i].nPos);
dataQ4.metrics[aInnerVtxSizeInd[i].nPos] = aInnerVtxSizeInd[i].dVal;
}
}
remesher_q4.settings.min_h = dMinH;
remesher_q4.settings.max_h = dMaxH;
remesher_q4.settings.initial_cleanup_flag = false;    // 입력되는 절점 순서를 지켜주기 위해서 필요하다.
remesher_q4.settings.final_optimization_flag = bMerge?true:false;
remesher_q4.settings.merge_tolerance = .05; // merge_tolerance는 상대값
remesher_q4.settings.patch_angle_tolerance = dPatchAng;
remesher_q4.settings.pure_quad_flag = (nMode == 1);
remesher_q4.settings.interrupt_hdl = interrupt_hdl;

clock_t  StartTime = clock();
remesher_q4.run (dataQ4);
if (CProduct::GetTestEnvValue(_T("ProfilingMesh")) == _T("yes")) 
{
CString sTimeMsg;
double dTime = (double(clock() - StartTime) / CLOCKS_PER_SEC);
sTimeMsg.Format(_T("[PROFILING] cm2::Remesh-Tria Core Time : %.5g[sec]"), dTime);
CMECDocBase::GetCurDoc()->OutputMessage(sTimeMsg, 0);	
}

#ifdef _DEBUG
CString sDebugFile;
sDebugFile.Format(_T("c:/temp/cm2_RemeshQuad%s.txt"),__TIME__);
sDebugFile.Remove(_T(':'));
sDebugFile.Insert(1,_T(':'));
dataQ4.save(LPCTSTR(sDebugFile));
#endif

//////////////////////////////////////////////////////////////////////////  
CString sErrMsg;
if (dataQ4.warning_code == surf_remesh_q4::mesher::data_type::CM2_INTERRUPTION) 
{
sErrMsg = _LS(IDS_MESH_ERR_INTERRUPTED);// interrupted state
}
else if (dataQ4.error_code != surf_remesh_q4::mesher::data_type::CM2_NO_ERROR)
{
if(dataQ4.error_code == surf_remesh_q4::mesher::data_type::CM2_REMESHING_ERROR)
{
sErrMsg = _LS(IDS_BASE_FAIL_REMESH_TRIA);
}
else
{
sErrMsg.Format(_T("1. %s\n2. %s"), dataQ4.msg1, dataQ4.msg2);
}
}
if(!sErrMsg.IsEmpty())
{
AfxMessageBox(sErrMsg);
return false;
}
//////////////////////////////////////////////////////////////////////////

// 아래의 cm2::meshtools::simplify함수를 이용하면 앞부분에 놓여질 지켜야할 hard-edge절점의 index를 지키지 못하므로
// 별도의 simplify함수를 만들어야 한다.
//cm2::meshtools::simplify(dataQ4.pos, dataQ4.connectM); 
CMap<int, int, UINT, UINT> mNewNodeKey;
mNewNodeKey.InitHashTable(GF_GetHashSizeByCount(nAllNodeCount));
T_NODE_D NewNodeD;
T_ELEM_D NewElemD;

raNewNodeD.SetSize(dataQ4.pos.cols());
int nCount = 0;
for(i=0; i<nHardFlag; i++) 
{
mNewNodeKey.SetAt(i, i);
NewNodeD.Initialize();
NewNodeD.x = dataQ4.pos(0,i);
NewNodeD.y = dataQ4.pos(1,i);
NewNodeD.z = dataQ4.pos(2,i);
NewNodeD.bHard = true;
raNewNodeD.SetAt(nCount++, NewNodeD);
}

int nNewElem = dataQ4.connectM.cols();
raNewElemD.SetSize(nNewElem);
NewElemD.Initialize();
for(i=0; i<nNewElem; i++)
{
NewElemD.vtktyp = VTK_QUAD;
NewElemD.aNodeK.SetSize(4);
for(j=0; j<4; j++) 
{
int NodeInd = dataQ4.connectM(j,i);
if(NodeInd < 0)
{
if(j!=3) { ASSERT(0); return false; }
NewElemD.vtktyp = VTK_TRIANGLE;
NewElemD.aNodeK.SetSize(3);
break;
}
UINT NewNodeInd;
if(!mNewNodeKey.Lookup(NodeInd, NewNodeInd))
{
NewNodeInd = nCount;
NewNodeD.x = dataQ4.pos(0,NodeInd);
NewNodeD.y = dataQ4.pos(1,NodeInd);
NewNodeD.z = dataQ4.pos(2,NodeInd);
raNewNodeD.SetAt(nCount++, NewNodeD);
mNewNodeKey.SetAt(NodeInd, NewNodeInd);
}
NewElemD.aNodeK.SetAt(j, NewNodeInd);
}
raNewElemD.SetAt(i, NewElemD);
}
raNewNodeD.SetSize(nCount);
raNewNodeD.FreeExtra();

return true;
}
*/

//============================================================================
// AdaptiveMeshTetra
// HJK : Because background mesh of CM2 remesh is not supported yet, hard node
// option is use to control mesh size
//============================================================================
/*
int CBaseMeshTool::AdaptiveMeshTria(const CArray<T_NODE_D,T_NODE_D&>& aNodes,
const CArray<T_ELEM_D,T_ELEM_D&>& aElems,
//const CArray<T_HARDEDGE, T_HARDEDGE&>& aHardEdge,
const double dSize,
const CMap<UINT, UINT, double, double>& mSize,
CArray<T_NODE_D,T_NODE_D&>& raNewNodeD,
CArray<T_ELEM_D,T_ELEM_D&>& raNewElemD,
const double dRefineFactor,
bool (*interrupt_hdl)(),
CArray<UINT, UINT>& raFailedInd,
const bool bMerge)
{
M_PROFILE(_T("bool CBaseMeshTool::AdaptiveMeshTria()"));
int i, j;
int nNode        = aNodes.GetSize();
int nElem         = aElems.GetSize();
cm2::DoubleMat pos    = cm2::DoubleMat(3, nNode);
cm2::UIntMat connectB = cm2::UIntMat(3, nElem);
cm2::UIntMat indicies;
FloatVec sizes;
surf_remesh_t3::mesher   remesher_t3;
surf_remesh_t3::mesher::data_type dataT3;
surf_remesh_t3::registration (_T("Licensed to MIDAS_IT."), _T("BEDAEDAF68CC"));
cm2::UIntMat hard_nodes;

// set nodes and hard nodes
dataT3.metrics.resize(nNode, 5.0);
hard_nodes.resize(nNode);
for(i=0;i<nNode;i++)
{
pos(0, i) = aNodes[i].x;
pos(1, i) = aNodes[i].y;
pos(2, i) = aNodes[i].z;

double dTmpSize;
if (mSize.Lookup(i, dTmpSize))
{
dataT3.metrics[i] = dTmpSize;
}
else
{
TRACE(_T("no size hard node:%d\n"), i+1);
}
hard_nodes[i] = i;
}

// set connectivities and hard nodes
CMapEx<UINT, UINT&, int, int&> mHardNode;
for (i=0;i<nElem;i++)
{
for (j=0;j<3;j++)
{
connectB(j, i) = aElems[i].aNodeK[j];
}
}

// remesh bondary
dataT3.pos = pos;                   // The initial coordinates matrix.
dataT3.connectM.copy(connectB);        // The connectivity matrix of the initial mesh.
dataT3.hard_nodes = hard_nodes;
remesher_t3.settings.refine_flag = true;
remesher_t3.settings.interrupt_hdl = interrupt_hdl;
remesher_t3.run(dataT3);

CString sErrMsg;
if (dataT3.warning_code == surf_remesh_t3::mesher::data_type::CM2_INTERRUPTION) 
{
sErrMsg = _LS(IDS_MESH_ERR_INTERRUPTED);// interrupted state
}
else if (dataT3.error_code != surf_remesh_t3::mesher::data_type::CM2_NO_ERROR)
{
if(dataT3.error_code == surf_remesh_t3::mesher::data_type::CM2_REMESHING_ERROR)
{
sErrMsg = _LS(IDS_BASE_FAIL_REMESH_TRIA);
}
else
{
sErrMsg.Format(_T("1. %s\n2. %s"), dataT3.msg1, dataT3.msg2);
}
}
if(!sErrMsg.IsEmpty())
{
AfxMessageBox(sErrMsg);
return false;
}
//////////////////////////////////////////////////////////////////////////

// 아래의 cm2::meshtools::simplify함수를 이용하면 앞부분에 놓여질 지켜야할 hard-edge절점의 index를 지키지 못하므로
// 별도의 simplify함수를 만들어야 한다.
//cm2::meshtools::simplify(dataT3.pos, dataT3.connectM); 
CMap<UINT, UINT, UINT, UINT> mNewNodeKey;
mNewNodeKey.InitHashTable(GF_GetHashSizeByCount(nNode));
T_NODE_D NewNodeD;
T_ELEM_D NewElemD;

nNode = dataT3.pos.cols();
raNewNodeD.SetSize(nNode);
int nCount = 0;
for(i=0; i<nNode; i++) 
{
mNewNodeKey.SetAt(i, i);
NewNodeD.Initialize();
NewNodeD.x = dataT3.pos(0,i);
NewNodeD.y = dataT3.pos(1,i);
NewNodeD.z = dataT3.pos(2,i);
raNewNodeD.SetAt(nCount++, NewNodeD);
}

int nNewElem = dataT3.connectM.cols();
raNewElemD.SetSize(nNewElem);
NewElemD.Initialize();
NewElemD.vtktyp = VTK_TRIANGLE;
NewElemD.aNodeK.SetSize(3);
for(i=0; i<nNewElem; i++)
{
for(j=0; j<3; j++) 
{
NewElemD.aNodeK.SetAt(j, dataT3.connectM(j,i));
}
raNewElemD.SetAt(i, NewElemD);
}
raNewNodeD.SetSize(nCount);
raNewNodeD.FreeExtra();

return 1;
}
*/

//============================================================================
// AdaptiveMeshTetra
// HJK
//============================================================================
/*
int CBaseMeshTool::AdaptiveMeshTetra(const CArray<T_NODE_D,T_NODE_D&>& aNodes,
const CArray<T_ELEM_D,T_ELEM_D&>& aBGMElems,
const CArray<T_ELEM_D,T_ELEM_D&>& aBndrElems,
const double dSize,
const CMap<UINT, UINT, double, double>& mBGMSize,
CArray<T_NODE_D,T_NODE_D&>& raNewNodeD,
CArray<T_ELEM_D,T_ELEM_D&>& raNewElemD,
const double dRefineFactor,
bool (*interrupt_hdl)(),
CArray<UINT, UINT>& raFailedInd,
const bool bMerge,
const bool bNoClamped)
{
M_PROFILE(_T("bool CBaseMeshTool::AdaptiveMeshTetra()"));
int i, j;
int nNode        = aNodes.GetSize();
int nBndr        = aBndrElems.GetSize();
int nBGM         = aBGMElems.GetSize();
cm2::DoubleMat pos    = cm2::DoubleMat(3, nNode);
cm2::UIntMat connectB = cm2::UIntMat(3, nBndr);
cm2::UIntMat connectB1;
cm2::UIntMat BGM      = cm2::UIntMat(4, nBGM);
cm2::UIntMat indicies;
FloatVec sizes;
surf_remesh_t3::mesher   remesher_t3;
surf_remesh_t3::mesher::data_type dataT3;
surf_remesh_t3::registration (_T("Licensed to MIDAS_IT."), _T("BEDAEDAF68CC"));
cm2::UIntMat hard_nodes;

for(i=0;i<nNode;i++)
{
pos(0, i) = aNodes[i].x;
pos(1, i) = aNodes[i].y;
pos(2, i) = aNodes[i].z;
}

for (i=0;i<nBndr;i++)
{
for (j=0;j<3;j++) connectB(j, i) = aBndrElems[i].aNodeK[j];
}

for (i=0;i<nBGM;i++) // CM2 is inverted to our tetra element numbering
{
BGM(0, i) = aBGMElems[i].aNodeK[0];
BGM(1, i) = aBGMElems[i].aNodeK[2];
BGM(2, i) = aBGMElems[i].aNodeK[1];
BGM(3, i) = aBGMElems[i].aNodeK[3];
}

// get hard nodes  
CMapEx<UINT, UINT&, int, int&> mHardNode;
for (i = 0; i < aBndrElems.GetSize(); i++)
{
for (j = 0; j < 3; j++)
{
int dummy = 0;
mHardNode.SetAt(aBndrElems[i].aNodeK[j], dummy);
}
}

hard_nodes.resize(mHardNode.GetCount());
int iNode = 0;
//dataT3.metrics.resize(pos.cols(), 0.0);
dataT3.metrics.resize(pos.cols(), 5.0);
POSITION posHardNode = mHardNode.GetStartPosition();

while (posHardNode)
{
UINT tmpNodeK;
mHardNode.GetNextKey(posHardNode, tmpNodeK);
hard_nodes[iNode++] = tmpNodeK;
double dSize;
if (mBGMSize.Lookup(tmpNodeK, dSize))
{
dataT3.metrics[tmpNodeK] = dSize;
}
else
{
TRACE(_T("no size hard node:%d\n"), tmpNodeK+1);
}
}

// remesh bondary
dataT3.pos = pos;                   // The initial coordinates matrix.
dataT3.connectM.copy(connectB);        // The connectivity matrix of the initial mesh.
dataT3.hard_nodes = hard_nodes;
remesher_t3.settings.refine_flag = true;
remesher_t3.settings.interrupt_hdl = interrupt_hdl;
remesher_t3.run (dataT3);
dataT3.extract(pos, connectB1);
//meshtools::vtk_output(_T("c:/vtk_win32/cm2/vtk/win32/AdaptiveMeshTetra1.vtk"), pos, connectB1, CM2_FACET3);
//meshtools::vtk_output(_T("c:/vtk_win32/cm2/vtk/win32/AdaptiveMeshTetra2.vtk"), pos, dataT3.connectM, CM2_FACET3);

// set hard edges
cm2::UIntMat   connectE;

tetramesh::mesher              tetramesher;
tetramesh::mesher::data_type   data;

// UNLOCK THE DLL.
tetramesh::registration (_LSX(Licensed to MIDAS_IT.), _T("46EF1F05BB3F"));
data.pos = pos;
data.connectB = connectB1;
//data.background_mesh = BGM;
//data.metrics = sizes;
tetramesher.settings.default_numbering_flag = false;
//tetramesher.settings.compute_Qh_flag = true;
//tetramesher.settings.optim_level = 3;

tetramesher.settings.no_clamped_tetra_flag = bNoClamped;
tetramesher.settings.interrupt_hdl        = interrupt_hdl;

tetramesher.run (data);
//meshtools::vtk_output(_T("c:/vtk_win32/cm2/vtk/win32/AdaptiveMeshTetra3.vtk"), data.pos, data.connectM, CM2_TETRA4);
//meshtools::vtk_output(_T("c:/vtk_win32/cm2/vtk/win32/AdaptiveMeshTetra4.vtk"), data.pos, data.background_mesh, CM2_TETRA4);

if (data.warning_code == tetramesh::mesher::data_type::CM2_INTERRUPTION) return 0;

if (data.error_code != tetramesh::mesher::data_type::CM2_NO_ERROR)
{
CString str;
switch(data.error_code)
{
case tetramesh::mesher::data_type::CM2_DATA_ERROR:
case tetramesh::mesher::data_type::CM2_INTERNAL_ERROR:
str = APPENDMSG_RECOMMEND_SIZE;
break;
case tetramesh::mesher::data_type::CM2_BOUNDARY_ERROR:
str = MSG_BOUNDARY3D_ERROR;
break;
case tetramesh::mesher::data_type::CM2_NODE_ERROR:
str = MSG_NODE_ERROR;
str += APPENDMSG_RECOMMEND_SIZE;
break;
case tetramesh::mesher::data_type::CM2_EDGE_ERROR:
str = MSG_EDGE_ERROR;
str += APPENDMSG_RECOMMEND_SIZE;
break;
case tetramesh::mesher::data_type::CM2_FACE_ERROR:
str = MSG_FACE_ERROR;
str += APPENDMSG_RECOMMEND_SIZE;
break;
default:
str.Format(_T("1. %s\n2. %s"), data.msg1, data.msg2);
break;
}
MessageBox(NULL, str, _LS(IDS_MESHGEN___MSG11), MB_ICONWARNING);

cm2::UIntMat vFailedBdy = data.unenforced_boundary_IDs;
raFailedInd.RemoveAll();
raFailedInd.SetSize(vFailedBdy.size());
for(i=0; i<vFailedBdy.size(); i++)
{
raFailedInd.SetAt(i, vFailedBdy[i]);
}
return -1;
}

////////////////////////////////////////////////////////////////////////////////////////
// simplify 함수에서 Tetra일때 connectivity를 update하는데에 문제있는것 같음.(20050211)
// 라이브러리가 수정되기 전까지 임시로 걸러내는 루틴을 만들어서 사용한다.
//int nChkSim = meshtools::simplify(data.pos, data.connectM); 
CMap<UINT, UINT, int, int> mNodeInd;
int nIndCount = 0;
int dum;

mNodeInd.InitHashTable(GF_GetHashSizeByCount(data.pos.cols()));
for(i=0; i<(int)data.connectM.cols(); i++)
{
for(j=0; j<4; j++)
{
UINT nNodeInd = data.connectM(j,i);
if(!mNodeInd.Lookup(nNodeInd, dum))
{
mNodeInd.SetAt(nNodeInd, 0);
nIndCount++;
}
}
}
int nPrevNodeCount = data.pos.cols();
if(nPrevNodeCount != nIndCount)
{
bool bChkBdyNodeShift = false;
bool bShift = true;
int nNewCount = 0;
cm2::DoubleMat TmpPos(3, nIndCount);
for(i=0; i<nPrevNodeCount; i++)
{
if(mNodeInd.Lookup(i, dum))
{
if(!bChkBdyNodeShift && i != nNewCount)
{
if(i<pos.cols()) { bShift = false; break; }
bChkBdyNodeShift = true;
}
TmpPos(0, nNewCount) = data.pos(0,i);
TmpPos(1, nNewCount) = data.pos(1,i);
TmpPos(2, nNewCount) = data.pos(2,i);
mNodeInd.SetAt(i, nNewCount);
nNewCount++;
}
}
if(bShift)
{
data.pos.clear_hard();
data.pos.copy(TmpPos);
for(i=0; i<data.connectM.cols(); i++)
{
data.connectM(0,i) = mNodeInd[data.connectM(0,i)];
data.connectM(1,i) = mNodeInd[data.connectM(1,i)];
data.connectM(2,i) = mNodeInd[data.connectM(2,i)];
data.connectM(3,i) = mNodeInd[data.connectM(3,i)];
}
}
}
////////////////////////////////////////////////////////////////////////////////////////

raNewNodeD.RemoveAll();
raNewElemD.RemoveAll();

raNewNodeD.SetSize(data.pos.cols());

T_NODE_D  node;
for(i=0; i<(int)data.pos.cols(); i++)
{
node.x = data.pos(0,i);
node.y = data.pos(1,i);
node.z = data.pos(2,i);
raNewNodeD.SetAt(i, node);
}

raNewElemD.SetSize(data.connectM.cols());
T_ELEM_D elem;
elem.Initialize();
elem.vtktyp = VTK_TETRA;
elem.aNodeK.SetSize(4);
for(i=0; i<(int)data.connectM.cols(); i++)
{
for(j=0; j<4; j++) elem.aNodeK.SetAt(j, data.connectM(j,i));
raNewElemD.SetAt(i, elem);
}

return 1;
}
*/

// Auxiliary function to generate the 6 faces of a cube centered at (0,0,0),
// with edge length equal to _T("L") and with _T("N") elements along each edge.
/*
void CBaseMeshTool::cube_boundary (double L, unsigned N, cm2::DoubleMat& pos, cm2::UIntMat& connectS)
{
cm2::UIntMat indices;
cm2::UIntMat connectE, connectS1, connectS2;
meshtools1D::extrude_translate (pos, DoubleVec3(-L/2,-L/2,+L/2),
DoubleVec3(L,0,0), N, indices);
meshtools1D::indices_to_connectE2 (indices, connectE);
meshtools2D::extrude_translate_T3 (pos, connectE, DoubleVec3(0,L,0), N,
true, connectS2);
meshtools::copy_mesh (pos, connectS1, connectS2);
meshtools::rotate (pos, DoubleVec3(0.0), DoubleVec3(0,M_PI,0),
connectS1);
meshtools::append (connectS2, connectS1);
meshtools::copy_mesh (pos, connectS1, connectS2);
meshtools::rotate (pos, DoubleVec3(0.0), DoubleVec3(M_PI/2,0,0),
connectS1);
meshtools::append (connectS2, connectS1);
meshtools::copy_mesh (pos, connectS1, connectS1);
meshtools::rotate (pos, DoubleVec3(0.0), DoubleVec3(0,0,M_PI/2),
connectS1);
meshtools::append (connectS2, connectS1);
meshtools::merge (pos, connectS2);
meshtools::append (connectS, connectS2);
}

void CBaseMeshTool::display_hdl(unsigned level, const char* msg) { TRACE(msg); }
*/
#endif // !defined(_NONE_MESH)
