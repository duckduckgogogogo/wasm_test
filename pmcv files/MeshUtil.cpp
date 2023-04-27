#include "stdafx.h"
#include "MeshUtil.h"

void T_MESH_NODE_D::Initialize()
{
    x = y = z = 0.0;
    CordK = 1;
    nNo = 0;
    bPost = false;
    bHard = false;
}

void T_MESH_ELEM_D::Initialize()
{
    vtktyp = 0;
    eltyp = 0;
    nKind = 0;
    elmat = 0;
    elpro = 0;
    nNo = 0;
    aNodeK.RemoveAll();
    aInfoK.RemoveAll();
    CordK = 1;
    bPost = false;
}

CMeshUtil::CMeshUtil()
{
}

CMeshUtil::~CMeshUtil()
{
}

bool CMeshUtil::GetMidNode(const T_MESH_NODE_D& Node1, const T_MESH_NODE_D& Node2, T_MESH_NODE_D& MidNode)
{
    MidNode.x = ( Node1.x + Node2.x ) * 0.5;
    MidNode.y = ( Node1.y + Node2.y ) * 0.5;
    MidNode.z = ( Node1.z + Node2.z ) * 0.5;
    return true;
}

bool CMeshUtil::checkESC_interruption(void* pass_thru, double progress)
{
    static MSG msg;
    UINT   message = WM_KEYDOWN;
    if ( ::PeekMessage(&msg, NULL, message, message, PM_NOREMOVE) )
    {
        ::GetMessage(&msg, NULL, message, message);
        if ( ( int ) msg.wParam == VK_ESCAPE ) return true;
    }
    return false;
}

double CMeshUtil::CalcLength(const T_MESH_NODE_D& N1, const T_MESH_NODE_D& N2)
{
    double dx = N2.x - N1.x;
    double dy = N2.y - N1.y;
    double dz = N2.z - N1.z;
    return sqrt(dx*dx + dy*dy + dz*dz);
}