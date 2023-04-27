// mesh_packet_2d.h: data packet for 2D auto-meshers and domain-maker

//////////////////////////////////////////////////////////////////////////////////
// Shin Dae-Seock (dsshin@midasIT.com)
// Don't believe in magic;
// understand what your libraries do, how they do it, and at what cost they do it.
//////////////////////////////////////////////////////////////////////////////////

#ifndef MESH_PACKET_2D_H_
#define MESH_PACKET_2D_H_

#include <vector>

/////////////////////////////////////////////////////////////////////
namespace mesh_packet_2d
{
    using namespace std;
    /////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////
    struct Mesh_Packet_2D
    {
        enum MESH_TYPE { MT_TRIA = 1, MT_COMB, MT_QUAD };
        /////////////////////////////////////////
        struct I_Node
        {
            bool   b_hard_flag;
            bool   b_island_flag;
            double d_mesh_size;
            double a_uv[2];	
            /////////////////////////////////////
            I_Node()
            {
                b_hard_flag   = false;
                b_island_flag = false;
                d_mesh_size   = 0.0;
                a_uv[0] = a_uv[1] = 0.0;
            }
        }; // end: struct Mesh_Packet_2D::I_Node
        /////////////////////////////////////////
        struct I_Element
        {
            vector<int> v_node_index;
        }; // end: struct Mesh_Packet_2D::I_Element 
        /////////////////////////////////////////
        struct I_Node_Chain
        {
            enum CHAIN_TYPE { CT_OUTER = 1, CT_INNER, CT_OPEN };
            /////////////////////////////////////
            int         n_tag;
            vector<int> v_tag;
            CHAIN_TYPE  t_type;
            vector<int> v_node_index;
        }; // end: struct Mesh_Packet_2D::I_Node_Chain
        /////////////////////////////////////////
        struct I_Domain
        {
            vector<int> v_node_chain_index;
            vector<int> v_island_node_index;
        }; // end: struct Mesh_Packet_2D::I_Domain
        /////////////////////////////////////////
        Mesh_Packet_2D()
        {
            d_mesh_size        = 0.0;
            t_mesh_type        = MT_QUAD;
            d_scale_factor     = 1.0;
            b_offset_flag      = true;
            b_relaxation_flag  = true;
            d_refinement_ratio = 1.0;
        }
        /////////////////////////////////////////
        double               d_mesh_size;
        MESH_TYPE            t_mesh_type;
        double               d_scale_factor;
        bool                 b_offset_flag;
        bool                 b_relaxation_flag;
        double               d_refinement_ratio;
        vector<I_Node>       v_node;
        vector<I_Element>    v_element;
        vector<I_Node_Chain> v_node_chain;
        vector<I_Domain>     v_domain;
        // data only for ortho_mesher
        vector<double>       v_uv_grid     [2];
        vector<double>       v_uv_mesh_size[2];
    }; // end: struct Mesh_Packet_2D
} // end: namespace mesh_packet_2d
/////////////////////////////////////////////////////////////////////

#endif // MESH_PACKET_2D_H_