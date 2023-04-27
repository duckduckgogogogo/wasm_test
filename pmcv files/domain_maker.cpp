// domain_maker.cpp: domain-maker for 2D domain makers

//////////////////////////////////////////////////////////////////////////////////
// Shin Dae-Seock (dsshin@midasIT.com)
// Don't believe in magic;
// understand what your libraries do, how they do it, and at what cost they do it.
//////////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"

#include "domain_maker.h"
#include "Math_Lib.h"
#include "Mathfunc.h"

#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <functional>

#pragma warning ( disable : 4996 )

/////////////////////////////////////////////////////////////////////
using namespace domain_maker;
using namespace std::placeholders;
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
#ifdef _DEBUG
int Edge::_n_edge_id = 1;
#endif
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
bool Domain::is_internal_point(const double* a_uv) const
{
    using namespace mesher_math_lib;
    /////////////////////////////////////////////
    for (vector<Wire*>::const_iterator itr_wr=_v_p_wire.begin(); itr_wr!=_v_p_wire.end(); ++itr_wr)
    {
        Wire* p_wire = *itr_wr;
        /////////////////////////////////////////
        if (p_wire->get_type() == Wire::WT_OPEN) continue;
        /////////////////////////////////////////
        const vector<Node*>& v_p_wire_node = p_wire->get_node();
        double               d_angle       = 0.0;
        for (int i=0; i<(int)v_p_wire_node.size()-1; ++i)
        {				
            Node* p_node_1 = v_p_wire_node[i  ];
            Node* p_node_2 = v_p_wire_node[i+1];
            d_angle       += angle_2d(a_uv, p_node_1->get_uv(), p_node_2->get_uv());
        }
        d_angle = fabs(d_angle);
        /////////////////////////////////////////
        if (p_wire->get_type() == Wire::WT_OUTER && d_angle < 350.0) return false;
        if (p_wire->get_type() == Wire::WT_INNER && d_angle >  10.0) return false;
    }
    /////////////////////////////////////////////
    return true;
} // end: Domain::is_internal_point()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
Domain_Maker::Domain_Maker()
{
    using namespace mesher_math_lib;
    /////////////////////////////////////////////
    _d_tolerance = 0.0;
} // class Domain_Maker - constructor
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
Domain_Maker::~Domain_Maker()
{
    for (vector<Node*>  ::const_iterator itr_nd=_v_p_node  .begin(); itr_nd!=_v_p_node  .end(); ++itr_nd) delete *itr_nd;
    for (vector<Edge*>  ::const_iterator itr_eg=_v_p_edge  .begin(); itr_eg!=_v_p_edge  .end(); ++itr_eg) delete *itr_eg;
    for (vector<Wire*>  ::const_iterator itr_wr=_v_p_wire  .begin(); itr_wr!=_v_p_wire  .end(); ++itr_wr) delete *itr_wr;
    for (vector<Domain*>::const_iterator itr_dm=_v_p_domain.begin(); itr_dm!=_v_p_domain.end(); ++itr_dm) delete *itr_dm;
    /////////////////////////////////////////////
    _v_p_node         .clear();
    _v_p_edge         .clear();
    _v_p_wire         .clear();
    _v_p_domain       .clear();
    _v_p_island_node  .clear();
    _m_node_index     .clear();
    _m_node_chain_node.clear();
} // end: class Domain_Maker::destructor
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Domain_Maker::make_domain(Mesh_Packet_2D& mesh_packet, vector<int>& v_unused_node_chain_index, vector<int>& v_unused_node_index, bool b_keep_inner_domain_flag, bool b_skip_inner_edge_flag, bool (*p_interrupt_handler)()) throw(exception)
{
    _p_interrupt_handler = p_interrupt_handler;
    /////////////////////////////////////////////
#ifdef _DEBUG
    dump_input_file(_T("c:/temp/domain.inp"), mesh_packet);
#endif
    /////////////////////////////////////////////
    try
    {
        set_input_data(mesh_packet);
        make_single_edge_domain();
        /////////////////////////////////////////
        find_open_edge();
        find_irreversible_edge();
        duplicate_reversible_edge();
        /////////////////////////////////////////
        Edge* p_edge = 0;
        while (p_edge = pick_wire_head_edge()) make_multi_edge_domain(p_edge);
        /////////////////////////////////////////
        if (_v_p_domain.empty())
        {
            throw exception("No area can be defined from the selected edge(s).");
        }
        /////////////////////////////////////////
        set_inner_domain(b_keep_inner_domain_flag);
        /////////////////////////////////////////
        if (!b_skip_inner_edge_flag) insert_internal_wire_node(v_unused_node_index);
        finalize_domain_making(mesh_packet, v_unused_node_chain_index);
    }
    catch (...)
    {
        throw;
    }
} // end: Domain_Maker::make_domain()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Domain_Maker::set_input_data(const Mesh_Packet_2D& mesh_packet) throw(exception)
{
    using namespace mesher_math_lib;
    /////////////////////////////////////////////
    bool b_scale_factor = fabs(mesh_packet.d_scale_factor-1.0) > 1.0e-6;
    /////////////////////////////////////////////
    map<int, Node*> m_p_node;
    double          a_uv[2];
    for (int i=0; i<(int)mesh_packet.v_node.size(); ++i)
    {
        const Mesh_Packet_2D::I_Node& i_node = mesh_packet.v_node[i];
        copy(i_node.a_uv, i_node.a_uv+2, a_uv);
        if (b_scale_factor) transform(a_uv, a_uv+2, a_uv, std::bind(multiplies<double>(), _1, mesh_packet.d_scale_factor));
        Node*                         p_node = new Node(a_uv, i+1, i_node.b_hard_flag, i_node.d_mesh_size);
        _v_p_node.push_back(p_node);
        if (i_node.b_island_flag) _v_p_island_node.push_back(p_node);
        m_p_node[i]                         = p_node;
        _m_node_index[p_node]                = i;
    }
    for (int i=0; i<(int)mesh_packet.v_node_chain.size(); ++i)
    {
        const Mesh_Packet_2D::I_Node_Chain& i_node_chain = mesh_packet.v_node_chain[i];
        Edge*                               p_edge       = new Edge(i_node_chain.n_tag);
        vector<Node*>                       v_p_edge_node;
        for (int j=0; j<(int)i_node_chain.v_node_index.size(); ++j) v_p_edge_node.push_back(m_p_node[i_node_chain.v_node_index[j]]);
        p_edge->set_node(v_p_edge_node);
        _v_p_edge.push_back(p_edge);
        /////////////////////////////////////////
        _m_node_chain_node[i] = v_p_edge_node;
    }
    /////////////////////////////////////////////
    try
    {
        check_input_data();
    }
    catch (...)
    {
        throw;
    }
} // end: Domain_Maker::set_input_data()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Domain_Maker::check_input_data() throw(exception)
{
    using namespace mesher_math_lib;
    /////////////////////////////////////////////
    // 1.e-7보다 작고 1.e-10보다 큰 간격이 있는 것을 처리할 수 있도록
    // ZERO(1.e-10)을 쓰지 않고, 상대적인 거리 Tolerance를 이용하기 위하여
    // dMinDist, dRegardSameTol 이용
    double dMinDist = DBL_MAX;
    vector<Edge*>::const_iterator itr_eg;
    for (itr_eg=_v_p_edge.begin(); itr_eg!=_v_p_edge.end(); ++itr_eg)
    {
        if (check_interrupt()) throw interrupt_exception();
        ///////////////////////////////////////////
        Edge*               p_edge        = *itr_eg;
        const vector<Node*> v_p_edge_node = p_edge->get_node();
        for (int i=0; i<(int)v_p_edge_node.size()-1; ++i)
        {
            Node*  p_node_1   = v_p_edge_node[i  ];
            Node*  p_node_2   = v_p_edge_node[i+1];
            double d_distance = distance_2d(p_node_1->get_uv(), p_node_2->get_uv());
            /////////////////////////////////////
            if (d_distance < ZERO)
            {
                throw exception("Zero-length seed segment exists in the selected edge(s).\nCheck the edges and/or their seeds.");
            }
            dMinDist = __min(dMinDist, d_distance);
        }
    }
    double dRegardSameTol = dMinDist * 1.e-4; //상대적인 거리를 위하여. 1.e-4는 경험상 수치.
    /////////////////////////////////////////////
    _d_tolerance = DBL_MAX;
    for (int i=0; i<(int)_v_p_node.size()-1; ++i)
    {
        Node* p_node_1 = _v_p_node[i];
        for (int j=i+1; j<(int)_v_p_node.size(); ++j)
        {
            if (check_interrupt()) throw interrupt_exception();
            /////////////////////////////////////////
            Node*  p_node_2   = _v_p_node[j];
            if (p_node_1 == p_node_2) continue;
            double d_distance = distance_2d(p_node_1->get_uv(), p_node_2->get_uv());
            if (d_distance > dRegardSameTol) _d_tolerance = __min(_d_tolerance, d_distance);
        }
    }

    // [2009-10-14] Kim, Geun Young (Tel: 2042, gykim@midasit.com) : 단위계 연동하여 Tol 적용
    //_d_tolerance /= 100.0;
    //double dFactor = CDBDoc::GetDocPoint()->m_pUnitCtrl->ConvertTgtUnitData2CurUnit(0, D_UNITSYS_LENGTH_INDEX_M, D_UNITSYS_BASE_LENGTH, 1);
    double dFactor = 1.0;
    if(dFactor <= 0) {ASSERT(0); dFactor=1;}
    _d_tolerance = _d_tolerance / (100*dFactor);

    /////////////////////////////////////////////
    map<Node*, Node*> m_node_pair;
    for (int i=0; i<(int)_v_p_node.size(); ++i) m_node_pair[_v_p_node[i]] = _v_p_node[i];
    /////////////////////////////////////////////
    for (int i=0; i<(int)_v_p_node.size()-1; ++i)
    {
        Node* p_node_1 = _v_p_node[i];
        if (m_node_pair[p_node_1] != p_node_1) continue;
        for (int j=i+1; j<(int)_v_p_node.size(); ++j)
        {
            if (check_interrupt()) throw interrupt_exception();
            /////////////////////////////////////////
            Node* p_node_2 = _v_p_node[j];
            if (m_node_pair[p_node_2] != p_node_2) continue;
            if (distance_2d(p_node_1->get_uv(), p_node_2->get_uv()) < _d_tolerance) m_node_pair[p_node_2] = p_node_1;
        }
    }
    for (itr_eg=_v_p_edge.begin(); itr_eg!=_v_p_edge.end(); ++itr_eg)
    {
        Edge*         p_edge        = *itr_eg;
        vector<Node*> v_p_edge_node = p_edge->get_node();
        for (vector<Node*>::iterator itr_nd=v_p_edge_node.begin(); itr_nd!=v_p_edge_node.end(); ++itr_nd) *itr_nd = m_node_pair[*itr_nd];
        p_edge->set_node(v_p_edge_node);
    }
    for (map< int, vector<Node*> >::iterator itr_cn=_m_node_chain_node.begin(); itr_cn!=_m_node_chain_node.end(); ++itr_cn)
    {
        vector<Node*>& v_p_node = itr_cn->second;
        for (vector<Node*>::iterator itr_nd=v_p_node.begin(); itr_nd!=v_p_node.end(); ++itr_nd) *itr_nd = m_node_pair[*itr_nd];
    }
    for (map<Node*, Node*>::const_iterator itr_mv=m_node_pair.begin(); itr_mv!=m_node_pair.end(); ++itr_mv)
    {
        if (itr_mv->first != itr_mv->second)
        {
            Node* p_node = itr_mv->first;
            _v_p_node.erase(find(_v_p_node.begin(), _v_p_node.end(), p_node));
            delete p_node;
        }
    }
    double line1[2][3], line2[2][3], dCorss[3];
    double dDist=0;
    /////////////////////////////////////////////
    for (int i=0; i<(int)_v_p_edge.size()-1; ++i)
    {
        Edge*               p_edge_1        = _v_p_edge[i];
        const vector<Node*> v_p_edge_node_1 = p_edge_1->get_node();
        for (int j=0; j<(int)v_p_edge_node_1.size()-1; ++j)
        {
            const double* a_uv_1 = v_p_edge_node_1[j  ]->get_uv();
            const double* a_uv_2 = v_p_edge_node_1[j+1]->get_uv();
            for (int k=i+1; k<(int)_v_p_edge.size(); ++k)
            {
                if (check_interrupt()) throw interrupt_exception();
                ///////////////////////////////////////
                Edge*               p_edge_2        = _v_p_edge[k];
                const vector<Node*> v_p_edge_node_2 = p_edge_2->get_node();
                for (int l=0; l<(int)v_p_edge_node_2.size()-1; ++l)
                {
                    const double* a_uv_3 = v_p_edge_node_2[l  ]->get_uv();
                    const double* a_uv_4 = v_p_edge_node_2[l+1]->get_uv();
                    if ((distance_2d(a_uv_1, a_uv_3) < _d_tolerance && distance_2d(a_uv_2, a_uv_4) < _d_tolerance) ||
                        (distance_2d(a_uv_1, a_uv_4) < _d_tolerance && distance_2d(a_uv_2, a_uv_3) < _d_tolerance))
                    {
                        throw exception("Duplicate seed segments.\nCheck the duplicate edges and/or their seeds");
                    }
                    // [2009-10-14] Kim, Geun Young (Tel: 2042, gykim@midasit.com)
                    //if (line_cross_2d(a_uv_1, a_uv_2, a_uv_3, a_uv_4))          
                    line1[0][0] = a_uv_1[0]; line1[0][1] = a_uv_1[1]; line1[0][2] = 0;
                    line1[1][0] = a_uv_2[0]; line1[1][1] = a_uv_2[1]; line1[1][2] = 0;

                    line2[0][0] = a_uv_3[0]; line2[0][1] = a_uv_3[1]; line2[0][2] = 0;
                    line2[1][0] = a_uv_4[0]; line2[1][1] = a_uv_4[1]; line2[1][2] = 0;
                    if (CDgnMathFunc::mathIntersectLine2(line1[0], line1[1], line2[0], line2[1], _d_tolerance, dDist, dCorss) == 1)					
                    {
                        // 한점을 공유하고 있을 경우
                        if( (distance_2d(a_uv_1, a_uv_3) < _d_tolerance) ||
                            (distance_2d(a_uv_1, a_uv_4) < _d_tolerance) ||
                            (distance_2d(a_uv_2, a_uv_3) < _d_tolerance) ||
                            (distance_2d(a_uv_2, a_uv_4) < _d_tolerance) )
                        {
                            continue;
                        }                         
                        throw exception("Crossing edges exist.\nSplit all edges and try again!");
                    }
                }
            }
        }
    }
} // end: Domain_Maker::check_input_data()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Domain_Maker::find_open_edge() const
{
    for (int i=0; i<(int)_v_p_edge.size()-1; ++i)
    {
        Edge*                p_edge_1        = _v_p_edge[i];
        const vector<Node*>& v_p_edge_node_1 = p_edge_1->get_node();
        Node*                p_node_1        = v_p_edge_node_1.front();
        Node*                p_node_2        = v_p_edge_node_1.back();
        for (int j=i+1; j<(int)_v_p_edge.size(); ++j)
        {
            Edge*               p_edge_2        = _v_p_edge[j];
            const vector<Node*> v_p_edge_node_2 = p_edge_2->get_node();
            Node*               p_node_3        = v_p_edge_node_2.front();
            Node*               p_node_4        = v_p_edge_node_2.back();
            /////////////////////////////////////
            if (p_node_3 == p_node_1 || p_node_4 == p_node_1) p_edge_1->add_preceder(p_edge_2);
            if (p_node_3 == p_node_2 || p_node_4 == p_node_2) p_edge_1->add_follower(p_edge_2);
            if (p_node_1 == p_node_3 || p_node_2 == p_node_3) p_edge_2->add_preceder(p_edge_1);
            if (p_node_1 == p_node_4 || p_node_2 == p_node_4) p_edge_2->add_follower(p_edge_1);
        }
    }
    /////////////////////////////////////////////
    Edge* p_edge = 0;
    while (p_edge = pick_leaf_edge())
    {
        for (int i=0; i<p_edge->get_preceder_size(); ++i) (p_edge->get_preceder(i))->remove_traversal(p_edge);
        for (int i=0; i<p_edge->get_follower_size(); ++i) (p_edge->get_follower(i))->remove_traversal(p_edge);
        /////////////////////////////////////////
        p_edge->set_open_edge();
        p_edge->set_irreversible_edge();
    }
    /////////////////////////////////////////////
#ifdef _VERBOSE
    {
        vector<Edge*> v_p_open_edge;
        for (int i=0; i<(int)_v_p_edge.size(); ++i)
        {
            if (_v_p_edge[i]->is_open_edge()) v_p_open_edge.push_back(_v_p_edge[i]);
        }
        if (!v_p_open_edge.empty())
        {
            cout << _T("* open edge [N=") << (int)v_p_open_edge.size() << "] : ";
            for (i=0; i<(int)v_p_open_edge.size(); ++i) cout << v_p_open_edge[i]->get_id() << ' ';
            cout << endl;
        }
    }
#endif
} // end: Domain_Maker::find_open_edge()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Domain_Maker::find_irreversible_edge() const throw(exception)
{
    try
    {
        vector<Edge*> v_p_bulk_edge;
        for (vector<Edge*>::const_iterator itr_eg=_v_p_edge.begin(); itr_eg!=_v_p_edge.end(); ++itr_eg)
        {
            Edge* p_edge = *itr_eg;
            if (!p_edge->is_open_edge()) v_p_bulk_edge.push_back(p_edge);
        }
        sort(v_p_bulk_edge.begin(), v_p_bulk_edge.end());
        /////////////////////////////////////////
        while (!v_p_bulk_edge.empty())
        {
            if (check_interrupt()) throw interrupt_exception();
            //////////////////////////////////////
            vector<Edge*> v_p_outermost_edge;
            /////////////////////////////////////
            find_outermost_edge(v_p_bulk_edge, v_p_outermost_edge);
            /////////////////////////////////////
            set<Edge*> s_p_duplicate_edge;
            vector<Edge*>::iterator itr_eg;
            for (itr_eg=v_p_outermost_edge.begin(); itr_eg!=v_p_outermost_edge.end(); ++itr_eg)
            {
                Edge* p_edge = *itr_eg;
                if (find(itr_eg+1, v_p_outermost_edge.end(), p_edge) != v_p_outermost_edge.end()) s_p_duplicate_edge.insert(p_edge);
            }
            if (!s_p_duplicate_edge.empty())
            {
                for (set<Edge*>::iterator itr_eg=s_p_duplicate_edge.begin(); itr_eg!=s_p_duplicate_edge.end(); ++itr_eg)
                {
                    Edge* p_edge = *itr_eg;
                    p_edge->set_open_edge();
                    p_edge->set_irreversible_edge();
                    for (int i=0; i<p_edge->get_preceder_size(); ++i) (p_edge->get_preceder(i))->remove_traversal(p_edge);
                    for (int i=0; i<p_edge->get_follower_size(); ++i) (p_edge->get_follower(i))->remove_traversal(p_edge);
                    v_p_bulk_edge.erase(find(v_p_bulk_edge.begin(), v_p_bulk_edge.end(), p_edge));
                }
                continue;
            }
            /////////////////////////////////////
            for (itr_eg=v_p_outermost_edge.begin(); itr_eg!=v_p_outermost_edge.end(); ++itr_eg) (*itr_eg)->set_irreversible_edge();
            /////////////////////////////////////
            vector<Edge*> v_p_chain_edge = v_p_outermost_edge;
            for (int i=0; i<(int)v_p_chain_edge.size(); ++i)
            {
                Edge* p_edge_1 = v_p_chain_edge[i];
                for (int j=0; j<p_edge_1->get_preceder_size(); ++j)
                {
                    Edge* p_edge_2 = p_edge_1->get_preceder(j);
                    if (find(v_p_chain_edge.begin(), v_p_chain_edge.end(), p_edge_2) == v_p_chain_edge.end()) v_p_chain_edge.push_back(p_edge_2);
                }
                for (int j=0; j<p_edge_1->get_follower_size(); ++j)
                {
                    Edge* p_edge_2 = p_edge_1->get_follower(j);
                    if (find(v_p_chain_edge.begin(), v_p_chain_edge.end(), p_edge_2) == v_p_chain_edge.end()) v_p_chain_edge.push_back(p_edge_2);
                }
            }
            sort(v_p_chain_edge.begin(), v_p_chain_edge.end());
            v_p_chain_edge.erase(unique(v_p_chain_edge.begin(), v_p_chain_edge.end()), v_p_chain_edge.end());
            /////////////////////////////////////
            v_p_bulk_edge.erase(set_difference(v_p_bulk_edge.begin(), v_p_bulk_edge.end(), v_p_chain_edge.begin(), v_p_chain_edge.end(), v_p_bulk_edge.begin()), v_p_bulk_edge.end());
        }
    }
    catch (...)
    {
        throw;
    }
} // end: Domain_Maker::find_irreversible_edge()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Domain_Maker::find_outermost_edge(const vector<Edge*>& v_p_bulk_edge, vector<Edge*>& v_p_outermost_edge) const throw(exception)
{
    using namespace mesher_math_lib;
    /////////////////////////////////////////////
    double a_max_uv[2] = { -DBL_MAX, -DBL_MAX };
    double a_min_uv[2] = {  DBL_MAX,  DBL_MAX };
    for (vector<Edge*>::const_iterator itr_eg=v_p_bulk_edge.begin(); itr_eg!=v_p_bulk_edge.end(); ++itr_eg)
    {
        Edge*                p_edge        = *itr_eg;
        const vector<Node*>& v_p_edge_node = p_edge->get_node();
        for (int i=0; i<(int)v_p_edge_node.size(); ++i)
        {
            const double* a_uv = v_p_edge_node[i]->get_uv();
            for (int j=0; j<2; ++j)
            {
                a_max_uv[j] = __max(a_max_uv[j], a_uv[j]);
                a_min_uv[j] = __min(a_min_uv[j], a_uv[j]);
            }
        }
    }
    double a_uv_1[2];
    mean_2<double>(a_min_uv, a_max_uv, a_uv_1, 2);
    double d_length  = distance_2d(a_min_uv, a_max_uv);
    double a_uv_2[2] = { a_uv_1[0] + d_length, a_uv_1[1] };
    /////////////////////////////////////////////
    Edge*        p_head_edge     = 0;
    const double a_angle_step[3] = { 12.34567890, 6.789012345, 3.456789012 };
    for (int i=0; i<3; ++i)
    {
        const int n_time = (int)(360.0 / a_angle_step[i]);
        for (int j=0; j<=n_time; ++j)
        {
            if (j > 0) rotation_2d(a_uv_2, a_uv_1, a_angle_step[i]);
            /////////////////////////////////////////
            double d_min_distance = DBL_MAX;
            Node*  p_node_ip_1    = 0;
            Node*  p_node_ip_2    = 0;
            double a_ip_uv[2];
            for (vector<Edge*>::const_iterator itr_eg=v_p_bulk_edge.begin(); itr_eg!=v_p_bulk_edge.end(); ++itr_eg)
            {
                Edge* p_edge = *itr_eg;
                /////////////////////////////////////
                if (p_edge->is_open_edge()) continue;
                /////////////////////////////////////
                const vector<Node*>& v_p_edge_node = p_edge->get_node();
                for (int k=0; k<(int)v_p_edge_node.size()-1; ++k)
                {
                    if (check_interrupt()) throw interrupt_exception();
                    /////////////////////////////////////
                    Node*  p_node_1 = v_p_edge_node[k  ];
                    Node*  p_node_2 = v_p_edge_node[k+1];
                    double a_uv[2];
                    if (line_intersection_2d(a_uv_1, a_uv_2, p_node_1->get_uv(), p_node_2->get_uv(), a_uv))
                    {
                        double d_distance = distance_2d(a_uv, a_uv_2);
                        if (d_distance < d_min_distance)
                        {
                            p_head_edge    = p_edge;
                            p_node_ip_1    = p_node_1;
                            p_node_ip_2    = p_node_2;
                            d_min_distance = d_distance;
                            copy(a_uv, a_uv+2, a_ip_uv);
                        }
                    }
                }
            }
            if (!p_head_edge) continue;
            /////////////////////////////////////////
            const double* a_uv_3 = (p_head_edge->get_start_node())->get_uv();
            const double* a_uv_4 = (p_head_edge->get_end_node()  )->get_uv();
            if (distance_2d(a_uv_3, a_ip_uv) < _d_tolerance || distance_2d(a_uv_4, a_ip_uv) < _d_tolerance)
            {
                p_head_edge = 0;
                continue;
            }
            double a_vector_12[2], a_vector_34[2];
            transform(a_uv_2,                a_uv_2+2,                a_uv_1,                a_vector_12, minus<double>());
            transform(p_node_ip_2->get_uv(), p_node_ip_2->get_uv()+2, p_node_ip_1->get_uv(), a_vector_34, minus<double>());
            if (!is_ccw_rotation_2d(a_vector_12, a_vector_34)) p_head_edge->reverse_direction();
            /////////////////////////////////////////
            if (p_head_edge) break;
        }
        if (p_head_edge) break;
    }
    /////////////////////////////////////////////
    if (!p_head_edge)
    {
        throw exception("Failed to construct planar area due to internal error!");
    }
    /////////////////////////////////////////////
    v_p_outermost_edge.push_back(p_head_edge);
    while (true)
    {
        Edge*                p_preceder        = v_p_outermost_edge.back();
        const vector<Node*>& v_p_preceder_node = p_preceder->get_node();
        Node*                p_node_1          = v_p_preceder_node[(int)v_p_preceder_node.size()-1];
        Node*                p_node_2          = v_p_preceder_node[(int)v_p_preceder_node.size()-2];
        Edge*                p_follower        = 0;
        double               d_max_angle       = 0.0;
        bool                 b_reverse_flag    = false;
        for (int i=0; i<p_preceder->get_follower_size(); ++i)
        {
            Edge*                p_edge        = p_preceder->get_follower(i);
            const vector<Node*>& v_p_edge_node = p_edge->get_node();
            Node*                p_node_3      = v_p_edge_node[0];
            Node*                p_node_4      = v_p_edge_node[1];
            if (p_node_3 == p_node_1)
            {
                double d_angle = positive_angle_2d(p_node_1->get_uv(), p_node_4->get_uv(), p_node_2->get_uv());
                if (d_angle > d_max_angle)
                {
                    p_follower     = p_edge;
                    d_max_angle    = d_angle;
                    b_reverse_flag = false;
                }
            }
            else
            {
                double d_angle = positive_angle_2d(p_node_1->get_uv(), p_node_3->get_uv(), p_node_2->get_uv());
                if (d_angle > d_max_angle)
                {
                    p_follower     = p_edge;
                    d_max_angle    = d_angle;
                    b_reverse_flag = true;
                }
            }
        }
        /////////////////////////////////////////
        if (!p_follower || p_follower == v_p_outermost_edge.front()) break;
        /////////////////////////////////////////
        if (b_reverse_flag) p_follower->reverse_direction();
        /////////////////////////////////////////
        v_p_outermost_edge.push_back(p_follower);
    }
    /////////////////////////////////////////////
#ifdef _VERBOSE
    {
        cout << _T("* outermost edge [N=") << (int)v_p_outermost_edge.size() << "] : ";
        for (vector<Edge*>::const_iterator itr_eg=v_p_outermost_edge.begin(); itr_eg!=v_p_outermost_edge.end(); ++itr_eg) cout << (*itr_eg)->get_id() << ' ';
        cout << endl;
    }
#endif
} // end: Domain_Maker::find_outermost_edge()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
inline Edge* Domain_Maker::pick_leaf_edge() const
{
    for (vector<Edge*>::const_iterator itr_eg=_v_p_edge.begin(); itr_eg!=_v_p_edge.end(); ++itr_eg)
    {
        Edge* p_edge = *itr_eg;
        /////////////////////////////////////////
        if (p_edge->is_open_edge()) continue;
        /////////////////////////////////////////
        if (p_edge->get_preceder_size() == 0 || p_edge->get_follower_size() == 0) return p_edge;
    }
    return 0;
} // end: Domain_Maker::pick_leaf_edge()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Domain_Maker::duplicate_reversible_edge()
{
    map<Edge*, Edge*> m_mate_edge;
    /////////////////////////////////////////////
    vector<Edge*> v_p_edge = _v_p_edge;
    for (vector<Edge*>::const_iterator itr_eg=v_p_edge.begin(); itr_eg!=v_p_edge.end(); ++itr_eg)
    {
        Edge* p_edge_1 = *itr_eg;
        /////////////////////////////////////////
        p_edge_1->clear_traversal();
        /////////////////////////////////////////
        if (!p_edge_1->is_reversible_edge()) continue;
        /////////////////////////////////////////
        Edge* p_edge_2 = new Edge();
        p_edge_2->set_node(p_edge_1->get_node());
        p_edge_2->reverse_direction();
        _v_p_edge.push_back(p_edge_2);
        /////////////////////////////////////////
        m_mate_edge[p_edge_2] = p_edge_1;
        /////////////////////////////////////////
        p_edge_1->set_irreversible_edge();
        p_edge_2->set_irreversible_edge();
    }
    /////////////////////////////////////////////
    for (int i=0; i<(int)_v_p_edge.size()-1; ++i)
    {
        Edge* p_edge_1 = _v_p_edge[i];
        /////////////////////////////////////////
        if (p_edge_1->is_open_edge()) continue;
        /////////////////////////////////////////
        Node* p_node_1 = p_edge_1->get_start_node();
        Node* p_node_2 = p_edge_1->get_end_node();
        for (int j=i+1; j<(int)_v_p_edge.size(); ++j)
        {
            Edge* p_edge_2 = _v_p_edge[j];
            /////////////////////////////////////
            if (p_edge_2->is_open_edge() || m_mate_edge[p_edge_2] == p_edge_1) continue;
            /////////////////////////////////////
            Node* p_node_3 = p_edge_2->get_start_node();
            Node* p_node_4 = p_edge_2->get_end_node();
            /////////////////////////////////////
            if (p_node_3 == p_node_2)
            {
                p_edge_1->add_follower(p_edge_2);
                p_edge_2->add_preceder(p_edge_1);
            }
            if (p_node_4 == p_node_1)
            {
                p_edge_1->add_preceder(p_edge_2);
                p_edge_2->add_follower(p_edge_1);
            }
        }
    }
} // end: Domain_Maker::duplicate_reversible_edge()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Domain_Maker::make_single_edge_domain()
{
    using namespace mesher_math_lib;
    /////////////////////////////////////////////
    //for (vector<Edge*>::iterator itr_eg=_v_p_edge.begin(); itr_eg!=_v_p_edge.end(); ++itr_eg)
    vector<Edge*>::iterator itr_next=_v_p_edge.begin();
    for (vector<Edge*>::iterator itr_eg=_v_p_edge.begin(); itr_eg!=_v_p_edge.end(); itr_eg=itr_next ) // mylee
    {
        itr_next = itr_eg; ++itr_next;

        Edge* p_edge = *itr_eg;
        /////////////////////////////////////////
        if (p_edge->get_start_node() != p_edge->get_end_node()) { continue; } // mylee	
        /////////////////////////////////////////
        const vector<Node*>& v_p_edge_node = p_edge->get_node();
        /////////////////////////////////////////
        const double* a_uv_1 = v_p_edge_node[0]->get_uv();
        const double* a_uv_2 = v_p_edge_node[1]->get_uv();
        double        a_vector_12[2];
        transform(a_uv_2, a_uv_2+2, a_uv_1, a_vector_12, minus<double>());
        double        d_norm = norm_2d(a_vector_12);
        transform(a_vector_12, a_vector_12+2, a_vector_12, std::bind(divides<double>(), _1, d_norm));
        double        a_uv[2];
        mean_2<double>(a_uv_1, a_uv_2, a_uv, 2);
        a_uv[0]             -= a_vector_12[1] * _d_tolerance;
        a_uv[1]             += a_vector_12[0] * _d_tolerance;
        /////////////////////////////////////////
        double d_angle = 0.0;
        for (int i=0; i<(int)v_p_edge_node.size()-1; ++i)
        {				
            Node* p_node_1 = v_p_edge_node[i  ];
            Node* p_node_2 = v_p_edge_node[i+1];
            d_angle       += angle_2d(a_uv, p_node_1->get_uv(), p_node_2->get_uv());
        }
        if (fabs(d_angle) < 350.0) p_edge->reverse_direction();
        /////////////////////////////////////////
        Wire* p_wire = new Wire();
        p_wire->set_type(Wire::WT_OUTER);
        p_wire->set_node(p_edge->get_node());
        _v_p_wire.push_back(p_wire);
        /////////////////////////////////////////
        Domain* p_domain = new Domain();
        p_domain->add_wire(p_wire);
        _v_p_domain.push_back(p_domain);
        /////////////////////////////////////////
        delete p_edge;
        //_v_p_edge.erase(itr_eg--);
        itr_next = _v_p_edge.erase(itr_eg); // mylee
    }
} // end: Domain_Maker::make_single_edge_domain()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Domain_Maker::make_multi_edge_domain(Edge* p_head_edge) throw(exception)
{
    using namespace mesher_math_lib;
    /////////////////////////////////////////////
    vector<Edge*> v_p_wire_edge;
    v_p_wire_edge.push_back(p_head_edge);
    while (true)
    {
        if (check_interrupt()) throw interrupt_exception();
        ///////////////////////////////////////////
        Edge*                p_preceder        = v_p_wire_edge.back();
        const vector<Node*>& v_p_preceder_node = p_preceder->get_node();
        Node*                p_node_1          = v_p_preceder_node[(int)v_p_preceder_node.size()-1];
        Node*                p_node_2          = v_p_preceder_node[(int)v_p_preceder_node.size()-2];
        Edge*                p_follower        = 0;
        double               d_min_angle       = DBL_MAX;
        bool                 b_reverse_flag    = false;
        /////////////////////////////////////////
        if (p_preceder->get_follower_size() == 1) p_follower = p_preceder->get_follower(0);
        else
        {
            for (int i=0; i<p_preceder->get_follower_size(); ++i)
            {
                Edge*                p_edge        = p_preceder->get_follower(i);
                const vector<Node*>& v_p_edge_node = p_edge->get_node();
                Node*                p_node_3      = v_p_edge_node[0];
                Node*                p_node_4      = v_p_edge_node[1];
                double               d_angle       = positive_angle_2d(p_node_1->get_uv(), p_node_4->get_uv(), p_node_2->get_uv());
                if (d_angle < d_min_angle)
                {
                    p_follower  = p_edge;
                    d_min_angle = d_angle;
                }
            }
        }
        /////////////////////////////////////////
        if (!p_follower || p_follower == v_p_wire_edge.front()) break;
        /////////////////////////////////////////
        v_p_wire_edge.push_back(p_follower);
    }
    /////////////////////////////////////////////
#ifdef _VERBOSE
    {
        Edge*                p_edge        = v_p_wire_edge[0];
        const vector<Node*>& v_p_edge_node = p_edge->get_node();
        const double*        a_uv_1        = v_p_edge_node[0]->get_uv();
        const double*        a_uv_2        = v_p_edge_node[1]->get_uv();
        double               a_vector_12[2];
        transform(a_uv_2, a_uv_2+2, a_uv_1, a_vector_12, minus<double>());
        double               d_norm        = norm_2d(a_vector_12);
        transform(a_vector_12, a_vector_12+2, a_vector_12, std::bind(divides<double>(), _1, d_norm));
        double               a_uv[2];
        mean_2<double>(a_uv_1, a_uv_2, a_uv, 2);
        a_uv[0]                           -= a_vector_12[1] * _d_tolerance;
        a_uv[1]                           += a_vector_12[0] * _d_tolerance;
        /////////////////////////////////////////
        double d_angle = 0.0;
        for (int i=0; i<(int)v_p_wire_edge.size(); ++i)
        {
            const vector<Node*>& v_p_edge_node = v_p_wire_edge[i]->get_node();
            for (int j=0; j<(int)v_p_edge_node.size()-1; ++j)
            {				
                Node* p_node_1 = v_p_edge_node[j  ];
                Node* p_node_2 = v_p_edge_node[j+1];
                d_angle       += angle_2d(a_uv, p_node_1->get_uv(), p_node_2->get_uv());
            }
        }
        if (fabs(d_angle) < 350.0)
        {
            throw exception(_T("Failed to construct planar area due to internal error!"));
        }
    }
#endif
    /////////////////////////////////////////////
    vector<Node*> v_p_wire_node;
    vector<Edge*>::iterator itr_next=v_p_wire_edge.begin();
    for (vector<Edge*>::iterator itr_eg=v_p_wire_edge.begin(); itr_eg!=v_p_wire_edge.end(); itr_eg=itr_next) // mylee
    {
        itr_next = itr_eg; ++itr_next;

        Edge*                p_edge        = *itr_eg;
        const vector<Node*>& v_p_edge_node = p_edge->get_node();
        copy(v_p_edge_node.begin(), v_p_edge_node.end(), back_inserter(v_p_wire_node));
        /////////////////////////////////////////
        for (int i=0; i<p_edge->get_preceder_size(); ++i) p_edge->get_preceder(i)->remove_follower(p_edge);
        for (int i=0; i<p_edge->get_follower_size(); ++i) p_edge->get_follower(i)->remove_preceder(p_edge);
        /////////////////////////////////////////
        _v_p_edge.erase(find(_v_p_edge.begin(), _v_p_edge.end(), p_edge));
        delete p_edge;
        //v_p_wire_edge.erase(itr_eg--);
        itr_next = v_p_wire_edge.erase(itr_eg); // mylee
    }
    v_p_wire_node.erase(unique(v_p_wire_node.begin(), v_p_wire_node.end()), v_p_wire_node.end());
    /////////////////////////////////////////////
    Wire* p_wire = new Wire();
    p_wire->set_type(Wire::WT_OUTER);
    p_wire->set_node(v_p_wire_node);
    _v_p_wire.push_back(p_wire);
    /////////////////////////////////////////////
    Domain* p_domain = new Domain();
    p_domain->add_wire(p_wire);
    _v_p_domain.push_back(p_domain);
} // end: Domain_Maker::make_multi_edge_domain()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
inline Edge* Domain_Maker::pick_wire_head_edge() const
{
    for (vector<Edge*>::const_iterator itr_eg=_v_p_edge.begin(); itr_eg!=_v_p_edge.end(); ++itr_eg)
    {
        Edge* p_edge = *itr_eg;
        if (!p_edge->is_open_edge()) return p_edge;
    }
    return 0;
} // end: Domain_Maker::pick_wire_head_edge()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Domain_Maker::set_inner_domain(bool b_keep_inner_domain_flag) throw(exception)
{
    using namespace mesher_math_lib;
    /////////////////////////////////////////////
    map<Domain*, UV> m_domain_uv;
    /////////////////////////////////////////////
    vector<Domain*>::const_iterator itr_dm;
    for (itr_dm=_v_p_domain.begin(); itr_dm!=_v_p_domain.end(); ++itr_dm)
    {
        Domain* p_domain = *itr_dm;
        ///////////////////////////////////////////
        double               a_uv[2];
        double               d_min_deviation = 360.0;
        const vector<Node*>& v_p_wire_node   = (p_domain->get_wire(0))->get_node();
        for (int i=0; i<(int)v_p_wire_node.size()-1; ++i)
        {
            Node*  p_node_0 = i == 0                      ? v_p_wire_node[v_p_wire_node.size()-2] : v_p_wire_node[i-1];
            Node*  p_node_1 = v_p_wire_node[i];
            Node*  p_node_2 = i == v_p_wire_node.size()-2 ? v_p_wire_node[1]                      : v_p_wire_node[i+1];
            double d_angle  = angle_2d(p_node_1->get_uv(), p_node_2->get_uv(), p_node_0->get_uv());
            if (d_angle > 0.0 && d_angle < 179)
            {
                double d_deviation = fabs(d_angle - 90.0);
                if (d_deviation < d_min_deviation)
                {
                    d_min_deviation = d_deviation;
                    for (int j=0; j<2; ++j) 
                    {
                        a_uv[j] = (p_node_0->get_uv()[j] + p_node_1->get_uv()[j] + p_node_2->get_uv()[j]) / 3.0;
                        a_uv[j] =  p_node_1->get_uv()[j] + (a_uv[j] - p_node_1->get_uv()[j]) * 0.05;
                    }
                }
            }
        }
        ///////////////////////////////////////////
        m_domain_uv[p_domain] = UV(a_uv);
    }
    /////////////////////////////////////////////
    map< Domain*, vector<Domain*> > m_sub_domain;
    map< Domain*, int             > m_mother_count;
    /////////////////////////////////////////////
    // [2009-09-29] Kim, Geun Young (Tel: 2042, gykim@midasit.com)
    // inner domain이 2차 이상일 경우 1차까지만 고려하기 위해 수정
    bool bInnerDomain = true;
    for (int i=0; i<(int)_v_p_domain.size(); ++i) m_mother_count[_v_p_domain[i]] = 0;
    for (int i=0; i<(int)_v_p_domain.size()-1; ++i)
    {
        Domain* p_domain_1 = _v_p_domain[i];
        for (int j=i+1; j<(int)_v_p_domain.size(); ++j)
        {
            if (check_interrupt()) throw interrupt_exception();
            /////////////////////////////////////
            Domain* p_domain_2 = _v_p_domain[j];
            /////////////////////////////////////
            if (p_domain_1->is_internal_point(m_domain_uv[p_domain_2].a_uv))
            {
                if(m_mother_count[p_domain_2] < 1)
                {
                    m_sub_domain[p_domain_1].push_back(p_domain_2);
                    ++m_mother_count[p_domain_2];
                }
                else
                {
                    if(b_keep_inner_domain_flag)          
                    {
                        bInnerDomain = false;
                    }      
                }
            }
            if (p_domain_2->is_internal_point(m_domain_uv[p_domain_1].a_uv))
            {
                if(m_mother_count[p_domain_1] < 1)
                {
                    m_sub_domain[p_domain_2].push_back(p_domain_1);
                    ++m_mother_count[p_domain_1];
                }   
                else
                {
                    if(b_keep_inner_domain_flag)
                    {
                        bInnerDomain = false;
                    } 
                }
            }
        }
    }
    /////////////////////////////////////////////
    /*
    for (i=0; i<(int)_v_p_domain.size(); ++i) 
    {
    if (m_mother_count[_v_p_domain[i]] > 1)
    throw exception(_T("[Warning] Duplicate inner domains exist.\nPlease generate mesh domain by inner domain!");
    }
    */
    if(!bInnerDomain) 
        AfxMessageBox(_T("Warning] Duplicate inner domains exist.\nPlease generate mesh domain by inner domain!"));  
    /////////////////////////////////////////////
    vector<Domain*> v_p_inner_domain;
    /////////////////////////////////////////////
    for (vector<Domain*>::const_iterator itr_dm_1=_v_p_domain.begin(); itr_dm_1!=_v_p_domain.end(); ++itr_dm_1)
    {
        Domain*          p_domain_1     = *itr_dm_1;
        vector<Domain*>& v_p_sub_domain = m_sub_domain[p_domain_1];
        /////////////////////////////////////////
        if (v_p_sub_domain.empty()) continue;
        /////////////////////////////////////////
        vector<Domain*>::iterator itr_next;
        vector<Domain*>::iterator itr_dm_2;
        for (itr_dm_2=v_p_sub_domain.begin(); itr_dm_2!=v_p_sub_domain.end(); itr_dm_2=itr_next) // mylee
        {
            itr_next = itr_dm_2; ++itr_next;

            if (check_interrupt()) throw interrupt_exception();
            ///////////////////////////////////////
            Domain* p_domain_2 = *itr_dm_2;
            for (vector<Domain*>::const_iterator itr_dm_3=_v_p_domain.begin(); itr_dm_3!=_v_p_domain.end(); ++itr_dm_3)
            {
                //Domain* p_domain_3 = v_p_sub_domain[i]; // mylee
                if(itr_dm_1 == itr_dm_3) continue; // mylee bug
                Domain* p_domain_3 = *itr_dm_3; // mylee bug
                if (find(m_sub_domain[p_domain_3].begin(), m_sub_domain[p_domain_3].end(), p_domain_2) != m_sub_domain[p_domain_3].end())
                {
                    //v_p_sub_domain.erase(itr_dm_2--);
                    itr_next = v_p_sub_domain.erase(itr_dm_2); // mylee
                    break;
                }
            }      
        }
        for (itr_dm_2=v_p_sub_domain.begin(); itr_dm_2!=v_p_sub_domain.end(); ++itr_dm_2)
        {
            Domain*       p_domain_2    = *itr_dm_2;
            Wire*         p_wire_1      = p_domain_2->get_wire(0);
            vector<Node*> v_p_wire_node = p_wire_1->get_node();
            reverse(v_p_wire_node.begin(), v_p_wire_node.end());
            /////////////////////////////////////
            Wire* p_wire_2 = new Wire();
            p_wire_2->set_type(Wire::WT_INNER);
            p_wire_2->set_node(v_p_wire_node);
            _v_p_wire.push_back(p_wire_2);
            p_domain_1->add_wire(p_wire_2);
            /////////////////////////////////////
            v_p_inner_domain.push_back(p_domain_2);
        }
        /////////////////////////////////////////
        check_inner_wire_connection(p_domain_1);
        //	check_matryoshka_domain(p_domain_1);
    }
    if (b_keep_inner_domain_flag && bInnerDomain) return;
    /////////////////////////////////////////////
    for (itr_dm=v_p_inner_domain.begin(); itr_dm!=v_p_inner_domain.end(); ++itr_dm)	
    {
        Domain* p_domain = *itr_dm;
        /////////////////////////////////////////
        for (int i=0; i<p_domain->get_wire_size(); ++i)
        {
            Wire* p_wire = p_domain->get_wire(i);
            _v_p_wire.erase(find(_v_p_wire.begin(), _v_p_wire.end(), p_wire));
            delete p_wire;
        }
        /////////////////////////////////////////
        _v_p_domain.erase(find(_v_p_domain.begin(), _v_p_domain.end(), p_domain));
        delete p_domain;
    }
} // end: Domain_Maker::set_inner_domain()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Domain_Maker::check_inner_wire_connection(Domain* p_domain) throw(exception)
{
    using namespace mesher_math_lib;
    /////////////////////////////////////////////
    vector<Wire*> v_p_inner_wire;
    for (int i=0; i<p_domain->get_wire_size(); ++i)
    {
        Wire* p_wire = p_domain->get_wire(i);
        if (p_wire->get_type() == Wire::WT_INNER) v_p_inner_wire.push_back(p_wire);
    }
    if (v_p_inner_wire.size() == 1) return;
    /////////////////////////////////////////////
    map<int, bool> m_obsolete_wire_flag;
    for (int i=0; i<(int)v_p_inner_wire.size(); ++i) m_obsolete_wire_flag[i] = false;
    for (int i=0; i<(int)v_p_inner_wire.size(); ++i)
    {
        if (m_obsolete_wire_flag[i]) continue;
        ///////////////////////////////////////////
        Wire* p_wire_1 = v_p_inner_wire[i];
        for (int j=0; j<(int)v_p_inner_wire.size(); ++j)
        {
            if (i == j || m_obsolete_wire_flag[j]) continue;
            /////////////////////////////////////////
            Wire*                p_wire_2              = v_p_inner_wire[j];
            const vector<Node*>& v_p_wire_node_1       = p_wire_1->get_node();
            const vector<Node*>& v_p_wire_node_2       = p_wire_2->get_node();
            const int            n_wire_segment_size_1 = (int)v_p_wire_node_1.size()-1;
            const int            n_wire_segment_size_2 = (int)v_p_wire_node_2.size()-1;
            /////////////////////////////////////////
            vector<Segment*>    v_p_segment;
            map<Segment*, bool> m_obsolete_segment_flag;
            for (int k=0; k<n_wire_segment_size_1; ++k)
            {
                Segment* p_segment                 = new Segment(v_p_wire_node_1[k], v_p_wire_node_1[k+1]);
                m_obsolete_segment_flag[p_segment] = false;
                v_p_segment.push_back(p_segment);
            }
            for (int k=0; k<n_wire_segment_size_2; ++k)
            {
                Segment* p_segment                 = new Segment(v_p_wire_node_2[k], v_p_wire_node_2[k+1]);
                m_obsolete_segment_flag[p_segment] = false;
                v_p_segment.push_back(p_segment);
            }
            /////////////////////////////////////////
            for (vector<Segment*>::iterator itr_sg_1=v_p_segment.begin(); itr_sg_1!=v_p_segment.end()-1; ++itr_sg_1)
            {
                Segment* p_segment_1 = *itr_sg_1;
                /////////////////////////////////////
                if (m_obsolete_segment_flag[p_segment_1]) continue;
                /////////////////////////////////////
                for (vector<Segment*>::iterator itr_sg_2=itr_sg_1+1; itr_sg_2!=v_p_segment.end(); ++itr_sg_2)
                {
                    Segment* p_segment_2 = *itr_sg_2;
                    if (p_segment_1->get_node_1() == p_segment_2->get_node_2() && p_segment_1->get_node_2() == p_segment_2->get_node_1())
                    {
                        m_obsolete_segment_flag[p_segment_1] = true;
                        m_obsolete_segment_flag[p_segment_2] = true;
                    }
                }
            }
            vector<Segment*>::iterator itr_next;
            vector<Segment*>::iterator itr_sg;
            for (itr_sg=v_p_segment.begin(); itr_sg!=v_p_segment.end(); itr_sg=itr_next) // mylee
            {
                itr_next = itr_sg; ++itr_next;

                Segment* p_segment = *itr_sg;
                if (m_obsolete_segment_flag[p_segment])
                {
                    delete p_segment;
                    //v_p_segment.erase(itr_sg--);
                    itr_next = v_p_segment.erase(itr_sg);
                }
            }
            /////////////////////////////////////////
            if ((int)v_p_segment.size() == n_wire_segment_size_1+n_wire_segment_size_2) 
            {
                // [2009-07-30] Kim, Geun Young (Tel: 2042, gykim@midasit.com)
                // memory leak 수정함.
                for (itr_sg=v_p_segment.begin(); itr_sg!=v_p_segment.end(); ++itr_sg) delete *itr_sg;
                continue;
            }
            /////////////////////////////////////////
            for (int k=0; k<(int)v_p_segment.size()-1; ++k)
            {
                Segment* p_segment_1 = v_p_segment[k];
                Node*    p_node_1    = p_segment_1->get_node_1();
                Node*    p_node_2    = p_segment_1->get_node_2();
                for (int l=k+1; l<(int)v_p_segment.size(); ++l)
                {
                    Segment* p_segment_2 = v_p_segment[l];
                    Node*    p_node_3    = p_segment_2->get_node_1();
                    Node*    p_node_4    = p_segment_2->get_node_2();
                    /////////////////////////////////
                    if (p_node_3 == p_node_2) p_segment_1->set_follower(p_segment_2);
                    if (p_node_4 == p_node_1) p_segment_2->set_follower(p_segment_1);
                }
            }
            /////////////////////////////////////////
            if (v_p_segment.empty())
            {
                throw exception("The configuration of inner domains are ambiguous or too complex to mesh simultaneously.\nPlease mesh the domains domain by domain!");
            }
            /////////////////////////////////////////
            vector<Segment*> v_p_chain_segment;
            v_p_chain_segment.push_back(v_p_segment.front());
            while (true)
            {
                Segment* p_follower = (v_p_chain_segment.back())->get_follower();
                /////////////////////////////////////
                if (!p_follower || p_follower == v_p_chain_segment.front()) break;
                /////////////////////////////////////
                v_p_chain_segment.push_back(p_follower);
            }
            /////////////////////////////////////////
            vector<Node*> v_p_union_node;
            for (itr_sg=v_p_chain_segment.begin(); itr_sg!=v_p_chain_segment.end(); ++itr_sg) v_p_union_node.push_back((*itr_sg)->get_node_1());
            v_p_union_node.push_back(v_p_union_node.front());
            /////////////////////////////////////////
            for (itr_sg=v_p_segment.begin(); itr_sg!=v_p_segment.end(); ++itr_sg) delete *itr_sg;
            /////////////////////////////////////////
            p_wire_1->set_node(v_p_union_node);
            /////////////////////////////////////////
            m_obsolete_wire_flag[j] = true;
        }
    }
    /////////////////////////////////////////////
    for (int i=0; i<(int)v_p_inner_wire.size(); ++i)
    {
        if (m_obsolete_wire_flag[i]) p_domain->remove_wire(v_p_inner_wire[i]);
    }
} // end: Domain_Maker::check_inner_wire_connection()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Domain_Maker::check_matryoshka_domain(Domain* p_domain) throw(exception)
{
    using namespace mesher_math_lib;
    /////////////////////////////////////////////
    vector<Wire*> v_p_inner_wire;
    for (int i=0; i<p_domain->get_wire_size(); ++i)
    {
        Wire* p_wire = p_domain->get_wire(i);
        if (p_wire->get_type() == Wire::WT_INNER) v_p_inner_wire.push_back(p_wire);
    }
    if (v_p_inner_wire.size() == 1) return;
    /////////////////////////////////////////////
    map<Wire*, UV> m_wire_uv;
    for (vector<Wire*>::const_iterator itr_wr=v_p_inner_wire.begin(); itr_wr!=v_p_inner_wire.end(); ++itr_wr)
    {
        Wire* p_wire = *itr_wr;
        ///////////////////////////////////////////
        const vector<Node*>& v_p_wire_node   = p_wire->get_node();
        double               a_uv[2];
        double               d_min_deviation = 360.0;
        for (int i=0; i<(int)v_p_wire_node.size()-1; ++i)
        {
            Node*  p_node_0 = i == 0                      ? v_p_wire_node[v_p_wire_node.size()-2] : v_p_wire_node[i-1];
            Node*  p_node_1 = v_p_wire_node[i];
            Node*  p_node_2 = i == v_p_wire_node.size()-2 ? v_p_wire_node[1]                      : v_p_wire_node[i+1];
            double d_angle  = angle_2d(p_node_1->get_uv(), p_node_0->get_uv(), p_node_2->get_uv());
            if (d_angle > 0.0 && d_angle < 179)
            {
                double d_deviation = fabs(d_angle - 90.0);
                if (d_deviation < d_min_deviation)
                {
                    d_min_deviation = d_deviation;
                    for (int j=0; j<2; ++j) a_uv[j] = (p_node_0->get_uv()[j] + p_node_1->get_uv()[j] + p_node_2->get_uv()[j]) / 3.0;
                }
            }
        }
        ///////////////////////////////////////////
        m_wire_uv[p_wire] = UV(a_uv);
    }
    /////////////////////////////////////////////
    for (vector<Wire*>::const_iterator itr_wr_1=v_p_inner_wire.begin(); itr_wr_1!=v_p_inner_wire.end(); ++itr_wr_1)
    {
        Wire*         p_wire_1 = *itr_wr_1;
        const double* a_uv     = m_wire_uv[p_wire_1].a_uv;
        for (vector<Wire*>::const_iterator itr_wr_2=v_p_inner_wire.begin(); itr_wr_2!=v_p_inner_wire.end(); ++itr_wr_2)
        {
            Wire* p_wire_2 = *itr_wr_2;
            /////////////////////////////////////////
            if (p_wire_1 == p_wire_2) continue;
            /////////////////////////////////////////
            vector<Node*> v_p_wire_node = p_wire_2->get_node();
            double        d_angle       = 0.0;
            for (int i=0; i<(int)v_p_wire_node.size()-1; ++i)
            {				
                Node* p_node_1 = v_p_wire_node[i  ];
                Node* p_node_2 = v_p_wire_node[i+1];
                d_angle       += angle_2d(a_uv, p_node_1->get_uv(), p_node_2->get_uv());
            }
            if (fabs(d_angle) > 350.0)
            {
                throw exception("The configuration of inner domains are ambiguous or too complex to mesh simultaneously.\nPlease mesh the domains domain by domain!");
            }
        }
    }
} // end: Domain_Maker::check_matryoshka_domain()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Domain_Maker::insert_internal_wire_node(vector<int>& v_unused_node_index) throw(exception)
{
    using namespace mesher_math_lib;
    /////////////////////////////////////////////
    v_unused_node_index.clear();
    /////////////////////////////////////////////
    vector<Edge*>::iterator itr_next;
    for (vector<Edge*>::iterator itr_eg=_v_p_edge.begin(); itr_eg!=_v_p_edge.end(); itr_eg=itr_next) // mylee
    {
        itr_next = itr_eg; ++itr_next;

        Edge*                p_edge        = *itr_eg;
        const vector<Node*>& v_p_edge_node = p_edge->get_node();
        const double*        a_uv_1        = v_p_edge_node[0]->get_uv();
        const double*        a_uv_2        = v_p_edge_node[1]->get_uv();
        double               a_uv[2];
        mean_2<double>(a_uv_1, a_uv_2, a_uv, 2);
        for (vector<Domain*>::const_iterator itr_dm=_v_p_domain.begin(); itr_dm!=_v_p_domain.end(); ++itr_dm)
        {
            if (check_interrupt()) throw interrupt_exception();
            /////////////////////////////////////
            Domain* p_domain = *itr_dm;
            if (p_domain->is_internal_point(a_uv))
            {
                Wire* p_wire = new Wire();
                p_wire->set_type(Wire::WT_OPEN);
                p_wire->set_node(v_p_edge_node);
                _v_p_wire.push_back(p_wire);
                p_domain->add_wire(p_wire);
                /////////////////////////////////
                delete p_edge;
                //_v_p_edge.erase(itr_eg--);
                itr_next = _v_p_edge.erase(itr_eg); // mylee
                /////////////////////////////////
                break;
            }
        }
    }
    /////////////////////////////////////////////
    for (int i=0; i<(int)_v_p_island_node.size(); ++i)
    {
        Node* p_node = _v_p_island_node[i];
        /////////////////////////////////////////
        bool b_used_node_flag = false;
        for (vector<Domain*>::const_iterator itr_dm=_v_p_domain.begin(); itr_dm!=_v_p_domain.end(); ++itr_dm)
        {
            if (check_interrupt()) throw interrupt_exception();
            ///////////////////////////////////////
            Domain* p_domain = *itr_dm;
            if (p_domain->is_internal_point(p_node->get_uv()))
            {
                p_domain->add_island_node(p_node);
                b_used_node_flag = true;
                break;
            }
        }
        if (!b_used_node_flag) v_unused_node_index.push_back(_m_node_index[p_node]);
    }
} // end: Domain_Maker::insert_internal_wire_node()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Domain_Maker::finalize_domain_making(Mesh_Packet_2D& mesh_packet, vector<int>& v_unused_node_chain_index)
{
    for (int i=0; i<(int)mesh_packet.v_node_chain.size(); ++i) mesh_packet.v_node_chain[i].v_tag.clear();
    /////////////////////////////////////////////
    map< Wire*, vector<int> > m_wire_tag;
    map< Edge*, vector<int> > m_edge_tag;
    for (vector<Wire*>::const_iterator itr_wr=_v_p_wire.begin(); itr_wr!=_v_p_wire.end(); ++itr_wr)
    {
        Wire*                p_wire        = *itr_wr;
        const vector<Node*>& v_p_wire_node = p_wire->get_node();
        for (int i=0; i<(int)mesh_packet.v_node_chain.size(); ++i)
        {
            const vector<Node*>& v_p_chain_node = _m_node_chain_node[i];
            bool                 b_subset_flag  = true;
            for (vector<Node*>::const_iterator itr_nd=v_p_chain_node.begin(); itr_nd!=v_p_chain_node.end(); ++itr_nd)
            {
                Node* p_node = *itr_nd;
                if (find(v_p_wire_node.begin(), v_p_wire_node.end(), p_node) == v_p_wire_node.end())
                {
                    b_subset_flag = false;
                    break;
                }
            }
            if (b_subset_flag) m_wire_tag[p_wire].push_back(mesh_packet.v_node_chain[i].n_tag);
        }
    }
    vector<Edge*>::const_iterator itr_eg;
    for (itr_eg=_v_p_edge.begin(); itr_eg!=_v_p_edge.end(); ++itr_eg)
    {
        Edge*                p_edge        = *itr_eg;
        const vector<Node*>& v_p_edge_node = p_edge->get_node();
        for (int i=0; i<(int)mesh_packet.v_node_chain.size(); ++i)
        {
            const vector<Node*>& v_p_chain_node = _m_node_chain_node[i];
            bool                 b_subset_flag  = true;
            for (vector<Node*>::const_iterator itr_nd=v_p_chain_node.begin(); itr_nd!=v_p_chain_node.end(); ++itr_nd)
            {
                Node* p_node = *itr_nd;
                if (find(v_p_edge_node.begin(), v_p_edge_node.end(), p_node) == v_p_edge_node.end())
                {
                    b_subset_flag = false;
                    break;
                }
            }
            if (b_subset_flag) m_edge_tag[p_edge].push_back(mesh_packet.v_node_chain[i].n_tag);
        }
    }
    /////////////////////////////////////////////
    mesh_packet.v_node_chain.clear();
    mesh_packet.v_domain    .clear();
    /////////////////////////////////////////////
    map<Wire*, int> m_wire_index;
    for (int i=0; i<(int)_v_p_wire.size(); ++i)
    {
        Wire* p_wire         = _v_p_wire[i];
        m_wire_index[p_wire] = i;
        /////////////////////////////////////////
        const vector<Node*>&         v_p_wire_node = p_wire->get_node();
        Mesh_Packet_2D::I_Node_Chain i_node_chain;
        i_node_chain.t_type                        = Mesh_Packet_2D::I_Node_Chain::CHAIN_TYPE((int)(p_wire->get_type()));
        for (int j=0; j<(int)v_p_wire_node.size(); ++j) i_node_chain.v_node_index.push_back(_m_node_index[v_p_wire_node[j]]);
        i_node_chain.v_tag                         = m_wire_tag[p_wire];
        mesh_packet.v_node_chain.push_back(i_node_chain);
    }
    /////////////////////////////////////////////
    for (vector<Domain*>::const_iterator itr_dm=_v_p_domain.begin(); itr_dm!=_v_p_domain.end(); ++itr_dm)
    {
        Domain* p_domain = *itr_dm;
        Mesh_Packet_2D::I_Domain i_domain;
        for (int i=0; i<p_domain->get_wire_size();        ++i) i_domain.v_node_chain_index .push_back( m_wire_index[p_domain->get_wire(i)       ]);
        for (int i=0; i<p_domain->get_island_node_size(); ++i) i_domain.v_island_node_index.push_back(_m_node_index[p_domain->get_island_node(i)]);
        mesh_packet.v_domain.push_back(i_domain);
    }
    /////////////////////////////////////////////
    int n_index = (int)mesh_packet.v_node_chain.size();
    for (itr_eg=_v_p_edge.begin(); itr_eg!=_v_p_edge.end(); ++itr_eg)
    {
        Edge*                p_edge        = *itr_eg;
        const vector<Node*>& v_p_edge_node = p_edge->get_node();
        /////////////////////////////////////////
        Mesh_Packet_2D::I_Node_Chain i_node_chain;
        i_node_chain.n_tag = p_edge->get_tag();
        for (int i=0; i<(int)v_p_edge_node.size(); ++i) i_node_chain.v_node_index.push_back(_m_node_index[v_p_edge_node[i]]);
        i_node_chain.v_tag = m_edge_tag[p_edge];
        mesh_packet.v_node_chain.push_back(i_node_chain);
        /////////////////////////////////////////
        v_unused_node_chain_index.push_back(n_index++);
    }
} // end: Domain_Maker::finalize_domain_making()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Domain_Maker::dump_input_file(tstring s_file_name, const Mesh_Packet_2D& mesh_packet)
{
    ofstream f_dump_file(s_file_name.c_str());
    if (!f_dump_file) { cerr << _T("* cannot open dump file - dump_input_file()") << endl; return; }
    /////////////////////////////////////////
    f_dump_file << mesh_packet.d_mesh_size << endl;
    f_dump_file << endl;
    f_dump_file << (int)mesh_packet.v_node.size() << endl;
    for (int j=0; j<(int)mesh_packet.v_node.size(); ++j) f_dump_file << mesh_packet.v_node[j].a_uv[0] << _T(" ") << mesh_packet.v_node[j].a_uv[1] << _T(" ") << mesh_packet.v_node[j].b_hard_flag << " " << mesh_packet.v_node[j].b_island_flag << " " << mesh_packet.v_node[j].d_mesh_size << endl;
    f_dump_file << endl;
    f_dump_file << (int)mesh_packet.v_node_chain.size() << endl;
    for (int j=0; j<(int)mesh_packet.v_node_chain.size(); ++j)
    {
        f_dump_file << _T("2 ") << (int)mesh_packet.v_node_chain[j].v_node_index.size() << endl;
        for (int k=0; k<(int)mesh_packet.v_node_chain[j].v_node_index.size(); ++k) f_dump_file << mesh_packet.v_node_chain[j].v_node_index[k] << _T(" ");
        f_dump_file << endl;
    }
    /////////////////////////////////////////
    f_dump_file.close();
} // end: Domain_Maker::dump_input_file()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Domain_Maker::write_domain_file(tstring s_file_name, const Mesh_Packet_2D& mesh_packet)
{
    if (s_file_name.empty()) return;
    /////////////////////////////////////////////
    TCHAR a_extension[1024];
    _tsplitpath(s_file_name.c_str(), 0, 0, 0, a_extension);
    int  n_size = (int)s_file_name.size() - (int)_tcsclen(a_extension);
    s_file_name = s_file_name.substr(0, n_size);
    /////////////////////////////////////////////
    for (int i=0; i<(int)mesh_packet.v_domain.size(); ++i)
    {
        Mesh_Packet_2D::I_Domain i_domain = mesh_packet.v_domain[i];
        /////////////////////////////////////////
        vector<Mesh_Packet_2D::I_Node_Chain> v_node_chain;
        for (int j=0; j<(int)i_domain.v_node_chain_index.size();  ++j) v_node_chain.push_back(mesh_packet.v_node_chain[i_domain.v_node_chain_index[j]]);
        /////////////////////////////////////////
        set<int>                       s_node_index;
        vector<Mesh_Packet_2D::I_Node> v_node;
        map<int, int>                  m_node_index;
        int                            n_index = 0;
        for (int j=0; j<(int)v_node_chain.size(); ++j)
        {
            const Mesh_Packet_2D::I_Node_Chain& i_node_chain = v_node_chain[j];
            for (int k=0; k<(int)i_node_chain.v_node_index.size(); ++k) s_node_index.insert(i_node_chain.v_node_index[k]);
        }
        for (int j=0; j<(int)i_domain.v_island_node_index.size(); ++j) s_node_index.insert(i_domain.v_island_node_index[j]);
        for (set<int>::const_iterator itr_ni=s_node_index.begin(); itr_ni!=s_node_index.end(); ++itr_ni)
        {
            v_node.push_back(mesh_packet.v_node[*itr_ni]);
            m_node_index[*itr_ni] = n_index++;
        }
        /////////////////////////////////////////
        TCHAR a_file_name[1024];
        _stprintf_s(a_file_name, _T("%s_%d.dom"), s_file_name.c_str(), i+1);
        ofstream f_domain_file(a_file_name);
        if (!f_domain_file) { cerr << _T("* cannot open output file - write_domain_file()") << endl; return; }
        /////////////////////////////////////////
        f_domain_file << mesh_packet.d_mesh_size << endl;
        f_domain_file << endl;
        f_domain_file << (int)v_node.size() << endl;
        for (int j=0; j<(int)v_node.size(); ++j) f_domain_file << v_node[j].a_uv[0] << _T(" ") << v_node[j].a_uv[1] << _T(" ") << v_node[j].b_hard_flag << _T(" ") << v_node[j].b_island_flag << _T(" ") << v_node[j].d_mesh_size << endl;
        f_domain_file << endl;
        f_domain_file << (int)v_node_chain.size() << endl;
        for (int j=0; j<(int)v_node_chain.size(); ++j)
        {
            f_domain_file << (int)v_node_chain[j].t_type << _T(" ") << (int)v_node_chain[j].v_node_index.size() << endl;
            for (int k=0; k<(int)v_node_chain[j].v_node_index.size(); ++k) f_domain_file << m_node_index[v_node_chain[j].v_node_index[k]] << _T(" ");
            f_domain_file << endl;
        }
        /////////////////////////////////////////
        f_domain_file.close();
    }
} // end: Domain_Maker::write_domain_file()
/////////////////////////////////////////////////////////////////////
