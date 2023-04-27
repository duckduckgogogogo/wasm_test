// mesher_2d.cpp: 2D auto-mesh generator

//////////////////////////////////////////////////////////////////////////////////
// Shin Dae-Seock (dsshin@midasIT.com)
// Don't believe in magic;
// understand what your libraries do, how they do it, and at what cost they do it.
//////////////////////////////////////////////////////////////////////////////////

#include "StdAfx.h"

#include "mesher_2d.h"
#include "Math_Lib.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <functional>
#include <numeric>
#include <iterator>
#include <iomanip>

/////////////////////////////////////////////////////////////////////
using namespace mesher_2d;
using namespace std::placeholders;
/////////////////////////////////////////////////////////////////////

#pragma warning ( disable : 4018 )
#pragma warning ( disable : 4267 )
#pragma warning ( disable : 4996 )

/////////////////////////////////////////////////////////////////////
#ifdef _VERBOSE
ofstream Mesher_2D::_f_debug;
#endif
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
Mesher_2D::Mesher_2D()
{
    _t_mesh_type        = MT_QUAD;
    _d_mesh_size        = 0.0;
    _d_refinement_ratio = 1.0;
    _d_tolerance        = 0.0;
    /////////////////////////////////////////////
    _v_p_node   .reserve(10000);
    _v_p_element.reserve(10000);
    /////////////////////////////////////////////
    _p_interrupt_handler = 0;
    /////////////////////////////////////////////
#ifdef _VERBOSE
    {
        _f_debug.open(_T("c:/temp/mesher.dmp"));
    }
#endif
} // end: class Mesher_2D - constructor
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
Mesher_2D::~Mesher_2D()
{
    clean_up();
    /////////////////////////////////////////////
#ifdef _VERBOSE
    {
        _f_debug.close();
    }
#endif
} // end: class Mesher_2D - destructor
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Mesher_2D::clean_up()
{
    for (vector<Node*>   ::const_iterator itr_nd=_v_p_node   .begin(); itr_nd!=_v_p_node   .end(); ++itr_nd) delete *itr_nd;
    for (vector<Element*>::const_iterator itr_el=_v_p_element.begin(); itr_el!=_v_p_element.end(); ++itr_el) delete *itr_el;
    for (vector<Wire*>   ::const_iterator itr_wr=_v_p_wire   .begin(); itr_wr!=_v_p_wire   .end(); ++itr_wr) delete *itr_wr;
    /////////////////////////////////////////////
    _v_p_node         .clear();
    _v_p_element      .clear();
    _v_p_wire         .clear();
    _v_p_island_node  .clear();
    _m_input_wire_node.clear();
} // end: Mesher_2D::clean_up()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Mesher_2D::set_input_data(const Mesh_Packet_2D& mesh_packet) throw(exception)
{
    using namespace mesher_math_lib;
    /////////////////////////////////////////////
    bool b_scale_factor = fabs(mesh_packet.d_scale_factor-1.0) > 1.0e-6;
    /////////////////////////////////////////////
    _t_mesh_type        = MESH_TYPE(mesh_packet.t_mesh_type);
    _d_mesh_size        = mesh_packet.d_mesh_size;
    _d_refinement_ratio = mesh_packet.d_refinement_ratio > ZERO ? mesh_packet.d_refinement_ratio : 1.0;
    /////////////////////////////////////////////
    if (b_scale_factor) _d_mesh_size *= mesh_packet.d_scale_factor;
    /////////////////////////////////////////////
    map<int, Node*> m_p_node;
    double          a_uv[2];
    for (int i=0; i<(int)mesh_packet.v_node.size(); ++i)
    {
        const Mesh_Packet_2D::I_Node& i_node = mesh_packet.v_node[i];
        copy(i_node.a_uv, i_node.a_uv+2, a_uv);
        if (b_scale_factor) transform(a_uv, a_uv+2, a_uv, std::bind(multiplies<double>(), _1, mesh_packet.d_scale_factor));
        Node* p_node = new Node(a_uv, i+1, i_node.b_hard_flag, i_node.d_mesh_size);
        _v_p_node.push_back(p_node);
        if (i_node.b_island_flag) _v_p_island_node.push_back(p_node);
        m_p_node[i]  = p_node;
    }
    for (vector<Mesh_Packet_2D::I_Node_Chain>::const_iterator itr_ch=mesh_packet.v_node_chain.begin(); itr_ch!=mesh_packet.v_node_chain.end(); ++itr_ch)
    {
        const Mesh_Packet_2D::I_Node_Chain& i_node_chain = *itr_ch;
        Wire*                               p_wire       = new Wire();
        vector<Node*>                       v_p_wire_node;
        for (int i=0; i<(int)i_node_chain.v_node_index.size(); ++i)
        {
            map<int, Node*>::iterator itr_nd = m_p_node.find(i_node_chain.v_node_index[i]);

            //if (itr_nd == m_p_node.end()) throw exception(CDgnMSGFunction::ConvertW2A(_LS(IDS_BASE_INVALID_NODE_INDEX)));
            Node* p_node = itr_nd->second;
            p_node->set_hard_node();
            v_p_wire_node.push_back(p_node);
        }
        p_wire->set_type(Wire::WIRE_TYPE((int)(i_node_chain.t_type)));
        p_wire->set_node(v_p_wire_node);
        _v_p_wire.push_back(p_wire);
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
} // end: Mesher_2D::set_input_data()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Mesher_2D::check_input_data() throw(exception)
{
    using namespace mesher_math_lib;
    /////////////////////////////////////////////

    //if (_v_p_node.empty()) throw exception(CDgnMSGFunction::ConvertW2A(_LS(IDS_BASE_NO_INPUT_NODE)));
    //if (_v_p_wire.empty()) throw exception(CDgnMSGFunction::ConvertW2A(_LS(IDS_BASE_NO_INPUT_WIRE)));
    /////////////////////////////////////////////
    int n_outer_wire_count = 0;
    vector<Wire*>::const_iterator itr_wr;
    for (itr_wr=_v_p_wire.begin(); itr_wr!=_v_p_wire.end(); ++itr_wr)
    {
        if ((*itr_wr)->get_type() == Wire::WT_OUTER) ++n_outer_wire_count;
    }
    //if (n_outer_wire_count == 0) throw exception(CDgnMSGFunction::ConvertW2A(_LS(IDS_BASE_NO_OUTER_WIRE)));
    //if (n_outer_wire_count >  1) throw exception(CDgnMSGFunction::ConvertW2A(_LS(IDS_BASE_MULTI_OUTER_WIRE)));
    /////////////////////////////////////////////
    _d_tolerance = DBL_MAX;
    for (int i=0; i<(int)_v_p_node.size()-1; ++i)
    {
        Node* p_node_1 = _v_p_node[i];
        for (int j=i+1; j<(int)_v_p_node.size(); ++j)
        {
            Node*  p_node_2   = _v_p_node[j];
            double d_distance = distance_2d(p_node_1->get_uv(), p_node_2->get_uv());
            /////////////////////////////////////
            //if (d_distance < ZERO) throw exception(CDgnMSGFunction::ConvertW2A(_LS(IDS_BASE_OVERLAP_NODE)));
            /////////////////////////////////////
            _d_tolerance = __min(_d_tolerance, d_distance);
        }
    }
    _d_tolerance /= 100.0;
    /////////////////////////////////////////////
    for (itr_wr=_v_p_wire.begin(); itr_wr!=_v_p_wire.end(); ++itr_wr)
    {
        Wire* p_wire = *itr_wr;
        /////////////////////////////////////////
        if (p_wire->get_type() == Wire::WT_OPEN) continue;
        /////////////////////////////////////////
        vector<Node*> v_p_wire_node = p_wire->get_node();
        //if (v_p_wire_node.front() != v_p_wire_node.back()) throw exception(CDgnMSGFunction::ConvertW2A(_LS(IDS_BASE_OPEN_OUTER_INNER)));
        /////////////////////////////////////////
        for (int i=0; i<(int)v_p_wire_node.size()-1; ++i)
        {
            //if (v_p_wire_node[i] == v_p_wire_node[i+1]) throw exception(CDgnMSGFunction::ConvertW2A(_LS(IDS_BASE_WIRE_NODE_ERROR)));
        }
        /////////////////////////////////////////
        const double* a_uv_1 = v_p_wire_node[0]->get_uv();
        const double* a_uv_2 = v_p_wire_node[1]->get_uv();
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
        for (int i=0; i<(int)v_p_wire_node.size()-1; ++i)
        {				
            Node* p_node_1 = v_p_wire_node[i  ];
            Node* p_node_2 = v_p_wire_node[i+1];
            d_angle       += angle_2d(a_uv, p_node_1->get_uv(), p_node_2->get_uv());
        }
        d_angle = fabs(d_angle);
        if ((p_wire->get_type() == Wire::WT_OUTER && d_angle < 350.0) || (p_wire->get_type() == Wire::WT_INNER && d_angle > 10.0))
        {
            reverse(v_p_wire_node.begin(), v_p_wire_node.end());
            p_wire->set_node(v_p_wire_node);
        }
    }
    /////////////////////////////////////////////
    for (vector<Node*>::iterator itr_nd=_v_p_island_node.begin(); itr_nd!=_v_p_island_node.end(); ++itr_nd)
    {
        Node*         p_node = *itr_nd;
        const double* a_uv   = p_node->get_uv();
        for (vector<Wire*>::const_iterator itr_wr=_v_p_wire.begin(); itr_wr!=_v_p_wire.end(); ++itr_wr)
        {
            Wire* p_wire = *itr_wr;
            /////////////////////////////////////
            if (p_wire->get_node_index(p_node) != -1)
            {
                _v_p_island_node.erase(itr_nd--);
                break;
            }
            /////////////////////////////////////
            if (p_wire->get_type() == Wire::WT_OPEN) continue;
            /////////////////////////////////////
            const vector<Node*>& v_p_wire_node = p_wire->get_node();
            double               d_angle       = 0.0;
            for (int i=0; i<(int)v_p_wire_node.size()-1; ++i)
            {				
                Node* p_node_1 = v_p_wire_node[i  ];
                Node* p_node_2 = v_p_wire_node[i+1];
                d_angle       += angle_2d(a_uv, p_node_1->get_uv(), p_node_2->get_uv());
            }
            d_angle = fabs(d_angle);
            if ((p_wire->get_type() == Wire::WT_OUTER && d_angle < 350.0) || (p_wire->get_type() == Wire::WT_INNER && d_angle > 10.0))
            {
                _v_p_node.erase(find(_v_p_node.begin(), _v_p_node.end(), p_node));
                delete p_node;
                _v_p_island_node.erase(itr_nd--);
                break;
            }
        }
    }
} // end: Mesher_2D::check_input_data()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Mesher_2D::prepare_mesh() throw(exception)
{
    using namespace mesher_math_lib;
    /////////////////////////////////////////////
    for (vector<Wire*>::const_iterator itr_wr=_v_p_wire.begin(); itr_wr!=_v_p_wire.end(); ++itr_wr)
    {
        Wire* p_wire               = *itr_wr;
        _m_input_wire_node[p_wire] = p_wire->get_node();
    }
    /////////////////////////////////////////////
    if (_d_tolerance < ZERO)
    {
        _d_tolerance = DBL_MAX;
        for (int i=0; i<(int)_v_p_node.size()-1; ++i)
        {
            Node* p_node_1 = _v_p_node[i];
            for (int j=i+1; j<(int)_v_p_node.size(); ++j)
            {
                Node*  p_node_2   = _v_p_node[j];
                double d_distance = distance_2d(p_node_1->get_uv(), p_node_2->get_uv());
                /////////////////////////////////
                //if (d_distance < ZERO) throw exception(CDgnMSGFunction::ConvertW2A(_LS(IDS_BASE_OVERLAP_NODE)));
                /////////////////////////////////
                _d_tolerance = __min(_d_tolerance, d_distance);
            }
        }
        _d_tolerance /= 100.0;
    }
    /////////////////////////////////////////////
    if (_d_mesh_size < ZERO)
    {
        double d_length_sum    = 0.0;
        int    n_segment_count = 0;
        for (vector<Wire*>::const_iterator itr_wr=_v_p_wire.begin(); itr_wr!=_v_p_wire.end(); ++itr_wr)
        {
            Wire*                p_wire        = *itr_wr;
            const vector<Node*>& v_p_wire_node = p_wire->get_node();
            for (int i=0; i<(int)v_p_wire_node.size()-1; ++i)
            {
                Node* p_node_1 = v_p_wire_node[i  ];
                Node* p_node_2 = v_p_wire_node[i+1];
                d_length_sum  += distance_2d(p_node_1->get_uv(), p_node_2->get_uv());
            }
            n_segment_count += (int)v_p_wire_node.size() - 1;
        }
        _d_mesh_size = d_length_sum / n_segment_count;
    }
} // end: Mesher_2D::prepare_mesh()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Mesher_2D::generate_offset_mesh() throw(exception)
{
    using namespace mesher_math_lib;
    /////////////////////////////////////////////
    if (check_interrupt()) throw interrupt_exception();
    /////////////////////////////////////////////
    vector<Wire*>::const_iterator itr_wr;
    for (itr_wr=_v_p_wire.begin(); itr_wr!=_v_p_wire.end(); ++itr_wr)
    {
        Wire*                p_wire        = *itr_wr;
        const vector<Node*>& v_p_wire_node = p_wire->get_node();
        int n_duplicate_count = 0;
        for (int i=0; i<v_p_wire_node.size()-2; ++i)
        {
            for (int j=i+1; j<v_p_wire_node.size()-1; ++j)
            {
                if (v_p_wire_node[i] == v_p_wire_node[j]) return;
            }
        }
    }
    /////////////////////////////////////////////
    vector<Wire*> v_p_wire;
    Wire*         p_outer_wire = 0;
    for (itr_wr=_v_p_wire.begin(); itr_wr!=_v_p_wire.end(); ++itr_wr)
    {
        Wire* p_wire = *itr_wr;
        /////////////////////////////////////////
        switch (p_wire->get_type())
        {
        case Wire::WT_INNER: v_p_wire.push_back(p_wire); break;
        case Wire::WT_OUTER:
            {
                if (p_wire->get_node_size() > 9) p_outer_wire = p_wire;
                break;
            }
        }
    }
    if (p_outer_wire) v_p_wire.push_back(p_outer_wire);
    /////////////////////////////////////////////
    if (v_p_wire.empty()) return;
    /////////////////////////////////////////////
    for (itr_wr=v_p_wire.begin(); itr_wr!=v_p_wire.end(); ++itr_wr)
    {
        Wire* p_wire = *itr_wr;
        /////////////////////////////////////////
        const vector<Node*>& v_p_wire_node = p_wire->get_node();
        /////////////////////////////////////////
        map<Node*, Offset> m_offset;
        double             d_min_deviation = DBL_MAX;
        int                n_base_index    = -1;
        for (int i=0; i<(int)v_p_wire_node.size()-1; ++i)
        {
            Offset offset;
            Node*  p_node_0    = i == 0 ? v_p_wire_node[(int)v_p_wire_node.size()-2] : v_p_wire_node[i-1];
            Node*  p_node_1    = v_p_wire_node[i  ];
            Node*  p_node_2    = v_p_wire_node[i+1];
            offset.set_angle(positive_angle_2d(p_node_1->get_uv(), p_node_2->get_uv(), p_node_0->get_uv()));
            m_offset[p_node_1] = offset;
            /////////////////////////////////////
            if (n_base_index == -1 && offset.get_type() == Offset::OT_SIDE) n_base_index    = i;
        }
        /////////////////////////////////////////
        vector<Node*> v_p_node;
        for (int i=n_base_index; i<(int)v_p_wire_node.size()-1; ++i) v_p_node.push_back(v_p_wire_node[i]);
        for (int i=0;            i<=n_base_index;               ++i) v_p_node.push_back(v_p_wire_node[i]);
        /////////////////////////////////////////
        for (int i=0; i<(int)v_p_node.size()-1; ++i)
        {
            if (check_interrupt()) throw interrupt_exception();
            /////////////////////////////////////
            Node*         p_node_0 = i == 0 ? v_p_node[(int)v_p_node.size()-2] : v_p_node[i-1];
            Node*         p_node_1 = v_p_node[i  ];
            Node*         p_node_2 = v_p_node[i+1];
            const Offset& offset_0 = m_offset[p_node_0];
            Offset& offset_1 = m_offset[p_node_1];
            const Offset& offset_2 = m_offset[p_node_2];
            /////////////////////////////////////
            if (offset_1.get_type() == Offset::OT_END  || (offset_0.get_type() == Offset::OT_END  && offset_2.get_type() == Offset::OT_END))  continue;
            if (offset_1.get_type() != Offset::OT_SIDE && (offset_0.get_type() != Offset::OT_SIDE || offset_2.get_type() != Offset::OT_SIDE)) continue;
            /////////////////////////////////////
            const double  d_node_angle = offset_1.get_angle();
            const double* a_uv_0       = p_node_0->get_uv();
            const double* a_uv_1       = p_node_1->get_uv();
            const double* a_uv_2       = p_node_2->get_uv();
            double        d_length_10  = distance_2d(a_uv_1, a_uv_0);
            double        d_length_12  = distance_2d(a_uv_1, a_uv_2);
            double        d_ratio      = d_length_10 / d_length_12;
            if (d_ratio > 2.0 || d_ratio < 0.5) continue;
            const double  d_interval   = mean_2<double>(d_length_10, d_length_12);
            offset_1.set_interval(d_interval);
            /////////////////////////////////////
            double a_vector_12[2];
            transform(a_uv_2, a_uv_2+2, a_uv_1, a_vector_12, minus<double>());
            transform(a_vector_12, a_vector_12+2, a_vector_12, std::bind(divides<double>(), _1, d_length_12));
            /////////////////////////////////////
            switch (offset_1.get_type())
            {
            case Offset::OT_SIDE:
                {
                    double d_length = d_interval / sin(radian(d_node_angle/2.0));
                    double a_uv[2];
                    transform(a_vector_12, a_vector_12+2, a_uv, std::bind(multiplies<double>(), _1, d_length));
                    transform(a_uv, a_uv+2, a_uv_1, a_uv, plus<double>());
                    rotation_2d(a_uv, a_uv_1, d_node_angle/2.0);
                    offset_1.add_uv(a_uv);
                    break;
                }
            case Offset::OT_CORNER:
                {
                    double d_length_1  = d_interval / sin(radian(d_node_angle/3.0));
                    double d_length_2  = sqrt(2.0) * d_length_1;
                    double a_length[3] = { d_length_1,           d_length_2,       d_length_1       };
                    double a_angle [3] = { d_node_angle*2.0/3.0, d_node_angle/2.0, d_node_angle/3.0 };
                    for (int j=0; j<3; ++j)
                    {
                        double a_uv[2];
                        transform(a_vector_12, a_vector_12+2, a_uv, std::bind(multiplies<double>(), _1, a_length[j]));
                        transform(a_uv, a_uv+2, a_uv_1, a_uv, plus<double>());
                        rotation_2d(a_uv, a_uv_1, a_angle[j]);
                        offset_1.add_uv(a_uv);
                    }
                    break;
                }
            case Offset::OT_REVERSAL:
                {
                    double d_length_1  = d_interval / sin(radian(d_node_angle/4.0));
                    double d_length_2  = sqrt(2.0) * d_length_1;
                    double a_length[5] = { d_length_1,           d_length_2,           d_length_1,       d_length_2,           d_length_1       };
                    double a_angle [5] = { d_node_angle*3.0/4.0, d_node_angle*5.0/8.0, d_node_angle/2.0, d_node_angle*3.0/8.0, d_node_angle/4.0 };
                    for (int j=0; j<5; ++j)
                    {
                        double a_uv[2];
                        transform(a_vector_12, a_vector_12+2, a_uv, std::bind(multiplies<double>(), _1, a_length[j]));
                        transform(a_uv, a_uv+2, a_uv_1, a_uv, plus<double>());
                        rotation_2d(a_uv, a_uv_1, a_angle[j]);
                        offset_1.add_uv(a_uv);
                    }
                    break;
                }
            }
            /////////////////////////////////////
            bool b_self_check = ((offset_0.get_type() == Offset::OT_END) + (offset_2.get_type() == Offset::OT_END) != 1);
            if (!check_offset_node(p_node_1, p_wire, offset_1, m_offset, b_self_check)) offset_1.clear();
        }
        /////////////////////////////////////////
        for (int i=0; i<(int)v_p_wire_node.size()-1; ++i)
        {
            Node*   p_node_0 = i == 0 ? v_p_wire_node[(int)v_p_wire_node.size()-2] : v_p_wire_node[i-1];
            Node*   p_node_1 = v_p_wire_node[i  ];
            Node*   p_node_2 = v_p_wire_node[i+1];
            Offset& offset_0 = m_offset[p_node_0];
            Offset& offset_1 = m_offset[p_node_1];
            Offset& offset_2 = m_offset[p_node_2];
            /////////////////////////////////////
            if (offset_1.get_size() == 1 && offset_2.get_size() == 1)
            {
                if (angle_2d(offset_1.get_uv(0), p_node_1->get_uv(), offset_2.get_uv(0)) < 1.0)
                {
                    offset_1.clear();
                    offset_2.clear();
                    continue;
                }
            }
            /////////////////////////////////////
            if (offset_1.get_size() != 0 && (offset_0.get_size() == 0 && offset_2.get_size() == 0)) offset_1.clear();
            /////////////////////////////////////
            if (offset_1.get_type() != Offset::OT_END)      continue;
            if (offset_0.is_empty() || offset_2.is_empty()) continue;
            double a_uv[2];
            mean_2<double>(offset_0.get_uv(0), offset_2.get_uv(0), a_uv, 2);
            /////////////////////////////////////
            double d_length_10 = distance_2d(p_node_1->get_uv(), p_node_0->get_uv());
            double d_length_12 = distance_2d(p_node_1->get_uv(), p_node_2->get_uv());
            double d_ratio     = d_length_10 / d_length_12;
            if (d_ratio > 2.0 || d_ratio < 0.5)
            {
                offset_0.clear();
                offset_2.clear();
                continue;
            }
            /////////////////////////////////////
            offset_0.set_uv(0, a_uv);
            offset_1.add_uv(a_uv);
            offset_2.set_uv(0, a_uv);
        }
        /////////////////////////////////////////
        map< Node*, vector<Node*> > m_node_offset;
        vector<Node*>               v_p_offset_wire_node;
        for (int i=0; i<(int)v_p_wire_node.size()-1; ++i)
        {
            Node*         p_node = v_p_wire_node[i];
            const Offset& offset = m_offset[p_node];
            vector<Node*> v_p_offset_node;
            if (offset.is_empty())
            {
                m_node_offset[p_node] = v_p_offset_node;
                continue;
            }
            for (int j=0; j<offset.get_size(); ++j)
            {
                v_p_offset_node.push_back(get_unique_node(offset.get_uv(j)));
            }
            m_node_offset[p_node] = v_p_offset_node;
        }
        for (int i=0; i<(int)v_p_wire_node.size()-1; ++i)
        {
            Node*                p_node_0          = i == 0 ? v_p_wire_node[(int)v_p_wire_node.size()-2] : v_p_wire_node[i-1];
            Node*                p_node_1          = v_p_wire_node[i  ];
            Node*                p_node_2          = v_p_wire_node[i+1];
            const Offset&        offset_0          = m_offset[p_node_0];
            const Offset&        offset_1          = m_offset[p_node_1];
            const Offset&        offset_2          = m_offset[p_node_2];
            const vector<Node*>& v_p_offset_node_0 = m_node_offset[p_node_0];
            const vector<Node*>& v_p_offset_node_1 = m_node_offset[p_node_1];
            const vector<Node*>& v_p_offset_node_2 = m_node_offset[p_node_2];
            /////////////////////////////////////
            if (v_p_offset_node_1.empty())
            {
                v_p_offset_wire_node.push_back(p_node_1);
                continue;
            }
            /////////////////////////////////////
            // (OT_END:0)4 -- 3
            //           |    |
            //    (0) -- 1 -- 2
            /////////////////////////////////////
            switch (offset_1.get_type())
            {
            case Offset::OT_END:
                {
                    if (v_p_offset_node_0.empty() || v_p_offset_node_2.empty()) break;
                    Node* p_node_3 = v_p_offset_node_1.front();
                    _v_p_element.push_back(new Element(p_node_1, p_node_2, p_node_3, p_node_0));
                    break;
                }
            case Offset::OT_SIDE:
                {
                    if (v_p_offset_node_2.empty())
                    {
                        v_p_offset_wire_node.push_back(p_node_1);
                        break;
                    }
                    if (offset_2.get_type() == Offset::OT_END)
                    {
                        if (offset_0.is_empty())
                        {
                            v_p_offset_wire_node.push_back(p_node_1);
                            v_p_offset_wire_node.push_back(v_p_offset_node_1.front());
                        }
                        break;
                    }
                    Node* p_node_3 = v_p_offset_node_2.front();
                    Node* p_node_4 = v_p_offset_node_1.front();
                    if (offset_0.is_empty())
                    {
                        v_p_offset_wire_node.push_back(p_node_1);
                        v_p_offset_wire_node.push_back(p_node_4);
                    }
                    _v_p_element.push_back(new Element(p_node_1, p_node_2, p_node_3, p_node_4));
                    v_p_offset_wire_node.push_back(p_node_3);
                    break;
                }
            case Offset::OT_CORNER:
                {
                    if (v_p_offset_node_0.empty() && v_p_offset_node_2.empty())
                    {
                        v_p_offset_wire_node.push_back(p_node_1);
                        break;
                    }
                    if (v_p_offset_node_0.empty())
                    {
                        v_p_offset_wire_node.push_back(p_node_1);
                        v_p_offset_wire_node.push_back(v_p_offset_node_1[0]);
                    }
                    _v_p_element.push_back(new Element(p_node_1, v_p_offset_node_1[2], v_p_offset_node_1[1], v_p_offset_node_1[0]));
                    copy(v_p_offset_node_1.begin()+1, v_p_offset_node_1.end(), back_inserter(v_p_offset_wire_node));
                    if (v_p_offset_node_2.empty())
                    {
                        v_p_offset_wire_node.push_back(p_node_1);
                        break;
                    }
                    Node* p_node_3 = v_p_offset_node_2.front();
                    Node* p_node_4 = v_p_offset_node_1[2];
                    _v_p_element.push_back(new Element(p_node_1, p_node_2, p_node_3, p_node_4));
                    v_p_offset_wire_node.push_back(p_node_3);
                    break;
                }
            case Offset::OT_REVERSAL:
                {
                    if (v_p_offset_node_0.empty() && v_p_offset_node_2.empty())
                    {
                        v_p_offset_wire_node.push_back(p_node_1);
                        break;
                    }
                    if (v_p_offset_node_0.empty()) 
                    {
                        v_p_offset_wire_node.push_back(p_node_1);
                        v_p_offset_wire_node.push_back(v_p_offset_node_1[0]);
                    }
                    _v_p_element.push_back(new Element(p_node_1, v_p_offset_node_1[2], v_p_offset_node_1[1], v_p_offset_node_1[0]));
                    _v_p_element.push_back(new Element(p_node_1, v_p_offset_node_1[4], v_p_offset_node_1[3], v_p_offset_node_1[2]));
                    copy(v_p_offset_node_1.begin()+1, v_p_offset_node_1.end(), back_inserter(v_p_offset_wire_node));
                    if (v_p_offset_node_2.empty())
                    {
                        v_p_offset_wire_node.push_back(p_node_1);
                        break;
                    }
                    Node* p_node_3 = v_p_offset_node_2.front();
                    Node* p_node_4 = v_p_offset_node_1[4];
                    _v_p_element.push_back(new Element(p_node_1, p_node_2, p_node_3, p_node_4));
                    v_p_offset_wire_node.push_back(p_node_3);
                    break;
                }
            }
        }
        v_p_offset_wire_node.push_back(v_p_offset_wire_node.front());
        /////////////////////////////////////////
        p_wire->set_node(v_p_offset_wire_node);
        /////////////////////////////////////////
#ifdef _VERBOSE
        {
            debug_out() << _T("* N(Wire Node) = ") << (int)v_p_offset_wire_node.size() << endl << _T("  > ");
            for (int i=0; i<(int)v_p_offset_wire_node.size(); ++i) debug_out() << v_p_offset_wire_node[i]->get_id() << ' ';
            debug_out() << endl;
        }
#endif
    }
    /////////////////////////////////////////////
    for (int i=0; i<(int)_v_p_wire.size(); ++i)
    {
        const vector<Node*>& v_p_wire_node_1 = _v_p_wire[i]->get_node();
        /////////////////////////////////////////
        for (int j=0; j<(int)v_p_wire_node_1.size()-2; ++j)
        {
            const double* a_uv_1 = v_p_wire_node_1[j  ]->get_uv();
            const double* a_uv_2 = v_p_wire_node_1[j+1]->get_uv();
            for (int k=0; k<j-1; ++k)
            {
                const double* a_uv_3 = v_p_wire_node_1[k  ]->get_uv();
                const double* a_uv_4 = v_p_wire_node_1[k+1]->get_uv();

                //if (line_cross_2d(a_uv_1, a_uv_2, a_uv_3, a_uv_4)) throw exception(CDgnMSGFunction::ConvertW2A(_LS(IDS_BASE_OFFSET_MESH_ERROR)));
            }
        }
        /////////////////////////////////////////
        for (int j=0; j<(int)v_p_wire_node_1.size()-1; ++j)
        {
            const double* a_uv_1 = v_p_wire_node_1[j  ]->get_uv();
            const double* a_uv_2 = v_p_wire_node_1[j+1]->get_uv();
            for (int k=0; k<i; ++k)
            {
                const vector<Node*>& v_p_wire_node_2 = _v_p_wire[k]->get_node();
                for (int l=0; l<(int)v_p_wire_node_2.size()-1; ++l)
                {
                    const double* a_uv_3 = v_p_wire_node_2[l  ]->get_uv();
                    const double* a_uv_4 = v_p_wire_node_2[l+1]->get_uv();
                    //if (line_cross_2d(a_uv_1, a_uv_2, a_uv_3, a_uv_4)) throw exception(CDgnMSGFunction::ConvertW2A(_LS(IDS_BASE_OFFSET_MESH_ERROR)));
                }
            }
        }
    }
    /////////////////////////////////////////////
    vector<Node*>::iterator itr_nd;
    for (itr_nd=_v_p_island_node.begin(); itr_nd!=_v_p_island_node.end(); ++itr_nd)
    {
        const double* a_uv = (*itr_nd)->get_uv();
        for (vector<Element*>::const_iterator itr_el=_v_p_element.begin(); itr_el!=_v_p_element.end(); ++itr_el)
        {
            Element* p_element = *itr_el;
            double   d_angle   = 0.0;
            for (int i=0; i<4; ++i)
            {
                const double* a_uv_1 = (p_element->get_node( i     ))->get_uv();
                const double* a_uv_2 = (p_element->get_node((i+1)%4))->get_uv();
                d_angle             += angle_2d(a_uv, a_uv_1, a_uv_2);
            }
            //if (fabs(d_angle) > 350.0) throw exception(CDgnMSGFunction::ConvertW2A(_LS(IDS_BASE_OFFSET_MESH_ERROR)));
        }
    }
    /////////////////////////////////////////////
    map <Node*, vector<Element*> > m_node_element;
    for (vector<Element*>::const_iterator itr_el=_v_p_element.begin(); itr_el!=_v_p_element.end(); ++itr_el)
    {
        Element* p_element = *itr_el;
        for (int i=0; i<p_element->get_node_size(); ++i) m_node_element[p_element->get_node(i)].push_back(p_element);
    }
    for (itr_nd=_v_p_node.begin(); itr_nd!=_v_p_node.end(); ++itr_nd)
    {
        Node* p_node = *itr_nd;
        if (!p_node->is_hard_node() && m_node_element[p_node].empty())
        {
            delete p_node;
            _v_p_node.erase(itr_nd--);
        }
    }
    /////////////////////////////////////////////
#ifdef _VERBOSE
    {
        dump_mesh(_v_p_node, _v_p_element, _T("OFFSET MESH"));
    }
#endif
} // end: Mesher_2D::generate_offset_mesh()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Mesher_2D::finalize_mesh(Mesh_Packet_2D& mesh_packet) const
{
    map<Node*, bool> m_island_node_flag;
    vector<Node*>::const_iterator itr_nd;
    for (itr_nd=_v_p_node       .begin(); itr_nd!=_v_p_node       .end(); ++itr_nd) m_island_node_flag[*itr_nd] = false;
    for (itr_nd=_v_p_island_node.begin(); itr_nd!=_v_p_island_node.end(); ++itr_nd) m_island_node_flag[*itr_nd] = true;
    /////////////////////////////////////////////
    map<Node*, int> m_node_index;
    int             n_index = 0;
    bool            b_scale_factor = fabs(mesh_packet.d_scale_factor-1.0) > 1.0e-6;
    for (itr_nd=_v_p_node.begin(); itr_nd!=_v_p_node.end(); ++itr_nd)
    {
        Node* p_node         = *itr_nd;
        m_node_index[p_node] = n_index++;
        /////////////////////////////////////////
        if (p_node->is_hard_node() || m_island_node_flag[p_node]) continue;
        /////////////////////////////////////////
        Mesh_Packet_2D::I_Node i_node;
        copy(p_node->get_uv(), p_node->get_uv()+2, i_node.a_uv);
        if (b_scale_factor) transform(i_node.a_uv, i_node.a_uv+2, i_node.a_uv, std::bind(divides<double>(), _1, mesh_packet.d_scale_factor));
        i_node.b_hard_flag   = false;
        i_node.b_island_flag = false;
        i_node.d_mesh_size   = p_node->get_mesh_size();
        /////////////////////////////////////////
        mesh_packet.v_node.push_back(i_node);
    }
    /////////////////////////////////////////////
    for (vector<Element*>::const_iterator itr_el=_v_p_element.begin(); itr_el!=_v_p_element.end(); ++itr_el)
    {
        Element*                  p_element = *itr_el;
        Mesh_Packet_2D::I_Element i_element;
        for (int i=0; i<p_element->get_node_size(); ++i)
        {
            Node* p_node = p_element->get_node(i);
            i_element.v_node_index.push_back(m_node_index[p_node]);
        }
        mesh_packet.v_element.push_back(i_element);
    }
} // end: Mesher_2D::finalize_mesh()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Mesher_2D::check_mesh_shape()
{
    using namespace mesher_math_lib;
    /////////////////////////////////////////////
    vector<Element*> v_p_quadrilateral;
    vector<Element*> v_p_triangle;
    vector<Element*>::iterator itr_el;
    for (itr_el=_v_p_element.begin(); itr_el!=_v_p_element.end(); ++itr_el)
    {
        Element* p_element = *itr_el;
        switch (p_element->get_node_size())
        {
        case 3: v_p_triangle     .push_back(p_element); break;
        case 4: v_p_quadrilateral.push_back(p_element); break;
        }
    }
    /////////////////////////////////////////////
    for (itr_el=v_p_quadrilateral.begin(); itr_el!=v_p_quadrilateral.end(); ++itr_el)
    {
        Element* p_element = *itr_el;
        for (int i=0; i<4; ++i)
        {
            Node*  p_node_1 = p_element->get_node( i     );
            Node*  p_node_2 = p_element->get_node((i+1)%4);
            Node*  p_node_3 = p_element->get_node((i+2)%4);
            Node*  p_node_4 = p_element->get_node((i+3)%4);
            double d_angle  = positive_angle_2d(p_node_1->get_uv(), p_node_2->get_uv(), p_node_4->get_uv());
            if (d_angle > 165.0)
            {
                Element* p_triangle_1 = new Element(p_node_1, p_node_2, p_node_3);
                Element* p_triangle_2 = new Element(p_node_1, p_node_3, p_node_4);
                _v_p_element.push_back(p_triangle_1);
                _v_p_element.push_back(p_triangle_2);
                v_p_triangle.push_back(p_triangle_1);
                v_p_triangle.push_back(p_triangle_2);
                _v_p_element.erase(find(_v_p_element.begin(), _v_p_element.end(), p_element));
                delete p_element;
                break;
            }
        }
    }
    /////////////////////////////////////////////
    for (itr_el=v_p_triangle.begin(); itr_el!=v_p_triangle.end(); ++itr_el)
    {
        Element* p_triangle = *itr_el;
        if (p_triangle->get_distortion_metric() < ZERO)
        {
            _v_p_element.erase(find(_v_p_element.begin(), _v_p_element.end(), p_triangle));
            delete p_triangle;
        }
    }
} // end: Mesher_2D::check_mesh_shape()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Mesher_2D::generate_front_mesh(vector<Front> v_front, bool b_check_in_domain_flag) throw(exception)
{
    using namespace mesher_math_lib;
    /////////////////////////////////////////////
    map<Element*, UV> m_element_center;
    for (vector<Element*>::const_iterator itr_el=_v_p_element.begin(); itr_el!=_v_p_element.end(); ++itr_el)
    {
        Element*  p_element   = *itr_el;
        double    a_uv[2]     = { 0.0, 0.0 };
        const int n_node_size = p_element->get_node_size();
        for (int i=0; i<2; ++i)
        {
            for (int j=0; j<n_node_size; ++j) a_uv[i] += (p_element->get_node(j))->get_uv()[i];
            a_uv[i] /= n_node_size;
        }
        m_element_center[p_element] = a_uv;
    }
    /////////////////////////////////////////////
    const double d_bound              = 5.0 * _d_mesh_size;
    const int    n_initial_front_size = (int)v_front.size();
    for (int i=0; i<(int)v_front.size(); ++i)
    {
        if (check_interrupt()) throw interrupt_exception();
        /////////////////////////////////////////
        if (v_front[i].is_obsolete_front() || !v_front[i].get_node_2()) continue;
        Node*         p_node_1 = v_front[i].get_node_1();
        Node*         p_node_2 = v_front[i].get_node_2();
        const double* a_uv_1   = p_node_1->get_uv();
        const double* a_uv_2   = p_node_2->get_uv();
        double        a_middle_uv[2];
        mean_2<double>(a_uv_1, a_uv_2, a_middle_uv, 2);
        /////////////////////////////////////////
        vector<Node*> v_p_candidate_node;
        /////////////////////////////////////////
        for (int j=0; j<(int)v_front.size(); ++j)
        {
            if (i == j || v_front[j].is_obsolete_front()) continue;
            Node* p_node_3 = v_front[j].get_node_1();
            Node* p_node_4 = v_front[j].get_node_2();
            if (            distance_2d(p_node_3->get_uv(), a_middle_uv) < d_bound) v_p_candidate_node.push_back(p_node_3);
            if (p_node_4 && distance_2d(p_node_4->get_uv(), a_middle_uv) < d_bound) v_p_candidate_node.push_back(p_node_4);
        }
        /////////////////////////////////////////
        Node* p_candidate_node = 0;
        while (!p_candidate_node)
        {
            p_candidate_node   = 0;
            double d_max_angle = 0.1;
            for (vector<Node*>::const_iterator itr_nd=v_p_candidate_node.begin(); itr_nd!=v_p_candidate_node.end(); ++itr_nd)
            {
                Node*  p_node_3 = *itr_nd;
                double d_angle  = angle_2d(p_node_3->get_uv(), a_uv_1, a_uv_2);
                if (d_angle < 179.0 && d_angle > d_max_angle)
                {
                    d_max_angle      = d_angle;
                    p_candidate_node = p_node_3;
                }
            }
            /////////////////////////////////////
            if (!p_candidate_node) break;
            /////////////////////////////////////
            for (int j=0; j<(int)v_front.size(); ++j)
            {
                if (i == j) continue;
                Node* p_node_3 = v_front[j].get_node_1();
                Node* p_node_4 = v_front[j].get_node_2();
                /////////////////////////////////
                if (!p_node_4) continue;
                /////////////////////////////////
                if (distance_2d(v_front[j].get_centroid_uv(), p_candidate_node->get_uv()) > d_bound) continue;
                /////////////////////////////////
                if (line_cross_2d(p_candidate_node->get_uv(), a_uv_1, p_node_3->get_uv(), p_node_4->get_uv()))
                {
                    v_p_candidate_node.erase(find(v_p_candidate_node.begin(), v_p_candidate_node.end(), p_candidate_node));
                    p_candidate_node = 0;
                    break;
                }
                if (line_cross_2d(p_candidate_node->get_uv(), a_uv_2, p_node_3->get_uv(), p_node_4->get_uv()))
                {
                    v_p_candidate_node.erase(find(v_p_candidate_node.begin(), v_p_candidate_node.end(), p_candidate_node));
                    p_candidate_node = 0;
                    break;
                }
            }
            if (!p_candidate_node) continue;
            /////////////////////////////////////
            for (vector<Element*>::const_iterator itr_el=_v_p_element.begin(); itr_el!=_v_p_element.end(); ++itr_el)
            {
                Element* p_element = *itr_el;
                /////////////////////////////////
                if (distance_2d(m_element_center[p_element].a_uv, p_candidate_node->get_uv()) > d_bound) continue;
                /////////////////////////////////
                const int n_node_size = p_element->get_node_size();
                for (int j=0; j<n_node_size; ++j)
                {
                    Node* p_node_3 = p_element->get_node( j   %n_node_size);
                    Node* p_node_4 = p_element->get_node((j+1)%n_node_size);
                    if (line_cross_2d(p_candidate_node->get_uv(), a_uv_1, p_node_3->get_uv(), p_node_4->get_uv()))
                    {
                        v_p_candidate_node.erase(find(v_p_candidate_node.begin(), v_p_candidate_node.end(), p_candidate_node));
                        p_candidate_node = 0;
                        break;
                    }
                    if (line_cross_2d(p_candidate_node->get_uv(), a_uv_2, p_node_3->get_uv(), p_node_4->get_uv()))
                    {
                        v_p_candidate_node.erase(find(v_p_candidate_node.begin(), v_p_candidate_node.end(), p_candidate_node));
                        p_candidate_node = 0;
                        break;
                    }
                }
                if (!p_candidate_node)  break;
            }
        }
        if (!p_candidate_node) continue;
        /////////////////////////////////////////
        bool   b_valid_flag = true;
        double a_uv[2]      = { 0.0, 0.0 };
        for (int j=0; j<2; ++j) a_uv[j] += (p_node_1->get_uv()[j] + p_node_2->get_uv()[j] + p_candidate_node->get_uv()[j]) / 3.0;
        if (b_check_in_domain_flag && !check_in_domain_uv(a_uv)) b_valid_flag = false;
        if (b_valid_flag)
        {
            Element* p_element           = new Element(p_node_1, p_node_2, p_candidate_node);
            _v_p_element.push_back(p_element);
            m_element_center[p_element]  = a_uv;
        }
        /////////////////////////////////////////
        v_front[i].set_obsolete_front();
        /////////////////////////////////////////
        bool b_flag_1 = true;
        bool b_flag_2 = true;
        bool b_flag_3 = true;
        bool b_flag_4 = true;
        for (int j=0; j<(int)v_front.size(); ++j)
        {
            if (i == j) continue;
            Node* p_node_3 = v_front[j].get_node_1();
            Node* p_node_4 = v_front[j].get_node_2();
            if (p_node_4)
            {
                if (p_candidate_node == p_node_3 && p_node_2         == p_node_4) b_flag_3 = false;
                if (p_node_1         == p_node_3 && p_candidate_node == p_node_4) b_flag_4 = false;
                if (p_node_2         == p_node_3 && p_candidate_node == p_node_4)
                {
                    b_flag_1 = false;
                    v_front[j].set_obsolete_front();
                    b_flag_3 = b_flag_3 && (j >= n_initial_front_size);
                }
                if (p_candidate_node == p_node_3 && p_node_1         == p_node_4)
                {
                    b_flag_2 = false;
                    v_front[j].set_obsolete_front();
                    b_flag_4 = b_flag_4 && (j >= n_initial_front_size);
                }
            }
        }
        if (b_flag_1)
        {
            Front front(p_node_2, p_candidate_node);
            front.set_obsolete_front();
            v_front.push_back(front);
        }
        if (b_flag_2)
        {
            Front front(p_candidate_node, p_node_1);
            front.set_obsolete_front();
            v_front.push_back(front);
        }
        if (b_flag_3)
        {
            Front front(p_candidate_node, p_node_2);
            v_front.push_back(front);
        }
        if (b_flag_4)
        {
            Front front(p_node_1, p_candidate_node);
            v_front.push_back(front);
        }
    }
} // end: Mesher_2D::generate_front_mesh()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Mesher_2D::transform_into_triangle_mesh() throw(exception)
{
    using namespace mesher_math_lib;
    /////////////////////////////////////////////
    vector<Element*> v_p_quadrilateral;
    vector<Element*>::iterator itr_el;
    for (itr_el=_v_p_element.begin(); itr_el!=_v_p_element.end(); ++itr_el)
    {
        Element* p_element = *itr_el;
        if (p_element->get_node_size() == 4)
        {
            v_p_quadrilateral.push_back(p_element);
            _v_p_element.erase(itr_el--);
        }
    }
    /////////////////////////////////////////////
    for (itr_el=v_p_quadrilateral.begin(); itr_el!=v_p_quadrilateral.end(); ++itr_el)
    {
        if (check_interrupt()) throw interrupt_exception();
        /////////////////////////////////////////
        Element* p_element     = *itr_el;
        Node*    p_node_1      = p_element->get_node(0);
        Node*    p_node_2      = p_element->get_node(1);
        Node*    p_node_3      = p_element->get_node(2);
        Node*    p_node_4      = p_element->get_node(3);
        Element* p_element_1   = 0;
        Element* p_element_2   = 0;
        double   d_distance_13 = distance_2d(p_node_1->get_uv(), p_node_3->get_uv());
        double   d_distance_24 = distance_2d(p_node_2->get_uv(), p_node_4->get_uv());
        /////////////////////////////////////////
        if (d_distance_13 < d_distance_24)
        {
            p_element_1 = new Element(p_node_1, p_node_2, p_node_3);
            p_element_2 = new Element(p_node_1, p_node_3, p_node_4);
        }
        else
        {
            p_element_1 = new Element(p_node_1, p_node_2, p_node_4);
            p_element_2 = new Element(p_node_2, p_node_3, p_node_4);
        }
        /////////////////////////////////////////
        delete p_element;
        /////////////////////////////////////////
        _v_p_element.push_back(p_element_1);
        _v_p_element.push_back(p_element_2);
    }
} // end: Mesher_2D::transform_into_triangle_mesh()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Mesher_2D::transform_into_combined_mesh(bool b_delaunay_flag, bool b_weak_merge_flag) throw(exception)
{
    using namespace mesher_math_lib;
    /////////////////////////////////////////////
    map< Node*, vector<Element*> > m_node_element;
    vector<Element*>::const_iterator itr_el;
    for (itr_el=_v_p_element.begin(); itr_el!=_v_p_element.end(); ++itr_el)
    {
        Element* p_element = *itr_el;
        for (int i=0; i<p_element->get_node_size(); ++i) m_node_element[p_element->get_node(i)].push_back(p_element);
    }
    /////////////////////////////////////////////
    //  (1-2-3) (1-4-2)           ||  (1-4-3) (3-4-2)
    /////////////////////////////////////////////
    //  * -------- 4 --------- *  ||  * -------- 4 --------- * 
    //  |\       /   \       / |  ||  |\       / | \       / | 
    //  |  \   /  [2]  \   /   |  ||  |  \   /   |   \   /   | 
    //  |    2 --------- 1     |  ||  |    2  [2]|[1]  1     | 
    //  |  /   \  [1]  /   \   |  ||  |  /   \   |   /   \   | 
    //  |/       \   /       \ |  ||  |/       \ | /       \ | 
    //  * -------- 3 --------- *  ||  * -------- 3 --------- * 
    /////////////////////////////////////////////
    for (vector<Element*>::const_iterator itr_el_1=_v_p_element.begin(); itr_el_1!=_v_p_element.end(); ++itr_el_1)
    {
        if (check_interrupt()) throw interrupt_exception();
        /////////////////////////////////////////
        Element* p_element_1 = *itr_el_1;
        /////////////////////////////////////////
        if (p_element_1->get_node_size() != 3) continue;
        /////////////////////////////////////////
        for (int i=0; i<3; ++i)
        {
            Node* p_node_1 = p_element_1->get_node( i     );
            Node* p_node_2 = p_element_1->get_node((i+1)%3);
            /////////////////////////////////////
            if (p_node_1->is_hard_node() || p_node_2->is_hard_node()) continue;
            /////////////////////////////////////
            vector<Element*>& v_p_element_1 = m_node_element[p_node_1]; 
            vector<Element*>& v_p_element_2 = m_node_element[p_node_2];
            if ((int)v_p_element_1.size() != 5 || (int)v_p_element_2.size() != 5) continue;
            /////////////////////////////////////
            bool b_triangle_flag = true;
            for (int j=0; j<(int)v_p_element_1.size(); ++j)
            {
                if (v_p_element_1[j]->get_node_size() != 3)
                {
                    b_triangle_flag = false;
                    break;
                }
            }
            if (!b_triangle_flag) continue;
            for (int j=0; j<(int)v_p_element_2.size(); ++j)
            {
                if (v_p_element_2[j]->get_node_size() != 3)
                {
                    b_triangle_flag = false;
                    break;
                }
            }
            if (!b_triangle_flag) continue;
            /////////////////////////////////////
            Node*    p_node_3    = p_element_1->get_node((i+2)%3);
            Element* p_element_2 = 0;
            for (int j=0; j<(int)v_p_element_1.size(); ++j)
            {
                if (v_p_element_1[j] == p_element_1) continue;
                if (find(v_p_element_2.begin(), v_p_element_2.end(), v_p_element_1[j]) != v_p_element_2.end())
                {
                    p_element_2 = v_p_element_1[j];
                    break;
                }
            }
            if (!p_element_2) continue;
            int   n_index  = p_element_2->get_node_index(p_node_1);
            Node* p_node_4 = p_element_2->get_node((n_index+1)%3);
            p_element_1->change_node((i+1)%3, p_node_4);
            p_element_2->change_node(n_index, p_node_3);
            v_p_element_1.erase(find(v_p_element_1.begin(), v_p_element_1.end(), p_element_2));
            v_p_element_2.erase(find(v_p_element_2.begin(), v_p_element_2.end(), p_element_1));
            m_node_element[p_node_3].push_back(p_element_2);
            m_node_element[p_node_4].push_back(p_element_1);
            break;
        }
    }
    /////////////////////////////////////////////
    //  3 ------- 2  ||  3 ------- 2
    //  | \ [1] / |  ||  |         |
    //  |   \ /   |  ||  |         |
    //  |[2] 1 [4]|  ||  |         |
    //  |   / \   |  ||  |         |
    //  | / [3] \ |  ||  |         |
    //  4 ------- 5  ||  4 ------- 5
    /////////////////////////////////////////////
    for (vector<Node*>::iterator itr_nd=_v_p_node.begin(); itr_nd!=_v_p_node.end(); ++itr_nd)
    {
        if (check_interrupt()) throw interrupt_exception();
        /////////////////////////////////////////
        Node*                   p_node_1    = *itr_nd;
        const vector<Element*>& v_p_element = m_node_element[p_node_1];
        /////////////////////////////////////////
        if (p_node_1->is_hard_node() || (int)v_p_element.size() != 4) continue;
        /////////////////////////////////////////
        bool b_triangle_flag = true;
        for (int i=0; i<(int)v_p_element.size(); ++i)
        {
            if (v_p_element[i]->get_node_size() != 3)
            {
                b_triangle_flag = false;
                break;
            }
        }
        if (!b_triangle_flag) continue;
        /////////////////////////////////////////
        Element* p_element_1 = v_p_element[0];
        int      n_index     = p_element_1->get_node_index(p_node_1);
        Node*    p_node_2    = p_element_1->get_node((n_index+1)%3);
        Node*    p_node_3    = p_element_1->get_node((n_index+2)%3);
        Element* p_element_2 = 0;
        Element* p_element_3 = 0;
        Element* p_element_4 = 0;
        for (int i=0; i<(int)v_p_element.size(); ++i)
        {
            if (v_p_element[i]->get_node_index(p_node_2) == -1 && v_p_element[i]->get_node_index(p_node_3) == -1)
            {
                p_element_3 = v_p_element[i];
                break;
            }
        }
        if (!p_element_3) continue;
        n_index        = p_element_3->get_node_index(p_node_1);
        Node* p_node_4 = p_element_3->get_node((n_index+1)%3);
        Node* p_node_5 = p_element_3->get_node((n_index+2)%3);
        /////////////////////////////////////////
        double d_angle = angle_2d(p_node_2->get_uv(), p_node_3->get_uv(), p_node_5->get_uv());
        if (d_angle < 60.0 || d_angle > 120.0) continue;
        d_angle        = angle_2d(p_node_3->get_uv(), p_node_4->get_uv(), p_node_2->get_uv());
        if (d_angle < 60.0 || d_angle > 120.0) continue;
        /////////////////////////////////////////
        for (int i=0; i<(int)v_p_element.size(); ++i)
        {
            if (v_p_element[i]->get_node_index(p_node_1) != -1 && v_p_element[i]->get_node_index(p_node_3) != -1 && v_p_element[i]->get_node_index(p_node_4) != -1) p_element_2 = v_p_element[i];
            if (v_p_element[i]->get_node_index(p_node_1) != -1 && v_p_element[i]->get_node_index(p_node_5) != -1 && v_p_element[i]->get_node_index(p_node_2) != -1) p_element_4 = v_p_element[i];
        }
        if(!p_element_2 || !p_element_4) continue;
        /////////////////////////////////////////
        Element* p_element = new Element(p_node_2, p_node_3, p_node_4, p_node_5);
        _v_p_element.push_back(p_element);
        /////////////////////////////////////////
        m_node_element[p_node_2].erase(find(m_node_element[p_node_2].begin(), m_node_element[p_node_2].end(), p_element_1));
        m_node_element[p_node_2].erase(find(m_node_element[p_node_2].begin(), m_node_element[p_node_2].end(), p_element_4));
        m_node_element[p_node_3].erase(find(m_node_element[p_node_3].begin(), m_node_element[p_node_3].end(), p_element_1));
        m_node_element[p_node_3].erase(find(m_node_element[p_node_3].begin(), m_node_element[p_node_3].end(), p_element_2));
        m_node_element[p_node_4].erase(find(m_node_element[p_node_4].begin(), m_node_element[p_node_4].end(), p_element_2));
        m_node_element[p_node_4].erase(find(m_node_element[p_node_4].begin(), m_node_element[p_node_4].end(), p_element_3));
        m_node_element[p_node_5].erase(find(m_node_element[p_node_5].begin(), m_node_element[p_node_5].end(), p_element_3));
        m_node_element[p_node_5].erase(find(m_node_element[p_node_5].begin(), m_node_element[p_node_5].end(), p_element_4));
        m_node_element[p_node_2].push_back(p_element);
        m_node_element[p_node_3].push_back(p_element);
        m_node_element[p_node_4].push_back(p_element);
        m_node_element[p_node_5].push_back(p_element);
        /////////////////////////////////////////
        for (int i=0; i<(int)v_p_element.size(); ++i)
        {
            _v_p_element.erase(find(_v_p_element.begin(), _v_p_element.end(), v_p_element[i]));
            delete v_p_element[i];
        }
        delete p_node_1;
        _v_p_node.erase(itr_nd--);
    }
    /////////////////////////////////////////////
    //  2 ------- 4  |  2 ------- 4
    //  | \       |  |  |         |
    //  |   \     |  |  |         |
    //  |     \   |  |  |         |
    //  |       \ |  |  |         |
    //  3 ------- 1  |  3 ------- 1
    /////////////////////////////////////////////
    vector<Element*> v_p_triangle;
    for (itr_el=_v_p_element.begin(); itr_el!=_v_p_element.end(); ++itr_el)
    {
        Element* p_element = *itr_el;
        if (p_element->get_node_size() == 3) v_p_triangle.push_back(p_element);
    }
    set_neighbor_relation(v_p_triangle);
    /////////////////////////////////////////////
    if (b_delaunay_flag) make_delaunay_triangle(v_p_triangle);
    /////////////////////////////////////////////
    const int    n_angle_size                = 4;
    const double a_lower_angle[n_angle_size] = {  79.0,  59.0,  34.0,  19.0 };
    const double a_upper_angle[n_angle_size] = { 101.0, 121.0, 146.0, 161.0 };
    /////////////////////////////////////////////
    map<Element*, bool> m_obsolete_triangle_flag;
    for (itr_el=v_p_triangle.begin(); itr_el!=v_p_triangle.end(); ++itr_el)  m_obsolete_triangle_flag[*itr_el] = false;
    /////////////////////////////////////////////
    for (int i=0; i<n_angle_size; ++i)
    {
        if (!b_weak_merge_flag && (i == n_angle_size-1)) break;
        /////////////////////////////////////////
        for (vector<Element*>::const_iterator itr_el=v_p_triangle.begin(); itr_el!=v_p_triangle.end(); ++itr_el)
        {
            if (check_interrupt()) throw interrupt_exception();
            /////////////////////////////////////
            Element* p_element_1 = *itr_el;
            /////////////////////////////////////
            if (m_obsolete_triangle_flag[p_element_1]) continue;
            /////////////////////////////////////
            int    n_candidate_index     = -1;
            double d_min_angle_deviation = DBL_MAX;
            for (int j=0; j<3; ++j) 
            {
                Element* p_neighbor = p_element_1->get_neighbor(j);
                /////////////////////////////////
                if (!p_neighbor || p_neighbor->get_node_size() != 3 || m_obsolete_triangle_flag[p_neighbor]) continue;
                /////////////////////////////////
                Node* p_node_1 = p_element_1->get_node( j     );
                Node* p_node_2 = p_element_1->get_node((j+1)%3);
                Node* p_node_3 = p_element_1->get_node((j+2)%3);
                Node* p_node_4 = p_neighbor ->get_node((p_neighbor->get_node_index(p_node_1)+1)%3);
                /////////////////////////////////
                if ((i == n_angle_size-1) && (p_node_1->is_hard_node() || p_node_2->is_hard_node())) continue;
                /////////////////////////////////
                double d_min_angle       = 0.0;
                double d_max_angle       = 0.0;
                double d_angle_deviation = quadrilateral_angle_deviation_2d(p_node_1->get_uv(), p_node_4->get_uv(), p_node_2->get_uv(), p_node_3->get_uv(), d_min_angle, d_max_angle);
                if (d_min_angle > a_lower_angle[i] && d_max_angle < a_upper_angle[i] && d_angle_deviation < d_min_angle_deviation)
                {
                    n_candidate_index     = j;
                    d_min_angle_deviation = d_angle_deviation;
                }
            }
            if (n_candidate_index == -1) continue;
            /////////////////////////////////////
            Element* p_element_2 = p_element_1->get_neighbor(n_candidate_index);
            Node*    p_node_1    = p_element_1->get_node( n_candidate_index     );
            Node*    p_node_2    = p_element_1->get_node((n_candidate_index+1)%3);
            Node*    p_node_3    = p_element_1->get_node((n_candidate_index+2)%3);
            Node*    p_node_4    = p_element_2->get_node((p_element_2->get_node_index(p_node_1)+1)%3);
            _v_p_element.push_back(new Element(p_node_1, p_node_4, p_node_2, p_node_3));
            /////////////////////////////////
            for (int j=0; j<3; ++j)
            {
                Element* p_neighbor = p_element_1->get_neighbor(j);
                if (p_neighbor)       p_neighbor ->change_neighbor(p_neighbor->get_neighbor_index(p_element_1), 0);
                p_neighbor          = p_element_2->get_neighbor(j);
                if (p_neighbor)       p_neighbor ->change_neighbor(p_neighbor->get_neighbor_index(p_element_2), 0);
            }
            m_obsolete_triangle_flag[p_element_1] = true;
            m_obsolete_triangle_flag[p_element_2] = true;
        }
    }
    /////////////////////////////////////////////
    for (itr_el=v_p_triangle.begin(); itr_el!=v_p_triangle.end(); ++itr_el)
    {
        Element* p_element = *itr_el;
        if (m_obsolete_triangle_flag[p_element])
        {
            _v_p_element.erase(find(_v_p_element.begin(), _v_p_element.end(), p_element));
            delete p_element;
        }
    }
} // end: Mesher_2D::transform_into_combined_mesh()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Mesher_2D::make_delaunay_triangle(const vector<Element*>& v_p_triangle) const
{
    using namespace mesher_math_lib;
    /////////////////////////////////////////////
    for (int i=0; i<(int)v_p_triangle.size(); ++i)
    {
        Element*      p_element_1 = v_p_triangle[i];
        const double* a_uv_1      = (p_element_1->get_node(0))->get_uv();
        const double* a_uv_2      = (p_element_1->get_node(1))->get_uv();
        const double* a_uv_3      = (p_element_1->get_node(2))->get_uv();
        double        a_center_uv[2], d_radius;
        circumscribed_circle_2d(a_uv_1, a_uv_2, a_uv_3, a_center_uv, d_radius);
        for (int j=0; j<3; ++j)
        {
            Element* p_element_2 = p_element_1->get_neighbor(j);
            if (!p_element_2) continue;
            Node*    p_node_1    = p_element_1->get_node(j);
            int      n_index     = p_element_2->get_node_index(p_node_1);
            Node*    p_node_4    = p_element_2->get_node((n_index+1)%3);
            if (distance_2d(a_center_uv, p_node_4->get_uv())+_d_tolerance < d_radius)
            {
                Node* p_node_3 = p_element_1->get_node((j+2)%3);
                p_element_1->change_node(j,             p_node_4);
                p_element_2->change_node((n_index+2)%3, p_node_3);
                /////////////////////////////////
                Element* p_neighbor_23 = p_element_1->get_neighbor((j      +1)%3);
                Element* p_neighbor_31 = p_element_1->get_neighbor((j      +2)%3);
                Element* p_neighbor_14 = p_element_2->get_neighbor( n_index     );
                Element* p_neighbor_42 = p_element_2->get_neighbor((n_index+1)%3);
                p_element_1->change_neighbor( j,            p_neighbor_42);
                p_element_1->change_neighbor((j      +1)%3, p_neighbor_23);
                p_element_1->change_neighbor((j      +2)%3, p_element_2);
                p_element_2->change_neighbor( n_index,      p_neighbor_14);
                p_element_2->change_neighbor((n_index+1)%3, p_element_1);
                p_element_2->change_neighbor((n_index+2)%3, p_neighbor_31);
                if (p_neighbor_42) p_neighbor_42->change_neighbor(p_neighbor_42->get_neighbor_index(p_element_2), p_element_1);
                if (p_neighbor_31) p_neighbor_31->change_neighbor(p_neighbor_31->get_neighbor_index(p_element_1), p_element_2);
                /////////////////////////////////
                --i;
                break;
            }
        }
    }
} // Mesher_2D::make_delaunay_triangle()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
bool Mesher_2D::relax_combined_mesh()
{
    using namespace mesher_math_lib;
    /////////////////////////////////////////////
    bool b_return_flag = false;
    /////////////////////////////////////////////
    map< Node*, vector<Element*> > m_node_element;
    vector<Element*>::iterator itr_el;
    for (itr_el=_v_p_element.begin(); itr_el!=_v_p_element.end(); ++itr_el)
    {
        Element* p_element = *itr_el;
        for (int i=0; i<p_element->get_node_size(); ++i) m_node_element[p_element->get_node(i)].push_back(p_element);
    }
    set_neighbor_relation(_v_p_element);
    /////////////////////////////////////////////
    //  Valence(1) >= 6
    /////////////////////////////////////////////
    //  4 ---------- 1           ||  4 ---------- 1           
    //  |           / \          ||  |\            \          
    //  |  [2]    /     \        ||  |  \            \        
    //  |       /   [1]   \      ||  |    \    [1]     \      
    //  |     /             \    ||  | [2]  \            \    
    //  * - 2 --------------- 3  ||  * ----- 2 ----------- 3                                
    /////////////////////////////////////////////
    //           1 ---------- 4  ||           1 ---------- 4  
    //          / \           |  ||          /            /|  
    //        /     \    [2]  |  ||        /            /  |  
    //      /   [1]   \       |  ||      /     [1]    /    |  
    //    /             \     |  ||    /            /  [2] |  
    //  2 --------------- 3 - *  ||  2 ----------- 3 ----- *  
    /////////////////////////////////////////////
    for (vector<Node*>::const_iterator itr_nd=_v_p_node.begin(); itr_nd!=_v_p_node.end(); ++itr_nd)
    {
        if (check_interrupt()) throw interrupt_exception();
        /////////////////////////////////////////
        Node* p_node_1 = *itr_nd;
        if ((int)m_node_element[p_node_1].size() < 6) continue;
        /////////////////////////////////////////
        Element* p_element_1           = 0;
        int      n_candidate_index     = -1;
        double   d_min_angle_deviation = DBL_MAX;
        /////////////////////////////////////////
        for (vector<Element*>::const_iterator itr_el=m_node_element[p_node_1].begin(); itr_el!=m_node_element[p_node_1].end(); ++itr_el)
        {
            Element* p_element = *itr_el;
            if (p_element->get_node_size() != 3) continue;
            int      n_index   = p_element->get_node_index(p_node_1);
            Node*    p_node_2  = p_element->get_node((n_index+1)%3);
            Node*    p_node_3  = p_element->get_node((n_index+2)%3);
            for (int i=0; i<2; ++i)
            {
                double d_angle_deviation = 0.0;
                double d_min_angle       = 0.0;
                double d_max_angle       = 0.0;
                if (i == 0)
                {
                    if (check_input_wire_node_pair(p_node_1, p_node_2))                       continue;
                    Element* p_neighbor = p_element->get_neighbor(n_index);
                    if (!p_neighbor || p_neighbor->get_node_size() != 4)                      continue;
                    Node*    p_node_4   = p_neighbor ->get_node((p_neighbor->get_node_index(p_node_1)+1)%4);
                    if (p_node_4->is_hard_node() || (int)m_node_element[p_node_4].size() > 4) continue;
                    d_angle_deviation   = quadrilateral_angle_deviation_2d(p_node_1->get_uv(), p_node_4->get_uv(), p_node_2->get_uv(), p_node_3->get_uv(), d_min_angle, d_max_angle);
                }
                else
                {
                    if (check_input_wire_node_pair(p_node_1, p_node_3))                       continue;
                    Element* p_neighbor = p_element->get_neighbor((n_index+2)%3);
                    if (!p_neighbor || p_neighbor->get_node_size() != 4)                      continue;
                    Node*    p_node_4   = p_neighbor ->get_node((p_neighbor->get_node_index(p_node_1)+3)%4);
                    if (p_node_4->is_hard_node() || (int)m_node_element[p_node_4].size() > 4) continue;
                    d_angle_deviation   = quadrilateral_angle_deviation_2d(p_node_1->get_uv(), p_node_2->get_uv(), p_node_3->get_uv(), p_node_4->get_uv(), d_min_angle, d_max_angle);
                }
                if (d_min_angle < 15.0 || d_max_angle > 165.0) continue;
                if (d_angle_deviation < d_min_angle_deviation)
                {
                    p_element_1           = p_element;
                    n_candidate_index     = i;
                    d_min_angle_deviation = d_angle_deviation;
                }
            }
            if (n_candidate_index == -1) continue;
        }
        if (!p_element_1) continue;
        /////////////////////////////////////////
        int      n_index     = p_element_1->get_node_index(p_node_1);
        Node*    p_node_2    = p_element_1->get_node((n_index+1)%3);
        Node*    p_node_3    = p_element_1->get_node((n_index+2)%3);
        Element* p_element_2 = n_candidate_index == 0
            ? p_element_1->get_neighbor( n_index     )
            : p_element_1->get_neighbor((n_index+2)%3);
        Node*    p_node_4    = n_candidate_index == 0
            ? p_element_2->get_node((p_element_2->get_node_index(p_node_1)+1)%4)
            : p_element_2->get_node((p_element_2->get_node_index(p_node_1)+3)%4);
        /////////////////////////////////////////
        set<Element*> s_p_neighbor;
        s_p_neighbor.insert(p_element_1);
        s_p_neighbor.insert(p_element_2);			
        for (int i=0; i<p_element_1->get_node_size(); ++i)
        {
            Element* p_neighbor = p_element_1->get_neighbor(i);
            if (p_neighbor) s_p_neighbor.insert(p_neighbor); 
        }
        for (int i=0; i<p_element_2->get_node_size(); ++i)
        {
            Element* p_neighbor = p_element_2->get_neighbor(i);
            if (p_neighbor) s_p_neighbor.insert(p_neighbor);
        }
        /////////////////////////////////////////
        switch (n_candidate_index)
        {
        case 0:
            p_element_1->insert_node((n_index+1)%3, p_node_4);
            p_element_2->remove_node(p_node_1);
            m_node_element[p_node_1].erase(find(m_node_element[p_node_1].begin(), m_node_element[p_node_1].end(), p_element_2));
            m_node_element[p_node_4].push_back(p_element_1);
            break;
        case 1:
            p_element_1->insert_node((n_index+3)%3, p_node_4);
            p_element_2->remove_node(p_node_1);
            m_node_element[p_node_1].erase(find(m_node_element[p_node_1].begin(), m_node_element[p_node_1].end(), p_element_2));
            m_node_element[p_node_4].push_back(p_element_1);
            break;
        }
        /////////////////////////////////////////
        update_neighbor_relation(s_p_neighbor);
        /////////////////////////////////////////
        b_return_flag = true;
    }
    /////////////////////////////////////////////
    bool b_continue_flag = true;
    while (b_continue_flag)
    {
        b_continue_flag = false;
        /////////////////////////////////////////
        //  * -------- 1 -------- *  ||  * -------- 1 -------- *  
        //    \  [2]  / \  [5]  /    ||    \   [2]  |  [5]   /    
        //      \   /     \   /      ||      \      |      /      
        //        2   [1]   4        ||        2 -- 5 -- 4        
        //      /   \     /   \      ||      /      |      \      
        //    /  [3]  \ /  [4]  \    ||    /   [3]  |  [4]   \    
        //  * -------- 3 -------- *  ||  * -------- 3 -------- *  
        /////////////////////////////////////////
        for (vector<Element*>::iterator itr_el=_v_p_element.begin(); itr_el!=_v_p_element.end(); ++itr_el)
        {
            if (check_interrupt()) throw interrupt_exception();
            /////////////////////////////////////
            Element* p_element_1 = *itr_el;
            /////////////////////////////////////
            if (p_element_1->get_node_size() != 4) continue;
            /////////////////////////////////////
            int n_triangle_count = 0;
            for (int i=0; i<4; ++i)
            {
                Element* p_neighbor = p_element_1->get_neighbor(i);
                if (!p_neighbor)
                {
                    n_triangle_count = 0;
                    break;
                }
                if (p_neighbor->get_node_size() == 3) ++n_triangle_count;
            }
            if (n_triangle_count != 4) continue;
            /////////////////////////////////////
            bool b_free_side_flag = true;
            for (int i=0; i<4; ++i)
            {
                Node* p_node_1 = p_element_1->get_node( i     );
                Node* p_node_2 = p_element_1->get_node((i+1)%4);
                if (check_input_wire_node_pair(p_node_1, p_node_2)) 
                {
                    b_free_side_flag = false;
                    break;
                }
            }
            if (!b_free_side_flag) continue;
            /////////////////////////////////////
            Node*    a_p_node   [4];
            Element* a_p_element[4];
            double   a_uv_5     [2] = { 0.0, 0.0 };
            for (int i=0; i<4; ++i)
            {
                a_p_node[i]    = p_element_1->get_node(i);
                a_p_element[i] = p_element_1->get_neighbor(i);
                for (int j=0; j<2; ++j) a_uv_5[j] += a_p_node[i]->get_uv()[j];
            }
            for (int i=0; i<2; ++i) a_uv_5[i] /= 4.0;
            /////////////////////////////////////
            bool b_neighbor_flag = false;
            for (int i=0; i<4; ++i)
            {
                for (int j=0; j<a_p_element[i]->get_neighbor_size(); ++j)
                {
                    for (int k=0; k<4; ++k)
                    {
                        if (i == k) continue;
                        if (a_p_element[i]->get_neighbor(j) == a_p_element[k])
                        {
                            b_neighbor_flag = true;
                            break;
                        }
                    }
                    if (b_neighbor_flag) break;
                }
            }
            if (b_neighbor_flag) continue;
            /////////////////////////////////////
            double d_angle = angle_2d(a_p_node[0]->get_uv(), a_p_node[1]->get_uv(), a_p_node[3]->get_uv());
            if (d_angle < 60.0 || d_angle > 120.0) continue;
            d_angle        = angle_2d(a_p_node[1]->get_uv(), a_p_node[2]->get_uv(), a_p_node[0]->get_uv());
            if (d_angle < 60.0 || d_angle > 120.0) continue;
            /////////////////////////////////////
            Node* p_node_5 = new Node(a_uv_5);
            _v_p_node.push_back(p_node_5);
            /////////////////////////////////////
            set<Element*> s_p_neighbor;
            for (int i=0; i<4; ++i)
            {
                a_p_element[i]->insert_node(a_p_element[i]->get_node_index(a_p_node[i]), p_node_5);
                m_node_element[a_p_node[i]].erase(find(m_node_element[a_p_node[i]].begin(), m_node_element[a_p_node[i]].end(), p_element_1));
                m_node_element[p_node_5].push_back(a_p_element[i]);
                /////////////////////////////////
                s_p_neighbor.insert(a_p_element[i]);
                for (int j=0; j<3; ++j)
                {
                    Element* p_neighbor = a_p_element[i]->get_neighbor(j);
                    if (p_neighbor && p_neighbor != p_element_1) s_p_neighbor.insert(p_neighbor);
                }
            }
            update_neighbor_relation(s_p_neighbor);
            /////////////////////////////////////
            delete p_element_1;
            _v_p_element.erase(itr_el--);
            /////////////////////////////////////
            b_continue_flag = true;
        }
        /////////////////////////////////////////
        //  * -------- 1 -------- *  ||  * -------- 1 -------- *
        //  |    [2]  / \  [5]  /    ||  | [2]   /  |  [5]   /  
        //  |       /     \   /      ||  |     / [1]|      /    
        //  * --- 2   [1]   4        ||  * - 2 ---- 5 -- 4         
        //      /   \     /   \      ||     /       |      \    
        //    /  [3]  \ /  [4]  \    ||   /    [3]  |  [4]   \  
        //  * -------- 3 -------- *  ||  * -------- 3 -------- *
        /////////////////////////////////////////
        for (vector<Element*>::iterator itr_el=_v_p_element.begin(); itr_el!=_v_p_element.end(); ++itr_el)
        {
            if (check_interrupt()) throw interrupt_exception();
            /////////////////////////////////////
            Element* p_element_1 = *itr_el;
            /////////////////////////////////////
            if (p_element_1->get_node_size() != 4) continue;
            /////////////////////////////////////
            int n_triangle_count =  0;
            int n_index          = -1;
            for (int i=0; i<4; ++i)
            {
                Element* p_neighbor = p_element_1->get_neighbor(i);
                if (!p_neighbor)
                {
                    n_triangle_count = 0;
                    break;
                }
                if (p_neighbor->get_node_size() == 3) ++n_triangle_count;
                if (p_neighbor->get_node_size() == 4) n_index = i;
            }
            if (n_triangle_count != 3) continue;
            /////////////////////////////////////
            bool b_free_side_flag = true;
            for (int i=0; i<4; ++i)
            {
                Node* p_node_1 = p_element_1->get_node( i     );
                Node* p_node_2 = p_element_1->get_node((i+1)%4);
                if (check_input_wire_node_pair(p_node_1, p_node_2)) 
                {
                    b_free_side_flag = false;
                    break;
                }
            }
            if (!b_free_side_flag) continue;
            /////////////////////////////////////
            Node*    a_p_node   [4];
            Element* a_p_element[4];
            double   a_uv_5     [2] = { 0.0, 0.0 };
            for (int i=0; i<4; ++i)
            {
                a_p_node[i]    = p_element_1->get_node((n_index+i)%4);
                a_p_element[i] = p_element_1->get_neighbor((n_index+i)%4);
                for (int j=0; j<2; ++j) a_uv_5[j] += a_p_node[i]->get_uv()[j];
            }
            for (int i=0; i<2; ++i) a_uv_5[i] /= 4.0;
            /////////////////////////////////////
            bool b_neighbor_flag = false;
            for (int i=0; i<4; ++i)
            {
                for (int j=0; j<a_p_element[i]->get_neighbor_size(); ++j)
                {
                    for (int k=0; k<4; ++k)
                    {
                        if (i == k) continue;
                        if (a_p_element[i]->get_neighbor(j) == a_p_element[k])
                        {
                            b_neighbor_flag = true;
                            break;
                        }
                    }
                    if (b_neighbor_flag) break;
                }
            }
            if (b_neighbor_flag) continue;
            /////////////////////////////////////
            double d_angle = angle_2d(a_p_node[0]->get_uv(), a_p_node[1]->get_uv(), a_p_node[3]->get_uv());
            if (d_angle < 60.0 || d_angle > 120.0) continue;
            d_angle        = angle_2d(a_p_node[1]->get_uv(), a_p_node[2]->get_uv(), a_p_node[0]->get_uv());
            if (d_angle < 60.0 || d_angle > 120.0) continue;
            /////////////////////////////////////
            Node* p_node_5 = new Node(a_uv_5);
            _v_p_node.push_back(p_node_5);
            /////////////////////////////////////
            for (int i=1; i<4; ++i) a_p_element[i]->insert_node(a_p_element[i]->get_node_index(a_p_node[i]), p_node_5);
            for (int i=2; i<4; ++i) m_node_element[a_p_node[i]].erase(find(m_node_element[a_p_node[i]].begin(), m_node_element[a_p_node[i]].end(), p_element_1));
            for (int i=1; i<4; ++i) m_node_element[p_node_5].push_back(a_p_element[i]);
            m_node_element[p_node_5].push_back(p_element_1);
            p_element_1->change_node((n_index+2)%4, p_node_5);
            p_element_1->remove_node(a_p_node[3]);			
            /////////////////////////////////
            set<Element*> s_p_neighbor;
            s_p_neighbor.insert(p_element_1);
            for (int i=0; i<4; ++i)
            {
                s_p_neighbor.insert(a_p_element[i]);
                for (int j=0; j<a_p_element[i]->get_neighbor_size(); ++j)
                {
                    Element* p_neighbor = a_p_element[i]->get_neighbor(j);
                    if (p_neighbor) s_p_neighbor.insert(p_neighbor);
                }
            }
            /////////////////////////////////////
            update_neighbor_relation(s_p_neighbor);
            /////////////////////////////////////
            b_continue_flag = true;
        }
    }
    /////////////////////////////////////////////
    // (Edge Nodes: 1, 2, 3)     ||              6
    //             6             ||             / \
    //            / \            ||           /     \
    //          /     \          ||         4   [3]   5
    //        4   [3]   5        ||       /   \     /   \
    //      /   \     /   \      ||     /        7        \
    //    /  [1]  \ /  [2]  \    ||   /    [1]   |   [2]    \
    //  1 -------- 2 -------- 3  ||  1 --------- 2 ---------- 3
    /////////////////////////////////////////////
    vector<Wire*>::const_iterator itr_wr;
    for (itr_wr=_v_p_wire.begin(); itr_wr!=_v_p_wire.end(); ++itr_wr)
    {
        if (check_interrupt()) throw interrupt_exception();
        /////////////////////////////////////////
        Wire*                p_wire        = *itr_wr;
        const vector<Node*>& v_p_wire_node = p_wire->get_node();
        for (int i=1; i<(int)v_p_wire_node.size()-1; ++i)
        {
            Node* p_node_2 = v_p_wire_node[i  ];
            if ((int)m_node_element[p_node_2].size() != 3) continue;
            Node* p_node_1 = v_p_wire_node[i-1];
            Node* p_node_3 = v_p_wire_node[i+1];
            /////////////////////////////////////
            double d_angle = angle_2d(p_node_2->get_uv(), p_node_3->get_uv(), p_node_1->get_uv());
            if (d_angle < 165.0 || d_angle > 195.0) continue;
            /////////////////////////////////////
            Element* p_element_1 = 0;
            Element* p_element_2 = 0;
            for (int j=0; j<3; ++j)
            {
                Element* p_element = m_node_element[p_node_2][j];
                if (find(m_node_element[p_node_1].begin(), m_node_element[p_node_1].end(), p_element) != m_node_element[p_node_1].end()) p_element_1 = p_element;
                if (find(m_node_element[p_node_3].begin(), m_node_element[p_node_3].end(), p_element) != m_node_element[p_node_3].end()) p_element_2 = p_element;
            }
            if (!p_element_1 || p_element_1->get_node_size() != 3) continue;
            if (!p_element_2 || p_element_2->get_node_size() != 3) continue;
            Element* p_element_3 = p_element_1->get_neighbor(p_element_1->get_node_index(p_node_2));
            if (!p_element_3 || p_element_3->get_node_size() != 4) continue;
            /////////////////////////////////////
            int   n_index_2 = p_element_3->get_node_index(p_node_2);
            Node* p_node_4  = p_element_3->get_node(((n_index_2)+3)%4);
            Node* p_node_5  = p_element_3->get_node(((n_index_2)+1)%4);
            Node* p_node_6  = p_element_3->get_node(((n_index_2)+2)%4);
            /////////////////////////////////////
            if (check_input_wire_node_pair(p_node_2, p_node_4)) continue;
            if (check_input_wire_node_pair(p_node_2, p_node_5)) continue;
            /////////////////////////////////////
            d_angle = angle_2d(p_node_2->get_uv(), p_node_5->get_uv(), p_node_4->get_uv());
            if (d_angle < 45.0 || d_angle > 135.0)                                                                       continue;
            if (p_node_1->is_hard_node() && angle_2d(p_node_1->get_uv(), p_node_2->get_uv(), p_node_4->get_uv()) < 45.0) continue;
            if (p_node_3->is_hard_node() && angle_2d(p_node_3->get_uv(), p_node_5->get_uv(), p_node_2->get_uv()) < 45.0) continue;
            /////////////////////////////////////
            set<Element*> s_p_neighbor;
            s_p_neighbor.insert(p_element_1);
            s_p_neighbor.insert(p_element_2);
            for (int j=0; j<p_element_1->get_node_size(); ++j)
            {
                Element* p_neighbor = p_element_1->get_neighbor(j);
                if (p_neighbor) s_p_neighbor.insert(p_neighbor);
            }
            for (int j=0; j<p_element_2->get_node_size(); ++j)
            {
                Element* p_neighbor = p_element_2->get_neighbor(j);
                if (p_neighbor) s_p_neighbor.insert(p_neighbor);
            }
            /////////////////////////////////////
            double a_uv_7[2];
            mean_2<double>(p_node_4->get_uv(), p_node_5->get_uv(), a_uv_7, 2);
            mean_2<double>(a_uv_7,             p_node_2->get_uv(), a_uv_7, 2);
            Node*  p_node_7 = get_unique_node(a_uv_7);
            p_element_1->insert_node((p_element_1->get_node_index(p_node_2)+1) % 3, p_node_7);
            p_element_2->insert_node((p_element_2->get_node_index(p_node_2)+3) % 3, p_node_7);
            p_element_3->change_node( p_element_3->get_node_index(p_node_2),        p_node_7);
            /////////////////////////////////////
            m_node_element[p_node_2].erase(find(m_node_element[p_node_2].begin(), m_node_element[p_node_2].end(), p_element_3));
            /////////////////////////////////////
            update_neighbor_relation(s_p_neighbor);
        }
    }
    /////////////////////////////////////////////
    // (Edge Nodes: 1, 2)       ||
    //                *         ||                *          
    //              /   \       ||              /   \        
    //            /       \     ||            /  [2]  \      
    //          /           \   ||          /           \    
    //         3     [2]     4  ||         3 ----------- 4   
    //       /   \          /   ||       /              /    
    //     /  [1]  \      /     ||     /     [1]      /      
    //   /           \  /       ||   /              /        
    //  1 ----------- 2         ||  1 ----------- 2          
    /////////////////////////////////////////////
    //         *                ||         *               
    //       /   \              ||       /   \             
    //     /       \            ||     /  [2]  \           
    //   /           \          ||   /           \         
    //  4     [2]     3         ||  4 ----------- 3        
    //   \          /   \       ||   \              \      
    //     \      /  [1]  \     ||     \     [1]      \    
    //       \  /           \   ||       \              \  
    //         1 ----------- 2  ||         1 ----------- 2 
    /////////////////////////////////////////////
    for (itr_wr=_v_p_wire.begin(); itr_wr!=_v_p_wire.end(); ++itr_wr)
    {
        if (check_interrupt()) throw interrupt_exception();
        /////////////////////////////////////////
        Wire*                p_wire        = *itr_wr;
        const vector<Node*>& v_p_wire_node = p_wire->get_node();
        for (int i=0; i<(int)v_p_wire_node.size()-1; ++i)
        {
            Node*    p_node_1    = v_p_wire_node[i  ];
            Node*    p_node_2    = v_p_wire_node[i+1];
            Element* p_element_1 = 0;
            for (vector<Element*>::iterator itr_el=m_node_element[p_node_1].begin(); itr_el!=m_node_element[p_node_1].end(); ++itr_el)
            {
                Element* p_element = *itr_el;
                if (p_element->get_node_size() == 3) 
                {
                    if (find(m_node_element[p_node_2].begin(), m_node_element[p_node_2].end(), p_element) != m_node_element[p_node_2].end())
                    {
                        p_element_1 = p_element;
                        break;
                    }
                }
            }
            if (!p_element_1) continue;
            /////////////////////////////////////
            int    n_index               = p_element_1->get_node_index(p_node_1);
            Node*  p_node_3              = p_element_1->get_node((n_index+2)%3);
            double d_min_angle_deviation = DBL_MAX;
            int    n_candidate_index     = -1;
            double d_min_angle           = 0.0;
            double d_max_angle           = 0.0;
            double d_angle_deviation     = 0.0;
            for (int j=0; j<2; ++j)
            {
                if (j == 0)
                {
                    if (check_input_wire_node_pair(p_node_2, p_node_3))                              continue;
                    if (angle_2d(p_node_1->get_uv(), p_node_2->get_uv(), p_node_3->get_uv()) < 45.0) continue;
                }
                else
                {
                    if (check_input_wire_node_pair(p_node_1, p_node_3))                              continue;
                    if (angle_2d(p_node_2->get_uv(), p_node_3->get_uv(), p_node_1->get_uv()) < 45.0) continue;
                }
                /////////////////////////////////
                Element* p_element_2 = p_element_1->get_neighbor((n_index+j+1)%3);
                if (!p_element_2 || p_element_2->get_node_size() != 4) continue;
                Node*    p_node_4    = p_element_2->get_node((p_element_2->get_node_index(p_node_3)+2)%4);
                if (p_node_4->is_hard_node()) continue;
                d_angle_deviation    = j == 0
                    ? quadrilateral_angle_deviation_2d(p_node_1->get_uv(), p_node_2->get_uv(), p_node_4->get_uv(), p_node_3->get_uv(), d_min_angle, d_max_angle)
                    : quadrilateral_angle_deviation_2d(p_node_1->get_uv(), p_node_2->get_uv(), p_node_3->get_uv(), p_node_4->get_uv(), d_min_angle, d_max_angle);
                if (d_min_angle < 15.0 || d_max_angle > 165.0) continue;
                if (d_angle_deviation < d_min_angle_deviation)
                {
                    n_candidate_index     = j;
                    d_min_angle_deviation = d_angle_deviation;
                }
            }
            if (n_candidate_index == -1) continue;
            /////////////////////////////////////
            Element* p_element_2 = p_element_1->get_neighbor((n_index+n_candidate_index+1)%3);
            Node*    p_node_4    = p_element_2->get_node((p_element_2->get_node_index(p_node_3)+2)%4);
            /////////////////////////////////////
            set<Element*> s_p_neighbor;
            s_p_neighbor.insert(p_element_1);
            s_p_neighbor.insert(p_element_2);			
            for (int j=0; j<p_element_1->get_node_size(); ++j)
            {
                Element* p_neighbor = p_element_1->get_neighbor(j);
                if (p_neighbor) s_p_neighbor.insert(p_neighbor); 
            }
            for (int j=0; j<p_element_2->get_node_size(); ++j)
            {
                Element* p_neighbor = p_element_2->get_neighbor(j);
                if (p_neighbor) s_p_neighbor.insert(p_neighbor);
            }
            /////////////////////////////////////
            switch (n_candidate_index)
            {
            case 0:
                p_element_1->insert_node((n_index+2)%3, p_node_4);
                p_element_2->remove_node(p_node_2);
                m_node_element[p_node_2].erase(find(m_node_element[p_node_2].begin(), m_node_element[p_node_2].end(), p_element_2));
                m_node_element[p_node_4].push_back(p_element_1);
                break;
            case 1:
                p_element_1->insert_node((n_index+3)%3, p_node_4);
                p_element_2->remove_node(p_node_1);
                m_node_element[p_node_1].erase(find(m_node_element[p_node_1].begin(), m_node_element[p_node_1].end(), p_element_2));
                m_node_element[p_node_4].push_back(p_element_1);
                break;
            }
            /////////////////////////////////////
            update_neighbor_relation(s_p_neighbor);
            /////////////////////////////////////
            b_return_flag = true;
        }
    }
    /////////////////////////////////////////////
    // element elimination
    /////////////////////////////////////////////
    //  * - * - *  |  * - * - *
    //  | a | d |  |  |   |   |
    //  |   1   |  |  | a | d | 
    //  | /   \ |  |  |   |   | 
    //  2       4  |  2 - * - 4 
    //  | \   / |  |  |   |   | 
    //  |   3   |  |  | b | c | 
    //  | b | c |  |  |   |   | 
    //  * - * - *  |  * - * - * 
    /////////////////////////////////////////
    for (itr_el=_v_p_element.begin(); itr_el!=_v_p_element.end(); ++itr_el)
    {
        if (check_interrupt()) throw interrupt_exception();
        /////////////////////////////////////////
        Element* p_element = *itr_el;
        /////////////////////////////////////////
        if (p_element->get_node_size() != 4) continue;
        /////////////////////////////////////////
        for (int i=0; i<2; ++i)
        {
            Node* p_node_1 = p_element->get_node( i     );
            Node* p_node_3 = p_element->get_node((i+2)%4);
            /////////////////////////////////////
            if (p_node_1->is_hard_node() || (int)m_node_element[p_node_1].size() != 3) continue;
            if (p_node_3->is_hard_node() || (int)m_node_element[p_node_3].size() != 3) continue;
            /////////////////////////////////////
            Node* p_node_2 = p_element->get_node((i+1)%4);
            Node* p_node_4 = p_element->get_node((i+3)%4);
            if ((int)m_node_element[p_node_2].size() == 3)
            {
                if (!p_node_2->is_hard_node())                                                             continue;
                if (positive_angle_2d(p_node_2->get_uv(), p_node_3->get_uv(), p_node_1->get_uv()) > 100.0) continue;
            }
            if ((int)m_node_element[p_node_4].size() == 3)
            {
                if (!p_node_4->is_hard_node())                                                             continue;
                if (positive_angle_2d(p_node_4->get_uv(), p_node_1->get_uv(), p_node_3->get_uv()) > 100.0) continue;
            }
            /////////////////////////////////////
            Element* p_element_1 = p_element->get_neighbor( i     );
            Element* p_element_2 = p_element->get_neighbor((i+1)%4);
            Element* p_element_3 = p_element->get_neighbor((i+2)%4);
            Element* p_element_4 = p_element->get_neighbor((i+3)%4);
            if (!p_element_1 || p_element_1->get_node_size() != 4) continue;
            if (!p_element_2 || p_element_2->get_node_size() != 4) continue;
            if (!p_element_3 || p_element_3->get_node_size() != 4) continue;
            if (!p_element_4 || p_element_4->get_node_size() != 4) continue;
            /////////////////////////////////////
            double a_uv[2];
            mean_2<double>(p_node_1->get_uv(), p_node_3->get_uv(), a_uv, 2);
            Node*  p_node = new Node(a_uv);
            _v_p_node.push_back(p_node);
            /////////////////////////////////////
            p_element_1->change_node    (p_element_1->get_node_index(p_node_1),      p_node);
            p_element_1->change_neighbor(p_element_1->get_neighbor_index(p_element), p_element_2);
            p_element_2->change_node    (p_element_2->get_node_index(p_node_3),      p_node);
            p_element_2->change_neighbor(p_element_2->get_neighbor_index(p_element), p_element_1);
            p_element_3->change_node    (p_element_3->get_node_index(p_node_3),      p_node);
            p_element_3->change_neighbor(p_element_3->get_neighbor_index(p_element), p_element_4);
            p_element_4->change_node    (p_element_4->get_node_index(p_node_1),      p_node);
            p_element_4->change_neighbor(p_element_4->get_neighbor_index(p_element), p_element_3);
            /////////////////////////////////////
            m_node_element[p_node_2].erase(find(m_node_element[p_node_2].begin(), m_node_element[p_node_2].end(), p_element));
            m_node_element[p_node_4].erase(find(m_node_element[p_node_4].begin(), m_node_element[p_node_4].end(), p_element));
            m_node_element.erase(p_node_1);
            m_node_element.erase(p_node_3);
            m_node_element[p_node  ].push_back(p_element_1);
            m_node_element[p_node  ].push_back(p_element_2);
            m_node_element[p_node  ].push_back(p_element_3);
            m_node_element[p_node  ].push_back(p_element_4);
            /////////////////////////////////////
            _v_p_node.erase(find(_v_p_node.begin(), _v_p_node.end(), p_node_1));
            _v_p_node.erase(find(_v_p_node.begin(), _v_p_node.end(), p_node_3));
            delete p_node_1;
            delete p_node_3;
            delete p_element;
            _v_p_element.erase(itr_el--);
            /////////////////////////////////////
            b_return_flag = true;
            break;
        }
    }
    /////////////////////////////////////////////
    return b_return_flag;
} // end: Mesher_2D::relax_combined_mesh()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
bool Mesher_2D::relax_quadrilateral_mesh(bool b_node_elimination_flag, bool b_element_elimination_flag, bool b_diagonal_swapping_flag, bool b_side_elimination_flag)
{
    using namespace mesher_math_lib;
    /////////////////////////////////////////////
    bool b_return_flag = false;
    /////////////////////////////////////////////
    map< Node*, vector<Element*> > m_node_element;
    for (vector<Element*>::const_iterator itr_el=_v_p_element.begin(); itr_el!=_v_p_element.end(); ++itr_el)
    {
        Element* p_element = *itr_el;
        for (int i=0; i<p_element->get_node_size(); ++i) m_node_element[p_element->get_node(i)].push_back(p_element);
    }
    /////////////////////////////////////////////
    // node elimination
    /////////////////////////////////////////////
    //      1      |      1             
    //    / | \    |    /   \          
    //  2 a * b 4  |  2   a   4        
    //    \ | /    |    \   /          
    //      3      |      3            
    /////////////////////////////////////////////
    if (b_node_elimination_flag)
    {
        for (vector<Node*>::iterator itr_nd=_v_p_node.begin(); itr_nd!=_v_p_node.end(); ++itr_nd)
        {
            if (check_interrupt()) throw interrupt_exception();
            /////////////////////////////////////
            Node* p_node = *itr_nd;
            /////////////////////////////////////
            if (p_node->is_hard_node() || (int)m_node_element[p_node].size() != 2) continue;
            /////////////////////////////////////
            Element* p_element_1 = m_node_element[p_node][0];
            Element* p_element_2 = m_node_element[p_node][1];
            /////////////////////////////////////
            if (p_element_1->get_node_size() != 4 || p_element_2->get_node_size() != 4) continue;
            /////////////////////////////////////
            int   n_index_1 = p_element_1->get_node_index(p_node);
            int   n_index_2 = p_element_2->get_node_index(p_node);
            Node* p_node_1  = p_element_1->get_node((n_index_1+1)%4);
            Node* p_node_3  = p_element_1->get_node((n_index_1+3)%4);
            Node* p_node_4  = p_element_2->get_node((n_index_2+2)%4);
            p_element_1->change_node(n_index_1, p_node_4);
            /////////////////////////////////////
            m_node_element[p_node_1].erase(find(m_node_element[p_node_1].begin(), m_node_element[p_node_1].end(), p_element_2));
            m_node_element[p_node_3].erase(find(m_node_element[p_node_3].begin(), m_node_element[p_node_3].end(), p_element_2));
            *find(m_node_element[p_node_4].begin(), m_node_element[p_node_4].end(), p_element_2) = p_element_1;
            m_node_element.erase(p_node);
            /////////////////////////////////////
            _v_p_element.erase(find(_v_p_element.begin(), _v_p_element.end(), p_element_2));
            delete p_element_2;
            delete p_node;
            _v_p_node.erase(itr_nd--);
        }
    }
    /////////////////////////////////////////////
    bool b_continue_flag = b_element_elimination_flag || b_diagonal_swapping_flag || b_side_elimination_flag;
    if (!b_continue_flag) return false;
    /////////////////////////////////////////////
    set_neighbor_relation(_v_p_element);
    /////////////////////////////////////////////
    while (b_continue_flag)
    {
        b_continue_flag = false;
        /////////////////////////////////////////
        // element elimination
        /////////////////////////////////////////
        //  * - * - *  |  * - * - * 
        //  | a | d |  |  |   |   | 
        //  |   1   |  |  | a | d |  
        //  | /   \ |  |  |   |   |  
        //  2       4  |  2 - * - 4  
        //  | \   / |  |  |   |   |  
        //  |   3   |  |  | b | c |  
        //  | b | c |  |  |   |   |  
        //  * - * - *  |  * - * - *  
        /////////////////////////////////////////
        if (b_element_elimination_flag)
        {
            for (vector<Element*>::iterator itr_el=_v_p_element.begin(); itr_el!=_v_p_element.end(); ++itr_el)
            {
                if (check_interrupt()) throw interrupt_exception();
                /////////////////////////////////
                Element* p_element = *itr_el;
                /////////////////////////////////
                if (p_element->get_node_size() != 4) continue;
                /////////////////////////////////
                for (int i=0; i<2; ++i)
                {
                    Node* p_node_1 = p_element->get_node( i     );
                    Node* p_node_3 = p_element->get_node((i+2)%4);
                    /////////////////////////////
                    if (p_node_1->is_hard_node() || (int)m_node_element[p_node_1].size() != 3) continue;
                    if (p_node_3->is_hard_node() || (int)m_node_element[p_node_3].size() != 3) continue;
                    /////////////////////////////
                    Node* p_node_2 = p_element->get_node((i+1)%4);
                    Node* p_node_4 = p_element->get_node((i+3)%4);
                    if ((int)m_node_element[p_node_2].size() == 3)
                    {
                        if (!p_node_2->is_hard_node())                                                             continue;
                        if (positive_angle_2d(p_node_2->get_uv(), p_node_3->get_uv(), p_node_1->get_uv()) > 100.0) continue;
                    }
                    if ((int)m_node_element[p_node_4].size() == 3)
                    {
                        if (!p_node_4->is_hard_node())                                                             continue;
                        if (positive_angle_2d(p_node_4->get_uv(), p_node_1->get_uv(), p_node_3->get_uv()) > 100.0) continue;
                    }
                    /////////////////////////////
                    Element* p_element_1 = p_element->get_neighbor( i     );
                    Element* p_element_2 = p_element->get_neighbor((i+1)%4);
                    Element* p_element_3 = p_element->get_neighbor((i+2)%4);
                    Element* p_element_4 = p_element->get_neighbor((i+3)%4);
                    if (!p_element_1 || p_element_1->get_node_size() != 4) continue;
                    if (!p_element_2 || p_element_2->get_node_size() != 4) continue;
                    if (!p_element_3 || p_element_3->get_node_size() != 4) continue;
                    if (!p_element_4 || p_element_4->get_node_size() != 4) continue;
                    /////////////////////////////
                    double a_uv[2];
                    mean_2<double>(p_node_1->get_uv(), p_node_3->get_uv(), a_uv, 2);
                    Node*  p_node = new Node(a_uv);
                    _v_p_node.push_back(p_node);
                    /////////////////////////////
                    p_element_1->change_node    (p_element_1->get_node_index(p_node_1),      p_node);
                    p_element_1->change_neighbor(p_element_1->get_neighbor_index(p_element), p_element_2);
                    p_element_2->change_node    (p_element_2->get_node_index(p_node_3),      p_node);
                    p_element_2->change_neighbor(p_element_2->get_neighbor_index(p_element), p_element_1);
                    p_element_3->change_node    (p_element_3->get_node_index(p_node_3),      p_node);
                    p_element_3->change_neighbor(p_element_3->get_neighbor_index(p_element), p_element_4);
                    p_element_4->change_node    (p_element_4->get_node_index(p_node_1),      p_node);
                    p_element_4->change_neighbor(p_element_4->get_neighbor_index(p_element), p_element_3);
                    /////////////////////////////
                    m_node_element[p_node_2].erase(find(m_node_element[p_node_2].begin(), m_node_element[p_node_2].end(), p_element));
                    m_node_element[p_node_4].erase(find(m_node_element[p_node_4].begin(), m_node_element[p_node_4].end(), p_element));
                    m_node_element.erase(p_node_1);
                    m_node_element.erase(p_node_3);
                    m_node_element[p_node  ].push_back(p_element_1);
                    m_node_element[p_node  ].push_back(p_element_2);
                    m_node_element[p_node  ].push_back(p_element_3);
                    m_node_element[p_node  ].push_back(p_element_4);
                    /////////////////////////////
                    _v_p_node.erase(find(_v_p_node.begin(), _v_p_node.end(), p_node_1));
                    _v_p_node.erase(find(_v_p_node.begin(), _v_p_node.end(), p_node_3));
                    delete p_node_1;
                    delete p_node_3;
                    delete p_element;
                    _v_p_element.erase(itr_el--);
                    /////////////////////////////
                    b_continue_flag = true;
                    b_return_flag   = true;
                    break;
                }
            }
        }
        /////////////////////////////////////////
        // diagonal swapping
        /////////////////////////////////////////
        //      2      |[1]   2      |[2]   2     
        //    / | \    |    /   \    |    /   \   
        //  5   |   4  |  5  a  / 4  |  5 \  b  4 
        //  | a | b |  |  |   /   |  |  |   \   | 
        //  3   |   6  |  3 /  b  6  |  3  a  \ 6 
        //    \ | /    |    \   /    |    \   /   
        //      1      |      1      |      1     
        /////////////////////////////////////////
        if (b_diagonal_swapping_flag)
        {
            map<Element*, bool> m_passed_element_flag;
            vector<Element*>::iterator itr_el;
            for (itr_el=_v_p_element.begin(); itr_el!=_v_p_element.end(); ++itr_el) m_passed_element_flag[*itr_el] = false;
            /////////////////////////////////////
            for (itr_el=_v_p_element.begin(); itr_el!=_v_p_element.end(); ++itr_el)
            {
                if (check_interrupt()) throw interrupt_exception();
                /////////////////////////////////
                Element* p_element_1               = *itr_el;
                m_passed_element_flag[p_element_1] = true;
                /////////////////////////////////
                if (p_element_1->get_node_size() != 4) continue;
                /////////////////////////////////
                for (int i=0; i<4; ++i)
                {
                    Node* p_node_1     = p_element_1->get_node( i     );
                    Node* p_node_2     = p_element_1->get_node((i+1)%4);
                    int   n_valence_12 = (int)m_node_element[p_node_1].size() + (int)m_node_element[p_node_2].size();
                    if (p_node_1->is_hard_node() || p_node_2->is_hard_node() || n_valence_12 < 9) continue;
                    /////////////////////////////
                    Element* p_element_2 = p_element_1->get_neighbor(i);
                    if (!p_element_2 || p_element_2->get_node_size() != 4 || m_passed_element_flag[p_element_2]) continue;
                    /////////////////////////////
                    int   n_index  = p_element_2->get_node_index(p_node_1);
                    Node* p_node_3 = p_element_1->get_node((i      +3)%4);
                    Node* p_node_4 = p_element_2->get_node((n_index+2)%4);
                    Node* p_node_5 = p_element_1->get_node((i      +2)%4);
                    Node* p_node_6 = p_element_2->get_node((n_index+1)%4);
                    /////////////////////////////
                    if (p_node_3->is_hard_node() || p_node_4->is_hard_node() || p_node_5->is_hard_node() || p_node_6->is_hard_node()) continue;
                    /////////////////////////////
                    int n_valence_34 = (int)m_node_element[p_node_3].size() + (int)m_node_element[p_node_4].size();
                    int n_valence_56 = (int)m_node_element[p_node_5].size() + (int)m_node_element[p_node_6].size();
                    /////////////////////////////
                    if ((n_valence_12-n_valence_34 >= n_valence_12-n_valence_56) && (n_valence_12-n_valence_56 >= 3))
                    {
                        p_element_1->change_node( i,            p_node_4);
                        p_element_2->change_node((n_index+3)%4, p_node_3);
                        /////////////////////////
                        Element* p_neighbor_31 = p_element_1->get_neighbor((i      +3)%4);
                        Element* p_neighbor_42 = p_element_2->get_neighbor((n_index+2)%4);
                        p_element_1->change_neighbor( i,            p_neighbor_42);
                        p_element_1->change_neighbor((i      +3)%4, p_element_2);
                        p_element_2->change_neighbor((n_index+2)%4, p_element_1);
                        p_element_2->change_neighbor((n_index+3)%4, p_neighbor_31);
                        if (p_neighbor_31) p_neighbor_31->change_neighbor(p_neighbor_31->get_neighbor_index(p_element_1), p_element_2);
                        if (p_neighbor_42) p_neighbor_42->change_neighbor(p_neighbor_42->get_neighbor_index(p_element_2), p_element_1);
                        /////////////////////////
                        m_node_element[p_node_1].erase(find(m_node_element[p_node_1].begin(), m_node_element[p_node_1].end(), p_element_1));
                        m_node_element[p_node_2].erase(find(m_node_element[p_node_2].begin(), m_node_element[p_node_2].end(), p_element_2));
                        m_node_element[p_node_3].push_back(p_element_2);
                        m_node_element[p_node_4].push_back(p_element_1);
                        /////////////////////////
                        m_passed_element_flag[p_element_1] = m_passed_element_flag[p_element_2] = false;
                        /////////////////////////
                        b_continue_flag = true;
                        b_return_flag   = true;
                        break;
                    }
                    /////////////////////////////
                    if ((n_valence_12-n_valence_56 > n_valence_12-n_valence_34) && (n_valence_12-n_valence_34 >= 3))
                    {
                        p_element_1->change_node((i+1)%4,  p_node_6);
                        p_element_2->change_node( n_index, p_node_5);
                        /////////////////////////
                        Element* p_neighbor_25 = p_element_1->get_neighbor((i+1)%4);
                        Element* p_neighbor_16 = p_element_2->get_neighbor( n_index);
                        p_element_1->change_neighbor( i,            p_neighbor_16);
                        p_element_1->change_neighbor((i      +1)%4, p_element_2);
                        p_element_2->change_neighbor( n_index,      p_element_1);
                        p_element_2->change_neighbor((n_index+3)%4, p_neighbor_25);
                        if (p_neighbor_25) p_neighbor_25->change_neighbor(p_neighbor_25->get_neighbor_index(p_element_1), p_element_2);
                        if (p_neighbor_16) p_neighbor_16->change_neighbor(p_neighbor_16->get_neighbor_index(p_element_2), p_element_1);
                        /////////////////////////
                        m_node_element[p_node_1].erase(find(m_node_element[p_node_1].begin(), m_node_element[p_node_1].end(), p_element_2));
                        m_node_element[p_node_2].erase(find(m_node_element[p_node_2].begin(), m_node_element[p_node_2].end(), p_element_1));
                        m_node_element[p_node_5].push_back(p_element_2);
                        m_node_element[p_node_6].push_back(p_element_1);
                        /////////////////////////
                        m_passed_element_flag[p_element_1] = m_passed_element_flag[p_element_2] = false;
                        /////////////////////////
                        b_continue_flag = true;
                        b_return_flag   = true;
                        break;
                    }
                }
            }
        }
        /////////////////////////////////////////
        // side elimination
        /////////////////////////////////////////
        //     1 ----- 2     |[1] 1 ---- 2     |[2]  1 ---- 2  
        //   /  \  a  /  \   |   /  \      \   |   /      /  \ 
        //  * c  5 - 6  d *  |  *  c  \  d  *  |  *  c  /  d  *
        //   \  /  b  \  /   |   \      \  /   |   \  /      / 
        //     3 ----- 4     |    3 ---- 4     |     3 ---- 4  
        /////////////////////////////////////////
        if (b_side_elimination_flag)
        {
            map<Element*, bool> m_obsolete_element_flag;
            map<Element*, bool> m_passed_element_flag;
            vector<Element*>::iterator itr_el;
            for (itr_el=_v_p_element.begin(); itr_el!=_v_p_element.end(); ++itr_el)
            {
                m_passed_element_flag  [*itr_el] = false;
                m_obsolete_element_flag[*itr_el] = false;
            }
            /////////////////////////////////////
            for (itr_el=_v_p_element.begin(); itr_el!=_v_p_element.end(); ++itr_el)
            {
                if (check_interrupt()) throw interrupt_exception();
                /////////////////////////////////
                Element* p_element_1 = *itr_el;
                /////////////////////////////////
                if (m_obsolete_element_flag[p_element_1]) continue;
                /////////////////////////////////////
                m_passed_element_flag[p_element_1] = true;
                /////////////////////////////////
                if (p_element_1->get_node_size() != 4) continue;
                /////////////////////////////////
                for (int i=0; i<4; ++i)
                {
                    Node* p_node_5 = p_element_1->get_node( i     );
                    Node* p_node_6 = p_element_1->get_node((i+1)%4);
                    if (p_node_5->is_hard_node() || (int)m_node_element[p_node_5].size() != 3) continue;
                    if (p_node_6->is_hard_node() || (int)m_node_element[p_node_6].size() != 3) continue;
                    /////////////////////////////
                    Element* p_element_2 = p_element_1->get_neighbor(i);
                    if (!p_element_2 || p_element_2->get_node_size() != 4 || m_obsolete_element_flag[p_element_2] || m_passed_element_flag[p_element_2]) continue;
                    /////////////////////////////
                    int   n_index  = p_element_2->get_node_index(p_node_5);
                    Node* p_node_1 = p_element_1->get_node((i      +3)%4);
                    Node* p_node_2 = p_element_1->get_node((i      +2)%4);
                    Node* p_node_3 = p_element_2->get_node((n_index+1)%4);
                    Node* p_node_4 = p_element_2->get_node((n_index+2)%4);
                    /////////////////////////////
                    if (p_node_1->is_hard_node() || p_node_2->is_hard_node() || p_node_3->is_hard_node() || p_node_4->is_hard_node()) continue;
                    /////////////////////////////
                    Element* p_element_3   = p_element_1->get_neighbor((i+3)%4);
                    Element* p_element_4   = p_element_1->get_neighbor((i+1)%4);
                    if (!p_element_3 || p_element_3->get_node_size() != 4 || m_obsolete_element_flag[p_element_3]) continue;
                    if (!p_element_4 || p_element_4->get_node_size() != 4 || m_obsolete_element_flag[p_element_4]) continue;
                    Element* p_neighbor_21 = p_element_1->get_neighbor((i      +2)%4);
                    Element* p_neighbor_34 = p_element_2->get_neighbor((n_index+1)%4);
                    if (m_obsolete_element_flag[p_neighbor_21]) p_neighbor_21 = 0;
                    if (m_obsolete_element_flag[p_neighbor_34]) p_neighbor_34 = 0;
                    /////////////////////////////
                    if ((int)m_node_element[p_node_1].size()+(int)m_node_element[p_node_4].size() <= (int)m_node_element[p_node_2].size()+(int)m_node_element[p_node_3].size())
                    {
                        int n_index = p_element_3->get_node_index(p_node_3);
                        p_element_3->change_node    ((n_index+1)%4, p_node_4);
                        p_element_3->change_neighbor( n_index,      p_neighbor_34);
                        p_element_3->change_neighbor((n_index+1)%4, p_element_4);
                        if (p_neighbor_34) p_neighbor_34->change_neighbor(p_neighbor_34->get_neighbor_index(p_element_2), p_element_3);
                        /////////////////////////
                        n_index = p_element_4->get_node_index(p_node_2);
                        p_element_4->change_node    ((n_index+1)%4, p_node_1);
                        p_element_4->change_neighbor( n_index,      p_neighbor_21);
                        p_element_4->change_neighbor((n_index+1)%4, p_element_3);
                        if (p_neighbor_21) p_neighbor_21->change_neighbor(p_neighbor_21->get_neighbor_index(p_element_1), p_element_4);
                        /////////////////////////
                        m_node_element[p_node_1].erase(find(m_node_element[p_node_1].begin(), m_node_element[p_node_1].end(), p_element_1));
                        m_node_element[p_node_2].erase(find(m_node_element[p_node_2].begin(), m_node_element[p_node_2].end(), p_element_1));
                        m_node_element[p_node_3].erase(find(m_node_element[p_node_3].begin(), m_node_element[p_node_3].end(), p_element_2));
                        m_node_element[p_node_4].erase(find(m_node_element[p_node_4].begin(), m_node_element[p_node_4].end(), p_element_2));
                        m_node_element[p_node_1].push_back(p_element_4);
                        m_node_element[p_node_4].push_back(p_element_3);
                    }
                    else
                    {
                        int n_index = p_element_3->get_node_index(p_node_3);
                        p_element_3->change_node    ((n_index+1)%4, p_node_2);
                        p_element_3->change_neighbor( n_index,      p_element_4);
                        p_element_3->change_neighbor((n_index+1)%4, p_neighbor_21);
                        if (p_neighbor_21) p_neighbor_21->change_neighbor(p_neighbor_21->get_neighbor_index(p_element_1), p_element_3);
                        /////////////////////////
                        n_index = p_element_4->get_node_index(p_node_2);
                        p_element_4->change_node    ((n_index+1)%4, p_node_3);
                        p_element_4->change_neighbor( n_index,      p_element_3);
                        p_element_4->change_neighbor((n_index+1)%4, p_neighbor_34);
                        if (p_neighbor_34) p_neighbor_34->change_neighbor(p_neighbor_34->get_neighbor_index(p_element_2), p_element_4);
                        /////////////////////////
                        m_node_element[p_node_1].erase(find(m_node_element[p_node_1].begin(), m_node_element[p_node_1].end(), p_element_1));
                        m_node_element[p_node_2].erase(find(m_node_element[p_node_2].begin(), m_node_element[p_node_2].end(), p_element_1));
                        m_node_element[p_node_3].erase(find(m_node_element[p_node_3].begin(), m_node_element[p_node_3].end(), p_element_2));
                        m_node_element[p_node_4].erase(find(m_node_element[p_node_4].begin(), m_node_element[p_node_4].end(), p_element_2));
                        m_node_element[p_node_2].push_back(p_element_3);
                        m_node_element[p_node_3].push_back(p_element_4);
                    }
                    /////////////////////////////
                    _v_p_node     .erase(find(_v_p_node.begin(), _v_p_node.end(), p_node_5));  
                    _v_p_node     .erase(find(_v_p_node.begin(), _v_p_node.end(), p_node_6));  
                    m_node_element.erase(p_node_5);
                    m_node_element.erase(p_node_6);
                    delete p_node_5;
                    delete p_node_6;
                    /////////////////////////////
                    m_obsolete_element_flag[p_element_1] = m_obsolete_element_flag[p_element_2] = true;
                    m_passed_element_flag  [p_element_1] = m_passed_element_flag  [p_element_2] = true;
                    m_passed_element_flag  [p_element_3] = m_passed_element_flag  [p_element_4] = false;
                    /////////////////////////////
                    b_continue_flag = true; 
                    b_return_flag   = true;
                    break;
                }
            }	
            /////////////////////////////////////
            for (itr_el=_v_p_element.begin(); itr_el!=_v_p_element.end(); ++itr_el)
            {
                Element* p_element = *itr_el;
                if (m_obsolete_element_flag[p_element])
                {
                    delete p_element;
                    _v_p_element.erase(itr_el--);
                }
            }
        }
    }
    /////////////////////////////////////////////
    return b_return_flag;
} // end: Mesher_2D::relax_quadrilateral_mesh()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
bool Mesher_2D::check_input_wire_node_pair(Node* p_node_1, Node* p_node_2) const
{
    for (map< Wire*, vector<Node*> >::const_iterator itr_wn=_m_input_wire_node.begin(); itr_wn!=_m_input_wire_node.end(); ++itr_wn)
    {
        Wire*                p_wire        = itr_wn->first;
        const vector<Node*>& v_p_wire_node = itr_wn->second;
        for (int i=0; i<(int)v_p_wire_node.size()-1; ++i)
        {
            if (v_p_wire_node[i] == p_node_1 && v_p_wire_node[i+1] == p_node_2) return true;
            if (v_p_wire_node[i] == p_node_2 && v_p_wire_node[i+1] == p_node_1) return true;
        }
    }
    return false;
} // end: Mesher_2D::check_input_wire_node_pair()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Mesher_2D::set_neighbor_relation(vector<Element*> v_p_element) const
{
    map< Node*,    vector<Element*> > m_node_element;
    map< Element*, vector<Element*> > m_neighbor;
    vector<Element*>::const_iterator itr_el;
    for (itr_el=v_p_element.begin(); itr_el!=v_p_element.end(); ++itr_el)
    {
        Element*  p_element   = *itr_el;
        const int n_node_size = p_element->get_node_size();
        for (int i=0; i<n_node_size; ++i) m_node_element[p_element->get_node(i)].push_back(p_element);
        m_neighbor[p_element] = vector<Element*>(n_node_size, static_cast<Element*>(0));
    }
    /////////////////////////////////////////////
    for (vector<Element*>::const_iterator itr_el_1=v_p_element.begin(); itr_el_1!=v_p_element.end(); ++itr_el_1)
    {
        Element*  p_element_1   = *itr_el_1;
        const int n_node_size_1 = p_element_1->get_node_size();
        for (int i=0; i<n_node_size_1; ++i)
        {
            if (m_neighbor[p_element_1][i]) continue;
            /////////////////////////////////////
            Node* p_node_1 = p_element_1->get_node( i                 );
            Node* p_node_2 = p_element_1->get_node((i+1)%n_node_size_1);
            /////////////////////////////////////
            if (check_input_wire_node_pair(p_node_1, p_node_2)) continue;
            /////////////////////////////////////
            const vector<Element*>& v_p_element_1 = m_node_element[p_node_1];
            const vector<Element*>& v_p_element_2 = m_node_element[p_node_2];
            for (vector<Element*>::const_iterator itr_el_2=v_p_element_1.begin(); itr_el_2!=v_p_element_1.end(); ++itr_el_2)
            {
                Element* p_element_2 = *itr_el_2;
                if (p_element_2 == p_element_1) continue;
                if (find(v_p_element_2.begin(), v_p_element_2.end(), p_element_2) != v_p_element_2.end())
                {
                    m_neighbor[p_element_1][i                                    ] = p_element_2;
                    m_neighbor[p_element_2][p_element_2->get_node_index(p_node_2)] = p_element_1;
                    break;
                }
            }
        }
    }
    for (itr_el=v_p_element.begin(); itr_el!=v_p_element.end(); ++itr_el)
    {
        Element* p_element = *itr_el;
        p_element->set_neighbor(m_neighbor[p_element]);
    }
} // end: Mesher_2D::set_neighbor_relation()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Mesher_2D::update_neighbor_relation(set<Element*> v_p_element) const
{
    for (set<Element*>::const_iterator itr_el=v_p_element.begin(); itr_el!=v_p_element.end(); ++itr_el)
    {
        Element* p_element = *itr_el;
        if (p_element->get_neighbor_size() != p_element->get_node_size()) p_element->reset_neighbor();
    }
    /////////////////////////////////////////////
    for (set<Element*>::const_iterator itr_el_1=v_p_element.begin(); itr_el_1!=v_p_element.end(); ++itr_el_1)
    {
        Element*  p_element_1   = *itr_el_1;
        const int n_node_size_1 = p_element_1->get_node_size();
        for (int i=0; i<n_node_size_1; ++i)
        {
            Node* p_node_1 = p_element_1->get_node( i                 );
            Node* p_node_2 = p_element_1->get_node((i+1)%n_node_size_1);
            /////////////////////////////////////
            if (check_input_wire_node_pair(p_node_1, p_node_2)) continue;
            /////////////////////////////////////
            for (set<Element*>::const_iterator itr_el_2=v_p_element.begin(); itr_el_2!=v_p_element.end(); ++itr_el_2)
            {
                Element*  p_element_2   = *itr_el_2;
                if (p_element_2 == p_element_1) continue;
                const int n_node_size_2 = p_element_2->get_node_size();
                for (int j=0; j<n_node_size_2; ++j)
                {
                    Node* p_node_3 = p_element_2->get_node( j                 );
                    Node* p_node_4 = p_element_2->get_node((j+1)%n_node_size_2);
                    if (p_node_1 == p_node_4 && p_node_2 == p_node_3)
                    {
                        p_element_1->change_neighbor(i, p_element_2);
                        p_element_2->change_neighbor(j, p_element_1);
                        break;
                    }
                }
            }
        }
    }
} // end: Mesher_2D::update_neighbor_relation()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
Node* Mesher_2D::get_unique_node(const double* a_uv, bool b_hard_node_flag)
{
    using namespace mesher_math_lib;
    /////////////////////////////////////////////
    for (vector<Node*>::const_iterator itr_nd=_v_p_node.begin(); itr_nd!=_v_p_node.end(); ++itr_nd)
    {
        Node* p_node = *itr_nd;
        if (distance_2d(p_node->get_uv(), a_uv) < _d_tolerance)
        {
            if (b_hard_node_flag) p_node->set_hard_node();
            return p_node;
        }
    }
    /////////////////////////////////////////////
    Node* p_node = new Node(a_uv, 0, b_hard_node_flag);
    _v_p_node.push_back(p_node);
    return p_node;
} // end: Mesher_2D::get_unique_node()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Mesher_2D::seed_node(vector<const double*> v_a_uv, const vector<double>& v_xi, bool b_hard_node_flag, vector<Node*>& v_p_node, double d_length)
{
    using namespace mesher_math_lib;
    /////////////////////////////////////////////
    v_p_node.clear();
    /////////////////////////////////////////////
    if (fabs(d_length) < ZERO)
    {
        for (int i=0; i<(int)v_a_uv.size()-1; ++i)
        {
            for (int j=i+1; j<(int)v_a_uv.size(); ++j) d_length += distance_2d(v_a_uv[i], v_a_uv[j]);
        }
    }
    /////////////////////////////////////////////
    for (int i=0; i<(int)v_xi.size(); i++)
    {
        double d_residual = d_length * v_xi[i];
        for (int j=0; j<(int)v_a_uv.size()-1; j++)
        {
            const double* a_uv_1     = v_a_uv[j  ];
            const double* a_uv_2     = v_a_uv[j+1];
            double        d_interval = distance_2d(a_uv_1, a_uv_2);
            if (d_interval > d_residual || fabs(d_interval - d_residual) < _d_tolerance)
            {
                double a_uv[2];
                a_uv[0] = a_uv_1[0] + d_residual / d_interval * (a_uv_2[0] - a_uv_1[0]);
                a_uv[1] = a_uv_1[1] + d_residual / d_interval * (a_uv_2[1] - a_uv_1[1]);
                v_p_node.push_back(get_unique_node(a_uv, b_hard_node_flag));
                break;
            }
            d_residual -= d_interval;
        }
    }
} // end: Mesher_2D::seed_node()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Mesher_2D::estimate_seed_interval(vector<const double*> v_a_uv, double d_seed_size_1, double d_seed_size_2, int n_division, vector<double>& v_interval, double d_length) const
{
    using namespace mesher_math_lib;
    /////////////////////////////////////////////
    if (fabs(d_seed_size_1 - d_seed_size_2) < _d_tolerance)
    {
        v_interval.resize(n_division);
        fill(v_interval.begin(), v_interval.end(), d_seed_size_1);
        return;
    }
    /////////////////////////////////////////////
    if (n_division == 1)
    {
        if (fabs(d_length) < ZERO)
        {
            for (int i=0; i<(int)v_a_uv.size()-1; ++i)
            {
                for (int j=i+1; j<(int)v_a_uv.size(); ++j) d_length += distance_2d(v_a_uv[i], v_a_uv[j]);
            }
        }
        v_interval.push_back(d_length);
        return;
    }
    /////////////////////////////////////////////
    // interval = a i + b
    const double d_a = (d_seed_size_2 - d_seed_size_1) / (n_division - 1);
    const double d_b = d_seed_size_1 - d_a;
    v_interval.resize(n_division);
    for (int i=0; i<n_division; ++i) v_interval[i] = d_a * (i+1) + d_b;
} // end: Mesher_2D::estimate_seed_interval()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
double Mesher_2D::estimate_ideal_division(vector<const double*> v_a_uv, double d_seed_size_1, double d_seed_size_2, double d_length) const
{
    using namespace mesher_math_lib;
    /////////////////////////////////////////////
    if (fabs(d_length) < ZERO)
    {
        for (int i=0; i<(int)v_a_uv.size()-1; ++i)
        {
            for (int j=i+1; j<(int)v_a_uv.size(); ++j) d_length += distance_2d(v_a_uv[i], v_a_uv[j]);
        }
    }
    /////////////////////////////////////////////
    return fabs(d_seed_size_1 - d_seed_size_2) < _d_tolerance
        ? d_length / d_seed_size_1
        : d_length * log(d_seed_size_2 / d_seed_size_1) / (d_seed_size_2 - d_seed_size_1);
} // end: Mesher_2D::estimate_ideal_division()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
bool Mesher_2D::check_offset_node(Node* p_node_1, Wire* p_mother_wire, const Offset& offset_1, const map<Node*, Offset>& m_offset, bool b_self_check) const
{
    using namespace mesher_math_lib;
    /////////////////////////////////////////////
    double d_nearness = offset_1.get_interval() * 0.75;
    for (int i=0; i<offset_1.get_size(); ++i)
    {
        const double* a_uv = offset_1.get_uv(i);
        /////////////////////////////////////////
        vector<Wire*>::const_iterator itr_wr;
        for (itr_wr=_v_p_wire.begin(); itr_wr!=_v_p_wire.end(); ++itr_wr)
        {
            Wire*                p_wire        = *itr_wr;
            const vector<Node*>& v_p_wire_node = p_wire->get_node();
            /////////////////////////////////////
            if (p_wire->get_type() != Wire::WT_OPEN)
            {
                double d_angle = 0.0;
                for (int j=0; j<(int)v_p_wire_node.size()-1; ++j)
                {
                    const double* a_uv_1 = v_p_wire_node[j  ]->get_uv();
                    const double* a_uv_2 = v_p_wire_node[j+1]->get_uv();
                    d_angle             += angle_2d(a_uv, a_uv_1, a_uv_2);
                }
                d_angle = fabs(d_angle);
                if ((p_wire->get_type() == Wire::WT_OUTER && d_angle < 10.0) || (p_wire->get_type() == Wire::WT_INNER && d_angle > 350.0)) return false;
            }
            /////////////////////////////////////
            if (p_wire == p_mother_wire) continue;
            /////////////////////////////////////
            for (int j=0; j<(int)v_p_wire_node.size(); ++j)
            {
                const double* a_uv_2 = v_p_wire_node[j]->get_uv();
                if (distance_2d(a_uv, a_uv_2) < d_nearness) return false;
            }
            /////////////////////////////////////
            for (int j=0; j<(int)v_p_wire_node.size()-1; ++j)
            {
                const double* a_uv_1 = v_p_wire_node[j  ]->get_uv();
                const double* a_uv_2 = v_p_wire_node[j+1]->get_uv();
                /////////////////////////////////
                double d_determinant = -pow(a_uv_1[0]-a_uv_2[0], 2.0) - pow(a_uv_1[1]-a_uv_2[1], 2.0);
                if (fabs(d_determinant) > ZERO)
                {
                    double d_xi = ((a_uv[0]-a_uv_1[0]) * (a_uv_1[0]-a_uv_2[0]) + (a_uv[1]-a_uv_1[1]) * (a_uv_1[1]-a_uv_2[1])) / d_determinant;
                    if (d_xi > -ZERO && d_xi < 1.0 + ZERO)
                    {
                        double a_uv_3[2];
                        for (int k=0; k<2; ++k) a_uv_3[k] = a_uv_1[k] + d_xi * (a_uv_2[k] - a_uv_1[k]);
                        if (distance_2d(a_uv, a_uv_3) < d_nearness) return false;
                    }
                }
            }
        }
        /////////////////////////////////////////
        for (itr_wr=_v_p_wire.begin(); itr_wr!=_v_p_wire.end(); ++itr_wr)
        {
            Wire* p_wire = *itr_wr;
            /////////////////////////////////////
            if (p_wire->get_type() != Wire::WT_OPEN) continue;
            /////////////////////////////////////
            const vector<Node*>& v_p_wire_node = p_wire->get_node();
            for (int j=0; j<(int)v_p_wire_node.size(); ++j)
            {
                const double* a_uv_2 = v_p_wire_node[j]->get_uv();
                if (distance_2d(a_uv, a_uv_2) < d_nearness) return false;
            }
            /////////////////////////////////////
            for (int j=0; j<(int)v_p_wire_node.size()-1; ++j)
            {
                const double* a_uv_1 = v_p_wire_node[j  ]->get_uv();
                const double* a_uv_2 = v_p_wire_node[j+1]->get_uv();
                /////////////////////////////////
                double d_determinant = -pow(a_uv_1[0]-a_uv_2[0], 2.0) - pow(a_uv_1[1]-a_uv_2[1], 2.0);
                if (fabs(d_determinant) > ZERO)
                {
                    double d_xi = ((a_uv[0]-a_uv_1[0]) * (a_uv_1[0]-a_uv_2[0]) + (a_uv[1]-a_uv_1[1]) * (a_uv_1[1]-a_uv_2[1])) / d_determinant;
                    if (d_xi > -ZERO && d_xi < 1.0 + ZERO)
                    {
                        double a_uv_3[2];
                        for (int k=0; k<2; ++k) a_uv_3[k] = a_uv_1[k] + d_xi * (a_uv_2[k] - a_uv_1[k]);
                        if (distance_2d(a_uv, a_uv_3) < d_nearness) return false;
                    }
                }
            }
        }
        /////////////////////////////////////////
        for (vector<Node*>::const_iterator itr_nd=_v_p_island_node.begin(); itr_nd!=_v_p_island_node.end(); ++itr_nd)
        {
            Node* p_node = *itr_nd;
            if (distance_2d(p_node->get_uv(), a_uv) < d_nearness) return false;
        }
    }
    /////////////////////////////////////////////
    for (int i=0; i<offset_1.get_size(); ++i)
    {
        const double* a_uv = offset_1.get_uv(i);
        /////////////////////////////////////////
        const vector<Node*>& v_p_wire_node = p_mother_wire->get_node();
        for (int j=0; j<(int)v_p_wire_node.size(); ++j)
        {
            const double* a_uv_2 = v_p_wire_node[j]->get_uv();
            if (distance_2d(a_uv, a_uv_2) < d_nearness) return false;
        }
        /////////////////////////////////////
        for (int j=0; j<(int)v_p_wire_node.size()-1; ++j)
        {
            const double* a_uv_1 = v_p_wire_node[j  ]->get_uv();
            const double* a_uv_2 = v_p_wire_node[j+1]->get_uv();
            /////////////////////////////////
            double d_determinant = -pow(a_uv_1[0]-a_uv_2[0], 2.0) - pow(a_uv_1[1]-a_uv_2[1], 2.0);
            if (fabs(d_determinant) > ZERO) 
            {
                double d_xi = ((a_uv[0]-a_uv_1[0]) * (a_uv_1[0]-a_uv_2[0]) + (a_uv[1]-a_uv_1[1]) * (a_uv_1[1]-a_uv_2[1])) / d_determinant;
                if (d_xi > -ZERO && d_xi < 1.0 + ZERO)
                {
                    double a_uv_3[2];
                    for (int k=0; k<2; ++k) a_uv_3[k] = a_uv_1[k] + d_xi * (a_uv_2[k] - a_uv_1[k]);
                    if (distance_2d(a_uv, a_uv_3) < d_nearness) return false;
                }
            }
        }		
    }
    /////////////////////////////////////////////
    if (!b_self_check) return true;
    /////////////////////////////////////////////
    const int n_wire_node_size = p_mother_wire->get_node_size();
    /////////////////////////////////////////////
    const double* a_uv_1 = p_node_1->get_uv();
    for (int i=0; i<offset_1.get_size(); ++i)
    {
        const double* a_uv_2 = offset_1.get_uv(i);
        for (map<Node*, Offset>::const_iterator itr_mo=m_offset.begin(); itr_mo!=m_offset.end(); ++itr_mo)
        {
            Node* p_node_2 = itr_mo->first;
            /////////////////////////////////////
            if (p_node_2 == p_node_1) continue;
            /////////////////////////////////////
            int n_index_1    = p_mother_wire->get_node_index(p_node_1);
            int n_index_2    = p_mother_wire->get_node_index(p_node_2);
            int n_difference = abs(n_index_1 - n_index_2);
            if (n_difference == 1 || n_difference == n_wire_node_size-2)
            {
                const Offset& offset_2 = itr_mo->second;
                for (int j=0; j<offset_2.get_size(); ++j)
                {
                    const double* a_uv_3 = offset_2.get_uv(j);
                    if (distance_2d(a_uv_2, a_uv_3) < _d_tolerance) return false;
                }
                continue;
            }
            /////////////////////////////////////
            const double* a_uv_3   = p_node_2->get_uv();
            const Offset& offset_2 = itr_mo->second;
            double        d_margin = mean_2<double>(offset_1.get_interval(), offset_2.get_interval()) * 0.75;
            for (int j=0; j<offset_2.get_size(); ++j)
            {
                const double* a_uv_4 = offset_2.get_uv(j);
                if (distance_2d(a_uv_2, a_uv_4) < d_margin) return false;
                /////////////////////////////////
                double d_determinant = -pow(a_uv_3[0]-a_uv_4[0], 2.0) - pow(a_uv_3[1]-a_uv_4[1], 2.0);
                if (fabs(d_determinant) > ZERO)
                {
                    double d_xi = ((a_uv_2[0]-a_uv_3[0]) * (a_uv_3[0]-a_uv_4[0]) + (a_uv_2[1]-a_uv_3[1]) * (a_uv_3[1]-a_uv_4[1])) / d_determinant;
                    if (d_xi > -ZERO && d_xi < 1.0 + ZERO)
                    {
                        double a_uv[2];
                        for (int k=0; k<2; ++k) a_uv[k] = a_uv_3[k] + d_xi * (a_uv_4[k] - a_uv_3[k]);
                        if (distance_2d(a_uv_2, a_uv) < d_margin) return false;
                    }
                }
                d_determinant = -pow(a_uv_1[0]-a_uv_2[0], 2.0) - pow(a_uv_1[1]-a_uv_2[1], 2.0);
                if (fabs(d_determinant) > ZERO)
                {
                    double d_xi = ((a_uv_4[0]-a_uv_1[0]) * (a_uv_1[0]-a_uv_2[0]) + (a_uv_4[1]-a_uv_1[1]) * (a_uv_1[1]-a_uv_2[1])) / d_determinant;
                    if (d_xi > -ZERO && d_xi < 1.0 + ZERO)
                    {
                        double a_uv[2];
                        for (int k=0; k<2; ++k) a_uv[k] = a_uv_1[k] + d_xi * (a_uv_2[k] - a_uv_1[k]);
                        if (distance_2d(a_uv_4, a_uv) < d_margin) return false;
                    }
                }
                /////////////////////////////////
                if (line_cross_2d(a_uv_1, a_uv_2, a_uv_3, a_uv_4)) return false;
            }
        }
    }	
    /////////////////////////////////////////////
    return true;
} // end: Mesher_2D::check_offset_node()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
bool Mesher_2D::check_wire_cross(const double* a_uv_1, const double* a_uv_2) const
{
    using namespace mesher_math_lib;
    /////////////////////////////////////////////
    for (vector<Wire*>::const_iterator itr_wr=_v_p_wire.begin(); itr_wr!=_v_p_wire.end(); ++itr_wr)
    {
        Wire*                p_wire        = *itr_wr;
        const vector<Node*>& v_p_wire_node = p_wire->get_node();
        for (int i=0; i<(int)v_p_wire_node.size()-1; ++i)
        {
            const double* a_uv_3 = (v_p_wire_node[i  ])->get_uv();
            const double* a_uv_4 = (v_p_wire_node[i+1])->get_uv();
            if (line_cross_2d(a_uv_1, a_uv_2, a_uv_3, a_uv_4)) return true;
        }
    }
    /////////////////////////////////////////////
    return false;
}  // end: Mesher_2D::check_wire_cross()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
bool Mesher_2D::check_in_domain_uv(const double* a_uv)
{
    using namespace mesher_math_lib;
    ///////////////////////////////////////////////
    for (vector<Wire*>::const_iterator itr_wr=_v_p_wire.begin(); itr_wr!=_v_p_wire.end(); ++itr_wr)
    {
        Wire* p_wire = *itr_wr;
        /////////////////////////////////////
        if (p_wire->get_type() == Wire::WT_OPEN) continue;
        /////////////////////////////////////
        const vector<Node*>& v_p_wire_node = p_wire->get_node();
        double               d_angle       = 0.0;
        for (int i=0; i<(int)v_p_wire_node.size()-1; ++i)
        {
            const double* a_uv_1 = v_p_wire_node[i  ]->get_uv();
            const double* a_uv_2 = v_p_wire_node[i+1]->get_uv();
            d_angle             += angle_2d(a_uv, a_uv_1, a_uv_2);
        }
        d_angle = fabs(d_angle);
        if ((p_wire->get_type() == Wire::WT_OUTER && d_angle < 10.0) || (p_wire->get_type() == Wire::WT_INNER && d_angle > 350.0)) return false;
    }
    return true;
} // end: Mesher_2D::check_in_domain_uv()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Mesher_2D::smooth_mesh(int n_time) const throw(exception)
{
    try 
    {
        smooth_mesh_winslow(n_time);
        smooth_mesh_laplace(n_time);
    }
    catch (...)
    {
        throw;
    }
} // end: Mesher_2D::smooth_mesh()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Mesher_2D::smooth_mesh_laplace(int n_time) const throw(exception)
{
    using namespace mesher_math_lib;
    /////////////////////////////////////////////
    map<Node*, vector<Element*> > m_node_element;
    map<Node*, set<Node*> >       m_node_neighbor;
    map<Node*, bool>              m_node_skip_flag;
    /////////////////////////////////////////////
    vector<Node*>::const_iterator itr_nd;
    for (itr_nd=_v_p_node.begin(); itr_nd!=_v_p_node.end(); ++itr_nd)
    {
        Node* p_node             = *itr_nd;
        m_node_skip_flag[p_node] = p_node->is_hard_node();
    }
    /////////////////////////////////////////////
    for (vector<Element*>::const_iterator itr_el=_v_p_element.begin(); itr_el!=_v_p_element.end(); ++itr_el)
    {
        Element* p_element = *itr_el;
        for (int i=0; i<p_element->get_node_size(); ++i) m_node_element[p_element->get_node(i)].push_back(p_element);
        /////////////////////////////////////////
        if (p_element->get_node_size() == 4)
        {
            for (int i=0; i<p_element->get_node_size(); ++i)
            {

                Node*  p_node_1 = p_element->get_node( i     );
                Node*  p_node_2 = p_element->get_node((i+1)%4);
                Node*  p_node_3 = p_element->get_node((i+2)%4);
                Node*  p_node_4 = p_element->get_node((i+3)%4);
                double d_angle  = positive_angle_2d(p_node_1->get_uv(), p_node_2->get_uv(), p_node_4->get_uv());
                if (p_node_1->is_hard_node() && d_angle > 165.0)
                {
                    m_node_skip_flag[p_node_1] = true;
                    m_node_skip_flag[p_node_3] = true;
                    break;
                }
            }
        }
    }
    for (itr_nd=_v_p_node.begin(); itr_nd!=_v_p_node.end(); ++itr_nd)
    {
        Node* p_node = *itr_nd;
        /////////////////////////////////////////
        if (m_node_skip_flag[p_node]) continue;
        /////////////////////////////////////////
        const vector<Element*>& v_p_element = m_node_element[p_node];
        if (v_p_element.empty())
        {
            m_node_skip_flag[p_node] = true;
            continue;
        }
        /////////////////////////////////////
        for (vector<Element*>::const_iterator itr_el=v_p_element.begin(); itr_el!=v_p_element.end(); ++itr_el)
        {
            Element*  p_element   = *itr_el;
            int       n_index     = p_element->get_node_index(p_node);
            const int n_node_size = p_element->get_node_size();
            m_node_neighbor[p_node].insert(p_element->get_node((n_index+1            )%n_node_size));
            m_node_neighbor[p_node].insert(p_element->get_node((n_index+n_node_size-1)%n_node_size));
        }
    }
    /////////////////////////////////////////////
    for (int i=0; i<n_time; ++i)
    {
        if (check_interrupt()) throw interrupt_exception();
        /////////////////////////////////////////
        for (vector<Node*>::const_iterator itr_nd_1=_v_p_node.begin(); itr_nd_1!=_v_p_node.end(); ++itr_nd_1)
        {
            Node* p_node = *itr_nd_1;
            /////////////////////////////////////
            if (m_node_skip_flag[p_node]) continue;
            /////////////////////////////////////
            const set<Node*>& s_p_neighbor = m_node_neighbor[p_node];
            if (s_p_neighbor.empty()) continue;
            /////////////////////////////////////
            int    n_hard_neighbor_count = 0;
            double a_uv[2]               = { 0.0, 0.0 };
            for (set<Node*>::const_iterator itr_nd_2=s_p_neighbor.begin(); itr_nd_2!=s_p_neighbor.end(); ++itr_nd_2)
            {
                Node* p_neighbor = *itr_nd_2;
                transform(a_uv, a_uv+2, p_neighbor->get_uv(), a_uv, plus<double>());
                if (p_neighbor->is_hard_node()) ++n_hard_neighbor_count;
            }
            transform(a_uv, a_uv+2, a_uv, std::bind(divides<double>(), _1, (int)s_p_neighbor.size()));
            /////////////////////////////////////
            if (n_hard_neighbor_count > 0 && ((double)n_hard_neighbor_count)/(int)s_p_neighbor.size() >= 0.5)
            {
                const vector<Element*>& v_p_element = m_node_element[p_node];
                vector<double>          v_quality_1((int)v_p_element.size());
                vector<double>          v_quality_2((int)v_p_element.size());
                for (int j=0; j<(int)v_p_element.size(); ++j)
                {
                    Element* p_element = v_p_element[j];
                    /////////////////////////////
                    v_quality_1[j] = p_element->get_distortion_metric();
                    /////////////////////////////
                    int           n_index     = p_element->get_node_index(p_node);
                    const int     n_node_size = p_element->get_node_size();
                    const double* a_uv_2      = (p_element->get_node((n_index+1)%n_node_size))->get_uv();
                    const double* a_uv_3      = (p_element->get_node((n_index+2)%n_node_size))->get_uv();
                    switch (n_node_size)
                    {
                    case 3: v_quality_2[j] = distortion_metric_2d(a_uv, a_uv_2, a_uv_3);
                        break;
                    case 4:
                        {
                            const double* a_uv_4 = (p_element->get_node((n_index+3)%n_node_size))->get_uv();
                            v_quality_2[j]       = distortion_metric_2d(a_uv, a_uv_2, a_uv_3, a_uv_4);
                            break;
                        }
                    }
                }
                double d_minimum_1 = *min_element(v_quality_1.begin(), v_quality_1.end());
                double d_minimum_2 = *min_element(v_quality_2.begin(), v_quality_2.end());
                /////////////////////////////////
                if (d_minimum_2 < d_minimum_1) continue;
            }
            /////////////////////////////////////
            p_node->set_uv(a_uv);
        }
    }
} // end: Mesher_2D::smooth_mesh_laplace()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Mesher_2D::smooth_mesh_winslow(int n_time) const throw(exception)
{
    using namespace mesher_math_lib;
    /////////////////////////////////////////////
    map<Node*, vector<Element*> > m_node_element;
    map<Node*, bool>              m_node_skip_flag;
    /////////////////////////////////////////////
    vector<Node*>::const_iterator itr_nd;
    for (itr_nd=_v_p_node.begin(); itr_nd!=_v_p_node.end(); ++itr_nd)
    {
        Node* p_node             = *itr_nd;
        m_node_skip_flag[p_node] = p_node->is_hard_node();
    }
    /////////////////////////////////////////////
    for (vector<Element*>::const_iterator itr_el=_v_p_element.begin(); itr_el!=_v_p_element.end(); ++itr_el)
    {
        Element* p_element = *itr_el;
        for (int i=0; i<p_element->get_node_size(); ++i) m_node_element[p_element->get_node(i)].push_back(p_element);
        /////////////////////////////////////////
        if (p_element->get_node_size() == 4)
        {
            for (int i=0; i<p_element->get_node_size(); ++i)
            {

                Node*  p_node_1 = p_element->get_node( i     );
                Node*  p_node_2 = p_element->get_node((i+1)%4);
                Node*  p_node_3 = p_element->get_node((i+2)%4);
                Node*  p_node_4 = p_element->get_node((i+3)%4);
                double d_angle  = positive_angle_2d(p_node_1->get_uv(), p_node_2->get_uv(), p_node_4->get_uv());
                if (p_node_1->is_hard_node() && d_angle > 165.0)
                {
                    m_node_skip_flag[p_node_1] = true;
                    m_node_skip_flag[p_node_3] = true;
                    break;
                }
            }
        }
    }
    /////////////////////////////////////////////
    for (itr_nd=_v_p_node.begin(); itr_nd!=_v_p_node.end(); ++itr_nd)
    {
        Node* p_node = *itr_nd;
        if (m_node_element[p_node].empty()) m_node_skip_flag[p_node] = true;
    }
    /////////////////////////////////////////////
    for (int i=0; i<n_time; ++i)
    {
        if (check_interrupt()) throw interrupt_exception();
        /////////////////////////////////////////
        for (vector<Node*>::const_iterator itr_nd_1=_v_p_node.begin(); itr_nd_1!=_v_p_node.end(); ++itr_nd_1)
        {
            Node*                   p_node      = *itr_nd_1;
            const vector<Element*>& v_p_element = m_node_element[p_node];
            /////////////////////////////////////
            if (m_node_skip_flag[p_node]) continue;
            /////////////////////////////////////
            set<Node*> s_p_neighbor;
            double     a_uv[2]      = { 0.0, 0.0 };
            double     d_total_area = 0.0;
            for (vector<Element*>::const_iterator itr_el=v_p_element.begin(); itr_el!=v_p_element.end(); ++itr_el)
            {
                Element*      p_element   = *itr_el;
                const double* a_uv_1      = (p_element->get_node(0))->get_uv();
                const double* a_uv_2      = (p_element->get_node(1))->get_uv();
                const double* a_uv_3      = (p_element->get_node(2))->get_uv();
                int           n_index     = p_element->get_node_index(p_node);
                const int     n_node_size = p_element->get_node_size();
                switch (n_node_size)
                {
                case 3:
                    {
                        double d_area = triangle_area_2d(a_uv_1, a_uv_2, a_uv_3);
                        for (int j=0; j<2; ++j) a_uv[j] += (a_uv_1[j] + a_uv_2[j] + a_uv_3[j]) / 3.0 * d_area;
                        d_total_area += d_area;
                        break;
                    }
                case 4:
                    {
                        const double* a_uv_4 = (p_element->get_node(3))->get_uv();
                        double        d_area = triangle_area_2d(a_uv_1, a_uv_2, a_uv_3) + triangle_area_2d(a_uv_1, a_uv_3, a_uv_4);
                        for (int j=0; j<2; ++j) a_uv[j] += (a_uv_1[j] + a_uv_2[j] + a_uv_3[j] + a_uv_4[j]) / 4.0 * d_area;
                        d_total_area        += d_area;
                        break;
                    }
                }
                /////////////////////////////////
                s_p_neighbor.insert(p_element->get_node((n_index+1            )%n_node_size));
                s_p_neighbor.insert(p_element->get_node((n_index+n_node_size-1)%n_node_size));
            }
            transform(a_uv, a_uv+2, a_uv, std::bind(divides<double>(), _1, d_total_area));
            /////////////////////////////////////
            int n_hard_neighbor_count = 0;
            for (set<Node*>::const_iterator itr_nd_2=s_p_neighbor.begin(); itr_nd_2!=s_p_neighbor.end(); ++itr_nd_2)
            {
                Node* p_neighbor = *itr_nd_2;
                if (p_neighbor->is_hard_node()) ++n_hard_neighbor_count;
            }
            /////////////////////////////////////
            if (n_hard_neighbor_count > 0 && ((double)n_hard_neighbor_count)/(int)s_p_neighbor.size() >= 0.5)
            {
                const vector<Element*>& v_p_element = m_node_element[p_node];
                vector<double>          v_quality_1((int)v_p_element.size());
                vector<double>          v_quality_2((int)v_p_element.size());
                for (int j=0; j<(int)v_p_element.size(); ++j)
                {
                    Element* p_element = v_p_element[j];
                    /////////////////////////////////
                    v_quality_1[j] = p_element->get_distortion_metric();
                    /////////////////////////////////
                    int           n_index     = p_element->get_node_index(p_node);
                    const int     n_node_size = p_element->get_node_size();
                    const double* a_uv_2      = (p_element->get_node((n_index+1)%n_node_size))->get_uv();
                    const double* a_uv_3      = (p_element->get_node((n_index+2)%n_node_size))->get_uv();
                    switch (n_node_size)
                    {
                    case 3: v_quality_2[j] = distortion_metric_2d(a_uv, a_uv_2, a_uv_3);
                        break;
                    case 4:
                        {
                            const double* a_uv_4 = (p_element->get_node((n_index+3)%n_node_size))->get_uv();
                            v_quality_2[j]       = distortion_metric_2d(a_uv, a_uv_2, a_uv_3, a_uv_4);
                            break;
                        }
                    }
                }
                double d_minimum_1 = *min_element(v_quality_1.begin(), v_quality_1.end());
                double d_minimum_2 = *min_element(v_quality_2.begin(), v_quality_2.end());
                /////////////////////////////////
                if (d_minimum_2 < d_minimum_1) continue;
            }
            /////////////////////////////////////
            p_node->set_uv(a_uv);
        }
    }
} // end: Mesher_2D::smooth_mesh_winslow()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Mesher_2D::smooth_mesh_angle(int n_time) const throw(exception)
{
    using namespace mesher_math_lib;
    /////////////////////////////////////////////
    map< Node*, vector<Element*> > m_node_element;
    /////////////////////////////////////////////
    for (vector<Element*>::const_iterator itr_el=_v_p_element.begin(); itr_el!=_v_p_element.end(); ++itr_el)
    {
        Element* p_element = *itr_el;
        for (int i=0; i<p_element->get_node_size(); ++i) m_node_element[p_element->get_node(i)].push_back(p_element);
    }
    /////////////////////////////////////////////
    for (vector<Node*>::const_iterator itr_nd=_v_p_node.begin(); itr_nd!=_v_p_node.end(); ++itr_nd)
    {
        Node*             p_node      = *itr_nd;
        vector<Element*>& v_p_element = m_node_element[p_node];
        for (int i=0; i<(int)v_p_element.size(); ++i)
        {
            if (v_p_element[i]->get_node_size() != 4)
            {
                v_p_element.clear();
                break;
            }
        }
    }
    /////////////////////////////////////////////
    for (int i=0; i<n_time; ++i)
    {
        if (check_interrupt()) throw interrupt_exception();
        /////////////////////////////////////////
        for (vector<Node*>::const_iterator itr_nd=_v_p_node.begin(); itr_nd!=_v_p_node.end(); ++itr_nd)
        {
            Node*                   p_node      = *itr_nd;
            const vector<Element*>& v_p_element = m_node_element[p_node];
            /////////////////////////////////////
            if (p_node->is_hard_node() || v_p_element.empty()) continue;
            /////////////////////////////////////
            const double* a_uv_1 = p_node->get_uv();
            /////////////////////////////////////
            double a_uv[2]         = { 0.0, 0.0 };
            double d_r             = 0.0;
            bool   b_continue_flag = true;
            for (vector<Element*>::const_iterator itr_el=v_p_element.begin(); itr_el!=v_p_element.end(); ++itr_el)
            {
                Element*      p_element   = *itr_el;
                const int     n_node_size = p_element->get_node_size();
                int           n_index     = p_element->get_node_index(p_node);
                const double* a_uv_2      = (p_element->get_node((n_index+1)%4))->get_uv();
                const double* a_uv_3      = (p_element->get_node((n_index+2)%4))->get_uv();
                const double* a_uv_4      = (p_element->get_node((n_index+3)%4))->get_uv();
                /////////////////////////////////
                double a_vector_13[2];
                transform(a_uv_3, a_uv_3+2, a_uv_1, a_vector_13, minus<double>());
                double d_norm  = norm_2d(a_vector_13);
                d_r           += d_norm;
                transform(a_vector_13, a_vector_13+2, a_vector_13, std::bind(divides<double>(), _1, d_norm));
                double d_angle = positive_angle_2d(a_uv_1, a_uv_2, a_uv_4);
                double d_delta = 0.0;
                if (fabs(90.0 - d_angle) > 45.0)
                {
                    d_delta         = 90.0 - d_angle;
                    b_continue_flag = false;
                }
                transform(a_vector_13, a_vector_13+2, a_vector_13, std::bind(multiplies<double>(), _1, d_delta));
                transform(a_uv, a_uv+2, a_vector_13, a_uv, plus<double>());
            }
            if (b_continue_flag) continue;
            /////////////////////////////////////
            d_r /= (int)v_p_element.size();
            transform(a_uv, a_uv+2, a_uv, std::bind(multiplies<double>(), _1, d_r/360.0));
            transform(a_uv, a_uv+2, a_uv_1, a_uv, plus<double>());
            /////////////////////////////////////
            vector<double> v_quality((int)v_p_element.size());
            for (int j=0; j<(int)v_p_element.size(); ++j)
            {
                Element*      p_element   = v_p_element[j];
                int           n_index     = p_element->get_node_index(p_node);
                const int     n_node_size = p_element->get_node_size();
                const double* a_uv_2      = (p_element->get_node((n_index+1)%n_node_size))->get_uv();
                const double* a_uv_3      = (p_element->get_node((n_index+2)%n_node_size))->get_uv();
                const double* a_uv_4      = (p_element->get_node((n_index+3)%n_node_size))->get_uv();
                v_quality[j]              = distortion_metric_2d(a_uv, a_uv_2, a_uv_3, a_uv_4);
            }
            /////////////////////////////////////
            if (*min_element(v_quality.begin(), v_quality.end()) >= 0.36) p_node->set_uv(a_uv);
        }
    }
} // end: Mesher_2D::smooth_mesh_angle()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Mesher_2D::dump_mesh_information() const
{
    using namespace mesher_math_lib;
    /////////////////////////////////////////////
    cout << _T("* N(node)   : ") << (int)_v_p_node   .size() << endl;
    cout << _T("* N(element): ") << (int)_v_p_element.size() << endl;
    /////////////////////////////////////////////
    double d_mean_distortion = 0.0;
    double d_max_distortion  = 0.0;
    double d_min_distortion  = 1.0;
    for (vector<Element*>::const_iterator itr_el=_v_p_element.begin(); itr_el!=_v_p_element.end(); ++itr_el)
    {
        Element* p_element    = *itr_el;
        double   d_distortion = p_element->get_distortion_metric();
        d_mean_distortion    += d_distortion;
        d_max_distortion      = __max(d_max_distortion, d_distortion);
        d_min_distortion      = __min(d_min_distortion, d_distortion);
    }
    d_mean_distortion /= (int)_v_p_element.size();
    /////////////////////////////////////////////
    cout << _T("* distortion metric") << endl;
    cout << setprecision(3);
    cout << _T("  - mean: ") << d_mean_distortion << endl; 
    cout << _T("  - max : ") << d_max_distortion  << endl;
    cout << _T("  - min : ") << d_min_distortion  << endl;
} // end: Mesher_2D::dump_mesh_information()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Mesher_2D::write_mesh_file(tstring s_file_name, const Mesh_Packet_2D& mesh_packet)
{
    if (s_file_name.empty()) return;
    /////////////////////////////////////////////
    TCHAR a_extension[1024];
    _tsplitpath(s_file_name.c_str(), 0, 0, 0, a_extension);
    int  n_size = s_file_name.size() - _tcsclen(a_extension);
    s_file_name = s_file_name.substr(0, n_size);
    /////////////////////////////////////////////
    write_mesh_m2_file ((s_file_name+tstring(_T(".m2") )).c_str(), mesh_packet);
    write_mesh_mgt_file((s_file_name+tstring(_T(".mgt"))).c_str(), mesh_packet);
} // end: Mesher_2D::write_mesh_file
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Mesher_2D::write_mesh_m2_file(tstring s_file_name, const Mesh_Packet_2D& mesh_packet)
{
    if (s_file_name.empty()) return;
    /////////////////////////////////////////////
    ofstream f_m2_file(s_file_name.c_str());
    if (!f_m2_file)
    {
        cerr << _T("* cannot open output file - write_mesh_m2_file()") << endl;
        return;
    }
    /////////////////////////////////////////////
    f_m2_file << mesh_packet.v_node.size() << _T(" ") << mesh_packet.v_element.size() << endl;
    for (int i=0; i<mesh_packet.v_node.size(); ++i)
    {
        Mesh_Packet_2D::I_Node i_node = mesh_packet.v_node[i];
        f_m2_file << i+1 << _T(" ") << i_node.a_uv[0] << _T(" ") << i_node.a_uv[1] << endl;
    }
    /////////////////////////////////////////////
    for (int i=0; i<mesh_packet.v_element.size(); ++i)
    {
        Mesh_Packet_2D::I_Element i_element = mesh_packet.v_element[i];
        f_m2_file << i+1 << _T(" ") << i_element.v_node_index.size() << " ";
        for (int j=0; j<i_element.v_node_index.size(); ++j) f_m2_file << i_element.v_node_index[j]+1 << _T(" ");
        f_m2_file << endl;
    }
} // end: Mesher_2D::write_mesh_m2_file()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Mesher_2D::write_mesh_mgt_file(tstring s_file_name, const Mesh_Packet_2D& mesh_packet)
{
    if (s_file_name.empty()) return;
    /////////////////////////////////////////////
    ofstream f_mgt_file(s_file_name.c_str());
    if (!f_mgt_file)
    {
        cerr << _T("* cannot open output file - write_mesh_mgt_file()") << endl;
        return;
    }
    /////////////////////////////////////////////
    f_mgt_file << _T("*VERSION") << endl;
    f_mgt_file << _T("4.3.2")    << endl << endl;
    /////////////////////////////////////////////
    f_mgt_file << _T("*NODE")    << endl;
    for (int i=0; i<mesh_packet.v_node.size(); ++i)
    {
        Mesh_Packet_2D::I_Node i_node = mesh_packet.v_node[i];
        f_mgt_file << i+1 << _T(", ") << i_node.a_uv[0] << _T(", ") << i_node.a_uv[1] << _T(", 0.0") << endl;
    }
    f_mgt_file << endl;
    /////////////////////////////////////////////
    f_mgt_file << _T("*ELEMENT") << endl;
    for (int i=0; i<mesh_packet.v_element.size(); ++i)
    {
        Mesh_Packet_2D::I_Element i_element = mesh_packet.v_element[i];
        f_mgt_file << i+1 << _T(", PLATE, 1, 1, ");
        for (int j=0; j<i_element.v_node_index.size(); ++j) f_mgt_file << i_element.v_node_index[j]+1 << _T(", ");
        switch (i_element.v_node_index.size())
        {
        case 4: f_mgt_file << _T("1, 0")    << endl; break;
        case 3: f_mgt_file << _T("0, 1, 0") << endl; break;
        }
    }
} // end: Mesher_2D::write_mesh_mgt_file()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
#ifdef _VERBOSE
void Mesher_2D::dump_input_file(tstring s_file_name, const Mesh_Packet_2D& mesh_packet)
{
    ofstream f_dump_file(s_file_name.c_str());
    if (!f_dump_file) { cerr << _T("* cannot open dump file - dump_input_file()") << endl; return; }
    /////////////////////////////////////////
    f_dump_file << mesh_packet.d_mesh_size << endl;
    f_dump_file << endl;
    f_dump_file << (int)mesh_packet.v_node.size() << endl;
    for (int j=0; j<(int)mesh_packet.v_node.size(); ++j) f_dump_file << mesh_packet.v_node[j].a_uv[0] << _T(" ") << mesh_packet.v_node[j].a_uv[1] << _T(" ") << mesh_packet.v_node[j].b_hard_flag << _T(" ") << mesh_packet.v_node[j].b_island_flag << _T(" ") << mesh_packet.v_node[j].d_mesh_size << endl;
    f_dump_file << endl;
    f_dump_file << (int)mesh_packet.v_node_chain.size() << endl;
    for (j=0; j<(int)mesh_packet.v_node_chain.size(); ++j)
    {
        f_dump_file << _T("2 ") << (int)mesh_packet.v_node_chain[j].v_node_index.size() << endl;
        for (int k=0; k<(int)mesh_packet.v_node_chain[j].v_node_index.size(); ++k) f_dump_file << mesh_packet.v_node_chain[j].v_node_index[k] << _T(" ");
        f_dump_file << endl;
    }
    /////////////////////////////////////////
    f_dump_file.close();
} // end: Mesher_2D::dump_input_file()
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
void Mesher_2D::dump_mesh(const vector<Node*>& v_p_node, const vector<Element*>& v_p_element, tstring s_message) const
{
    debug_out() << _T("=================================================") << endl;
    debug_out() << _T("= BEGIN: ") << s_message                            << endl;
    debug_out() << _T("=================================================") << endl;
    /////////////////////////////////////////////
    debug_out() << _T("*VERSION") << endl;
    debug_out() << _T("4.3.2")    << endl << endl;
    /////////////////////////////////////////////
    debug_out() << _T("*NODE")    << endl;
    for (vector<Node*>::const_iterator itr_nd=v_p_node.begin(); itr_nd!=v_p_node.end(); ++itr_nd)
    {
        Node* p_node = *itr_nd;
        debug_out() << p_node->get_id() << _T(", ") << p_node->get_uv()[0] << _T(", ") << p_node->get_uv()[1] << _T(", 0.0") << endl;
    }
    debug_out() << endl;
    /////////////////////////////////////////////
    debug_out() << _T("*ELEMENT") << endl;
    for (vector<Element*>::const_iterator itr_el=v_p_element.begin(); itr_el!=v_p_element.end(); ++itr_el)
    {
        Element* p_element = *itr_el;
        debug_out() << p_element->get_id() << _T(", PLATE, 1, 1, ");
        for (int i=0; i<p_element->get_node_size(); ++i)
        {
            Node* p_node = p_element->get_node(i);
            debug_out() << p_node->get_id() << _T(", ");
        }
        switch (p_element->get_node_size())
        {
        case 4: debug_out() << _T("1, 0")    << endl; break;
        case 3: debug_out() << _T("0, 1, 0") << endl; break;
        }
    }
    /////////////////////////////////////////////
    debug_out() << _T("=================================================") << endl;
    debug_out() << _T("= END: ") << s_message                              << endl;
    debug_out() << _T("=================================================") << endl;
}
#endif
/////////////////////////////////////////////////////////////////////
