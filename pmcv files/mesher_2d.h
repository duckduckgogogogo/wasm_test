// mesher_2d.h: 2D auto-mesh generator

//////////////////////////////////////////////////////////////////////////////////
// Shin Dae-Seock (dsshin@midasIT.com)
// Don't believe in magic;
// understand what your libraries do, how they do it, and at what cost they do it.
//////////////////////////////////////////////////////////////////////////////////

#ifndef MESHER_2D_H_
#define MESHER_2D_H_

#pragma warning( disable : 4786 )
#pragma warning( disable : 4290 )

#include "data_2d.h"
#include "mesh_packet_2d.h"

#include <vector>
#include <map>
#include <set>
#include <string>
#include <algorithm>
#include <exception>
#include <functional>

/////////////////////////////////////////////////////////////////////
/*
#if defined(_DEBUG)
#define _VERBOSE
#endif
*/
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
namespace mesher_2d
{
    using namespace std;
    using namespace data_2d;
    using namespace mesh_packet_2d;
    /////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////
    class interrupt_exception : public exception {};
    /////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////
    class Offset
    {
    public:
        /////////////////////////////////////////
        enum OFFSET_TYPE { OT_END = 1, OT_SIDE, OT_CORNER, OT_REVERSAL }; 
        /////////////////////////////////////////
        void set_angle(double d_angle)
        {
            _d_angle = d_angle;
            if      (                     _d_angle < 120.0) _t_type = OT_END;
            else if (_d_angle >= 120.0 && _d_angle < 216.0) _t_type = OT_SIDE;
            else if (_d_angle >= 216.0 && _d_angle < 308.5) _t_type = OT_CORNER;
            else                                            _t_type = OT_REVERSAL;
        }
        void set_interval(double d_interval)         { _d_interval = d_interval;                }
        void add_uv(const double* a_uv)              { _v_uv.push_back(UV(a_uv));               }
        void set_uv(int n_index, const double* a_uv) { copy(a_uv, a_uv+2, _v_uv[n_index].a_uv); }
        void clear()                                 { _v_uv.clear();                           }
        /////////////////////////////////////////
        OFFSET_TYPE   get_type()          const { return _t_type;             }
        double        get_angle()         const { return _d_angle;            }
        double        get_interval()      const { return _d_interval;         }
        bool          is_empty()          const { return _v_uv.empty();       }
        int           get_size()          const { return (int)_v_uv.size();   }
        const double* get_uv(int n_index) const { return _v_uv[n_index].a_uv; }
        /////////////////////////////////////////
    private:
        OFFSET_TYPE _t_type;
        double      _d_angle;
        double      _d_interval;
        vector<UV>  _v_uv;
    }; // end: class Offset
    /////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////
    class Front
    {
    public:
        Front(Node* p_node_1, Node* p_node_2 = 0) : _p_node_1(p_node_1), _p_node_2(p_node_2), _b_obsolete_flag(false)
        {
            if (_p_node_2) for (int i=0; i<2; ++i) _a_centroid_uv[i] = (_p_node_1->get_uv()[i] + _p_node_2->get_uv()[i]) / 2.0;
            else           copy(_p_node_1->get_uv(), _p_node_1->get_uv()+2, _a_centroid_uv);
        }
        /////////////////////////////////////////
        void set_node_1(Node* p_node_1) { _p_node_1        = p_node_1; }
        void set_node_2(Node* p_node_2) { _p_node_2        = p_node_2; }
        void set_obsolete_front()       { _b_obsolete_flag = true;     }
        /////////////////////////////////////////
        Node*         get_node_1()        const { return _p_node_1;        }
        Node*         get_node_2()        const { return _p_node_2;        }
        const double* get_centroid_uv()   const { return _a_centroid_uv;   }
        bool          is_obsolete_front() const { return _b_obsolete_flag; }
        /////////////////////////////////////////
    private:
        Node*  _p_node_1;
        Node*  _p_node_2;
        double _a_centroid_uv[2];
        bool   _b_obsolete_flag;
    }; // end: class Front
    /////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////
    class Mesher_2D
    {
    public:
        enum MESH_TYPE { MT_TRIA = 1, MT_COMB, MT_QUAD };
        /////////////////////////////////////////
    public:
        Mesher_2D();
        ~Mesher_2D();
        /////////////////////////////////////////
        virtual void generate_mesh(Mesh_Packet_2D& mesh_packet, bool (*p_interrupt_handler)()=0) throw(exception) = 0;
        /////////////////////////////////////////
        void dump_mesh_information() const;
        /////////////////////////////////////////
        static void write_mesh_file(tstring s_file_name, const Mesh_Packet_2D& mesh_packet);
        static void write_mesh_m2_file(tstring s_file_name, const Mesh_Packet_2D& mesh_packet);
        static void write_mesh_mgt_file(tstring s_file_name, const Mesh_Packet_2D& mesh_packet);
        /////////////////////////////////////////
    protected:
        virtual void clean_up();
        /////////////////////////////////////////
        bool check_interrupt() const { return _p_interrupt_handler ? _p_interrupt_handler() : false; }
        /////////////////////////////////////////
        void set_input_data(const Mesh_Packet_2D& mesh_packet) throw(exception);
        void check_input_data() throw(exception);
        /////////////////////////////////////////
        virtual void prepare_mesh() throw(exception);
        virtual void generate_offset_mesh() throw(exception);
        virtual void finalize_mesh(Mesh_Packet_2D& mesh_packet) const;
        virtual void check_mesh_shape();
        /////////////////////////////////////////
        void generate_front_mesh(vector<Front> v_front, bool b_check_in_domain_flag=false) throw(exception);
        /////////////////////////////////////////
        void transform_into_triangle_mesh() throw(exception);
        void transform_into_combined_mesh(bool b_delaunay_flag = true, bool b_weak_merge_flag = false) throw(exception);
        void make_delaunay_triangle(const vector<Element*>& v_triangle) const;
        bool relax_combined_mesh() throw(exception);
        bool relax_quadrilateral_mesh(bool b_node_elimination_flag = true, bool b_element_elimination_flag = true, bool b_diagonal_swapping_flag = true, bool b_side_elimination_flag = true) throw(exception);
        /////////////////////////////////////////
        bool check_input_wire_node_pair(Node* p_node_1, Node* p_node_2) const;
        void set_neighbor_relation(vector<Element*> v_p_element) const;
        void update_neighbor_relation(set<Element*> v_p_element) const;
        /////////////////////////////////////////
        Node*  get_unique_node(const double* a_uv, bool b_hard_node_flag = false);
        void   seed_node(vector<const double*> v_a_uv, const vector<double>& v_xi, bool b_hard_node_flag, vector<Node*>& v_p_node, double d_length = 0.0);
        void   estimate_seed_interval(vector<const double*> v_a_uv, double d_seed_size_1, double d_seed_size_2, int n_division, vector<double>& v_interval, double d_length = 0.0) const;
        double estimate_ideal_division(vector<const double*> v_a_uv, double d_seed_size_1, double d_seed_size_2, double d_length = 0.0) const;
        /////////////////////////////////////////
        bool check_offset_node(Node* p_node, Wire* p_mother_wire, const Offset& offset, const map<Node*, Offset>& m_offset, bool b_self_check) const;
        bool check_wire_cross(const double* a_uv_1, const double* a_uv_2) const;
        bool check_in_domain_uv(const double* a_uv);
        /////////////////////////////////////////
        void smooth_mesh(int n_time = 2) const throw(exception);
        void smooth_mesh_laplace(int n_time = 3) const throw(exception);
        void smooth_mesh_winslow(int n_time = 3) const throw(exception);
        void smooth_mesh_angle(int n_time = 1) const throw(exception);
        /////////////////////////////////////////
#ifdef _VERBOSE
        void dump_input_file(tstring s_file_name, const Mesh_Packet_2D& mesh_packet);
        void dump_mesh(const vector<Node*>& v_p_node, const vector<Element*>& v_p_element, tstring s_message) const;
        /////////////////////////////////////////
        static ofstream& debug_out() { return _f_debug; }
#endif
        /////////////////////////////////////////
    protected:
        vector<Node*>               _v_p_node;
        vector<Element*>            _v_p_element;
        vector<Wire*>               _v_p_wire;
        vector<Node*>               _v_p_island_node;
        map< Wire*, vector<Node*> > _m_input_wire_node;
        /////////////////////////////////////////
        MESH_TYPE _t_mesh_type;
        double    _d_mesh_size;
        double    _d_refinement_ratio;
        double    _d_tolerance;
        /////////////////////////////////////////
        bool (*_p_interrupt_handler)();
        /////////////////////////////////////////
#ifdef _VERBOSE
        static ofstream _f_debug;
#endif
    }; // end: class Mesher_2D
    /////////////////////////////////////////////
} // end: namespace mesher_2d
/////////////////////////////////////////////////////////////////////

#endif // MESHER_2D_H_