// domain_maker.h: domain-maker for 2D domain makers

//////////////////////////////////////////////////////////////////////////////////
// Shin Dae-Seock (dsshin@midasIT.com)
// Don't believe in magic;
// understand what your libraries do, how they do it, and at what cost they do it.
//////////////////////////////////////////////////////////////////////////////////

#ifndef DOMAIN_MAKER_H_
#define DOMAIN_MAKER_H_

#pragma warning( disable : 4786 )
#pragma warning( disable : 4290 )

#include "data_2d.h"
#include "mesh_packet_2d.h"

#include <map>
#include <string>

/////////////////////////////////////////////////////////////////////
namespace domain_maker
{
    using namespace std;
    using namespace data_2d;
    using namespace mesh_packet_2d;
    /////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////
    class interrupt_exception : public exception {};
    /////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////
    class Segment
    {
    public:
        Segment(Node* p_node_1, Node* p_node_2) : _p_node_1(p_node_1), _p_node_2(p_node_2) { _p_preceder=0;  _p_follower=0;}
        /////////////////////////////////////////
        Node*    get_node_1()   const { return _p_node_1;   }
        Node*    get_node_2()   const { return _p_node_2;   }
        Segment* get_preceder() const { return _p_preceder; }
        Segment* get_follower() const { return _p_follower; }
        /////////////////////////////////////////
        void set_node(Node* p_node_1, Node* p_node_2) 
        {
            _p_node_1 = p_node_1;
            _p_node_2 = p_node_2;
        }
        void set_preceder(Segment* p_segment) { _p_preceder = p_segment; }
        void set_follower(Segment* p_segment) { _p_follower = p_segment; }
        /////////////////////////////////////////
    private:
        Node*    _p_node_1;
        Node*    _p_node_2;
        Segment* _p_preceder;
        Segment* _p_follower;
    }; // end: class Segment
    /////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////
    class Edge : public Wire
    {
    public:
        Edge(int n_tag = -999) : _n_tag(n_tag)
        {
            _b_open_flag       = false;
            _b_reversible_flag = true;
            _n_id              = 0;
#ifdef _DEBUG
            _n_id = _n_edge_id++;
#endif
        }
        ~Edge() { _v_p_preceder.clear(); _v_p_follower.clear(); }
        /////////////////////////////////////////
        int   get_id()         const { return _n_id;             }
        int   get_tag()        const { return _n_tag;            }
        Node* get_start_node() const { return _v_p_node.front(); }
        Node* get_end_node()   const { return _v_p_node.back();  }
        /////////////////////////////////////////
        void set_open_edge()            { _b_open_flag       = true;       }
        void set_irreversible_edge()    { _b_reversible_flag = false;      }
        void add_preceder(Edge* p_edge) { _v_p_preceder.push_back(p_edge); }
        void add_follower(Edge* p_edge) { _v_p_follower.push_back(p_edge); }
        void clear_preceder()           { _v_p_preceder.clear();           }
        void clear_follower()           { _v_p_follower.clear();           }
        void clear_traversal()
        {
            clear_preceder();
            clear_follower();
        }
        void remove_preceder(Edge* p_edge)
        {
            vector<Edge*>::iterator itr_eg=find(_v_p_preceder.begin(), _v_p_preceder.end(), p_edge);
            if (itr_eg != _v_p_preceder.end()) _v_p_preceder.erase(itr_eg);
        }
        void remove_follower(Edge* p_edge)
        {
            vector<Edge*>::iterator itr_eg=find(_v_p_follower.begin(), _v_p_follower.end(), p_edge);
            if (itr_eg != _v_p_follower.end()) _v_p_follower.erase(itr_eg);
        }
        void remove_traversal(Edge* p_edge)
        {
            remove_preceder(p_edge);
            remove_follower(p_edge);
        }
        void reverse_direction()
        {
            vector<Edge*> v_temp = _v_p_preceder;
            _v_p_preceder        = _v_p_follower;
            _v_p_follower        = v_temp;
            reverse(_v_p_node.begin(), _v_p_node.end());
        }
        /////////////////////////////////////////
        bool  is_open_edge()            const { return _b_open_flag;           }
        bool  is_reversible_edge()      const { return _b_reversible_flag;     } 
        int   get_preceder_size()       const { return (int)_v_p_preceder.size();   }
        Edge* get_preceder(int n_index) const { return _v_p_preceder[n_index]; }
        int   get_follower_size()       const { return (int)_v_p_follower.size();   }
        Edge* get_follower(int n_index) const { return _v_p_follower[n_index]; }
        /////////////////////////////////////////
    private:
        int           _n_id;
        int           _n_tag;
        bool          _b_open_flag;
        bool          _b_reversible_flag;
        vector<Edge*> _v_p_preceder;
        vector<Edge*> _v_p_follower;
#ifdef _DEBUG
        static int _n_edge_id;
#endif
    }; // end: class Edge
    /////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////
    class Domain
    {
    public:
        Domain() {}
        ~Domain() { _v_p_wire.clear(); _v_p_island_node.clear(); }
        /////////////////////////////////////////
        void set_wire(const vector<Wire*>& v_p_wire)        { _v_p_wire        = v_p_wire;                                       }
        void add_wire(Wire* p_wire)                         { _v_p_wire.push_back(p_wire);                                       }
        void set_island_node(const vector<Node*>& v_p_node) { _v_p_island_node = v_p_node;                                       }
        void add_island_node(Node* p_node)                  { _v_p_island_node.push_back(p_node);                                }
        void remove_wire(int n_index)                       { _v_p_wire.erase(_v_p_wire.begin()+n_index);                        }
        void remove_wire(Wire* p_wire)                      { _v_p_wire.erase(find(_v_p_wire.begin(), _v_p_wire.end(), p_wire)); }
        /////////////////////////////////////////
        int   get_wire_size()                       const { return (int)_v_p_wire.size();          }
        Wire* get_wire(int n_index)                 const { return _v_p_wire[n_index];        }
        int   get_island_node_size()                const { return (int)_v_p_island_node.size();   }
        Node* get_island_node(int n_index)          const { return _v_p_island_node[n_index]; }
        bool  is_internal_point(const double* a_uv) const;
        /////////////////////////////////////////
    private:
        vector<Wire*> _v_p_wire;
        vector<Node*> _v_p_island_node;
    }; // end: class Domain
    /////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////
    class Domain_Maker
    {
    public:
        Domain_Maker();
        ~Domain_Maker();
        /////////////////////////////////////////
        void make_domain(Mesh_Packet_2D& mesh_packet, vector<int>& v_unused_node_chain_index, vector<int>& v_unused_node_index, bool b_keep_inner_domain_flag, bool b_skip_inner_edge_flag=false, bool (*p_interrupt_handler)()=0) throw(exception);
        /////////////////////////////////////////
        static void dump_input_file(tstring s_file_name, const Mesh_Packet_2D& mesh_packet);
        static void write_domain_file(tstring s_file_name, const Mesh_Packet_2D& mesh_packet);
        /////////////////////////////////////////
    private:
        bool check_interrupt() const { return _p_interrupt_handler ? _p_interrupt_handler() : false; }
        void set_input_data(const Mesh_Packet_2D& mesh_packet) throw(exception);
        void check_input_data() throw(exception);
        /////////////////////////////////////////
        void  find_open_edge() const;
        void  find_irreversible_edge() const throw(exception);
        void  find_outermost_edge(const vector<Edge*>& v_p_bulk_edge, vector<Edge*>& v_p_outermost_edge) const throw(exception);
        Edge* pick_leaf_edge() const;
        void  duplicate_reversible_edge();
        void  make_single_edge_domain();
        void  make_multi_edge_domain(Edge* p_edge) throw(exception);
        Edge* pick_wire_head_edge() const;
        void  set_inner_domain(bool b_keep_inner_domain_flag) throw(exception);
        void  check_inner_wire_connection(Domain* p_domain) throw(exception);
        void  check_matryoshka_domain(Domain* p_domain) throw(exception);
        void  insert_internal_wire_node(vector<int>& v_unused_node_index) throw(exception);
        void  finalize_domain_making(Mesh_Packet_2D& mesh_packet, vector<int>& v_unused_node_chain_index);
        /////////////////////////////////////////
        vector<Node*>             _v_p_node;
        vector<Edge*>             _v_p_edge;
        vector<Wire*>             _v_p_wire;
        vector<Domain*>           _v_p_domain;
        vector<Node*>             _v_p_island_node;
        map<Node*, int>           _m_node_index;
        map< int, vector<Node*> > _m_node_chain_node;
        /////////////////////////////////////////
        double _d_tolerance;
        /////////////////////////////////////////
        bool (*_p_interrupt_handler)();
    }; // end: class Domain_Maker
} // end: namespace domain_maker
/////////////////////////////////////////////////////////////////////

#endif // DOMAIN_MAKER_H_