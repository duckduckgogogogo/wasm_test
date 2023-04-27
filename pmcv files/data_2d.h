// data_2d.h: data definition for 2D auto-meshers

//////////////////////////////////////////////////////////////////////////////////
// Shin Dae-Seock (dsshin@midasIT.com)
// Don't believe in magic;
// understand what your libraries do, how they do it, and at what cost they do it.
//////////////////////////////////////////////////////////////////////////////////

#ifndef DATA_2D_H_
#define DATA_2D_H_

#include <vector>
#include <algorithm>
#include <exception>

#ifdef _UNICODE
#define tstring wstring
#else
#define tstring string
#endif

/////////////////////////////////////////////////////////////////////
namespace data_2d
{
    using namespace std;
    /////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////
    class Node;
    class Element;
    class Wire;
    /////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////
    struct UV
    {
        UV() { a_uv[0] = a_uv[1] = 0.0; }
        UV(const double* uv) { copy(uv, uv+2, a_uv); }
        /////////////////////////////////////////
        double a_uv[2];
    }; // end: struct UV
    /////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////
    class Node
    {
    public:
        Node(const double* a_uv, int n_id = 0, bool b_hard_flag = false, double d_mesh_size = 0.0)
            : _n_id(n_id), _b_hard_flag(b_hard_flag), _d_mesh_size(d_mesh_size)
        {
            copy(a_uv, a_uv+2, _a_uv);
#ifdef _DEBUG
            _n_id = _n_node_id++;
#endif
        }
        virtual ~Node() {}
        /////////////////////////////////////////
        void set_id(int n_id)                  { _n_id          = n_id;      }
        void set_hard_node()                   { _b_hard_flag   = true;      }
        void reset_hard_node()                 { _b_hard_flag   = false;     }
        void set_uv(const double* a_uv)        { copy(a_uv, a_uv+2, _a_uv);  }
        void set_mesh_size(double d_mesh_size) { _d_mesh_size = d_mesh_size; }
        /////////////////////////////////////////
        int           get_id()         const { return _n_id;          }
        bool          is_hard_node()   const { return _b_hard_flag;   }
        const double* get_uv()         const { return _a_uv;          }
        double        get_mesh_size()  const { return _d_mesh_size;   }
        /////////////////////////////////////////
    protected:
        int    _n_id;
        bool   _b_hard_flag;
        double _d_mesh_size;
        double _a_uv[2];		
#ifdef _DEBUG
        static int _n_node_id;
#endif
    }; // end: class Node
    /////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////
    class Element
    {
    public:
        Element(Node* p_node_1, Node* p_node_2, Node* p_node_3)
        {
            _v_p_node.push_back(p_node_1);
            _v_p_node.push_back(p_node_2);
            _v_p_node.push_back(p_node_3);
            /////////////////////////////////////
#ifdef _DEBUG
            _n_id = _n_element_id++;
#endif
        }
        Element(Node* p_node_1, Node* p_node_2, Node* p_node_3, Node* p_node_4)
        {
            _v_p_node.push_back(p_node_1);
            _v_p_node.push_back(p_node_2);
            _v_p_node.push_back(p_node_3);
            _v_p_node.push_back(p_node_4);
            /////////////////////////////////////
#ifdef _DEBUG
            _n_id = _n_element_id++;
#endif
        }
        virtual ~Element() { _v_p_node.clear(); _v_p_neighbor.clear(); }
        /////////////////////////////////////////
        void set_id(int n_id)                                   { _n_id                  = n_id;                                     }
        void set_node(const vector<Node*>& v_p_node)            { _v_p_node              = v_p_node;                                 }
        void change_node(int n_index, Node* p_node)             { _v_p_node[n_index]     = p_node;                                   }
        void set_neighbor(const vector<Element*>& v_p_neighbor) { _v_p_neighbor          = v_p_neighbor;                             }
        void change_neighbor(int n_index, Element* p_neighbor)  { _v_p_neighbor[n_index] = p_neighbor;                               }
        void reset_neighbor()                                   { _v_p_neighbor          = vector<Element*>(_v_p_node.size(), 0);    }
        void remove_node(int n_index)                           { _v_p_node.erase(_v_p_node.begin()+n_index);                        }
        void remove_node(Node* p_node)                          { _v_p_node.erase(find(_v_p_node.begin(), _v_p_node.end(), p_node)); }
        void insert_node(int n_index, Node* p_node)
        {
            vector<Node*> v_p_node;
            for (int i=0;       i<n_index;               ++i) v_p_node.push_back(_v_p_node[i]);
            v_p_node.push_back(p_node);
            for (int i=n_index; i<(int)_v_p_node.size(); ++i) v_p_node.push_back(_v_p_node[i]);
            _v_p_node = v_p_node;
        }
        /////////////////////////////////////////
        int      get_id()                                const { return _n_id;                  }
        int      get_node_size()                         const { return (int)_v_p_node.size();  }
        Node*    get_node(int n_index)                   const { return _v_p_node[n_index];     }
        int      get_node_index(Node* p_node)            const
        {
            for (int i=0; i<(int)_v_p_node.size(); ++i)
            {
                if (_v_p_node[i] == p_node) return i;
            }
            return -1;
        }
        int      get_neighbor_size()                     const { return (int)_v_p_neighbor.size();   }
        Element* get_neighbor(int n_index)               const { return _v_p_neighbor[n_index]; }
        int      get_neighbor_index(Element* p_neighbor) const
        {
            for (int i=0; i<(int)_v_p_neighbor.size(); ++i)
            {
                if (_v_p_neighbor[i] == p_neighbor) return i;
            }
            return -1;
        }
        double   get_distortion_metric()                 const;
        /////////////////////////////////////////
    protected:
        int              _n_id;
        vector<Node*>    _v_p_node;
        vector<Element*> _v_p_neighbor;
        /////////////////////////////////////////
#ifdef _DEBUG
        static int _n_element_id;
#endif
    }; // end: class Element
    /////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////////////////////////
    class Wire
    {
    public:
        enum WIRE_TYPE { WT_OUTER = 1, WT_INNER, WT_OPEN };
        /////////////////////////////////////////
    public:
        Wire() { _t_type = WT_OUTER; }
        virtual ~Wire() { _v_p_node.clear(); }
        /////////////////////////////////////////
        void set_type(WIRE_TYPE t_type)              { _t_type   = t_type;   }
        void set_node(const vector<Node*>& v_p_node) { _v_p_node = v_p_node; }
        /////////////////////////////////////////
        WIRE_TYPE            get_type()                   const { return _t_type;          }
        int                  get_node_size()              const { return (int)_v_p_node.size(); }
        const vector<Node*>& get_node()                   const { return _v_p_node;        }
        int                  get_node_index(Node* p_node) const
        {
            for (int i=0; i<(int)_v_p_node.size(); ++i)
            {
                if (_v_p_node[i] == p_node) return i;
            }
            return -1;
        }
        /////////////////////////////////////////
    protected:
        WIRE_TYPE     _t_type;
        vector<Node*> _v_p_node;
    }; // end: class Wire
} // end: namespace data_2d
/////////////////////////////////////////////////////////////////////

#endif // DATA_2D_H_
