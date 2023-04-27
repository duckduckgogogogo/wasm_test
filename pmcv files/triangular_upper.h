/**
   \file       triangular_upper.h
   \brief      Definition of upper part of a symmetric template matrix.
   \copyright  (C)1999-2017, Computing Objects, France. info@computing-objects.com

   $Rev: 2963 $
   $Date: 2017-12-15 12:24:31 +0100 (ven., 15 déc. 2017) $ 
   */

/******************************************************************************************
   This program is not free software. It is the sole property of Computing Objects, France.
   You are allowed to modify it, to recompile it, to port it on several platforms,
   but not to redistribute it, even the modified version.
 ******************************************************************************************/

#ifndef __CM2_TRIANGULAR_UPPER_H__
#define __CM2_TRIANGULAR_UPPER_H__


 namespace cm2 { 

/**
   Upper triangular part of a symmetric matrix.

   This class is almost like its template parameter symmetric matrix class.  \n
   It has the same storage. \n
   The main difference is that triangular_upper is not considered as
   a symmetric matrix. \n
   Its mathematical shape is the same as its storage shape, i.e. a triangular shape.

   \see        symmetric_fixed, symmetric_full, triangular_lower.
   */
template <class SymMatrix>
class triangular_upper : public SymMatrix 
{

public:

/// self_type.
typedef triangular_upper<SymMatrix>                      self_type;
/// parent_type.
typedef SymMatrix                                        parent_type;

/// value_type
typedef typename parent_type::value_type                 value_type;
/// reference
typedef typename parent_type::reference                  reference;
/// const_reference
typedef typename parent_type::const_reference            const_reference;
/// pointer.
typedef typename parent_type::pointer                    pointer;
/// const_pointer.
typedef typename parent_type::const_pointer              const_pointer;
/// iterator
typedef typename parent_type::iterator                   iterator;
/// const_iterator
typedef typename parent_type::const_iterator             const_iterator;
/// size_type
typedef typename parent_type::size_type                  size_type;

/// shape_tag.
typedef cm2::non_symmetric_tag                           shape_tag;
/// orientation_tag
typedef typename cm2::orientation_as_upper<typename parent_type::storage_tag,
                                           typename parent_type::orientation_tag>::tag  orientation_tag;
/// storage_tag
typedef cm2::upper_right_storage_tag                     storage_tag;
/// fullness_tag
typedef cm2::non_full_tag                                fullness_tag;



/**
   Default constructor.
   Elements are uninitialized.
   */
triangular_upper() { }

/**
   Constructor with a common initializing value.
   All the elements are initialized to \p v.
   */
INLINE explicit triangular_upper (value_type v)           
   : parent_type(v) 
{ }

/**
   Constructor from a symmetric matrix.
   The lower part excluding the diagonal (considered null).
   */
triangular_upper (const parent_type& A) 
   : parent_type(A) 
{ }

/**
   Copy constructor.
   */
triangular_upper (const self_type& A) 
   : parent_type(A) 
{ }



/**@name Copy and assignment operators */
//@{
/// Scalar assignment.
INLINE value_type
operator= (value_type v) 
{ 
   return parent_type::operator=(v);
}
//@}


/**@name Value and reference access operators */
//@{
/**
   Returns the element (\p i, \p j) by reference.
   Range check performed in DEBUG mode.

   \pre \p j >= i. 
   */
template <class Index_>
reference      
operator() (Index_ i, Index_ j) 
{
   assert (j >= i);
   return parent_type::operator()(i, j); 
}

/**
   Returns the element (\p i, \p j) by value.
   Range check performed in DEBUG mode.
   */
INLINE value_type  
operator() (size_type i, size_type j) const 
{ 
   return (j >= i) ? parent_type::operator()(i, j) : value_type(0); 
}
//@}

};

} // end namespace cm2

#endif
