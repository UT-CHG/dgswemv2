// Substantial portitions of this file are copied out of the blaze library.
// Original License listed below

//=================================================================================================
/*!
//  \file blaze/math/serialization/MatrixSerializer.h
//  \brief Serialization of dense and sparse matrices
//
//  Copyright (C) 2012-2018 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. You can redistribute it and/or modify it under
//  the terms of the New (Revised) BSD License. Redistribution and use in source and binary
//  forms, with or without modification, are permitted provided that the following conditions
//  are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice, this list
//     of conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//  3. Neither the names of the Blaze development group nor the names of its contributors
//     may be used to endorse or promote products derived from this software without specific
//     prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
//  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
//  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
//  SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
//  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
//  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
//  DAMAGE.
*/
//=================================================================================================

#ifndef SERIALIZATION_BLAZE_MATRIX_HPP
#define SERIALIZATION_BLAZE_MATRIX_HPP

#include <blaze/math/serialization/MatrixSerializer.h>

namespace hpx {
namespace serialization {
namespace detail {
template<typename MT >
void serialize_blaze_matrix_header(output_archive& ar, const MT& mat ) {
    using ET = blaze::ElementType_t<MT>;

   ar << uint8_t ( 1U );
   ar << uint8_t ( blaze::MatrixSerializer::MatrixValueMapping<MT>::value );
   ar << uint8_t ( blaze::TypeValueMapping<ET>::value );
   ar << uint8_t ( sizeof( ET ) );
   ar << uint64_t( mat.rows() );
   ar << uint64_t( mat.columns() );
   ar << uint64_t( ( blaze::IsDenseMatrix_v<MT> ) ? ( mat.rows()*mat.columns() ) : ( mat.nonZeros() ) );
}

//Dense specialization
template< typename MT, bool SO >
void serialize_blaze_matrix(output_archive& ar, const blaze::DenseMatrix<MT,SO>& mat ) {
    if( blaze::IsRowMajorMatrix_v<MT> ) {
        for( size_t i=0UL; i<(~mat).rows(); ++i ) {
            for( size_t j=0UL; j<(~mat).columns(); ++j ) {
                ar << (~mat)(i,j);
            }
        }
    } else {
        for( size_t j=0UL; j<(~mat).columns(); ++j ) {
            for( size_t i=0UL; i<(~mat).rows(); ++i ) {
                ar << (~mat)(i,j);
            }
        }
    }
}

//Sparse specialization
template<typename MT, bool SO >
void serialize_blaze_matrix(output_archive& ar, const blaze::SparseMatrix<MT,SO>& mat )
{
    using ConstIterator = blaze::ConstIterator_t<MT>;

    if( blaze::IsRowMajorMatrix_v<MT> ) {
        for( size_t i=0UL; i<(~mat).rows(); ++i ) {
            ar << uint64_t( (~mat).nonZeros( i ) );
            for( ConstIterator element=(~mat).begin(i); element!=(~mat).end(i); ++element ) {
                ar << element->index() << element->value();
            }
        }
    } else {
        for( size_t j=0UL; j<(~mat).columns(); ++j ) {
            ar << uint64_t( (~mat).nonZeros( j ) );
            for( ConstIterator element=(~mat).begin(j); element!=(~mat).end(j); ++element ) {
                ar << element->index() << element->value();
            }
        }
    }
}
}

template <typename MT, bool SO>
void serialize(output_archive& ar, const blaze::Matrix<MT,SO>& mat, unsigned) {
    detail::serialize_blaze_matrix_header(ar, ~mat);
    detail::serialize_blaze_matrix(ar, ~mat);
}


namespace detail {

struct blaze_matrix_serializer_helper {
   uint8_t  version_;      //!< The version of the archive.
   uint8_t  type_;         //!< The type of the matrix.
   uint8_t  elementType_;  //!< The type of an element.
   uint8_t  elementSize_;  //!< The size in bytes of a single element of the matrix.
   uint64_t rows_;         //!< The number of rows of the matrix.
   uint64_t columns_;      //!< The number of columns of the matrix.
   uint64_t number_;       //!< The total number of elements contained in the matrix.

    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        ar & version_ & type_ & elementType_ & elementSize_ & rows_ & columns_ & number_;
    }

    template <typename MT>
    bool check_header(const MT& mat) {
        using ET = blaze::ElementType_t<MT>;
        bool error_found{false};
        if( version_ != 1UL ) {
            std::cerr << "Invalid version detected" << std::endl;
            error_found = true;
        } else if ( ( type_ & 1U ) != 1U || ( type_ & (~7U) ) != 0U ) {
            std::cerr << "Invalid matrix type detected" << std::endl;
            error_found = true;
        } else if ( elementType_ != blaze::TypeValueMapping<ET>::value ) {
            std::cerr << "Invalid element type detected" << std::endl;
            error_found = true;
        } else if ( elementSize_ != sizeof( ET ) ) {
            std::cerr << "Invalid element size detected" << std::endl;
            error_found = true;
        } else if ( !blaze::IsResizable_v<MT> && ( rows_ != mat.rows() || columns_ != mat.columns() ) ) {
            std::cerr << "Invalid matrix size detected" << std::endl;
            error_found = true;
        } else if( number_ > rows_*columns_ ) {
            std::cerr << "Invalid number of elements detected" << std::endl;
            error_found = true;
        }
        return !error_found;
    }


    template< typename MT, bool SO >
    blaze::DisableIf_t< blaze::IsResizable_v<MT> > prepare_matrix( blaze::DenseMatrix<MT,SO>& mat ) {
        reset( ~mat );
    }

    template< typename MT, bool SO >
    blaze::DisableIf_t< blaze::IsResizable_v<MT> > prepare_matrix( blaze::SparseMatrix<MT,SO>& mat ) {
        (~mat).reserve( number_ );
        reset( ~mat );
    }

    template< typename MT >
    blaze::EnableIf_t< blaze::IsResizable_v<MT> > prepare_matrix( MT& mat ) {
        mat.resize ( rows_, columns_, false );
        mat.reserve( number_ );
        reset( mat );
    }

    template<typename MT >
    void deserializeMatrix(input_archive& ar, MT& mat ) {
        if( type_ == 1U ) {
            deserializeDenseRowMatrix( ar, ~mat );
        } else if( type_ == 5UL ) {
            deserializeDenseColumnMatrix( ar, ~mat );
        } else if( type_ == 3UL ) {
            deserializeSparseRowMatrix( ar, ~mat );
        } else if( type_ == 7UL ) {
            deserializeSparseColumnMatrix( ar, ~mat );
        } else {
            BLAZE_INTERNAL_ASSERT( false, "Undefined type flag" );
        }
    }

    template<typename MT, bool SO >
    void deserializeDenseRowMatrix(input_archive& ar, blaze::DenseMatrix<MT,SO>& mat ) {
        using ET = blaze::ElementType_t<MT>;

        ET value{};

        for( size_t i=0UL; i<rows_; ++i ) {
            for ( size_t j = 0UL; j < columns_; ++j ) {
                ar >> value;
                (~mat)(i,j) = value;
            }
        }
    }

    template<typename MT, bool SO >
    blaze::EnableIf_t< blaze::IsNumeric_v< blaze::ElementType_t<MT> > >
    deserializeDenseRowMatrix(input_archive& ar, blaze::SparseMatrix<MT,SO>& mat ) {
        blaze::DynamicMatrix< blaze::ElementType_t<MT>, blaze::rowMajor > tmp( rows_, columns_ );
        deserializeDenseRowMatrix( ar, tmp );
        (~mat) = tmp;
    }

    template<typename MT, bool SO >
    blaze::DisableIf_t< blaze::IsNumeric_v< blaze::ElementType_t<MT> > >
    deserializeDenseRowMatrix(input_archive& ar, blaze::SparseMatrix<MT,SO>& mat ) {
        using ET = blaze::ElementType_t<MT>;

        ET value{};

        const size_t dim1( ( SO == blaze::rowMajor )?( rows_ ):( columns_ ) );
        const size_t dim2( ( SO != blaze::rowMajor )?( rows_ ):( columns_ ) );

        for( size_t i=0UL; i<dim1; ++i ) {
            (~mat).reserve( i, dim2 );
        }

        for( size_t i=0UL; i<rows_; ++i ) {
            for( size_t j=0UL; j < columns_; ++j ) {
                ar >> value;
                (~mat).append( i, j, value, false );
            }
        }
    }

    template<typename MT, bool SO >
    void deserializeDenseColumnMatrix(input_archive& ar, blaze::DenseMatrix<MT,SO>& mat )
    {
        using ET = blaze::ElementType_t<MT>;

        ET value{};

        for( size_t j=0UL; j<columns_; ++j ) {
            for( size_t i=0UL; i < rows_; ++i ) {
                ar >> value;
                (~mat)(i,j) = value;
            }
        }
    }

    template< typename MT, bool SO >
    blaze::EnableIf_t< blaze::IsNumeric_v< blaze::ElementType_t<MT> > >
    deserializeDenseColumnMatrix(input_archive& ar, blaze::SparseMatrix<MT,SO>& mat ) {
        blaze::DynamicMatrix< blaze::ElementType_t<MT>, blaze::columnMajor > tmp( rows_, columns_ );
        deserializeDenseColumnMatrix( ar, tmp );
        (~mat) = tmp;
    }

    template<typename MT, bool SO >
    blaze::DisableIf_t< blaze::IsNumeric_v< blaze::ElementType_t<MT> > >
    deserializeDenseColumnMatrix(input_archive& ar, blaze::SparseMatrix<MT,SO>& mat )
    {
        using ET = blaze::ElementType_t<MT>;

        ET value{};

        const size_t dim1( ( SO == blaze::rowMajor )?( rows_ ):( columns_ ) );
        const size_t dim2( ( SO != blaze::rowMajor )?( rows_ ):( columns_ ) );

        for( size_t i=0UL; i<dim1; ++i ) {
            (~mat).reserve( i, dim2 );
        }

        for( size_t j=0UL; j<columns_; ++j ) {
            for ( size_t i=0UL; i < rows_; ++i ) {
                ar >> value;
                (~mat).append( i, j, value, false );
            }
        }
    }

    template<typename MT, bool SO >
    void deserializeSparseRowMatrix(input_archive& ar, blaze::DenseMatrix<MT,SO>& mat ) {
        using ET = blaze::ElementType_t<MT>;

        uint64_t number( 0UL );
        size_t   index ( 0UL );
        ET       value {};

        for( size_t i=0UL; i<rows_; ++i ) {
            ar >> number;
            for ( size_t j=0UL; j < number; ++j )  {
                ar >> index >> value;
                (~mat)(i,index) = value;
            }
        }
    }

    template<typename MT >
    void deserializeSparseRowMatrix(input_archive& ar, blaze::SparseMatrix<MT,blaze::rowMajor>& mat )
    {
        using ET = blaze::ElementType_t<MT>;

        uint64_t number( 0UL );
        size_t   index ( 0UL );
        ET       value {};

        for( size_t i=0UL; i<rows_; ++i )
        {
            ar >> number;

            for ( size_t j = 0UL; j < number; ++j ) {
                ar >> index >> value;
                (~mat).append( i, index, value, false );
            }

            (~mat).finalize( i );
        }
    }

    template<typename MT >     // Type of the matrix
    void deserializeSparseRowMatrix(input_archive& ar, blaze::SparseMatrix<MT,blaze::columnMajor>& mat )
    {
        blaze::CompressedMatrix< blaze::ElementType_t<MT>, blaze::rowMajor > tmp( rows_, columns_, number_ );
        deserializeSparseRowMatrix( ar, tmp );
        (~mat) = tmp;
    }

    template<typename MT, bool SO >
    void deserializeSparseColumnMatrix(input_archive& ar, blaze::DenseMatrix<MT,SO>& mat ) {
        using ET = blaze::ElementType_t<MT>;

        uint64_t number( 0UL );
        size_t   index ( 0UL );
        ET       value {};

        for( size_t j=0UL; j<columns_; ++j ) {
            ar >> number;
            for ( size_t i=0UL; i < number; ++i ) {
                ar >> index >> value;
                (~mat)(index,j) = value;
            }
        }
    }

    template<typename MT >
    void deserializeSparseColumnMatrix(input_archive& ar, blaze::SparseMatrix<MT,blaze::rowMajor>& mat ) {
        blaze::CompressedMatrix< blaze::ElementType_t<MT>, blaze::columnMajor > tmp( rows_, columns_, number_ );
        deserializeSparseColumnMatrix( ar, tmp );
        (~mat) = tmp;
    }

    template<typename MT >
    void deserializeSparseColumnMatrix(input_archive& ar, blaze::SparseMatrix<MT,blaze::columnMajor>& mat ) {
        using ET = blaze::ElementType_t<MT>;

        uint64_t number( 0UL );
        size_t   index ( 0UL );
        ET       value {};

        for( size_t j=0UL; j<columns_; ++j )
        {
            ar >> number;

            for ( size_t i=0UL; i < number; ++i ) {
                ar >> index >> value;
                (~mat).append( index, j, value, false );
            }

            (~mat).finalize( j );
        }
    }
};
}

template<typename MT, bool SO >
void serialize(input_archive& ar, blaze::Matrix<MT,SO>& mat, unsigned)
{
    detail::blaze_matrix_serializer_helper ms;

    ar >> ms;
    assert(ms.check_header(~mat));

    ms.prepare_matrix( ~mat );
    ms.deserializeMatrix( ar, ~mat );
}

/*
template <typename Type, bool SO>
void serialize(output_archive& ar, const blaze::IdentityMatrix<Type,SO>& id, unsigned) {
    ar << id.capacity();
}

template <typename Type, bool SO>
void serialize(input_archive& ar, blaze::IdentityMatrix<Type,SO>& id, unsigned) {
    size_t n;
    ar >> n;
    (~id).reset();
    (~id) = blaze::IdentityMatrix<Type,SO>(n);
    }*/
}
}
#endif
