// Substantial portitions of this file are copied out of the blaze library.
// Original License listed below

//=================================================================================================
/*!
//  \FILE67 blaze/math/serialization/VectorSerializer.h
//  \brief Serialization of dense and sparse vectors
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

#ifndef SERIALIZATION_BLAZE_VECTOR_HPP
#define SERIALIZATION_BLAZE_VECTOR_HPP

#include <blaze/math/serialization/VectorSerializer.h>

namespace hpx {
namespace serialization {
namespace detail {
template <typename VT>
void serialize_blaze_header(output_archive& ar, const VT& vec) {
    using ET = blaze::ElementType_t<VT>;

    ar << uint8_t(1U);
    ar << uint8_t(blaze::VectorSerializer::VectorValueMapping<VT>::value);
    ar << uint8_t(blaze::TypeValueMapping<ET>::value);
    ar << uint8_t(sizeof(ET));
    ar << uint64_t(vec.size());
    ar << uint64_t(blaze::IsDenseVector_v<VT> ? vec.size() : vec.nonZeros());
}

// Dense specialization
template <typename VT, bool TF>
void serialize_blaze_vector(output_archive& ar, const blaze::DenseVector<VT, TF>& vec) {
    for (size_t i = 0Ul; i < (~vec).size(); ++i) {
        ar << (~vec)[i];
    }
}

// Sparse specialization
template <typename VT, bool TF>
void serialize_blaze_vector(output_archive& ar, const blaze::SparseVector<VT, TF>& vec) {
    using ConstIterator = blaze::ConstIterator_t<VT>;

    for (ConstIterator element = (~vec).begin(); element != (~vec).end(); ++element) {
        ar << element->index() << element->value();
    }
}
}

template <typename VT, bool TF>
void serialize(output_archive& ar, const blaze::Vector<VT, TF>& vec, unsigned) {
    detail::serialize_blaze_header(ar, ~vec);
    detail::serialize_blaze_vector(ar, ~vec);
}

namespace detail {

struct blaze_vector_serializer_helper {
    uint8_t version_;      //!< The version of the archive.
    uint8_t type_;         //!< The type of the vector.
    uint8_t elementType_;  //!< The type of an element.
    uint8_t elementSize_;  //!< The size in bytes of a single element of the vector.
    uint64_t size_;        //!< The size of the vector.
    uint64_t number_;      //!< The total number of elements contained in the vector.

    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        ar& version_& type_& elementType_& elementSize_& size_& number_;
    }

    template <typename VT>
    bool check_header(const VT& vec) {
        using ET = blaze::ElementType_t<VT>;
        bool error_found{false};
        if (version_ != 1UL) {
            std::cerr << "Invalid version detected" << std::endl;
            error_found = true;
        } else if ((type_ & 1U) != 0U || (type_ & (~3U)) != 0U) {
            std::cerr << "Invalid vector type detected" << std::endl;
            error_found = true;
        } else if (elementType_ != blaze::TypeValueMapping<ET>::value) {
            std::cerr << "Invalid element type detected" << std::endl;
            error_found = true;
        } else if (elementSize_ != sizeof(ET)) {
            std::cerr << "Invalid element size detected" << std::endl;
            error_found = true;
        } else if (!blaze::IsResizable_v<VT> && size_ != vec.size()) {
            std::cerr << "Invalid vector size detected" << std::endl;
            error_found = true;
        } else if (number_ > size_) {
            std::cerr << "Invalid number of elements detected" << std::endl;
            error_found = true;
        }
        return !error_found;
    }

    template <typename VT, bool TF>
    blaze::DisableIf_t<blaze::IsResizable_v<VT>> prepare_vector(blaze::DenseVector<VT, TF>& vec) {
        reset(~vec);
    }

    template <typename VT, bool TF>
    blaze::DisableIf_t<blaze::IsResizable_v<VT>> prepare_vector(blaze::SparseVector<VT, TF>& vec) {
        (~vec).reserve(number_);
        reset(~vec);
    }

    template <typename VT>
    blaze::EnableIf_t<blaze::IsResizable_v<VT>> prepare_vector(VT& vec) {
        vec.resize(size_, false);
        vec.reserve(number_);
        reset(vec);
    }

    template <typename VT>
    void deserializeVector(input_archive& ar, VT& vec) {
        if (type_ == 0U) {
            deserializeDenseVector(ar, vec);
        } else if (type_ == 2U) {
            deserializeSparseVector(ar, vec);
        } else {
            BLAZE_INTERNAL_ASSERT(false, "Undefined type flag");
        }
    }

    template <typename VT, bool TF>
    void deserializeDenseVector(input_archive& ar, blaze::DenseVector<VT, TF>& vec) {
        using ET = blaze::ElementType_t<VT>;
        ET value{};

        for (size_t i = 0UL; i < size_; ++i) {
            ar >> value;
            (~vec)[i] = value;
        }
    }

    template <typename VT, bool TF>
    void deserializeDenseVector(input_archive& ar, blaze::SparseVector<VT, TF>& vec) {
        using ET = blaze::ElementType_t<VT>;
        ET value{};

        for (size_t i = 0UL; i < size_; ++i) {
            ar >> value;
            (~vec)[i] = value;
        }
    }

    template <typename VT, bool TF>
    void deserializeSparseVector(input_archive& ar, blaze::DenseVector<VT, TF>& vec) {
        using ET = blaze::ElementType_t<VT>;

        size_t index{0UL};
        ET value{};

        for (size_t i = 0UL; i < number_; ++i) {
            ar >> index >> value;
            (~vec)[index] = value;
        }
    }

    template <typename VT, bool TF>
    void deserializeSparseVector(input_archive& ar, blaze::SparseVector<VT, TF>& vec) {
        using ET = blaze::ElementType_t<VT>;

        size_t index{0UL};
        ET value{};

        for (size_t i = 0; i < number_; ++i) {
            ar >> index >> value;
            (~vec).append(index, value, false);
        }
    }
};
}

template <typename VT, bool TF>
void serialize(input_archive& ar, blaze::Vector<VT, TF>& vec, unsigned) {
    detail::blaze_vector_serializer_helper vs;

    ar >> vs;
    assert(vs.check_header(~vec));

    vs.prepare_vector(~vec);
    vs.deserializeVector(ar, ~vec);
}
}
}
#endif
