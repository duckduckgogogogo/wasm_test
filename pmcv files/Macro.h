#pragma once

#include <afxtempl.h>
#include <array>
#include "..\..\include\IDGN_frx\AllIDgnInclude.h"

#define MACRO_T1 template <typename T1>
#define MACRO_T1_T2 template <typename T1, typename T2>
#define MACRO_T1_SIZE template <typename T1, std::size_t _size>

#define MACRO_T1_F1  template <typename T1, class _MonadicOperation>
#define MACRO_T1_FITA template <typename T1, template <typename> class FITArray>
#define MACRO_T1_FITA_F1  template <typename T1, template <typename> class FITArray, class _MonadicOperation>
#define MACRO_T2_FITA_F1  template <typename T1, template <typename> class FITArray, typename T2, class _MonadicOperation>

#define MACRO_T1_SEQ_CONT template < \
    typename T1, \
    template <typename T, typename Alloc = std::allocator<T>> class SeqCont \
    >

#define MACRO_T1_SEQ_CONT_FM template < \
    typename T1, \
    template <typename T, typename Alloc = std::allocator<T>> class SeqCont, \
    class Function \
    >

#define MACRO_T2_SEQ_CONT_FM template < \
    typename T1, \
    typename T2, \
    template <typename T, typename Alloc = std::allocator<T>> class SeqCont, \
    class Function \
    >

namespace Macro
{
    // Sequence Container를 모두 지원하고 싶은 원대한 꿈이 있었으나...
    // 쉽지 않다 ㄷㄷ
    
    MACRO_T1_SEQ_CONT int GetSize(const SeqCont<T1>& crsT);
    MACRO_T1_SEQ_CONT int Copy(const SeqCont<T1>& crsTA, SeqCont<T1>& rsTB);
    MACRO_T1_SEQ_CONT void Append(SeqCont<T1>& rsT, const SeqCont<T1>& crsT);
    MACRO_T1_SEQ_CONT bool IsSeqValid(const SeqCont<T1>& crsT);
    MACRO_T1_SEQ_CONT bool IsPSeqValid(const SeqCont<T1*>& crspT);

    MACRO_T1_SEQ_CONT bool IsExist(const T1& crT, const SeqCont<T1>& crsT);
    MACRO_T1_SEQ_CONT bool IsInBound(const T1& crT, const SeqCont<T1>& crsT);
    MACRO_T1_SEQ_CONT_FM bool IsInBound(const T1& crT, const SeqCont<T1>& crsT, Function LT);

    MACRO_T1_SEQ_CONT int UniqueSort(SeqCont<T1>& rsT);
    MACRO_T1_SEQ_CONT int UniqueSort(const SeqCont<T1>& crsT, SeqCont<T1>& rsT);

    MACRO_T1_SEQ_CONT_FM int Filter(SeqCont<T1>& rsT, Function Func);
    MACRO_T1_SEQ_CONT_FM int Filter(const SeqCont<T1>& crsT, SeqCont<T1>& rsT, Function Func);
    MACRO_T1_SEQ_CONT_FM void Split(const SeqCont<T1>& crsT, SeqCont<T1>& rsTrue, SeqCont<T1>& rsFalse, Function Func);

    MACRO_T2_SEQ_CONT_FM void Convert(const SeqCont<T1>& crsT1, SeqCont<T2>& rsT2, Function Func);
    
    MACRO_T1_SEQ_CONT T1    Fold(const SeqCont<T1>& crsT1, T1 Init);
    MACRO_T2_SEQ_CONT_FM T2 Fold(const SeqCont<T1>& crsT1, T2 Init, Function Func);

    MACRO_T1_SEQ_CONT    bool Equal(const SeqCont<T1>& crsT1);
    MACRO_T1_SEQ_CONT_FM bool Equal(const SeqCont<T1>& crsT1, Function Func);

    MACRO_T1_SEQ_CONT int GetUnionSet(const SeqCont<T1>& crsTA, const SeqCont<T1>& crsTB, SeqCont<T1>& rsT);
    MACRO_T1_SEQ_CONT int GetIntersectionSet(const SeqCont<T1>& crsTA, const SeqCont<T1>& crsTB, SeqCont<T1>& rsT);
    MACRO_T1_SEQ_CONT int GetComplementSet(const SeqCont<T1>& crsTA, const SeqCont<T1>& crsTB, SeqCont<T1>& rsT);

    const int cnErrorIndex = -1;
    MACRO_T1 bool IsErrorIndex(const T1& nIndex);

    MACRO_T1_SEQ_CONT int GetMinElementIndex(const SeqCont<T1>& crsT);
    MACRO_T1_SEQ_CONT int GetMaxElementIndex(const SeqCont<T1>& crsT);
    MACRO_T1_SEQ_CONT_FM int GetMinElementIndex(const SeqCont<T1>& crsT, Function Func);
    MACRO_T1_SEQ_CONT_FM int GetMaxElementIndex(const SeqCont<T1>& crsT, Function Func);

    template <typename TEnum>
    struct TNumber
    {
        enum : int
        {
            Number = static_cast<int>(TEnum::Number),
        };
    };
    MACRO_T1_SIZE UINT FindIdx(const std::array<T1, _size>& crA, const T1& crT);
    MACRO_T1 typename std::underlying_type<T1>::type Ordinal(const T1& cEnum);

    MACRO_T1_T2 T2* SafeDownCast(T1* pI);

    // IDGN_FRX
    MACRO_T1_FITA std::vector<T1> FITArray2Vec(const FITArray<T1>& craT);
    MACRO_T1_T2 std::vector<T2> CArray2Vec(const CArray<T1, T1>& craT);
    MACRO_T1 int FTArrayUnique(idgn::FTArray<T1>& raT1);
    MACRO_T1 int FArrayUnique(idgn::FArray<T1>& raT1);

    MACRO_T1_FITA void static_assert_FITArray(const FITArray<T1>& raT1);
    MACRO_T1_FITA_F1 void FITArrayConvert(FITArray<T1>& raT1, _MonadicOperation Func);
    MACRO_T2_FITA_F1 void FITArrayConvert(const FITArray<T1>& craT1, idgn::FArray<T2>& raT2, _MonadicOperation Func);
    MACRO_T1_FITA_F1 int  FITArrayFilter(const FITArray<T1>& craT1, idgn::FArray<T1>& raT1, _MonadicOperation Func);
    MACRO_T1_F1      int  FArrayFilter(idgn::FArray<T1>& raT1, _MonadicOperation Func);

    MACRO_T1_FITA    int GetFITArrayMaxIdx(const FITArray<T1>& craT1);
    MACRO_T1_FITA_F1 int GetFITArrayMaxIdx(const FITArray<T1>& craT1, _MonadicOperation Func);
    MACRO_T1_FITA    int GetFITArrayMinIdx(const FITArray<T1>& craT1);
    MACRO_T1_FITA_F1 int GetFITArrayMinIdx(const FITArray<T1>& craT1, _MonadicOperation Func);
}

#include "Macro.inl"

#undef MACRO_T1
#undef MACRO_T1_T2
#undef MACRO_T1_SIZE

#undef MACRO_T1_F1
#undef MACRO_T1_FITA
#undef MACRO_T1_FITA_F1
#undef MACRO_T2_FITA_F1

#undef MACRO_T1_SEQ_CONT
#undef MACRO_T1_SEQ_CONT_FM
#undef MACRO_T2_SEQ_CONT_FM