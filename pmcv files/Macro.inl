namespace Macro
{
    MACRO_T1 void assert_is_floating(const T1& crT)
    {
        static_assert(!std::is_floating_point<T1>::value,
            "Macro : This function is not support floating point type");
    }

    MACRO_T1 void assert_is_not_floating(const T1& crT)
    {
        static_assert(std::is_floating_point<T1>::value,
            "Macro : This function must be floating point type");
    }

    MACRO_T1 void assert_is_not_enumclass(const T1& crT)
    {
        typedef std::integral_constant<bool, 
            std::is_enum<T1>::value && !std::is_convertible<T1, int>::value> is_enum_class;

        static_assert(is_enum_class::value,
            "Macro : This function must be enum class type");
    }

    MACRO_T1_T2 void assert_is_not_base_of()
    {
        static_assert(std::is_base_of<T1, T2>::value,
            "Macro : This function must be inheritance relationship");
    }

    MACRO_T1_SEQ_CONT int GetSize(const SeqCont<T1>& crsT)
    {
        return static_cast<int>(crsT.size());
    }

    MACRO_T1_SEQ_CONT int Copy(const SeqCont<T1>& crsTA, SeqCont<T1>& rsTB)
    {
        rsTB.resize(crsTA.size());
        std::copy(crsTA.begin(), crsTA.end(), rsTB.begin());

        return GetSize(rsTB);
    }

    MACRO_T1_SEQ_CONT void Append(SeqCont<T1>& rsT, const SeqCont<T1>& crsT)
    {
        rsT.insert(rsT.end(), crsT.begin(), crsT.end());

        return;
    }

    MACRO_T1_SEQ_CONT bool IsSeqValid(const SeqCont<T1>& crsT)
    {
        for ( const auto& crT : crsT )
        {
            if ( !crT.IsValid() )
            {
                ASSERT(0); return false;
            }
        }

        return true;
    }

    MACRO_T1_SEQ_CONT bool IsPSeqValid(const SeqCont<T1*>& crspT)
    {
        for ( const auto& crpT : crspT )
        {
            if ( !crpT->IsValid() )
            {
                ASSERT(0); return false;
            }
        }

        return true;
    }

    MACRO_T1_SEQ_CONT bool IsExist(const T1& crT, const SeqCont<T1>& crsT)
    {
        assert_is_floating(crT);

        return std::find(crsT.begin(), crsT.end(), crT) != crsT.end();
    }

    MACRO_T1_SEQ_CONT bool IsInBound(const T1& crT, const SeqCont<T1>& crsT)
    {
        assert_is_floating(crT);

        const auto max = *std::max_element(crsT.begin(), crsT.end());
        const auto min = *std::min_element(crsT.begin(), crsT.end());

        return (crT <= max) && (crT >= min);
    }

    MACRO_T1_SEQ_CONT_FM bool IsInBound(const T1& crT, const SeqCont<T1>& crsT, Function LT)
    {
        const auto max = *std::max_element(crsT.begin(), crsT.end(), LT);
        const auto min = *std::min_element(crsT.begin(), crsT.end(), LT);

        return ( !LT(max, crT) ) && ( !LT(min, crT) );
    }

    MACRO_T1_SEQ_CONT int UniqueSort(SeqCont<T1>& rsT)
    {
        assert_is_floating(T1());

        std::sort(rsT.begin(), rsT.end());
        rsT.erase(std::unique(rsT.begin(), rsT.end()), rsT.end());

        return GetSize(rsT);
    }

    MACRO_T1_SEQ_CONT int UniqueSort(const SeqCont<T1>& crsT, SeqCont<T1>& rsT)
    {
        rsT.clear();

        std::copy(crsT.begin(), crsT.end(),
            std::back_inserter(rsT));

        return UniqueSort(rsT);
    }

    MACRO_T1_SEQ_CONT_FM int Filter(SeqCont<T1>& rsT, Function Func)
    {
        assert_is_floating(T1());

        rsT.erase(std::remove_if(rsT.begin(), rsT.end(), [&Func] (const T1& crT)
        {
            return !Func(crT);
        }), rsT.end());

        return GetSize(rsT);
    }

    MACRO_T1_SEQ_CONT_FM int Filter(const SeqCont<T1>& crsT, SeqCont<T1>& rsT, Function Func)
    {
        assert_is_floating(T1());
        
        rsT.resize(crsT.size());

        auto it = std::copy_if(crsT.begin(), crsT.end(), rsT.begin(), Func);
        rsT.resize(std::distance(rsT.begin(), it));

        return GetSize(rsT);
    }

    MACRO_T1_SEQ_CONT_FM void Split(const SeqCont<T1>& crsT, SeqCont<T1>& rsTrue, SeqCont<T1>& rsFalse, Function Func)
    {
        Filter(crsT, rsTrue, Func);
        Filter(crsT, rsFalse, [&Func] (const T1& crT)
        {
            return !Func(crT);
        });

        return;
    }

    MACRO_T2_SEQ_CONT_FM void Convert(const SeqCont<T1>& crsT1, SeqCont<T2>& rsT2, Function Func)
    {
        rsT2.clear();

        std::transform(crsT1.begin(), crsT1.end(),
            std::back_inserter(rsT2), Func);

        return;
    }

    MACRO_T1_SEQ_CONT T1 Fold(const SeqCont<T1>& crsT1, T1 Init)
    {
        for ( const auto& crT1 : crsT1 )
        {
            Init += crT1;
        }

        return Init;
    }

    MACRO_T2_SEQ_CONT_FM T2 Fold(const SeqCont<T1>& crsT1, T2 Init, Function Func)
    {
        for ( const auto& crT1 : crsT1 )
        {
            Init += Func(crT1);
        }

        return Init;
    }

    MACRO_T1_SEQ_CONT bool Equal(const SeqCont<T1>& crsT1)
    {
        assert_is_floating(T1());

        if ( crsT1.empty() ) return true;

        return std::all_of(crsT1.begin() + 1, crsT1.end(),
            [&crsT1] (const T1& Val)
        {
            return (crsT1[0] == Val);
        });
    }

    MACRO_T1_SEQ_CONT_FM bool Equal(const SeqCont<T1>& crsT1, Function Func)
    {
        //TODO 부동소수점에 대해서는 지금 DGNCompare에 존재하는 EQ만 사용 가능.
        assert_is_not_floating(T1());

        if ( crsT1.empty() ) return true;

        return std::all_of(crsT1.begin() + 1, crsT1.end(),
            [&Func, &crsT1] (const T1& Val)
        {
            return Func(crsT1[0], Val, static_cast<T1>(dgn::D_ZERO_LIMIT));
        });
    }

    MACRO_T1_SEQ_CONT_FM int GetSetByFunc(const SeqCont<T1>& crsTA, const SeqCont<T1>& crsTB, SeqCont<T1>& rsT, Function Func)
    {
        rsT.clear();

        SeqCont<T1> sTA;
        SeqCont<T1> sTB;

        UniqueSort(crsTA, sTA);
        UniqueSort(crsTB, sTB);

        Func(sTA, sTB, rsT);

        return GetSize(rsT);
    }

    MACRO_T1_SEQ_CONT int GetUnionSet(const SeqCont<T1>& crsTA, const SeqCont<T1>& crsTB, SeqCont<T1>& rsT)
    {
        assert_is_floating(T1());

        return GetSetByFunc(crsTA, crsTB, rsT, [] (const SeqCont<T1>& crsTA, const SeqCont<T1>& crsTB, SeqCont<T1>& rsT)
        {
            return std::set_union(crsTA.begin(), crsTA.end(), crsTB.begin(), crsTB.end(), std::back_inserter(rsT));
        });
    }

    MACRO_T1_SEQ_CONT int GetIntersectionSet(const SeqCont<T1>& crsTA, const SeqCont<T1>& crsTB, SeqCont<T1>& rsT)
    {
        assert_is_floating(T1());

        return GetSetByFunc(crsTA, crsTB, rsT, [] (const SeqCont<T1>& crsTA, const SeqCont<T1>& crsTB, SeqCont<T1>& rsT)
        {
            return std::set_intersection(crsTA.begin(), crsTA.end(), crsTB.begin(), crsTB.end(), std::back_inserter(rsT));
        });
    }

    MACRO_T1_SEQ_CONT int GetComplementSet(const SeqCont<T1>& crsTA, const SeqCont<T1>& crsTB, SeqCont<T1>& rsT)
    {
        assert_is_floating(T1());

        return GetSetByFunc(crsTA, crsTB, rsT, [] (const SeqCont<T1>& crsTA, const SeqCont<T1>& crsTB, SeqCont<T1>& rsT)
        {
            return std::set_difference(crsTA.begin(), crsTA.end(), crsTB.begin(), crsTB.end(), std::back_inserter(rsT));
        });
    }

    MACRO_T1 bool IsErrorIndex(const T1& nIndex)
    {
        static_assert(std::is_integral<T1>::value,
            "Macro::IsErrorIndex : The type of nIndex must be int");

        return nIndex == cnErrorIndex;
    }

    MACRO_T1_SEQ_CONT int GetMinElementIndex(const SeqCont<T1>& crsT)
    {
        if ( crsT.empty() )
        {
            ASSERT(0); return cnErrorIndex;
        }

        const auto& itr = std::min_element(crsT.begin(), crsT.end());

        return static_cast<int>(std::distance(crsT.begin(), itr));
    }

    MACRO_T1_SEQ_CONT int GetMaxElementIndex(const SeqCont<T1>& crsT)
    {
        if ( crsT.empty() )
        {
            ASSERT(0); return cnErrorIndex;
        }

        const auto& itr = std::max_element(crsT.begin(), crsT.end());

        return static_cast<int>(std::distance(crsT.begin(), itr));
    }

    MACRO_T1_SEQ_CONT_FM int GetMinElementIndex(const SeqCont<T1>& crsT, Function Func)
    {
        if ( crsT.empty() )
        {
            ASSERT(0); return cnErrorIndex;
        }

        SeqCont<decltype(Func(T1()))> sT;
        Convert(crsT, sT, Func);

        return GetMinElementIndex(sT);
    }

    MACRO_T1_SEQ_CONT_FM int GetMaxElementIndex(const SeqCont<T1>& crsT, Function Func)
    {
        if ( crsT.empty() )
        {
            ASSERT(0); return cnErrorIndex;
        }

        SeqCont<decltype(Func(T1()))> sT;
        Convert(crsT, sT, Func);

        return GetMaxElementIndex(sT);
    }

    MACRO_T1_SIZE UINT FindIdx(const std::array<T1, _size>& crA, const T1& crT)
    {
        auto itr = std::find(crA.begin(), crA.end(), crT);
        if ( itr == crA.end() )
        {
            ASSERT(0); // 없다니까 첫번째 값으로 ㄷㄷ
            itr = crA.begin();
        }

        return static_cast<UINT>(std::distance(crA.begin(), itr));
    }
    
    MACRO_T1 typename std::underlying_type<T1>::type Ordinal(const T1& cEnum)
    {
        assert_is_not_enumclass(T1());

        return static_cast<std::underlying_type<T1>::type>(cEnum);
    }

    MACRO_T1_T2 T2* SafeDownCast(T1* pT1)
    {
        assert_is_not_base_of<T1, T2>();

        auto* pT2 = dynamic_cast<T2*>(pT1);
        if ( !pT2 )
        {
            ASSERT(0); return nullptr;
        }

        return pT2;
    }

    MACRO_T1_FITA std::vector<T1> FITArray2Vec(const FITArray<T1>& craT)
    {
        auto nSize = craT.GetSize();

        std::vector<T1> vT;
        vT.resize(nSize);

        for ( decltype(nSize) iIdx = 0; iIdx < nSize; ++iIdx )
        {
            vT[iIdx] = craT[iIdx];
        }

        return vT;
    }

    MACRO_T1_T2 std::vector<T2> CArray2Vec(const CArray<T1, T1>& craT)
    {
        auto nSize = craT.GetSize();

        std::vector<T2> vT;
        vT.resize(nSize);

        for ( decltype(nSize) iIdx = 0; iIdx < nSize; ++iIdx )
        {
            vT[iIdx] = craT[iIdx];
        }

        return vT;
    }

    MACRO_T1 int FTArrayUnique(idgn::FTArray<T1>& raT1)
    {
        std::vector<T1> vT = FITArray2Vec(raT1);

        std::sort(vT.begin(), vT.end());
        vT.erase(std::unique(vT.begin(), vT.end()), vT.end());

        raT1.RemoveAll();
        raT1.SetSize(vT.size());
        for ( const auto& cT : vT )
        {
            raT1.PushBack(cT);
        }

        return raT1.GetSize();
    }

    MACRO_T1 int FArrayUnique(idgn::FArray<T1>& raT1)
    {
        std::vector<T1> vT = FITArray2Vec(raT1);

        std::sort(vT.begin(), vT.end());
        vT.erase(std::unique(vT.begin(), vT.end()), vT.end());

        raT1.RemoveAll();
        for ( const auto& cT : vT )
        {
            raT1.Add(cT);
        }

        return raT1.GetSize();
    }

    MACRO_T1_FITA inline void static_assert_FITArray(const FITArray<T1>& raT1)
    {
        // GetSize땜에 빌드가 안되긴 해도, 붙여주는게 나을듯
        static_assert(
            std::is_same<FITArray<T1>, FArray<T1>>::value ||
            std::is_same<FITArray<T1>, FTArray<T1>>::value,
            "Macro : Input argument must be FArray or FTArray");
    }

    MACRO_T1_FITA_F1 void FITArrayConvert(FITArray<T1>& raT1, _MonadicOperation Func)
    {
        static_assert_FITArray(raT1);

        auto Size = raT.GetSize();
        for ( decltype(Size) Idx = 0; Idx < Size; ++Idx )
        {
            raT[Idx] = Func(raT[Idx]);
        }
    }

    MACRO_T2_FITA_F1 void FITArrayConvert(const FITArray<T1>& craT1, idgn::FArray<T2>& raT2, _MonadicOperation Func)
    {
        static_assert_FITArray(craT1);

        raT2.RemoveAll();
        if ( craT1.IsEmpty() ) { return; }

        auto Size = craT1.GetSize();
        raT2.SetSize(Size);

        for ( decltype(Size) Idx = 0; Idx < Size; ++Idx )
        {
            raT2[Idx] = Func(craT1[Idx]);
        }
    }

    MACRO_T1_FITA_F1 int FITArrayFilter(const FITArray<T1>& craT1, idgn::FArray<T1>& raT1, _MonadicOperation Func)
    {
        static_assert_FITArray(craT1);

        raT1.RemoveAll();
        if ( craT1.IsEmpty() ) { return 0; }

        auto Size = craT1.GetSize();
        raT1.Reserve(Size);

        for ( decltype(Size) Idx = 0; Idx < Size; ++Idx )
        {
            const T1& crT = craT1[Idx];
            if ( Func(crT) )
            {
                raT1.Add(crT);
            }
        }

        return raT1.GetSize();
    }

    MACRO_T1_F1 int FArrayFilter(idgn::FArray<T1>& raT1, _MonadicOperation Func)
    {
        int nSize = raT1.GetSize();
        if ( nSize == 0 )
        {
            return 0;
        }

        idgn::FArrayi aiDel;
        aiDel.Reserve(nSize);
        for ( int nIdx = 0; nIdx < nSize; ++nIdx )
        {
            const T1& crT = raT1[nIdx];
            if ( !Func(crT) )
            {
                aiDel.Add(nIdx);
            }
        }

        std::reverse(aiDel.begin(), aiDel.end());
        for ( const auto& iDel :aiDel )
        {
            raT1.RemoveAt(iDel);
        }

        return raT1.GetSize();
    }

    MACRO_T1_FITA int GetFITArrayMaxIdx(const FITArray<T1>& craT1)
    {
        static_assert_FITArray(craT1);

        if ( craT1.IsEmpty() )
        {
            ASSERT(0); return -1;
        }

        const auto& itr = std::max_element(craT1.begin(), craT1.end());

        return static_cast<int>(std::distance(craT1.begin(), itr));
    }

    MACRO_T1_FITA_F1 int GetFITArrayMaxIdx(const FITArray<T1>& craT1, _MonadicOperation Func)
    {
        static_assert_FITArray(craT1);

        if ( craT1.IsEmpty() )
        {
            ASSERT(0); return -1;
        }

        idgn::FArray<decltype(Func(craT1[0]))> aTReturn;
        FITArrayConvert(craT1, aTReturn, Func);

        return GetFITArrayMaxIdx(aTReturn);
    }

    MACRO_T1_FITA int GetFITArrayMinIdx(const FITArray<T1>& craT1)
    {
        static_assert_FITArray(craT1);

        if ( craT1.IsEmpty() )
        {
            ASSERT(0); return -1;
        }

        const auto& itr = std::min_element(craT1.begin(), craT1.end());

        return static_cast<int>(std::distance(craT1.begin(), itr));
    }

    MACRO_T1_FITA_F1 int GetFITArrayMinIdx(const FITArray<T1>& craT1, _MonadicOperation Func)
    {
        static_assert_FITArray(craT1);

        if ( craT1.IsEmpty() )
        {
            ASSERT(0); return -1;
        }

        idgn::FArray<decltype(Func(craT1[0]))> aTReturn;
        FITArrayConvert(craT1, aTReturn, Func);

        return GetFITArrayMinIdx(aTReturn);
    }
}