#include <emscripten.h>
#include <iostream>
#include "MathFunc.h"

extern "C"
{
    EMSCRIPTEN_KEEPALIVE
    int main()
    {
        return CDgnMathFunc::mathTan(5);   
        return 2;
    }
}


