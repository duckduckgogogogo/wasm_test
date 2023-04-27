#if defined(__HEADERPRE_H__)
#error Export header file error. check following causes. (ask jbseon)
// 원인:
// 1. HeaderPre.h와 HeaderPost.h 사이에 다른 Export Header 파일이 포함된 경우
// 2. 다른 Export Header 파일을 포함한 경우 그 파일에서 HeaderPost.h를 Include하지 않은 경우
//    찾는 방법
//    Find에서 Find What에 "HeaderP.*" 를 입력하고 Regular Expression Check 켠다.
//    검색된 것 중에서 "HeaderPre.h 만 있는 것을 찾으면 된다.
#endif

#define __HEADERPRE_H__

//#if defined( __USE_DEF_FILE__) // Module Define File을 사용 할 경우
//  #error Module definition file can't be used this project (ask jbseon)

//#else   // Module Define File을 사용하지 않을 경우

#if defined(__PROJ_DGN_BASE__)  // 현재 Project에 포함되는 경우
#define       AFX_EXT_CLASS
#define __MY_EXT_API__        AFX_EXT_API
#define __MY_EXT_DATA__       AFX_EXT_DATA
#else
#if defined(_WINDLL)      // 다른 DLL Project에 포함되는 경우
#define       __declspec(dllimport)
#define __MY_EXT_API__        __declspec(dllimport)
#define __MY_EXT_DATA__       __declspec(dllimport)
#else                     // 다른 EXE Project에 포함되는 경우
#define __MY_EXT_API__        AFX_EXT_API
#define __MY_EXT_DATA__       AFX_EXT_DATA
#endif
#endif
//#endif

