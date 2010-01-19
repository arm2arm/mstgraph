#ifndef _MYRENDER_
#define _MYRENDER_
////////////////////////////////////////////////////
#ifdef WITH_CIMG

#define cimg_plugin  "plugins/draw_galaxypoints.h"
#define cimg_plugin1 "plugins/LUT_ALL.rgb.h"
#define cimg_plugin2 "plugins/draw_gradient.h"
#include "CImg.h"
using namespace cimg_library;
#undef min
#undef max
////////////////////////////////////////////////////
/***********************************************/
#include "global_typedefs.h"
#include "readers.h"
#include <utility>
#include <algorithm>
#include <map>
#include <set>
#include <iostream>
#include <vector>
#include <cstring>
#include "functions.h"
using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;

class CRender {
public:
    CRender(void);

    CRender(int _visu) : m_visu(_visu) {
    };
    ~CRender(void);
    void DoJob();
    void setup_cic_image();
    void PlotData(std::vector<MyFloat> &vec);
private:

    CImg<float> img, dest, res;
    int m_visu;
protected:

};
#endif
#endif


