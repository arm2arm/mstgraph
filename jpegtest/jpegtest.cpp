#define cimg_plugin "plugins/draw_galaxypoints.h"
#include "CImg.h"
using namespace cimg_library;
#undef min
#undef max
#include <vector>
#include "readers.h"
// Main procedure for drawing 
//----------------
int DrawVec(std::vector<float> xy, const float xrange[2], const float yrange[2]) {
  
  // Read command line arguments.
  const unsigned int s_nb = xy.size()/2;
  float s_xmin = float(xrange[0]), s_xmax = float(xrange[1]);
  
    
  CImg<> samples(s_nb,2);
  cimg_forX(samples,i) {
    const float
      x = (float)xy[i],
      y = (float)(xy[i+s_nb]);
    samples(i,0) = x;
    samples(i,1) = y;
  }

  // Fit Gaussian function on the sample points and display curve iterations.
  CImgDisplay disp(640,640,"Galaxy points",0);
 
 
  for (;!disp.is_closed() && !disp.is_keyQ() && !disp.is_keyESC();) {    
    // Display points.
	  const unsigned char black[] = { 0,0,0 }, gray[] = { 228,228,228 };
	CImg<unsigned char>(disp.width(),disp.height(),1,3,255).draw_galaxypoints(samples,xrange, yrange).
   //   draw_rectangle(5,7,150,100,gray,0.9f).draw_rectangle(5,7,150,100,black,1,~0U).
      display(disp.resize(false).wait(1));
  }

  return 0;
}

void read_gadget(std::string fname)
	{
	All.ONLY_TYPE=4;
	read_ic12(fname.c_str());
	}
int main() {
	std::vector<float> xy;
	read_gadget("C:\\arm2arm\\DATA\\MODEL7\\MODELS\\MODEL7\\RUNG2\\SNAPS\\snap_gal_sfr_0450");
	xy.resize(All.NumPart*2);
	for(int i=0;i<All.NumPart;i++)
		{
		xy[i]=Part[i].Pos[0];xy[i+All.NumPart]=Part[i].Pos[1];
		}
	float xrange[]={-10,10};
	DrawVec(xy,  xrange, xrange);
	return 0;
	}
