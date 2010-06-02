#include "OAFHelper.h"

/*COAFHelper::COAFHelper(void)
	{
	}

COAFHelper::~COAFHelper(void)
	{
	}
*/
bool COAFHelper::is_id_inselection(int id)
	{
	if(idu_sel.find(id)==idu_sel.end())
		return false;
	return true;
	}
bool COAFHelper::insert(int id, int snap, unsigned char type, float* pos, float *vel,
			float mass, float hsml,float rho, float u, float *z, float sfr)
	{
	bool is_id_there=is_id_inselection(id);
	if(is_id_there)
		{
		
		}
	return is_id_there;
	}