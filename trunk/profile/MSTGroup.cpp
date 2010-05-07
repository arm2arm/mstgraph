#include "MSTGroup.h"
//int CMSTGroup::ID = 0; // definition outside class declaration
bool operator<(const CMSTGroup& left, const CMSTGroup& right){
	return left.Ntotal > right.Ntotal;
};
