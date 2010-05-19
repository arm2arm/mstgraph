#include "MSTree.h"
#include "utils.h"

CMSTree::~CMSTree(void)
	{
	save();
	delete tree;

	}


