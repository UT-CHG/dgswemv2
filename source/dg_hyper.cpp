#define _CRTDBG_MAP_ALLOC
#include <cstdlib>
#include <crtdbg.h>

#include <iostream>
#include <fstream>
#include <iomanip>

#include "class_mesh_v2.h"

int main(int argc, const char* argv[])
{
	MESH* mesh = new MESH(2, 0);
	
	delete mesh;

	_CrtDumpMemoryLeaks();
}