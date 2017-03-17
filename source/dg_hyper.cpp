#define _CRTDBG_MAP_ALLOC
#include <cstdlib>
#include <crtdbg.h>

#include <iostream>
#include <fstream>
#include <iomanip>

#include "general_definitions.h"
#include "class_problem.h"

int main(int argc, const char* argv[])
{
	PROBLEM* problem = new PROBLEM();

	delete problem;

	_CrtDumpMemoryLeaks();
}