#define _CRTDBG_MAP_ALLOC
#include <cstdlib>
#include <crtdbg.h>

#include <iostream>
#include <fstream>
#include <iomanip>

#include "class_problem_SWE_2D_LINEAR.h"

int main(int argc, const char* argv[])
{
    PROBLEM* problem = new PROBLEM();

	problem->Solve(EEULER, 1, 432000.0, 300.0);


    delete problem;

    _CrtDumpMemoryLeaks();
}