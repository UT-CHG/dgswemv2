#define _CRTDBG_MAP_ALLOC
#include <cstdlib>
#include <crtdbg.h>

#include <iostream>
#include <fstream>
#include <iomanip>

#include "general_definitions.h"
#include "problem\SWE_2D\class_problem_SWE_2D.h"

int main(int argc, const char* argv[])
{
    PROBLEM* problem = new PROBLEM();

	problem->Solve(EEULER, 0.005, 1.0, 0.01);


    delete problem;

    _CrtDumpMemoryLeaks();
}