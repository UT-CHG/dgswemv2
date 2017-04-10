#define _CRTDBG_MAP_ALLOC
#include <cstdlib>
#include <crtdbg.h>

#include <iostream>
#include <fstream>
#include <iomanip>

#include "general_definitions.h"
#include "problem\class_problem_SWE_2D.h"

int main(int argc, const char* argv[])
{
    PROBLEM* problem = new PROBLEM();

    problem->EETimeStepper(10000);
    problem->WriteDataVTK();

    delete problem;

    _CrtDumpMemoryLeaks();
}