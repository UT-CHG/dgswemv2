#include <iostream>
#include <fstream>
#include <iomanip>

#include "class_problem_SWE_2D_LINEAR.h"
#include "../SWE/swe_data.hpp"

int main(int argc, const char* argv[])
{  
    PROBLEM* problem = new PROBLEM();

	problem->Solve(EEULER, 1, 172800.0, 300.0);

    delete problem;
}