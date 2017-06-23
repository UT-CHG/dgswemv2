#ifndef CLASS_MESH_HPP
#define CLASS_MESH_HPP

#include "problem/SWE/swe_data.hpp"
#include "problem/SWE/swe_kernels.hpp"
#include "problem/SWE/swe_boundary_conditions.hpp"
#include "problem/SWE/swe_definitions.hpp"

#include "general_definitions.hpp"
#include "stepper.hpp"

#include "class_element.hpp"
#include "class_boundary.hpp"

class MESH {
private:
    uint p;
	
	Master::Triangle<Basis::Dubiner_2D, Integration::Dunavant_2D>* triangle; 
	Shape::StraightTriangle* shape;

    std::map<uint, Element<2, Master::Triangle<Basis::Dubiner_2D, Integration::Dunavant_2D>, 
		Shape::StraightTriangle, SWE::Data>*> elements;

	std::map<unsigned char, std::vector<Boundary<>*>> boundaries;
	std::vector<Interface<>*> interfaces;

public:
	MESH(uint p) : p(p), snapshot(0) {}
	~MESH();

	void RectangularDomainTest(double, double, uint, uint, uint);
	void Solve();
	void WriteDataVTK();
	uint snapshot;

private:
	void InitializeBoundariesInterfaces();
    void InitializeVTK();
};

#endif