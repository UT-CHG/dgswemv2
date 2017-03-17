#ifndef CLASS_INTERFACE_H
#define CLASS_INTERFACE_H

class INTERFACE {
	friend class PROBLEM;

protected:
	int number_gp;

	double* normal_x;
	double* normal_y;

	double** u_boundary_in;
	double** u_boundary_ex;

public:
	INTERFACE(int number_gp, double** u_boundary)
		: number_gp(number_gp), u_boundary_in(u_boundary) {};
    virtual ~INTERFACE() = default;

	void SetPointerEX(double** u_boundary) { this->u_boundary_ex = u_boundary; };
};

#endif