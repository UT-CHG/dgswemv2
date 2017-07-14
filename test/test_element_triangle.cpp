#include "general_definitions.hpp"
#include "utilities/almost_equal.hpp"
#include "geometry/mesh_definitions.hpp"
#include "problem/SWE/swe_problem.hpp"

const std::vector<double> IntegrationPhi_true = 
{
1.563315270799324e-03,
2.390314963233519e-03,
0.000000000000000e-00,
2.312867350286335e-03,
0.000000000000000e-00,
4.22441525166454e-06,
1.640762883746507e-03,
0.000000000000000e-00,
-2.253021467554421e-06,
0.000000000000000e-00,
8.98955565554214e-04,
0.000000000000000e-00,
9.85696892055059e-07,
0.000000000000000e-00,
1.971393784110118e-06,
3.801600983975879e-04,
0.000000000000000e-00,
-3.47893020725315e-07,
0.000000000000000e-00,
-6.9578604145063e-07,
0.000000000000000e-00,
1.233570669321846e-04,
0.000000000000000e-00,
9.66369502014764e-08,
0.000000000000000e-00,
1.932739004029528e-07,
0.000000000000000e-00,
6.18476481289449e-07,
2.96471990907582e-05,
0.000000000000000e-00,
-2.034462109504767e-08,
0.000000000000000e-00,
-4.068924219009533e-08,
0.000000000000000e-00,
-1.30205575008305e-07,
0.000000000000000e-00,
5.007828482545982e-06,
0.000000000000000e-00,
3.051693164257149e-09
};

int main() {
  using Utilities::almost_equal;
  bool error_found = false;

  using MasterType  = Master::Triangle<Basis::Dubiner_2D, Integration::Dunavant_2D>;
  using ShapeType   = Shape::StraightTriangle;
  using ElementType = Geometry::Element<2, MasterType, ShapeType, SWE::Data>;

  //make an equilateral triangle
  std::vector<Point<2> > vrtxs(3);
  vrtxs[0] = { -0.5, 0. };
  vrtxs[1] = { 0.5, 0. };
  vrtxs[2] = { 0, std::sqrt(3.) / 2. };

  MasterType master(10);
  ShapeType shape(vrtxs);

  ElementType triangle(0, master, vrtxs, std::vector<uint>(0), std::vector<unsigned char>(0));
  
  Integration::Dunavant_2D integ;
	std::vector<Point<2>> gp = integ.GetRule(20).second;

  std::vector<double> x = shape.InterpolateNodalValues({-0.5,0.5,0}, gp); 
  std::vector<double> y = shape.InterpolateNodalValues({0 ,0 ,std::sqrt(3.) / 2.}, gp);;
  
  std::vector<double> f_vals(triangle.data.get_ngp_internal());

  for(uint gp=0;gp<triangle.data.get_ngp_internal();gp++){
    f_vals[gp] = std::pow(x[gp],10)+std::pow(y[gp],10);
  }

  for(uint dof=0;dof<40;dof++){
    if(!almost_equal(IntegrationPhi_true[dof],triangle.IntegrationPhi(dof,f_vals)))
    {
      printf("%d, %.15e\n",dof, triangle.IntegrationPhi(dof,f_vals));
    }
  }
/*ShapeType shape(vrtxs);

  ElementType elt(0, m_tri, shape);

  //for this example we will consider the function f(x,y) = y. We will project it to the triangle and
  //ensurer that all integrals get evaulated correctly

  //check the function is being evaluated correctly at Gauss Points
  {
    std::function<double(Point)> x_ = [](Point pt) {
      return pt[0];
    };

    std::function<double(Point)> y_ = [](Point pt) {
      return pt[1];
    };

    std::vector<double> f_at_gp(elt.get_num_area_gauss_points() );

    Geometry::Quadrature::Dunavant quadr(2);
    const std::vector<Point>& gps = quadr.get_gp();

    elt.fill_area_gauss_points(x_, f_at_gp);
    for ( uint i = 0; i < gps.size(); ++i ) {
      if ( !almost_equal(f_at_gp[i], shape.master2elt(gps[i])[0] ) ) {
        std::cerr << "Error in evaluating function\n";
        std::cerr << "Got: " << f_at_gp[i] << " Should be: " << shape.master2elt(gps[i])[0] << "\n";
      }
    }

    elt.fill_area_gauss_points(y_, f_at_gp);
    for ( uint i = 0; i < gps.size(); ++i ) {
      if ( !almost_equal(f_at_gp[i], shape.master2elt(gps[i])[1] ) ) {
        std::cerr << "Error in evaluating function\n";
        std::cerr << "Got: " << f_at_gp[i] << " Should be: " << shape.master2elt(gps[i])[0] << "\n";
      }
    }
  }

  //check if Jacobian is correct by integrating 1 over the element
  {
    std::vector<double> f_at_gp(elt.get_num_area_gauss_points(), 1);
    double area = elt.integrate_area_phi(0,f_at_gp);

    if ( !almost_equal(area, std::sqrt(3.)/4.) ) {
      std::cerr << "Error in integrating 1 over the triangle\n";
      std::cerr << "Got: " << area << " Should be: " << std::sqrt(3.)/4. << "\n\n";
    }
  }

  //check normals of edges
  {
    std::array<double,2> normal;
    normal = {-std::sqrt(3.)*0.5, 0.5};

    std::array<double,2> computed = elt.get_normal(0,0);
    if ( !almost_equal(computed[0],normal[0]) || !almost_equal(computed[1],normal[1]) ) {
      std::cerr << "Error in computing normals\n";
      std::cerr << " Got: { " << computed[0] << " , " << computed[1] << " } Should be: { "
                << normal[0] << " , " << normal[1] << " }\n";
      error_found = true;
    }

    normal = {0, -1};

    computed = elt.get_normal(0,1);
    if ( !almost_equal(computed[0],normal[0]) || !almost_equal(computed[1],normal[1]) ) {
      std::cerr << "Error in computing normals\n";
      std::cerr << " Got: { " << computed[0] << " , " << computed[1] << " } Should be: { "
                << normal[0] << " , " << normal[1] << " }\n";
    }

    normal = {std::sqrt(3.)*0.5, 0.5};

    computed = elt.get_normal(0,2);
    if ( !almost_equal(computed[0],normal[0]) || !almost_equal(computed[1],normal[1]) ) {
      std::cerr << "Error in computing normals\n";
      std::cerr << " Got: { " << computed[0] << " , " << computed[1] << " } Should be: { "
                << normal[0] << " , " << normal[1] << " }\n\n";
      error_found = true;
    }
  }

  //check surface jacobians
  for ( uint edge = 0; edge < 3; ++edge ) {
    double ej = elt.get_edge_jacobian(0, edge);
    if ( !almost_equal(ej, 0.5) ) {
      std::cerr << "Error Edge jacobian for an equilateral\n";
      std::cerr << " For Edge( " << edge << " ): Got:  " << ej
                << "Should be: " << 0.5 << "\n";
      error_found = true;
    }
  }

  //check some integral evaluations
  {
    std::function<double(Point)> x_ = [](Point pt) {
      return pt[0];
    };

    std::vector<double> f_at_gp(elt.get_num_area_gauss_points() );
    elt.fill_area_gauss_points(x_, f_at_gp);
    double eval = elt.integrate_area_phi(0, f_at_gp);

    if ( !almost_equal(eval,0.) ) {
      std::cerr<< "Error in integrating x\n";
      std::cerr<< "Got: " << eval << " Should be: " << 0 << "\n";
      error_found = true;
    }

    std::function<double(Point)> y_ = [](Point pt) {
      return pt[1];
    };

    elt.fill_area_gauss_points(y_, f_at_gp);
    eval = elt.integrate_area_phi(0, f_at_gp);

    if ( !almost_equal(eval,1./8.) ) {
      std::cerr<< "Error in integrating y\n";
      std::cerr<< "Got: " << eval << " Should be: " << 1./8. << "\n";
      error_found = true;
    }

  }

  {//check get_master2elt
    std::array<double,2> master_pt = {-1,-1};
    std::array<double,2> elt_pt = elt.get_master2elt(master_pt);

    if ( !almost_equal(elt_pt[0], vrtxs[0][0]) || !almost_equal(elt_pt[1], vrtxs[0][1]) ) {
      std::cerr << "Error in get_master2elt\n";
      std::cerr << "Got: ( " << elt_pt[0] << ", " << elt_pt[1] << ") Should be: ( "
                << vrtxs[0][0] << ", " << vrtxs[0][1] << ")\n";
      error_found = true;
    }

    master_pt = { 1,-1};
    elt_pt = elt.get_master2elt(master_pt);

    if ( !almost_equal(elt_pt[0], vrtxs[1][0]) || !almost_equal(elt_pt[1], vrtxs[1][1]) ) {
      std::cerr << "Error in get_master2elt\n";
      std::cerr << "Got: ( " << elt_pt[0] << ", " << elt_pt[1] << ") Should be: ( "
                << vrtxs[1][0] << ", " << vrtxs[1][1] << ")\n";
      error_found = true;
    }

    master_pt = {-1, 1};
    elt_pt = elt.get_master2elt(master_pt);

    if ( !almost_equal(elt_pt[0], vrtxs[2][0]) || !almost_equal(elt_pt[1], vrtxs[2][1]) ) {
      std::cerr << "Error in get_master2elt\n";
      std::cerr << "Got: ( " << elt_pt[0] << ", " << elt_pt[1] << ") Should be: ( "
                << vrtxs[2][0] << ", " << vrtxs[2][1] << ")\n";
      error_found = true;
    }
  }

  if (error_found) {
    return 1;
  }

  return 0;
  */
}
