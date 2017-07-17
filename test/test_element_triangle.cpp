#include "general_definitions.hpp"
#include "utilities/almost_equal.hpp"
#include "geometry/mesh_definitions.hpp"
#include "problem/SWE/swe_problem.hpp"

const std::vector<double> IntegrationPhi_true = {
  6.881941874331419e-01,
  -8.89156081756484e-02,
  7.216878364870322e-02,
  1.20281306081172e-02,
  0.00000000000000e+00,
  4.811252243246882e-03,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00,
  0.00000000000000e+00
};

const std::vector<double> IntegrationDPhiDX_true = {
  0.00000000000000e+00,
  0.00000000000000e+00,
  1.376388374866284,
  0.00000000000000e+00,
  6.212068893253612e-01,
  4.330127018922193e-01,
  0.00000000000000e+00,
  4.541823898251592e-01,
  1.732050807568877e-01,
  8.09844822527753e-01,
  0.00000000000000e+00,
  3.593266739736606e-01,
  8.66025403784439e-02,
  4.752989329117831e-01,
  2.886751345948129e-01,
  0.00000000000000e+00,
  2.976844657418943e-01,
  4.948716593053934e-02,
  3.477531724792759e-01,
  1.649572197684645e-01,
  5.677806551742285e-01,
  0.00000000000000e+00,
  2.542695064773368e-01,
  3.092947870658709e-02,
  2.799779589772887e-01,
  1.030982623552903e-01,
  3.839882003961648e-01,
  2.165063509461096e-01,
  0.00000000000000e+00,
  2.219874007710344e-01,
  2.061965247105806e-02,
  2.368352966075515e-01,
  6.873217490352688e-02,
  2.924290225727635e-01,
  1.443375672974064e-01,
  4.382557445972097e-01,
  0.00000000000000e+00,
  1.970208683922866e-01,
  1.443375672974064e-02,
  2.063521709566122e-01,
  4.811252243246881e-02,
  2.388944320048427e-01,
  1.010362971081845e-01,
  3.216406511354539e-01,
  1.732050807568877e-01,
  0.00000000000000e+00,
  1.771259393312325e-01,
  1.049727762162956e-02,
  1.833727084180284e-01,
  3.499092540543186e-02,
  2.038020106931869e-01,
  7.348094335140692e-02,
  2.541597625891523e-01,
  1.259673314595547e-01,
  3.571900656194553e-01,
  0.00000000000000e+00,
  1.6089440460229e-01,
  7.87295821622217e-03,
  1.652850866122932e-01,
  2.624319405407389e-02,
  1.788406944880994e-01,
  5.511070751355518e-02,
  2.112510229573667e-01,
  9.4475498594666e-02,
  2.765124421822678e-01,
  1.443375672974064e-01
};

const std::vector<double> IntegrationDPhiDY_true = {
  0.00000000000000e+00,
  2.383974596215561,
  0.00000000000000e+00,
  -1.744016935856293,
  4.166666666666668e-01,
  5.639957660359269e-01,
  2.337820631441114,
  -3.00000000000000e-01,
  3.443590455313351e-01,
  2.00000000000000e-01,
  -1.879743494847109,
  3.50000000000000e-01,
  1.913828550551446e-01,
  1.00000000000000e-01,
  3.851221159381804e-01,
  2.252991532071854,
  -3.238095238095238e-01,
  2.131929837166789e-01,
  5.714285714285716e-02,
  2.453297042212707e-01,
  1.428571428571429e-01,
  -1.937912020128888,
  3.392857142857143e-01,
  1.178190458792872e-01,
  3.571428571428573e-02,
  1.83043385083546e-01,
  8.92857142857143e-02,
  2.855387752616669e-01,
  2.210576982387224,
  -3.293650793650794e-01,
  1.641331427947897e-01,
  2.380952380952382e-02,
  1.483567763383696e-01,
  5.952380952380954e-02,
  2.021367719068524e-01,
  1.111111111111111e-01,
  -1.970227867507653,
  3.361111111111112e-01,
  8.21387896549881e-02,
  1.666666666666667e-02,
  1.26044828084989e-01,
  4.166666666666668e-02,
  1.569447932368951e-01,
  7.777777777777779e-02,
  2.272196374916589e-01,
  2.185128252576446,
  -3.313131313131313e-01,
  1.367871463596913e-01,
  1.212121212121213e-02,
  1.102706931333926e-01,
  3.030303030303032e-02,
  1.293162666349232e-01,
  5.656565656565658e-02,
  1.71701644642061e-01,
  9.09090909090909e-02,
  -1.990792497657777,
  3.348484848484849e-01,
  6.041767889566378e-02,
  9.0909090909091e-03
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
    f_vals[gp] = std::pow(x[gp]+1.,2.)+std::pow(y[gp]-1.,2.);
  }

  for(uint dof=0;dof<66;dof++){
    if(!almost_equal(IntegrationPhi_true[dof],triangle.IntegrationPhi(dof,f_vals),1.e+04)) {
            error_found = true;
            std::cerr << "Error found in Tringle element in IntegrationPhi" << 
                         " - integration true value: " << IntegrationPhi_true[dof] <<
                         ", integration computed value: " << triangle.IntegrationPhi(dof,f_vals) <<
            std::endl;
    }
  }

  for(uint dof=0;dof<66;dof++){
    if(!almost_equal(IntegrationDPhiDX_true[dof],triangle.IntegrationDPhi(GlobalCoord::x,dof,f_vals),1.e+04)) {
            error_found = true;
            std::cerr << "Error found in Tringle element in IntegrationDPhi in x direction" << 
                         " - integration true value: " << IntegrationDPhiDX_true[dof] <<
                         ", integration computed value: " << triangle.IntegrationDPhi(GlobalCoord::x,dof,f_vals) <<
            std::endl;
    }
  }

  for(uint dof=0;dof<59;dof++){
    if(!almost_equal(IntegrationDPhiDY_true[dof],triangle.IntegrationDPhi(GlobalCoord::y,dof,f_vals),1.e+04)) {
            error_found = true;
            std::cerr << "Error found in Tringle element in IntegrationDPhi in y direction" << 
                         " - integration true value: " << IntegrationDPhiDY_true[dof] <<
                         ", integration computed value: " << triangle.IntegrationDPhi(GlobalCoord::y,dof,f_vals) <<
            std::endl;
    }
  }

  std::vector<double> mod_vals(triangle.data.get_ndof());
  std::vector<double> gp_vals(triangle.data.get_ngp_internal());

  for(uint dof=0;dof<66;dof++){
    std::fill(mod_vals.begin(), mod_vals.end(),0.0);
    mod_vals[dof] = 1.0;

    triangle.ComputeUgp(mod_vals, gp_vals);
    
    for(uint doff=0;doff<66;doff++){
      if(dof == doff){
        printf("%f\n", triangle.IntegrationPhi(doff, gp_vals));
      }
      else {
        if(!almost_equal(triangle.IntegrationPhi(doff, gp_vals), 0.0)) {
          printf("%f\n", triangle.IntegrationPhi(doff, gp_vals));       
        }
      }
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
