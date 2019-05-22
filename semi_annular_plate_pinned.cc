//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//    Version 1.0; svn revision $LastChangedRevision: 1097 $
//LIC//
//LIC// $LastChangedDate: 2015-12-17 11:53:17 +0000 (Thu, 17 Dec 2015) $
//LIC// 
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================
#include <fenv.h> 
//Generic routines
#include "generic.h" 

// The equations
#include "C1_linear_plate_bending.h"

// The mesh
#include "meshes/triangle_mesh.h"

using namespace std;
using namespace oomph;
using MathematicalConstants::Pi;

namespace TestSoln
{
//Shape of the domain
double A = 1.0;
double B = 1.0;
double inner_radius=0.5;
double max_angle=Pi/2.;
double A_h = inner_radius;
double B_h = inner_radius;
// The coupling of the stretching energy
double eta = 0;
double p_mag = 1; 
const double nu = 0.0;

/*                     PARAMETRIC BOUNDARY DEFINITIONS                        */
// Here we create the geom objects for the Parametric Boundary Definition 
CurvilineCircleTop parametric_curve_top;
CurvilineCircleTop parametric_curve_bottom(inner_radius,true);

// Assigns the value of pressure depending on the position (x,y)
void get_pressure(const Vector<double>& xi, double& pressure)
{
 const double x=xi[0],y=xi[1];
 //const double theta =(x<0 ? atan2(-x,y)+Pi/2 : atan2(y,x)), r = sqrt(x*x+y*y);
 const double theta =(atan2(y,x)), r = sqrt(x*x+y*y);
 // Pressure
// pressure = (756 - 15*r*(495 - 560*r + 288*pow(r,3)))/
//    (70.*std::exp(sqrt(11)*theta)*pow(r,3));
//
 pressure = (756 - 15*r*(495 - 560*r + 288*pow(r,3)))/
   (70.*std::exp(sqrt(11)*theta)*pow(r,3));

}

// Assigns the value of pressure depending on the position (x,y)
void get_pressure(const Vector<double>& xi,Vector<double>& pressure)
{
 get_pressure(xi,pressure[0]);
}

// The normal and tangential directions.
void get_normal_and_tangent(const Vector<double>& x, Vector<double>& t, 
 Vector<double>& n, DenseMatrix<double>& Dt, DenseMatrix<double>& Dn)
{
 // Fill in the normal and derivatives of the normal
 n[0] = x[0]/sqrt(x[0]*x[0]+x[1]*x[1]);
 n[1] = x[1]/sqrt(x[0]*x[0]+x[1]*x[1]);

 // The derivatives of the x and y components
 Dn(0,0) = x[1]*x[1] * pow(x[0]*x[0]+x[1]*x[1],-1.5);
 Dn(1,0) =-x[1]*x[0] * pow(x[0]*x[0]+x[1]*x[1],-1.5);
 Dn(0,1) =-x[0]*x[1] * pow(x[0]*x[0]+x[1]*x[1],-1.5);
 Dn(1,1) = x[0]*x[0] * pow(x[0]*x[0]+x[1]*x[1],-1.5);

  // Fill in the tangent and derivatives of the tangent
 t[0] =-x[1]/sqrt(x[0]*x[0]+x[1]*x[1]);
 t[1] = x[0]/sqrt(x[0]*x[0]+x[1]*x[1]);

 Dt(0,0) =-Dn(1,0);
 Dt(1,0) = Dn(0,0); 
 Dt(0,1) =-Dn(1,1);
 Dt(1,1) = Dn(0,1);

  // Fill in the tangent and derivatives of the tangent
 t[0] =-x[1]/(x[0]*x[0]+x[1]*x[1]);
 t[1] = x[0]/(x[0]*x[0]+x[1]*x[1]);

 // The derivatives of the x and y components
 Dt(0,0)=(+2*x[0]*x[1])/pow(x[0]*x[0] + x[1]*x[1],2);
 Dt(0,1)=(-x[0]*x[0] + x[1]*x[1])/pow(x[0]*x[0] + x[1]*x[1],2);
 Dt(1,0)=(-x[0]*x[0] + x[1]*x[1])/pow(x[0]*x[0] + x[1]*x[1],2);
 Dt(1,1)=(-2*x[0]*x[1])/pow(x[0]*x[0] + x[1]*x[1],2);
}


 //Exact solution for constant pressure, circular domain and resting boundary conditions
 void get_exact_polar_w(const Vector<double>& ri, Vector<double>& w)
 {
  // Solution : w =  r exp[-\theta]
  const double phi = ri[1], r = ri[0];
  using namespace std;
  w[0] = -exp(-sqrt(11)*phi) *r*(-21. + 180.*r - 140.*r*r +24.*pow(r,4))/280.;
  w[1] = -exp(-sqrt(11)*phi) *(-21 + 360*r - 420*r*r +120*pow(r,4))/280.;
  w[2] =  r*exp(-sqrt(11)*phi)*sqrt(11)*(-21 + 180*r - 140*r*r + 24*pow(r,4))/280.;
  w[3] = -exp(-sqrt(11)*phi) *(360 - 840*r +480*pow(r,3))/280.;
  w[4] =  exp(-sqrt(11)*phi)*sqrt(11)*(-21 + 360*r - 420*r*r +120*pow(r,4))/280.;
  w[5] = -11*exp(-sqrt(11)*phi)*r*(-21 + 180*r - 140*r*r +24*pow(r,4))/(280.);
 }


//Exact solution for constant pressure, circular domain and resting boundary conditions
void get_exact_w(const Vector<double>& xi, Vector<double>& w)
{
 const double x=xi[0],y=xi[1];
 Vector<double> ri (2),wp(6);
 ri[1]=(atan2(y,x)); 
 ri[0] = sqrt(x*x+y*y);
 w.resize(6);
 get_exact_polar_w(ri,wp);
 // Now dtheta / dx etc.
 const double r2 = x*x+ y*y,  dtheta_dx = - y / r2 , dtheta_dy = x / r2, 
   d2theta_dx2 = 2*x*y / (r2*r2), d2theta_dxdy = (y*y-x*x) / (r2*r2),
   d2theta_dy2 =-2*x*y / (r2*r2); 

 const double r = sqrt(r2),  dr_dx =  x / r , dr_dy = y / r, 
   d2r_dx2 = y*y / (r2*r), d2r_dxdy = -x*y / (r2*r),
   d2r_dy2 = x*x / (r2*r); 
 // Now construct
 w[0] = wp[0];
 w[1] = wp[1]* dr_dx + wp[2]* dtheta_dx;
 w[2] = wp[1]* dr_dy + wp[2]* dtheta_dy;
 w[3] = wp[3]* pow(dr_dx,2) + 2*wp[4]* dr_dx * dtheta_dx + wp[5] * pow(dr_dy,2)
      + wp[1] * d2r_dx2 + wp[2] * d2theta_dx2;
 w[4] = wp[3]* dr_dx*dr_dy + wp[4]* (dr_dy * dtheta_dx + dr_dx * dtheta_dy) 
      + wp[5] * dtheta_dx * dtheta_dy + wp[1] * d2r_dxdy + wp[2] * d2theta_dxdy;
 w[5] = wp[3]* pow(dr_dy,2) + 2*wp[4]* dr_dy * dtheta_dy + wp[5] * pow(dr_dy,2)
      + wp[1] * d2r_dy2 + wp[2] * d2theta_dy2;
}

//Exact solution for constant pressure, circular domain and resting boundary conditions
void get_exact_w_on_straight_edge(const double r, Vector<double>& w, const unsigned& ibound)
{
 w.resize(6);
 if(ibound == 2)
   {   
   // Get the polar coodinate          
   Vector<double> ri(2);
   ri[1] = Pi/2, ri[0] = r;
   // Get the radial solution
   get_exact_polar_w(ri,w);
   // Tmp storage
   const Vector<double> wr = w ;
   // dw / dt =  -dw /  dx
   w[1] = -wr[2];
   // dw / dr =   dw /  dy
   w[2] =  wr[1];
   // d2w / dr2 =   d2w /  dy2
   w[3] = -wr[5]/r; /* we don't know the value of this! HERE */
   // d2w / dr dt =  -d2w /  dx dy
   w[4] = -wr[4]; 
   // d2w / dr2 =   d2w /  dy2
   w[5] = wr[3];
   }
 else if (ibound ==3) 
   {
   // Get the polar coodinate          
   Vector<double> ri(2);
   ri[1] = 0, ri[0] = r;
   // Get the radial solution
   get_exact_polar_w(ri,w);
   w[5] = w[5] / r;
   }
 else
  {
    throw OomphLibError(
     "I have encountered a boundary number that I wasn't expecting. This function\
 is only defined for straight edges.",
     "TestSoln:get_exact_w_on_straight_edge(...)",
     OOMPH_EXCEPTION_LOCATION);
  }
 // The last second derivative is not the second y derivative, so remove it
 //w.pop_back();
}

}

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////


//==start_of_problem_class============================================
/// Class definition
//====================================================================
template<class ELEMENT>
class UnstructuredFvKProblem : public virtual Problem
{

public:

/// Constructor
UnstructuredFvKProblem(double element_area = 0.1);

/// Destructor
~UnstructuredFvKProblem()
{
};

/// Update after solve (empty)
void actions_after_newton_solve()
{
}

/// Update the problem specs before solve: Re-apply boundary conditions
/// Empty as the boundary conditions stay fixed
void actions_before_newton_solve()
{
apply_boundary_conditions();
}

/// Doc the solution
void doc_solution(const std::string& comment="");

/// \short Overloaded version of the problem's access function to
/// the mesh. Recasts the pointer to the base Mesh object to
/// the actual mesh type.
TriangleMesh<ELEMENT>* mesh_pt()
{
return dynamic_cast<TriangleMesh<ELEMENT>*> (Problem::mesh_pt()); 
}

/// Doc info object for labeling output
DocInfo Doc_info;

private:


/// Helper function to apply boundary conditions
void apply_boundary_conditions();

/// \short Helper function to (re-)set boundary condition
/// and complete the build of  all elements
void complete_problem_setup();

/// Trace file to document norm of solution
ofstream Trace_file;

// Keep track of boundary ids
enum
{
 Outer_boundary0 = 0,
 Inner_boundary0 = 1,
 Straight_boundary0 =2,
 Straight_boundary1 =3,
 Internal_boundary0 =4,
 Internal_boundary1 =5,
 Internal_boundary2 =6,
 Internal_boundary3 =7
};

double Element_area;

void upgrade_edge_elements_to_curve(const unsigned &b, Mesh* const & 
 bulk_mesh_pt);

void rotate_edge_degrees_of_freedom(Mesh* const &bulk_mesh_pt);



/// Pointer to "bulk" mesh
TriangleMesh<ELEMENT>* Bulk_mesh_pt;

/// Pointer to "surface" mesh
Mesh* Surface_mesh_pt;

}; // end_of_problem_class


template<class ELEMENT>
UnstructuredFvKProblem<ELEMENT>::UnstructuredFvKProblem(double element_area)
:
Element_area(element_area)
{
Vector<double> zeta(1);
Vector<double> posn(2);

//Outer boundary
//--------------

Ellipse* outer_boundary_ellipse_pt = new Ellipse(TestSoln::A, TestSoln::B);
Ellipse* inner_boundary_ellipse_pt = new Ellipse(TestSoln::A_h, TestSoln::B_h);

// Create storage for curve sections
Vector<TriangleMeshCurveSection*> outer_curvilinear_boundary_pt(2);

//First bit
double zeta_start = 0.0;
double zeta_end = MathematicalConstants::Pi/2.;
unsigned nsegment = (int)(MathematicalConstants::Pi/(2*sqrt(element_area)));
outer_curvilinear_boundary_pt[0] =
new TriangleMeshCurviLine(outer_boundary_ellipse_pt, zeta_start,
zeta_end, nsegment, Outer_boundary0);

// Second bit
Vector<Vector<double> > vertices(3,Vector<double>(2,0.0));

vertices[0][1] = 1.0;

vertices[1][1] = 0.75;

vertices[2][1] = 0.5;

TriangleMeshPolyLine *boundary1_pt =
  new TriangleMeshPolyLine(vertices, Straight_boundary0);
 

// Create storage for curve sections
Vector<TriangleMeshCurveSection*> inner_curvilinear_boundary_pt(2);
//Third bit
zeta_start = MathematicalConstants::Pi/2.;
zeta_end = 0.0;
nsegment = (int)(TestSoln::inner_radius*MathematicalConstants::Pi/(2*sqrt(element_area)));
inner_curvilinear_boundary_pt[0] =
new TriangleMeshCurviLine(inner_boundary_ellipse_pt, zeta_start,
zeta_end, nsegment, Inner_boundary0);

// Fourth bit 
vertices[0][1] = 0.0;
vertices[0][0] = 0.5;

vertices[1][1] = 0.0;
vertices[1][0] = 0.75;

vertices[2][1] = 0.0;
vertices[2][0] = 1.0;

TriangleMeshPolyLine *boundary3_pt =
  new TriangleMeshPolyLine(vertices, Straight_boundary1);

 // Setting up the domain with PolyLines -> polymorphism
 Vector<TriangleMeshCurveSection*> outer_boundary_polyline_pt(4);
 outer_boundary_polyline_pt[0] = outer_curvilinear_boundary_pt[0];
 outer_boundary_polyline_pt[1] = boundary1_pt;
 outer_boundary_polyline_pt[2] = inner_curvilinear_boundary_pt[0];
 outer_boundary_polyline_pt[3] = boundary3_pt;

 vertices.resize(2);
 for (unsigned i = 0; i < 2; i++)
   {
     vertices[i].resize(2);
   }
 // First internal boundary
 vertices[0][0] = 1.0;
 vertices[0][1] = 0.0;

 vertices[1][0] =  0.75;
 vertices[1][1] =  0.25;

 TriangleMeshPolyLine *inner_boundary0_pt =
   new TriangleMeshPolyLine(vertices, Internal_boundary0);

 // Second internal boundary
 vertices[0][0] = 0.5;
 vertices[0][1] = 0.0;

 vertices[1][0] = 0.75;
 vertices[1][1] = 0.25;

 TriangleMeshPolyLine *inner_boundary1_pt =
   new TriangleMeshPolyLine(vertices, Internal_boundary1);

 // Third internal boundary
 vertices[0][0] = 0.0;
 vertices[0][1] = 0.5;

 vertices[1][0] = 0.25;
 vertices[1][1] = 0.75;

 TriangleMeshPolyLine *inner_boundary2_pt =
   new TriangleMeshPolyLine(vertices, Internal_boundary2);

 // Fourth internal boundary
 vertices[0][0] = 0.0;
 vertices[0][1] = 1.0;

 vertices[1][0] = 0.25;
 vertices[1][1] = 0.75;

 TriangleMeshPolyLine *inner_boundary3_pt =
   new TriangleMeshPolyLine(vertices, Internal_boundary3);

 // Now join them up
 inner_boundary0_pt->connect_initial_vertex_to_polyline(boundary3_pt,2);
 inner_boundary0_pt->connect_initial_vertex_to_curviline(
  dynamic_cast <TriangleMeshCurviLine*>( outer_curvilinear_boundary_pt[0]),0.0);
 inner_boundary1_pt->connect_initial_vertex_to_polyline(boundary3_pt,0);
 inner_boundary2_pt->connect_initial_vertex_to_polyline(boundary1_pt,0);
 inner_boundary3_pt->connect_initial_vertex_to_polyline(boundary1_pt,2);

 inner_boundary1_pt->connect_final_vertex_to_polyline(inner_boundary0_pt,1);
 inner_boundary3_pt->connect_final_vertex_to_polyline(inner_boundary2_pt,1);

 // The inner lines
 Vector<TriangleMeshCurveSection*> inner_boundary_polyline_pt(1);
 inner_boundary_polyline_pt[0] = inner_boundary0_pt;

 TriangleMeshOpenCurve* inner_open_curve_0_pt = 
  new TriangleMeshOpenCurve(inner_boundary_polyline_pt);

 inner_boundary_polyline_pt[0] = inner_boundary1_pt;

 TriangleMeshOpenCurve* inner_open_curve_1_pt = 
  new TriangleMeshOpenCurve(inner_boundary_polyline_pt);

 inner_boundary_polyline_pt[0] = inner_boundary2_pt;

 TriangleMeshOpenCurve* inner_open_curve_2_pt = 
  new TriangleMeshOpenCurve(inner_boundary_polyline_pt);

 inner_boundary_polyline_pt[0] = inner_boundary3_pt;

 TriangleMeshOpenCurve* inner_open_curve_3_pt = 
  new TriangleMeshOpenCurve(inner_boundary_polyline_pt);

 Vector<TriangleMeshOpenCurve*> internal_open_boundaries(4);
 internal_open_boundaries[0] = inner_open_curve_0_pt;
 internal_open_boundaries[1] = inner_open_curve_1_pt;
 internal_open_boundaries[2] = inner_open_curve_2_pt;
 internal_open_boundaries[3] = inner_open_curve_3_pt;

 // The outer shape
 TriangleMeshClosedCurve* outer_boundary_pt = 
  new TriangleMeshClosedCurve(outer_boundary_polyline_pt);

//Create the mesh
//---------------
//Create mesh parameters object
TriangleMeshParameters mesh_parameters(outer_boundary_pt);

mesh_parameters.internal_open_curves_pt() = internal_open_boundaries;
mesh_parameters.element_area() = element_area;

// Build an assign bulk mesh
Bulk_mesh_pt=new TriangleMesh<ELEMENT>(mesh_parameters);

// Create "surface mesh" that will contain only the prescribed-traction
// elements. The constructor creates the mesh without adding any nodes
// elements etc.
Surface_mesh_pt =  new Mesh;

//Add two submeshes to problem
add_sub_mesh(Bulk_mesh_pt);
add_sub_mesh(Surface_mesh_pt);

// Combine submeshes into a single Mesh
build_global_mesh();

// Curved Edge upgrade
upgrade_edge_elements_to_curve(0,Bulk_mesh_pt);
upgrade_edge_elements_to_curve(1,Bulk_mesh_pt);
 
// Rotate degrees of freedom
//rotate_edge_degrees_of_freedom(Bulk_mesh_pt);

// Store number of bulk elements
complete_problem_setup();

char filename[100];
sprintf(filename, "RESLT/trace.dat");
Trace_file.open(filename);

oomph_info << "Number of equations: "
        << assign_eqn_numbers() << '\n';
}



//==start_of_complete======================================================
/// Set boundary condition exactly, and complete the build of 
/// all elements
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::complete_problem_setup()
{   
// Complete the build of all elements so they are fully functional
unsigned n_element = Bulk_mesh_pt->nelement();
for(unsigned e=0;e<n_element;e++)
{
// Upcast from GeneralisedElement to the present element
ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

//Set the pressure function pointers and the physical constants
el_pt->pressure_fct_pt() = &TestSoln::get_pressure;
el_pt->nu_pt() = &TestSoln::nu;
}

// Re-apply Dirichlet boundary conditions (projection ignores
// boundary conditions!)
apply_boundary_conditions();
}

//==start_of_apply_bc=====================================================
/// Helper function to apply boundary conditions
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::apply_boundary_conditions()
{
unsigned nbound = Straight_boundary1+1;
// Set the boundary conditions for problem: All nodes are
// free by default -- just pin the ones that have Dirichlet conditions
// here. 
//Just loop over outer boundary since inner boundary doesn't have boundary
//conditions
for(unsigned ibound= Straight_boundary0;ibound<nbound;ibound++)
{
oomph_info<<"Pinning boundary: "<<ibound <<"\n";
unsigned num_nod=Bulk_mesh_pt->nboundary_node(ibound);
oomph_info<<"Pinning : "<<num_nod <<" nodes\n";
for (unsigned inod=0;inod<num_nod;inod++)
{
 // Get nod
 Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
 Vector<double> x(2,0.0),w(6,0.0);
 x[0]=nod_pt->x(0);
 x[1]=nod_pt->x(1);
 oomph_info << "Pinning node at x="<< x <<"\n";
 //const double r = (ibound ==2 ? x[1] : x[0]) ;
 // TestSoln::get_exact_w_on_straight_edge(r,w, ibound);
 TestSoln::get_exact_w(x,w);
 oomph_info << "To values x="<< w <<"\n";
// // Pin unknown values (everything except for the second normal derivative)
 for(unsigned i =0;i<6;++i)
  {
   if(ibound == 2 && i==3)
    { /* Do nothing */ }
   else if(ibound == 3 && i==5)
    { /* Do nothing */ }
   else 
    {
     nod_pt->pin(i);
     nod_pt->set_value(i,w[i]);
    }
  }
 }
} // end loop over boundaries 

} // end set bc


template <class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::
upgrade_edge_elements_to_curve(const unsigned &b, Mesh* const &bulk_mesh_pt) 
{
 // How many bulk elements adjacent to boundary b
 unsigned n_element = bulk_mesh_pt-> nboundary_element(b);
 
 // These depend on the boundary we are on
 CurvilineGeomObject* parametric_curve_pt; 
 
// Define the functions for each part of the boundary
 switch (b)
  {
   // Upper boundary
   case 0:
    parametric_curve_pt = &TestSoln::parametric_curve_top;
   break;

   // Lower boundary
   case 1:
    parametric_curve_pt = &TestSoln::parametric_curve_bottom;
   break;

   default:
    throw OomphLibError(
     "I have encountered a boundary number that I wasn't expecting. This is very\
 peculiar.",
     "UnstructuredFvKProblem::upgrade_edge_elements_to_curve(...)",
     OOMPH_EXCEPTION_LOCATION);
   break;
  }
 
 // Loop over the bulk elements adjacent to boundary b
 for(unsigned e=0; e<n_element; e++)
  {
   // Get pointer to bulk element adjacent to b
   ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(
    bulk_mesh_pt->boundary_element_pt(b,e));
   
   // Loop over nodes
   unsigned nnode=bulk_el_pt->nnode();
   unsigned index_of_interior_node=3;

   // Enum for the curved edge
   MyC1CurvedElements::Edge edge(MyC1CurvedElements::none);

   // Vertices positions
   Vector<Vector<double> > xn(3,Vector<double>(2,0.0));
 
   // Get vertices for debugging
   Vector<Vector<double> > verts(3,Vector<double>(2,0.0));
   // Loop nodes
   for(unsigned n=0;n<nnode;++n)
    {
     // If it is on boundary
     Node* nod_pt = bulk_el_pt->node_pt(n);
     verts[n][0]=nod_pt->x(0);
     verts[n][1]=nod_pt->x(1);
     if(nod_pt->is_on_boundary(b))
      {
       xn[n][0]=nod_pt->x(0);
       xn[n][1]=nod_pt->x(1);
      }
     // The edge is denoted by the index of the  opposite (interior) node
     else {index_of_interior_node = n;}
    }
   // Initialise s_ubar s_obar
   double s_ubar, s_obar;

   // s at the next (cyclic) node after interior
   s_ubar = parametric_curve_pt->get_zeta(xn[(index_of_interior_node+1) % 3]);
   // s at the previous (cyclic) node before interior
   s_obar = parametric_curve_pt->get_zeta(xn[(index_of_interior_node+2) % 3]);

   // Assign edge case
   switch(index_of_interior_node)
    {
     case 0: edge= MyC1CurvedElements::zero; 
      break;
     case 1: edge= MyC1CurvedElements::one; 
      break;
     case 2: edge= MyC1CurvedElements::two; 
      break;
     // Should break it here HERE
     default: edge= MyC1CurvedElements::none; 
      throw OomphLibError(
       "The edge number has been set to a value greater than two: either we have\
 quadrilateral elements or more likely the index_of_interior_node was never set\
 and remains at its default value.",
       "UnstructuredFvKProblem::upgrade_edge_elements_to_curve(...)",
       OOMPH_EXCEPTION_LOCATION);
      break;
     }
   if (s_ubar>s_obar)
    {
     oomph_info <<"Apparent clockwise direction of parametric coordinate."
                <<"This will probably result in an inverted element."
                <<"s_start="<<s_ubar<<"; s_end ="<<s_obar<<std::endl;
     throw OomphLibError(
       "The Edge coordinate appears to be decreasing from s_start to s_end. \
Either the parametric boundary is defined to be clockwise (a no-no) or \
the mesh has returned an inverted element (less likely)",
       "UnstructuredFvKProblem::upgrade_edge_elements_to_curve(...)",
       OOMPH_EXCEPTION_LOCATION);
    }

   // Upgrade it
   bulk_el_pt->upgrade_element_to_curved(edge,s_ubar,s_obar,
    parametric_curve_pt,5);     
    
  // // Get vertices for debugging
  // Vector<Vector<double> > lverts(3,Vector<double>(2,0.0));
  // lverts[0][0]=1.0;
  // lverts[1][1]=1.0;
  // Vector<Vector<double> > fkverts(3,Vector<double>(2,0.0));
  // bulk_el_pt->interpolated_x(lverts[0],fkverts[0]);
  // bulk_el_pt->interpolated_x(lverts[1],fkverts[1]);
  // bulk_el_pt->interpolated_x(lverts[2],fkverts[2]);

  }
}// end upgrade elements

template <class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::
rotate_edge_degrees_of_freedom( Mesh* const &bulk_mesh_pt)
{
 // How many bulk elements
 unsigned n_element = bulk_mesh_pt-> nelement();
 
 // Loop over the bulk elements adjacent to boundary b
 for(unsigned e=0; e<n_element; e++)
  {
   // Get pointer to bulk element adjacent to b
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
 
   // Loop nodes  
   unsigned nnode = el_pt->nnode();
   unsigned nbnode=0 ;
   // Count the number of boundary nodes
   for (unsigned n=0; n<nnode;++n)
     {nbnode+=unsigned(el_pt->node_pt(n)->is_on_boundary(Straight_boundary0) ||
       el_pt->node_pt(n)->is_on_boundary(Straight_boundary1));}

   // Now if we have nodes on boundary 
   if(nbnode>0)
    {
     // Set up vector
     Vector<unsigned> bnode (nbnode,0);
     unsigned inode(0);

     // Fill in the bnode Vector
     for (unsigned n=0; n<nnode;++n)
      {
       // If it is on the boundary
       if(el_pt->node_pt(n)->is_on_boundary(Straight_boundary0) ||
       el_pt->node_pt(n)->is_on_boundary(Straight_boundary1))
        {
         // Set up the Vector
         bnode[inode]=n;
         ++inode;
        }
      }
    // Output that we have found element HERE
    std::cout<<"Element "<<e<<" has "<<bnode<< " nodes on the boundary.\n";
    std::cout<<"at:\n";
    for(unsigned i=0;i<nbnode;++i)
     std::cout<<"("<<el_pt->node_pt(bnode[i])->x(0)<<","<<el_pt->node_pt(bnode[i])->x(1)<<")\n";
    el_pt->set_up_rotated_dofs(nbnode,bnode,&TestSoln::get_normal_and_tangent);
   // Now rotate the nodes
   }
 }
}// end create traction elements

//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::doc_solution(const 
                                                    std::string& comment)
{ 
ofstream some_file;
char filename[100];

// Number of plot points
unsigned npts = 2;

sprintf(filename,"RESLT/coarse_soln%i-%f.dat",Doc_info.number(),Element_area);
some_file.open(filename);
Bulk_mesh_pt->output(some_file,npts); 
some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \"" 
       << comment << "\"\n";
some_file.close();

// Number of plot points
npts = 25;

sprintf(filename,"RESLT/soln%i-%f.dat",Doc_info.number(),Element_area);
some_file.open(filename);
Bulk_mesh_pt->output(some_file,npts); 
some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \"" 
       << comment << "\"\n";
some_file.close();

//  Output exact solution
sprintf(filename,"%s/exact_interpolated_soln%i-%f.dat","RESLT",Doc_info.number(),Element_area);
some_file.open(filename);
Bulk_mesh_pt->output_fct(some_file,npts,TestSoln::get_exact_w); 
some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \"" 
       << comment << "\"\n";
some_file.close();

//  Output exact solution
sprintf(filename,"%s/pressure_soln%i-%f.dat","RESLT",Doc_info.number(),Element_area);
some_file.open(filename);
Bulk_mesh_pt->output_fct(some_file,npts,TestSoln::get_pressure); 
some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \"" 
       << comment << "\"\n";
some_file.close();

// Output boundaries
//------------------
sprintf(filename,"RESLT/boundaries%i-%f.dat",Doc_info.number(),Element_area);
some_file.open(filename);
Bulk_mesh_pt->output_boundaries(some_file);
some_file.close();

// Output regions
unsigned n_region = Bulk_mesh_pt->nregion();
if (n_region > 1)
{
for (unsigned r = 0; r < n_region; r++)
{
 //Attempt to output elements in different regions
 sprintf(filename,"RESLT/region%i%i-%f.dat",r,Doc_info.number(),
   Element_area);
 some_file.open(filename);
 unsigned nel = Bulk_mesh_pt->nregion_element(r);
 for (unsigned e = 0; e < nel; e++)
  {
   Bulk_mesh_pt->region_element_pt(r,e)->output(some_file,npts);
  }
 some_file.close();
}
}

// // Doc error and return of the square of the L2 error
// //---------------------------------------------------
// //double error,norm,dummy_error,zero_norm;
  double dummy_error,zero_norm;
 sprintf(filename,"RESLT/error%i-%f.dat",Doc_info.number(),Element_area);
 some_file.open(filename);
 
 Bulk_mesh_pt->compute_error(some_file,TestSoln::get_exact_w,
                        dummy_error,zero_norm);
 some_file.close();
 
 // Doc L2 error and norm of solution
 oomph_info << "Norm of computed solution: " << sqrt(dummy_error)<< std::endl;
 
 Trace_file << TestSoln::p_mag << " " << "\n ";

// Doc error and return of the square of the L2 error
//---------------------------------------------------
sprintf(filename,"RESLT/L2-norm%i-%f.dat",
        Doc_info.number(),
        Element_area);
some_file.open(filename);

some_file<<"### L2 Norm\n";
some_file<<"##  Format: err^2 norm^2 log(err/norm) \n";
// Print error in prescribed format
some_file<< dummy_error <<" "<< zero_norm <<" ";

// Only divide by norm if its nonzero
some_file<<0.5*(log10(fabs(dummy_error))-log10(zero_norm))<<"\n";
some_file.close();

// Increment the doc_info number
Doc_info.number()++;

} // end of doc


//=======start_of_main========================================
///Driver code for demo of inline triangle mesh generation
//============================================================
int main(int argc, char **argv)
{

 feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define possible command line arguments and parse the ones that
 // were actually specified

 // Validation?
 CommandLineArgs::specify_command_line_flag("--validation");

 // Directory for solution
 string output_dir="RSLT";
 CommandLineArgs::specify_command_line_flag("--dir", &output_dir);

 // Applied Pressure
 CommandLineArgs::specify_command_line_flag("--p", &TestSoln::p_mag);
 
 // P_step
 double p_step=10;
 CommandLineArgs::specify_command_line_flag("--dp", &p_step);
 
 // Element Area
 double element_area=0.2;
 CommandLineArgs::specify_command_line_flag("--element_area", &element_area);

 // Parse command line
 CommandLineArgs::parse_and_assign(); 

 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Problem instance
 UnstructuredFvKProblem<KirchhoffPlateBendingC1CurvedBellElement<2,2,5> >problem(element_area);
 problem.doc_solution();
  problem.self_test();
  problem.max_newton_iterations()=1;
// problem.newton_solver_tolerance()=1e12;
  problem.max_residuals()=1e3;
  oomph_info<<"Solving for p=" << TestSoln::p_mag<<"\n";
//  problem.doc_solution();
  problem.newton_solve();

  // Document
  problem.doc_solution();
  oomph_info << std::endl;
  oomph_info << "---------------------------------------------" << std::endl;
  oomph_info << " Pcos (" << TestSoln::p_mag << ")" << std::endl;
  oomph_info << "Current dp  (" << p_step << ")" << std::endl;
  oomph_info << "Poisson ratio (" << TestSoln::nu << ")" << std::endl;
  oomph_info << "Solution number (" <<problem.Doc_info.number()-1 << ")" << std::endl;
  oomph_info << "---------------------------------------------" << std::endl;
  oomph_info << std::endl;

  // Dump the data 
  char filename[100];
  std::ofstream filestream;
  filestream.precision(15);
  sprintf(filename,"%s/fvk_circle_data%i-%f.dump",
          problem.Doc_info.directory().c_str(),
          problem.Doc_info.number(),
          TestSoln::p_mag
         );
  filestream.open(filename);
  problem.dump(filestream);
  filestream.close();

} //End of main

