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

// Include the library elements
#include "C1_linear_plate_bending.h"

// The mesh
#include "meshes/triangle_mesh.h"

// Analytic Solutions for linear plate bending problems
#include "kirchhoff_plate_bending_analytic_solutions.h"

using namespace std;
using namespace oomph;
using MathematicalConstants::Pi;

namespace TestSoln
{
//Shape of the domain
double A = 1.0;
double B = 1.0;
// The coupling of the stretching energy
double eta = 0;
double p_mag = 1; 
double nu = 1/3.;

enum Boundary_case {Resting = 0,  Clamped = 1, Free  = 2, FreeSinusoidalLoad =3 };
Boundary_case boundary_case = Resting;

/*                     PARAMETRIC BOUNDARY DEFINITIONS                        */
// Here we create the geom objects for the Parametric Boundary Definition 
CurvilineCircleTop parametric_curve_top;
CurvilineCircleBottom parametric_curve_bottom;

// Assigns the value of pressure depending on the position (x,y)
void get_pressure(const Vector<double>& x, double& pressure)
{
 switch (boundary_case)
  {  
  // Uniform Pressure
  case Resting: case Clamped:
   pressure = p_mag;
   break;
  case Free:
   { 
   double r2 = x[0]*x[0]+x[1]*x[1];
   pressure = p_mag*(r2-0.5);
   break;
   }
  case FreeSinusoidalLoad: 
   { 
   const double theta = atan2(x[1],x[0]);
   pressure = p_mag*sin(3*theta);
   break;
   }
  default:
   /* Scream */
  break;
  }
}

// The normal and tangential directions.
void get_normal_and_tangent(const Vector<double>& x, Vector<double>& n, 
 Vector<double>& t, DenseMatrix<double>& Dn, DenseMatrix<double>& Dt)
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
}

//Exact solution for constant pressure, circular domain and resting boundary conditions
void get_exact_w(const Vector<double>& x, Vector<double>& w)
{
// Exact w depends on boundary case
switch (boundary_case)
  {
  case Resting:
   // See Timoshenko
   w[0]= p_mag*(1-x[0]*x[0]-x[1]*x[1])*((5.+nu)/(1.+nu)-x[0]*x[0]-x[1]*x[1])/64.;
  break;
  case Clamped:
   // See Timoshenko
   w[0]= p_mag*pow(1-x[0]*x[0]-x[1]*x[1],2)/64.;
  break;
  case Free:
   {
   const double r2 = x[0]*x[0] +x[1]*x[1];
   // Solved by DGR 2019 (almost certainly exists elsewhere) 
   w[0]= (p_mag*r2*(-9*r2 + 2*pow(r2,2) + 
       (12*(2 + nu))/(1 + nu)))/1152.;
   }
  case FreeSinusoidalLoad:
   {
   // Define r and theta
   const double r = sqrt(x[0]*x[0] +x[1]*x[1]), theta = atan2(x[1],x[0]);
   // Alias
   const double& Nu = nu;
   double (*Power)(double base, int exponent) = &std::pow;
   double (*Sin)(double arg) = &std::sin;
   
   // The solution for k =3 nu=?
   w[0] = p_mag*(Power(r,3)*(-202 + 2*Nu*(55 + 18*Nu) - 72*(-1 + Nu)*(3 + Nu)*r + 
       3*(-1 + Nu)*(23 + 12*Nu)*Power(r,2))*Sin(3*theta))/
   (2520.*(-1 + Nu)*(3 + Nu));
   /* 
   double (*Cos)(double arg) = &std::cos;
   // r - theta derivatives for k=3 nu = 1/3
   w[0] = (Power(r,3)*(242 + 3*r*(-80 + 27*r))*Sin(3*theta))/8400.; 
   w[1] = (Power(r,2)*(242 - 320*r + 135*Power(r,2))*Sin(3*theta))/2800.;
   w[2] = (Power(r,3)*(242 + 3*r*(-80 + 27*r))*Cos(3*theta))/2800.;
   w[3] = (r*(121 - 240*r + 135*Power(r,2))*Sin(3*theta))/700.;
   w[4] = (3*Power(r,2)*(242 - 320*r + 135*Power(r,2))*Cos(3*theta))/2800;
   w[5] = (-3*Power(r,3)*(242 + 3*r*(-80 + 27*r))*Sin(3*theta))/2800.;
   */
   }
  break;
  default:
    throw OomphLibError(
     "I have encountered a boundary case that I wasn't expecting..",
     "UnstructuredFvKProblem::get_exact_w(...)",
     OOMPH_EXCEPTION_LOCATION);
  break;
 }
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
 delete (Surface_mesh_pt);
 delete (Bulk_mesh_pt);
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
 Outer_boundary1 = 1,
 Inner_boundary0 = 2
};

double Element_area;

// The extra members required for flux type boundary conditions.
/// \short Number of "bulk" elements (We're attaching flux elements
/// to bulk mesh --> only the first Nkirchhoff_elements elements in
/// the mesh are bulk elements!)
// unsigned Nkirchhoff_elements;

/// \short Create bending moment elements on the b-th boundary of the
/// problems mesh 
void create_traction_elements(const unsigned &b, Mesh* const & bulk_mesh_py,
                            Mesh* const &surface_mesh_pt);

void upgrade_edge_elements_to_curve(const unsigned &b, Mesh* const & 
 bulk_mesh_pt);

void rotate_edge_degrees_of_freedom(Mesh* const &bulk_mesh_pt);

/// Helper to pin the displacement at the centre
void pin_displacement_at_centre_node();

/// Helper to pin the azimuthal dof on the edge
void pin_rotation_on_edge_node();

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

double A = 1.0;
double B = 1.0;
Ellipse* outer_boundary_ellipse_pt = new Ellipse(A, B);

TriangleMeshClosedCurve* outer_boundary_pt = 0;

Vector<TriangleMeshCurveSection*> outer_curvilinear_boundary_pt(2);

//First bit
double zeta_start = 0.0;
double zeta_end = MathematicalConstants::Pi;
unsigned nsegment = (int)(MathematicalConstants::Pi/sqrt(element_area));
outer_curvilinear_boundary_pt[0] =
new TriangleMeshCurviLine(outer_boundary_ellipse_pt, zeta_start,
zeta_end, nsegment, Outer_boundary0);

//Second bit
zeta_start = MathematicalConstants::Pi;
zeta_end = 2.0*MathematicalConstants::Pi;
nsegment = (int)(MathematicalConstants::Pi/sqrt(element_area));
outer_curvilinear_boundary_pt[1] =
new TriangleMeshCurviLine(outer_boundary_ellipse_pt, zeta_start,
zeta_end, nsegment, Outer_boundary1);

outer_boundary_pt =
new TriangleMeshClosedCurve(outer_curvilinear_boundary_pt);

// Internal bit - this means we can have a boundary which is just the centre
// We start by creating the internal boundaries
// The boundary 2 is defined by its two vertices
// Open curve 1
Vector<Vector<double> > vertices(2,Vector<double>(2,0.0));
vertices[0][0] = 1.0;

unsigned boundary_id = Inner_boundary0;

TriangleMeshPolyLine *boundary2_pt =
  new TriangleMeshPolyLine(vertices, boundary_id);

// Total number of open curves in the domain
unsigned n_open_curves = 1;
// We want internal open curves
Vector<TriangleMeshOpenCurve *> inner_open_boundaries_pt(n_open_curves);
// Connect it
boundary2_pt -> connect_initial_vertex_to_curviline(dynamic_cast<TriangleMeshCurviLine*>(outer_curvilinear_boundary_pt[0]),0.0);

// Each internal open curve is defined by a vector of
// TriangleMeshCurveSection,
// on this example we only need one curve section for each internal boundary
Vector<TriangleMeshCurveSection *> internal_curve_section1_pt(1);
internal_curve_section1_pt[0] = boundary2_pt;

// The open curve that define this boundary is composed of just one
// curve section
 inner_open_boundaries_pt[0] =
    new TriangleMeshOpenCurve(internal_curve_section1_pt);

//Create the mesh
//---------------
//Create mesh parameters object
TriangleMeshParameters mesh_parameters(outer_boundary_pt);

mesh_parameters.element_area() = element_area;

// Specify the internal open boundaries
mesh_parameters.internal_open_curves_pt() = inner_open_boundaries_pt;

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
rotate_edge_degrees_of_freedom(Bulk_mesh_pt);

// Store number of bulk elements
complete_problem_setup();

char filename[100];
sprintf(filename, "RESLT/trace.dat");
Trace_file.open(filename);

oomph_info << "Number of equations: "
        << assign_eqn_numbers() << '\n';

// Clean up memory - this may create dangling pointers but we don't need to
// So move the relevant pointers to private data members and clean this in
// destructor
delete outer_boundary_pt;
delete outer_boundary_ellipse_pt;
delete outer_curvilinear_boundary_pt[0];
delete outer_curvilinear_boundary_pt[1];
delete inner_open_boundaries_pt[0];
delete boundary2_pt;
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

// Apply Dirichlet boundary conditions (projection ignores
apply_boundary_conditions();
}


template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::pin_displacement_at_centre_node()
{
// Pin the node that is at the centre in the domain!
// Get the num of nods on internal_boundary 0
unsigned num_int_nod=Bulk_mesh_pt->nboundary_node(2);
bool found_centre_node = false;
for (unsigned inod=0;inod<num_int_nod;inod++)
{
 // Get node point
 Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(2,inod);
 // If the node is on the other internal boundary too
 if(nod_pt->x(0)==0.0 && nod_pt->x(1) == 0.0)
 {
  found_centre_node = true;
  // Be verbose and let everyone know we have found the centre
  oomph_info<<"Found centre point\n";
  oomph_info<<"x = ("<<nod_pt->x(0)<<" "<<nod_pt->x(1)<<") \n";
  // Pin it! It's the centre of the domain!
  nod_pt->pin(0);
  nod_pt->set_value(0,0.0);
  nod_pt->pin(1);
  nod_pt->set_value(1,0.0);
  nod_pt->pin(2);
  nod_pt->set_value(2,0.0);
  return;
 }
}
// This is bad
if(!found_centre_node)
 {
    throw OomphLibError(
     "I couldn't find the centre node to pin it! Exiting.",
     "UnstructuredFvKProblem::pin_displacement_on_centre_node(...)",
     OOMPH_EXCEPTION_LOCATION);
 }
}

template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::pin_rotation_on_edge_node()
{
// Pin the node that is at the centre in the domain!
// Get the num of nods on internal_boundary 0
unsigned num_int_nod=Bulk_mesh_pt->nboundary_node(1);
bool found_edge_node = false;
for (unsigned inod=0;inod<num_int_nod;inod++)
{
 // Get node point
 Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(1,inod);
 // If the node is on the internal boundary too
 if(nod_pt->is_on_boundary(2))
 {
  found_edge_node = true;
  // Be verbose and let everyone know we have found the centre
  oomph_info<<"Found edge point\n";
  oomph_info<<"x = ("<<nod_pt->x(0)<<" "<<nod_pt->x(1)<<") \n";
  // Pin the theta derivative (i.e y derivative at x = 1, y =0)
  nod_pt->pin(2);
  nod_pt->set_value(2,0.0);
  return;
 }
}
// This is bad
if(!found_edge_node)
 {
    throw OomphLibError(
     "I couldn't find the edge node to pin it! Exiting.",
     "UnstructuredFvKProblem::pin_rotation_on_edge_node(...)",
     OOMPH_EXCEPTION_LOCATION);
 }
}

//==start_of_apply_bc=====================================================
/// Helper function to apply boundary conditions
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::apply_boundary_conditions()
{
unsigned nbound = Outer_boundary1 + 1;
// Set the boundary conditions for problem: All nodes are
// free by default -- just pin the ones that have Dirichlet conditions
// here. 
//Just loop over outer boundary since inner boundary doesn't have boundary
//conditions
for(unsigned ibound=0;ibound<nbound;ibound++)
{
unsigned num_nod=Bulk_mesh_pt->nboundary_node(ibound);
for (unsigned inod=0;inod<num_nod;inod++)
{
 // Get nod
 Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
 /* Get x here if needed for b/c*/

 // Bell and Curved Bell have six dof types on each vertex node 
 // They are by default [ w w,x w,y w,xx w,xy w,yy ]
 // BUT if we rotate the Hermite degrees of freedom we replace them with the
 // specified normal / tangent degrees of freedom such that they are
 // [ w w,n w,t w,nn w,nt w,tt ]
 const unsigned ndof_types =6;
 for(unsigned idof=0; idof<ndof_types;++idof)
  {
   // Trace dofs are value, tangent deriv and second tangent deriv
   const bool w_trace_dof = (idof==0 || idof==2 || idof ==5);
   const bool dwdn_trace_dof = (idof==1 || idof==4);
   // If we are clamped do w trace dofs and dwdn trace dofs 
   if(TestSoln::boundary_case==TestSoln::Clamped && (w_trace_dof || dwdn_trace_dof))
    {
     nod_pt->pin(idof);
     nod_pt->set_value(idof,0.0);
    }
   // If we are resting only do w trace dofs 
   else if(TestSoln::boundary_case==TestSoln::Resting && (w_trace_dof))
    {
     nod_pt->pin(idof);
     nod_pt->set_value(idof,0.0);
    }
   else {/* Pin Nothing */}
 }
 }//Loop nodes
} // end loop over boundaries 

if( TestSoln::boundary_case == TestSoln::Free || TestSoln::boundary_case == TestSoln::FreeSinusoidalLoad )
 {
  // Pin the displacements
  pin_displacement_at_centre_node();
//  pin_rotation_on_edge_node();
 }
} // end set bc
/// A function that upgrades straight sided elements to be curved. This involves
// Setting up the parametric boundary, F(s) and the first derivative F'(s)
// We also need to set the edge number of the upgraded element and the positions
// of the nodes j and k (defined below) and set which edge (k) is to be exterior
/*            @ k                                                             */
/*           /(                                                               */
/*          /. \                                                              */
/*         /._._)                                                             */
/*      i @     @ j                                                           */
// For RESTING or FREE boundaries we need to have a C2 CONTINUOUS boundary
// representation. That is we need to have a continuous 2nd derivative defined 
// too. This is well discussed in by [Zenisek 1981] (Aplikace matematiky , 
// Vol. 26 (1981), No. 2, 121--141). This results in the necessity for F''(s) 
// as well.
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
   const unsigned nnode=3;
   unsigned index_of_interior_node=3;

   // The edge that is curved
   MyC1CurvedElements::Edge edge;

   // Vertices positions
   Vector<Vector<double> > xn(3,Vector<double>(2,0.0));
 
   // Loop nodes
   for(unsigned n=0;n<nnode;++n)
    {
     // If it is on boundary
     Node* nod_pt = bulk_el_pt->node_pt(n);
     if(nod_pt->is_on_boundary(0) || nod_pt->is_on_boundary(1))
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
   bulk_el_pt->upgrade_to_curved_element(edge,s_ubar,s_obar,
    parametric_curve_pt);     
  }
}// end upgrade elements

template <class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::
rotate_edge_degrees_of_freedom( Mesh* const &bulk_mesh_pt)
{
 // How many bulk elements
 unsigned n_element = bulk_mesh_pt-> nelement();
 
 // Loop over all the bulk elements
 for(unsigned e=0; e<n_element; e++)
  {
   // Get pointer to bulk element    
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
 
   // Loop nodes  
   unsigned nnode = el_pt->nnode();
   unsigned nbnode=0 ;
   // Count the number of boundary nodes that element e has on boundaries 0 and 1
   for (unsigned n=0; n<nnode;++n)
     {
      nbnode+=unsigned(el_pt->node_pt(n)->is_on_boundary(0)||
                       el_pt->node_pt(n)->is_on_boundary(1));
      }

   // Now, if we have nodes on boundary upgrade it
   if(nbnode>0)
    {
     // Set up vector
     Vector<unsigned> bnode (nbnode,0);
     unsigned inode(0);

     // Fill in the bnode Vector
     for (unsigned n=0; n<nnode;++n)
      {
       // If it is on the boundary
       if((el_pt->node_pt(n)->is_on_boundary(0)||
                       el_pt->node_pt(n)->is_on_boundary(1)))
        {
         // Set up the Vector
         bnode[inode]=n;
         ++inode;
        }
      }
    // Now rotate Hermite degrees of freedom
    el_pt->set_up_rotated_dofs(nbnode,bnode,&TestSoln::get_normal_and_tangent);
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
npts = 6;

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
 oomph_info << "L2 Norm of computed solution: " << sqrt(dummy_error)<< std::endl;
 
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

 // Poisson Ratio
 CommandLineArgs::specify_command_line_flag("--nu", &TestSoln::nu);
 
 // Applied Pressure
 CommandLineArgs::specify_command_line_flag("--p", &TestSoln::p_mag);
 
 // P_step
 double p_step=10;
 CommandLineArgs::specify_command_line_flag("--dp", &p_step);
 
 // Element Area
 double element_area=0.2;
 CommandLineArgs::specify_command_line_flag("--element_area", &element_area);

 // The `usage' flag
 CommandLineArgs::specify_command_line_flag("--usage");

 // Boundary case as an unsigned
 // Flag for bad user input
 bool invalid_input=false;
 int boundary_case_num=0;
 CommandLineArgs::specify_command_line_flag("--case", &boundary_case_num);

 // Parse command line
 CommandLineArgs::parse_and_assign(); 

 // Assign the boundary case
 if (boundary_case_num>=0 && boundary_case_num <=3)
  {
  // Cast int to enum
  TestSoln::boundary_case=(TestSoln::Boundary_case)(boundary_case_num);
  }
 else // Default to what is set in TestSoln
  { 
   oomph_info<<"Boundary case \""<<boundary_case_num<<"\" not recognised.\n";
   invalid_input=true;
  }

 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Print usage 
 if(CommandLineArgs::command_line_flag_has_been_set("--usage") || invalid_input)
 {
  oomph_info << "KIRCHOFF PLATE BENDING (TRI MESH)\n"<<
   "-------------------------------------------------------------------------\n";
  oomph_info<<"\nUsage of this driver code:\n\n";
  oomph_info<<"Flags \n\n";
  oomph_info<<"--dir DIR    Specifies output directory \n";
  oomph_info<<"--p NUM              Specifies applied (constant) pressure \n";
  oomph_info<<"--dp NUM             Specifies pressure step\n";
  oomph_info<<"--element_area A     Specifies typical element size for Triangle"
            <<" mesh generator. \n";
  oomph_info<<"--n_step NUM         Specifies number of steps to take in "
            <<"pressure \n";
  oomph_info<<"--nu NUM             Specifies poisson ratio \n";
  oomph_info<<"--case CASE          Specifies boundary case. Enumeration as "
            <<"follows: \n";
  oomph_info<<"                      0 - resting, uniform pressure\n";
  oomph_info<<"                      1 - clamped, uniform pressure\n";
  oomph_info<<"                      2 - free, quadratic loading\n";
  oomph_info<<"                      3 - free, sinusoidal loading\n";
                 
  // Terminate here
  return(0);
 }

 // Problem instance
 UnstructuredFvKProblem<KirchhoffPlateBendingC1CurvedBellElement<2,2,5> >problem(element_area);
 problem.max_newton_iterations()=1;
 oomph_info<<"Solving for p=" << TestSoln::p_mag<<"\n";
 problem.newton_solve();

 // Document
 problem.doc_solution();
 oomph_info << std::endl;
 oomph_info << "---------------------------------------------" << std::endl;
 oomph_info << "Pcos (" << TestSoln::p_mag << ")" << std::endl;
 oomph_info << "Current dp  (" << p_step << ")" << std::endl;
 oomph_info << "Poisson ratio (" << TestSoln::nu << ")" << std::endl;
 oomph_info << "Solution number (" <<problem.Doc_info.number()-1 << ")" << std::endl;
 oomph_info << "---------------------------------------------" << std::endl;
 oomph_info << std::endl;
} //End of main

