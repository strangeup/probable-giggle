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
#include "C1_basis/SubparametricTElement.h"

// The mesh
#include "meshes/triangle_mesh.h"

// Analytic Solutions for linear plate bending problems
#include "kirchhoff_plate_bending_analytic_solutions.h"

using namespace std;
using namespace oomph;
using MathematicalConstants::Pi;

// Global variables in this driver code
namespace TestSoln
{
//Shape of the domain
double A = 1.0;
double B = 1.0;
// The coupling of the stretching energy
double eta = 0;
double p_mag = 1; 
double nu = 0.3;

/*                     PARAMETRIC BOUNDARY DEFINITIONS                        */
// Here we create the geom objects for the Parametric Boundary Definition 
CurvilineCircleTop parametric_curve_top;
CurvilineCircleBottom parametric_curve_bottom;

// Assigns the value of pressure depending on the position (x,y)
void get_pressure(const Vector<double>& x, double& pressure)
{
 pressure = p_mag; //constant pressure
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
void get_exact_radial_w(const Vector<double>& xi, Vector<double>& w)
{
// //solution (r^4-br^2+c)/64 for w=0 and w'=0 or mrt=0
Vector<double> polar_coord(2);
polar_coord[0] = sqrt(xi[0]*xi[0]+xi[1]*xi[1]);
// Rotate the solution by -Pi/6
polar_coord[1] = atan2(xi[0],xi[1])-Pi/6.;
// Get analytic solution
Vector<double> wpolar(6,0.0);
// Compute the solution up to fifteenth term (worry about round-off / underflow
// here - but the solution is pap so we need a few terms...)
const unsigned nterms = 15 , npins = 3;
KirchhoffPlateBendingAnalyticSolutions::circle_pinned_symmetrically(polar_coord, npins, nu, wpolar, nterms);
w[0] = wpolar[0];
w[1] = wpolar[1];
w[2] = wpolar[2];
w[3] = wpolar[3];
w[4] = wpolar[4];
w[5] = wpolar[5];

}

// Output in High resolution
bool High_resolution =false ;
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
void actions_before_newton_solve()
{
 complete_problem_setup();
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
 // 0 to Pi
 Outer_boundary0 = 0,
 Outer_boundary1 = 1,
 Outer_boundary2 = 2,
 //Pi to 2Pi
 Outer_boundary3 = 3,
 Outer_boundary4 = 4,
 Outer_boundary5 = 5,
 // Internal boundary to make location of central node possible
 Inner_boundary0 = 6,
};

/// Element area
double Element_area;

/// Central node pointer
Node* Centre_node_pt;

/// Helper to upgrade edge elements on specified boundary to be curved
void upgrade_edge_elements_to_curve(const unsigned &b, Mesh* const & 
 bulk_mesh_pt);

/// Helper to rotate edge degrees of freedom
void rotate_edge_degrees_of_freedom(Mesh* const &bulk_mesh_pt);

public:

/// Get the central value of the deflection
void get_deflection_at_centre(double& deflection)
{
// IF we have already found it
if(Centre_node_pt!=0)
{
 deflection =  Centre_node_pt->raw_value(0);
 return;
}

// Loop over nodes 
unsigned num_int_nod=Bulk_mesh_pt->nboundary_node(Inner_boundary0);
for (unsigned inod=0;inod<num_int_nod;inod++)
{
 // Get node point
 Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(Inner_boundary0,inod);
 // IF it is at the centre
 if(nod_pt->x(0)==0.0 && nod_pt->x(1)==0.0)
  {
   // Get the central value
   Vector<double> w_exact(6);
   TestSoln::get_exact_radial_w(Vector<double>(2,0.0), w_exact);
   Centre_node_pt = nod_pt;

   // Found the central value
   oomph_info<<"Found central node in mesh at x=("
            <<nod_pt->x(0)<<" "<<nod_pt->x(1)<<")\n";
   oomph_info<<"The calculated central value is: " 
             << nod_pt->raw_value(0)<<std::endl;
   oomph_info<<"The exact central value is: " << w_exact[0] <<std::endl;
   deflection = nod_pt->raw_value(0);
   return;
  }
 }
 // IF we get here we couldn't find it
 oomph_info<<"Could not find central node in mesh.\n";
} 

private:
/// Pointer to "bulk" mesh
TriangleMesh<ELEMENT>* Bulk_mesh_pt;

/// Pointer to "surface" mesh
Mesh* Surface_mesh_pt;

}; // end_of_problem_class

/// Constructor
template<class ELEMENT>
UnstructuredFvKProblem<ELEMENT>::UnstructuredFvKProblem(double element_area)
: Element_area(element_area) , Centre_node_pt(0)
{
Vector<double> zeta(1);

//Outer boundary
//--------------

double A = 1.0;
double B = 1.0;
Ellipse* outer_boundary_ellipse_pt = new Ellipse(A, B);

TriangleMeshClosedCurve* outer_boundary_pt = 0;

// Here we split the boundary into six parts, so that we have the three pins
// accessible and we can use the existing parametric curves which define
// the top and bottom of a circle
Vector<TriangleMeshCurveSection*> outer_curvilinear_boundary_pt(6);
//First bit
double zeta_start = 0.0;
double zeta_end = MathematicalConstants::Pi/3.;
unsigned nsegment = (int)(MathematicalConstants::Pi/(3.0*sqrt(element_area)));
outer_curvilinear_boundary_pt[0] =
new TriangleMeshCurviLine(outer_boundary_ellipse_pt, zeta_start,
zeta_end, nsegment, Outer_boundary0);

//Second bit
zeta_start = MathematicalConstants::Pi/3.;
zeta_end = 2.*MathematicalConstants::Pi/3.;
nsegment = (int)(MathematicalConstants::Pi/(3.0*sqrt(element_area)));
outer_curvilinear_boundary_pt[1] =
new TriangleMeshCurviLine(outer_boundary_ellipse_pt, zeta_start,
zeta_end, nsegment, Outer_boundary1);

//Third bit
zeta_start = 2.*MathematicalConstants::Pi/3.;
zeta_end = MathematicalConstants::Pi;
nsegment = (int)(MathematicalConstants::Pi/(3.0*sqrt(element_area)));
outer_curvilinear_boundary_pt[2] =
new TriangleMeshCurviLine(outer_boundary_ellipse_pt, zeta_start,
zeta_end, nsegment, Outer_boundary2);

//Fourth bit
zeta_start = MathematicalConstants::Pi;
zeta_end = 4.0*MathematicalConstants::Pi/3.;
nsegment = (int)(MathematicalConstants::Pi/(3.0*sqrt(element_area)));
outer_curvilinear_boundary_pt[3] =
new TriangleMeshCurviLine(outer_boundary_ellipse_pt, zeta_start,
zeta_end, nsegment, Outer_boundary3);

//Fifth bit
zeta_start = 4.0*MathematicalConstants::Pi/3.;
zeta_end = 5.*MathematicalConstants::Pi/3.;
nsegment = (int)(MathematicalConstants::Pi/(3.0*sqrt(element_area)));
outer_curvilinear_boundary_pt[4] =
new TriangleMeshCurviLine(outer_boundary_ellipse_pt, zeta_start,
zeta_end, nsegment, Outer_boundary4);

//Sixth bit
zeta_start = 5.*MathematicalConstants::Pi/3.;
zeta_end = 2.0*MathematicalConstants::Pi;
nsegment = (int)(MathematicalConstants::Pi/(3.0*sqrt(element_area)));
outer_curvilinear_boundary_pt[5] =
new TriangleMeshCurviLine(outer_boundary_ellipse_pt, zeta_start,
zeta_end, nsegment, Outer_boundary5);
outer_boundary_pt =
new TriangleMeshClosedCurve(outer_curvilinear_boundary_pt);

// Internal bit - this means we can have a boundary which is just the centre
// We start by creating the internal boundaries
// The boundary 2 is defined by its two vertices
// Open curve 1
Vector<Vector<double> > vertices(2,Vector<double>(2,0.0));
vertices[0][0] = 0.5;
vertices[0][1] = 0.0;

vertices[1][0] = 0.0;
vertices[1][1] = 0.0;
unsigned boundary_id = Inner_boundary0;

TriangleMeshPolyLine *boundary2_pt =
  new TriangleMeshPolyLine(vertices, boundary_id);

// Total number of open curves in the domain
unsigned n_open_curves = 1;
// We want internal open curves
Vector<TriangleMeshOpenCurve *> inner_open_boundaries_pt(n_open_curves);

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
upgrade_edge_elements_to_curve(2,Bulk_mesh_pt);
upgrade_edge_elements_to_curve(3,Bulk_mesh_pt);
upgrade_edge_elements_to_curve(4,Bulk_mesh_pt);
upgrade_edge_elements_to_curve(5,Bulk_mesh_pt);
 
// Rotate degrees of freedom
// This is unecessary when we are only pinning single points
// rotate_edge_degrees_of_freedom(Bulk_mesh_pt);

// Store number of bulk elements
complete_problem_setup();

char filename[100];
sprintf(filename, "RESLT/trace.dat");
Trace_file.open(filename);

oomph_info << "Number of equations: "
        << assign_eqn_numbers() << '\n';

// Clean up memory
delete outer_boundary_pt;
outer_boundary_pt = 0;
delete outer_boundary_ellipse_pt;
outer_boundary_ellipse_pt = 0;
const unsigned n_bounds = 6;
// Loop over boundary_pt
for(unsigned ibound = 0; ibound<n_bounds;++ibound)
 {
 delete outer_curvilinear_boundary_pt[ibound];
 outer_curvilinear_boundary_pt[ibound] = 0;
 }
delete inner_open_boundaries_pt[0];
inner_open_boundaries_pt[0] = 0;
delete boundary2_pt;
boundary2_pt = 0;
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
 oomph_info<<"Applying boundary conditions."<<std::endl;
// Set the boundary conditions for problem: All nodes are
// free by default -- just pin the ones that have Dirichlet conditions
// here. 
// Just loop over outer boundary since inner boundary doesn't have boundary
// conditions
// We want to pin the two points that are on both boundaries (1,0) and (-1,0) 
// Loop over boundaries
const unsigned nbound = 6;
for(unsigned ibound =0; ibound<nbound;ibound+=2)
{
 // Now find the i+1 boundary
 unsigned jbound = (ibound + 1) % (nbound);
 unsigned num_int_nod=Bulk_mesh_pt->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_int_nod;inod++)
 {
  // Get node point
  Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
  // If the node is on the other internal boundary too
  if( nod_pt->is_on_boundary(jbound))
  {
   oomph_info<<"Found boundary node that belongs to two boundaries "
             <<"at x = ("<<nod_pt->x(0)<<" "<<nod_pt->x(1)<<") \n";
   oomph_info<<"Pinning it.\n";
   // Pin it - if it's on both boundaries
   nod_pt->pin(0);
   nod_pt->set_value(0,0.0);
  }
 } //end loop over boundary nodes
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
void UnstructuredFvKProblem<ELEMENT >::
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
  case 0: case 1: case 2:
   parametric_curve_pt = &TestSoln::parametric_curve_top;
  break;

  // Lower boundary
  case 3: case 4: case 5:
   parametric_curve_pt = &TestSoln::parametric_curve_bottom;
  break;

  default:
   throw OomphLibError(
    "I have encountered a boundary number that I wasn't expecting. Please fill \
me in if you want additional curved boundaries..",
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
   
   // Loop over (vertex) nodes
   const unsigned nnode=3; //This should always be = 3 for triangles
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

     // Check if it is on the outer boundaries
     bool is_on_boundary = false;
     const unsigned n_outer_boundaries=6;
     for(unsigned i=0; i<n_outer_boundaries;++i)
       {
        // On a boundary if on any boundary
        is_on_boundary = is_on_boundary || nod_pt->is_on_boundary(i);
       }
     if(is_on_boundary) 
      {
       xn[n][0]=nod_pt->x(0);
       xn[n][1]=nod_pt->x(1);
      }
     // The edge is denoted by the index of the  opposite (interior) node
     else {index_of_interior_node = n;}
    }
   // Initialise s_ubar s_obar (start and end respectively)
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
   // Check for inverted elements HERE
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
  }
}// end upgrade elements

template <class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::
rotate_edge_degrees_of_freedom(Mesh* const &bulk_mesh_pt)
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
     {nbnode+=unsigned(el_pt->node_pt(n)->is_on_boundary());}

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
       if(el_pt->node_pt(n)->is_on_boundary())
        {
         // Set up the Vector
         bnode[inode]=n;
         ++inode;
        }
      }
    // Output that we have found element HERE
    oomph_info<<"Element "<<e<<" has "<<bnode<< " nodes on the boundary.\n";

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
npts = (TestSoln::High_resolution ? 25 : 6 );

sprintf(filename,"RESLT/soln%i-%f.dat",Doc_info.number(),Element_area);
some_file.open(filename);
Bulk_mesh_pt->output(some_file,npts); 
some_file << "TEXT X = 22, Y = 92, CS=FRAME T = \"" 
       << comment << "\"\n";
some_file.close();

//  Output exact solution
sprintf(filename,"%s/exact_interpolated_soln%i-%f.dat","RESLT",Doc_info.number(),Element_area);
some_file.open(filename);
Bulk_mesh_pt->output_fct(some_file,npts,TestSoln::get_exact_radial_w); 
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

// Output the centre deflection
double centre_deflection(0);
get_deflection_at_centre(centre_deflection);

 // Doc error and return of the square of the L2 error
 //---------------------------------------------------
 //double error,norm,dummy_error,zero_norm;
  double dummy_error,zero_norm;
 sprintf(filename,"RESLT/error%i-%f.dat",Doc_info.number(),Element_area);
 some_file.open(filename);
 
 Bulk_mesh_pt->compute_error(some_file,TestSoln::get_exact_radial_w,
                        dummy_error,zero_norm);
 some_file.close();

 
 // Doc L2 error and norm of solution
 oomph_info << "L2 Norm of computed solution: " << sqrt(dummy_error/zero_norm)<< std::endl;

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

 // Validation
 CommandLineArgs::specify_command_line_flag("--validation");

 // High Resolution
 CommandLineArgs::specify_command_line_flag("--high_resolution");

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

 // Parse command line
 CommandLineArgs::parse_and_assign(); 

 // Output in high resolution
 TestSoln::High_resolution = CommandLineArgs::
   command_line_flag_has_been_set("--high_resolution");

 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Problem instance
 UnstructuredFvKProblem<KirchhoffPlateBendingC1CurvedBellElement<2,2,5> >problem(element_area);
 oomph_info<<"Solving for p=" << TestSoln::p_mag<<"\n";
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
} //End of main

