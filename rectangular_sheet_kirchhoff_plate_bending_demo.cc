//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented,
//LIC// multi-physics finite-element library, available
//LIC// at http://www.oomph-lib.org.
//LIC//
//LIC//           Version 0.85. June 9, 2008.
//LIC//
//LIC// Copyright (C) 2006-2008 Matthias Heil and Andrew Hazel
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
//LIC// The authors may be contacted at /oomph-lib@maths.man.ac.uk.
//LIC//
//LIC//====================================================================
#include <fenv.h>

//Generic routines
#include "generic.h"

// Include the library elements
#include "C1_linear_plate_bending.h"

// The mesh
#include "meshes/triangle_mesh.h"

using namespace std;
using namespace oomph;
using MathematicalConstants::Pi;

// =========================================================================
// The boundaries are enumerated as follows:
//
// 0 :  x = - 1/2, y = [ -L/2 , L/2], Normal : -(1,0)
// 1 :  y = + L/2, x = [ -1/2 , 1/2], Normal : +(0,1)
// 2 :  x = + 1/2, y = [ -L/2 , L/2], Normal : +(1,0)
// 3 :  y = - L/2, x = [ -1/2 , 1/2], Normal : -(0,1)
//
//           1
//        _______
//       |       |
//    0  |       |  2
//       |_______|
//
//           3

namespace TestSoln
{
 /// Nondim length of strip (width is lengthscale and therefore 1)
 // length is in the y-direction, width=1 in y
 double Length=1.0;

 /// Poisson's ratio
 double Nu = 0.3;

  /// Nondimensional pressure scale
 double Pressure=1.0;

 // Flag for high resolution output
 bool High_resolution = false;

 /// Bool to determine which boundaries are at zero slope
 bool is_boundary_held_at_zero_slope(unsigned& ibound)
 {
    // Free/Pinned/Sliding/Clamped
    return (ibound == 2 || ibound == 3);
 }

 /// Bool to determine which boundaries to pin
 bool is_boundary_pinned(unsigned& ibound)
 {
    // Free/Pinned/Sliding/Clamped
    return (ibound == 1 || ibound == 3);
 }

 // Assigns the value of pressure depending on the position (x,y)
 void get_pressure(const Vector<double>& coord, double& pressure)
 {
    pressure = Pressure;
 }
} // end of TestSoln 

//==start_of_problem_class============================================
/// Class definition
//====================================================================
template<class ELEMENT>
class UnstructuredFvKProblem : public virtual Problem
{

public:

  /// Constructor
  UnstructuredFvKProblem(double element_area);

 /// Destructor
 ~UnstructuredFvKProblem()
  {
   // Close the trace
   Trace_file.close();
   // Delete the Surface and Bulk mesh
   delete mesh_pt();
  };

  /// Update after solve (empty)
  void actions_after_newton_solve()
   {
   }

  /// Update the problem specs before solve: Re-apply boundary conditons
  void actions_before_newton_solve()
  {
    // This will also reapply boundary conditions
    complete_problem_setup();
    // Reassign the equation numbers
    assign_equation_numbers();
  }

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

 /// Unpin all of the dofs
 void unpin_all_dofs()
   {
    // Get total number of nodes
    unsigned n_node=mesh_pt()->nnode();
    // Loop over nodes
    for(unsigned inod=0;inod<n_node;++inod)
     {
      Node* nod_pt=mesh_pt()->node_pt(inod);
      const unsigned ndof_type=6;
      for(unsigned i=0;i<ndof_type;++i)
       { nod_pt->unpin(i); }
     }
   }

private:

  /// Trace file to document norm of solution
  ofstream Trace_file;

  /// Helper function to apply boundary conditions
  void apply_boundary_conditions();

  /// \short Helper function to (re-)set boundary condition
  /// and complete the build of  all elements
  void complete_problem_setup();

  /// Helper function to (re-)set equation numbers for problem
  void assign_equation_numbers()
   {
   // Assign equations numbers and document
   oomph_info << "Number of equations: "
              << this->assign_eqn_numbers() << '\n';
   oomph_info << "Number of elements: " <<mesh_pt()->nelement() << std::endl;
   }

  /// Create the mesh
  TriangleMeshParameters* set_up_triangle_mesh_parameters_for_rectangle(const double& element_area);

private:
  // The initial (and maximum) element area
  double Element_area;
  const unsigned Number_of_boundaries = 4;
}; // end_of_problem_class

// Free function helper in anonymous namepace (use a lambda in c++11)
namespace {
 // Helper to make three vertice line from two vertices
 Vector<Vector<double>> make_line(const Vector<double>& x1,
    const Vector<double>& x2)
  {
    const unsigned number_of_vertices = 3, DIM = x1.size();
    Vector<Vector<double>> vertices(number_of_vertices, Vector<double>(DIM));
    vertices[0] = x1;
    // Compute middle co-ordinate
    for(unsigned idim = 0; idim <DIM; ++idim)
      { vertices[1][idim] = (x1[idim] + x2[idim])/2.0; }
    vertices[2] = x2;
    return vertices;
  }
}

/// Set-up the rectangular mesh for the problem
template <class ELEMENT>
TriangleMeshParameters* UnstructuredFvKProblem<ELEMENT>::set_up_triangle_mesh_parameters_for_rectangle(
 const double& element_area)
{
 // Local copies
 const double length = TestSoln::Length/2.0;
 const double width = 1.0/2.0;

 // A Rectangular sheet has four corners
 const unsigned num_vertices = 4, DIM = 2, num_boundaries = 4;
 Vector<Vector <double> > corners(num_vertices, Vector<double>(DIM));
 corners[0][0] = -width; corners[0][1] = -length;
 corners[1][0] = -width; corners[1][1] =  length;
 corners[2][0] =  width; corners[2][1] =  length;
 corners[3][0] =  width; corners[3][1] = -length;

 // >> Setting up the domain with PolyLines
 auto outer_boundary_polyline_pt = Vector<TriangleMeshCurveSection*>();
 outer_boundary_polyline_pt.resize(num_boundaries);

 // Labelling from Boundary 0 (left boundary), fill create boundary lines
 for(unsigned boundary_id = 0; boundary_id < num_boundaries; ++boundary_id)
  {
   const unsigned icorner1 = boundary_id % 4,
                  icorner2 = (boundary_id + 1) % 4;
   outer_boundary_polyline_pt[boundary_id] = new TriangleMeshPolyLine(
      ::make_line(corners[icorner1], corners[icorner2]),
      boundary_id);
  }

 // The outer polygon
 auto outer_boundary_pt =
  new TriangleMeshClosedCurve(outer_boundary_polyline_pt);

 //Create mesh parameters object
 auto triangle_mesh_parameters_pt = new TriangleMeshParameters(outer_boundary_pt);

 // Set the element area into the mesh arguments
 triangle_mesh_parameters_pt->element_area() = element_area;

 return triangle_mesh_parameters_pt;
} // end of set_up_triangle_mesh_parameters_for_rectangle

/// Constructor for the FvK problem
template<class ELEMENT>
UnstructuredFvKProblem<ELEMENT>::UnstructuredFvKProblem(double element_area)
  : Element_area(element_area) 
{
 // Set the maximum residuals for Newton iterations
 Problem::Max_residuals = 1000;
 // Creates a triangle_mesh_parameters instance
 auto triangle_mesh_parameters_pt = set_up_triangle_mesh_parameters_for_rectangle(element_area);
 // Set up the mesh
 mesh_pt() = new TriangleMesh<ELEMENT>(*triangle_mesh_parameters_pt);

 oomph_info<<"Number Elements "<<mesh_pt()->nelement()<<std::endl;
 oomph_info<<"Number Nodes "<<mesh_pt()->nnode()<<std::endl;

 // Open the trace and start recording
 char filename[100];
 sprintf(filename, "RESLT/trace.dat");
 Trace_file.open(filename);

 // Complete the problem setup
 complete_problem_setup();

 // Assign the equation numbers
 assign_equation_numbers();

 // Clean up assigned memory
 auto outer_boundary_pt = triangle_mesh_parameters_pt->outer_boundary_pt()[0];
 const unsigned num_bounds = outer_boundary_pt->max_boundary_id();
 for(unsigned ibound = 0; ibound< num_bounds;  ++ibound)
  { 
   delete outer_boundary_pt->curve_section_pt(ibound);
  }
 delete outer_boundary_pt; 
 delete triangle_mesh_parameters_pt;
} // end_constructor

//==start_of_complete======================================================
 /// Set boundary condition exactly, and complete the build of
 /// all elements
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::complete_problem_setup()
{
 // -------------------------------------------------------------------
 // Complete the build of all elements so they are fully functional
 const unsigned n_element = mesh_pt()->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));

   //Set the pressure function pointers and the physical constants
   el_pt->nu_pt() = &TestSoln::Nu;
   el_pt->pressure_fct_pt() = &TestSoln::get_pressure;
  }
 // Apply boundary conditions
 apply_boundary_conditions();
} //end_of_complete

//==start_of_apply_bc=====================================================
/// Helper function to apply boundary conditions
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::apply_boundary_conditions()
{
 // First reset any boundary conditions
 unpin_all_dofs();

 // Print helpful message
 oomph_info<<"Applying boundary conditions to the sheet."<<std::endl;

 // loop over boundaries
 for(unsigned ibound=0;ibound<Number_of_boundaries;ibound++)
  {
   const unsigned num_nod=mesh_pt()->nboundary_node(ibound);

   for (unsigned inod=0;inod<num_nod;inod++)
    {
     Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);

     // Dofs (by default) at a node X = (x0, x1) are:
     // #0 - w(x1,x2)
     // #1 - d/dx w (x0,x1) 
     // #2 - d/dy w (x0,x1) 
     // #3 - d2/dx2 w (x0,x1) 
     // #4 - d2/dxdy w (x0,x1) 
     // #5 - d2/dy2 w(x0,x1) 

     // If the boundary is to be pinned:
     // Set w(x1, x2), dw/dt (x1, x2) and d2w/dt2 (x1, x2) for all nodes on the
     // boundary.
     if(TestSoln::is_boundary_pinned(ibound))
      {
       // Solution and tangential derivatives on the boundary {w, dw/dt, d2wdt2}
       const Vector<double> solution_on_boundary(3,0.0);
       // Is this boundary x-axis aligned or y-axis aligned
       const bool is_y_aligned = ibound % 2 == 0;
       // Tangent direction: y-axis on even and x-axis on odd boundaries
       // Dof number for value dof, derivative wrt. tangent dof and 2nd derivative
       // wrt tangent dof
       const unsigned iw = 0;
       const unsigned idwdt = (is_y_aligned ? 2 : 1);
       const unsigned id2wdt2 = (is_y_aligned  % 2 == 0 ? 5 : 3);

       // Pin value w(x,y) at boundary
       nod_pt->pin(iw);
       nod_pt->set_value(iw, solution_on_boundary[0]);

       // Pin tangent derivative d/dt w(x,y)
       nod_pt->pin(idwdt);
       nod_pt->set_value(idwdt, solution_on_boundary[1]);
       // Pin second tangent derivative d2/dt^2 w(x,y)
       nod_pt->pin(id2wdt2);
       nod_pt->set_value(id2wdt2, solution_on_boundary[2]);
      } // for (inod<num_nod)

     // If an angle is to be imposed:
     // Set dw/dn(x1, x2) and d2w/dndt (x1, x2) for all nodes.
     if(TestSoln::is_boundary_held_at_zero_slope(ibound))
      {
       // Normal and cross derivative on the boundary {dw/dn, dw2/dwdn}
       const Vector<double> normal_derivative_on_boundary(2,0.0);
       // Get node
       Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);

       // Normal direction: x-axis on even and y-axis on odd boundaries
       const unsigned idwdn = (ibound % 2 == 0 ? 1 : 2),
                      id2wdndt = 4;
       // Pin the normal derivatives, d/dn w(x,y)
       nod_pt->pin(idwdn);
       nod_pt->set_value(idwdn, normal_derivative_on_boundary[0]);
       // Pin the cross derivative, d2/dtdn w(x,y)
       nod_pt->pin(id2wdndt);
       nod_pt->set_value(id2wdndt, normal_derivative_on_boundary[1]);
      } // end if
    } // end for
  } // for (ibound<Number_of_boundaries)
} // end set bc


// Free function helper in anonymous namepace (use lambda in c++11)
namespace {
  void output_solution(const DocInfo& doc_info, Mesh* mesh_pt,
     const std::string& basename = "soln",
     const unsigned npts = 5)
    {
     char filename[100];
     ofstream some_file;
     sprintf(filename,"%s/%s%i.dat",
             basename.c_str(),
             doc_info.directory().c_str(),
             doc_info.number());
     some_file.open(filename);
     mesh_pt->output(some_file,npts);
     some_file.close();
   }
}

//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{
  oomph_info << "==============================================" << std::endl;
  oomph_info << "========= Documenting solution ===============" << std::endl;
  oomph_info << "==============================================" << std::endl;
  oomph_info << std::endl;

  // Number of plot points
  unsigned const npts_coarse= 2;
  ::output_solution(doc_info, this->mesh_pt(), "coarse_soln", npts_coarse);

  // Output with very high resolution
  unsigned const npts = TestSoln::High_resolution ? 25 : 5;
  ::output_solution(doc_info, this->mesh_pt(), "soln", npts);
  ::output_solution(doc_info, this->mesh_pt(), "bc", npts);

 // Find the solution at x = y = 0
 MeshAsGeomObject* Mesh_as_geom_obj_pt=
  new MeshAsGeomObject(mesh_pt());
 Vector<double> s(2);
 GeomObject* geom_obj_pt=0;
 Vector<double> r(2,0.0);
 Mesh_as_geom_obj_pt->locate_zeta(r,geom_obj_pt,s);

 // Get interpolated displacement at the centre
 const unsigned number_of_output_fields = 6;
 Vector<double> u_0(number_of_output_fields,0.0);
 dynamic_cast<ELEMENT*>(geom_obj_pt)->interpolated_u_biharmonic(s,u_0);
 oomph_info << "w in the middle: " << u_0[0] << std::endl;

 // Get the precision
 const unsigned default_precision = Trace_file.precision(),
                temporary_precision = 15;
 // Print to trace file
 Trace_file.precision(temporary_precision);
 Trace_file <<Element_area <<" "<< u_0[0] << std::endl;
 Trace_file.precision(default_precision);

 // Increment the doc_info number
 doc_info.number()++;

 delete Mesh_as_geom_obj_pt;
} // end of doc

//=======start_of_main========================================
///Driver code for demo of inline triangle mesh generation
//============================================================
int main(int argc, char **argv)
{
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Swith timings on
 Global_timings::Doc_comprehensive_timings = true;

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define possible command line arguments and parse the ones that
 // were actually specified
 string output_dir="RESLT";
 double element_area=0.1;
 CommandLineArgs::specify_command_line_flag("--dir", &output_dir);
 CommandLineArgs::specify_command_line_flag("--element_area", &element_area);
 CommandLineArgs::specify_command_line_flag("--high_resolution");
 CommandLineArgs::specify_command_line_flag("--p", &TestSoln::Pressure);
 CommandLineArgs::specify_command_line_flag("--nu", &TestSoln::Nu);
 CommandLineArgs::specify_command_line_flag("--length", &TestSoln::Length);
 CommandLineArgs::specify_command_line_flag("--validate");

 // Parse command line
 CommandLineArgs::parse_and_assign();

 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Output in high resolution
 TestSoln::High_resolution = CommandLineArgs::
   command_line_flag_has_been_set("--high_resolution");

 if(CommandLineArgs::command_line_flag_has_been_set("--validate"))
 {
   oomph_info<<"Validate set, ignoring other command line arguments and setting"
             <<" defaults."<<std::endl;
   element_area = 0.1;
   TestSoln::High_resolution = false;
   TestSoln::Pressure = 1.0;
   TestSoln::Nu = 0.3;
   TestSoln::Length = 1.0;
 };

 // Label for output
 DocInfo doc_info;

 // Output directory
 doc_info.set_directory(output_dir);

 // Create Problem instance, using (not upgraded) curved bell elements
 UnstructuredFvKProblem<KirchhoffPlateBendingC1CurvedBellElement>
   problem(element_area);
 problem.max_newton_iterations() = 1;

 // Newton Solve
 problem.newton_solve();

 //Output solution
 problem.doc_solution(doc_info);

 oomph_info << std::endl;
 oomph_info << "---------------------------------------------" << std::endl;
 oomph_info << "Pressure (" << TestSoln::Pressure << ")" << std::endl;
 oomph_info << "Poisson ratio (" << TestSoln::Nu << ")" << std::endl;
 oomph_info << "Solution number (" << doc_info.number()-1 << ")" << std::endl;
 oomph_info << "---------------------------------------------" << std::endl;
 oomph_info << std::endl;

} //End of main
