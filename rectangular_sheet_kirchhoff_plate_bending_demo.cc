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

// The mesh
#include "meshes/triangle_mesh.h"

// Analytic Solutions for linear plate bending problems
#include "kirchhoff_plate_bending_analytic_solutions.h"

using namespace std;
using namespace oomph;
using MathematicalConstants::Pi;

// =========================================================================
namespace TestSoln
{
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

 /// Nondim length of strip (width is lengthscale and therefore 1)
 // length is in the y-direction, width=1 in y
 double Length_of_strip=1.0;

 /// Poisson's ratio
 double nu = 0.3; /// hierher FvK parameter (Poisson's ratio)

 // Flag for high resolution output
 bool High_resolution = false;

  /// Perturbation
 double Pressure=1.0;

 /// Bool to determine which boundaries require traction elements
 bool is_boundary_moment_free(unsigned& ibound)
 {
    // Free/Pinned/Sliding/Clamped
    return ibound == 0 || ibound == 1;
 }

 bool is_boundary_held_at_zero_slope(unsigned& ibound)
 {
    // Free/Pinned/Sliding/Clamped
    return ibound == 2 || ibound == 3;
 }

 /// Bool to determine which boundaries to pin
 bool is_boundary_pinned(unsigned& ibound)
 {
    // Free/Pinned/Sliding/Clamped
    return ibound == 1 || ibound == 3;
 }

 // Assigns the value of pressure depending on the position (x,y)
 void get_pressure(const Vector<double>& coord, double& pressure)
 {
    pressure = Pressure;
 }

 // Get the exact solution
 void get_exact_w(const Vector<double>& x, Vector<double>& exact_w)
  {
     exact_w = Vector<double>(exact_w.size(),0.0);
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
  UnstructuredFvKProblem(double element_area);

 /// Destructor
 ~UnstructuredFvKProblem()
  {
   // Clean up boundary data
   delete Outer_boundary_pt;
   delete Boundary0_pt;
   delete Boundary1_pt;
   delete Boundary2_pt;
   delete Boundary3_pt;
   // Delete the parameters
   delete Triangle_mesh_parameters_pt;
   // Close the trace
   Trace_file.close();
   // Delete the Surface and Bulk mesh
   delete this->Surface_mesh_pt;
   delete Bulk_mesh_pt;
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

 /// \short Overloaded version of the problem's access function to
 /// the mesh. Recasts the pointer to the base Mesh object to
 /// the actual mesh type.
 TriangleMesh<ELEMENT>* mesh_pt()
  {
   return dynamic_cast<RefineableTriangleMesh<ELEMENT>*> (Problem::mesh_pt());
  }

 // Unpin all of the dofs
  void unpin_all_dofs()
    {
     // Get total number of nodes
     unsigned n_node=Bulk_mesh_pt->nnode();
     // Loop over nodes
     for(unsigned inod=0;inod<n_node;++inod)
      {
       Node* nod_pt=Bulk_mesh_pt->node_pt(inod);
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
  /// and //complete the build of  all elements
  void complete_problem_setup();

  /// Helper function to (re-)set equation numbers for problem
  void assign_equation_numbers()
   {
   // Assign equations numbers
   oomph_info << "Number of equations: "
              << this->assign_eqn_numbers() << '\n';

   // Document the number of elements in the mesh
   oomph_info << "Number of elements: " << Bulk_mesh_pt->nelement() << std::endl;
   }

  /// Create the mesh
  void set_up_rectangular_mesh(const double& element_area);

  /// Pointers to specific mesh
  TriangleMesh<ELEMENT>* Bulk_mesh_pt;

  /// Pointer to the "surface" mesh
  Mesh* Surface_mesh_pt;

protected:

  // The initial (and maximum) element area
  double Element_area;
  // The mesh parameters
  TriangleMeshParameters* Triangle_mesh_parameters_pt;
  TriangleMeshClosedCurve* Outer_boundary_pt;
  TriangleMeshPolyLine* Boundary0_pt;
  TriangleMeshPolyLine* Boundary1_pt;
  TriangleMeshPolyLine* Boundary2_pt;
  TriangleMeshPolyLine* Boundary3_pt;

private:
  const unsigned Number_of_boundaries = 4;
}; // end_of_problem_class

/// Set-up the rectangular mesh for the problem
template <class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::set_up_rectangular_mesh(
 const double& element_area)
{
 // Local copies
 const double length = TestSoln::Length_of_strip/2.0;
 const double width = 1.0/2.0;

 // Create the polylines
 // ---------------------------------------------------------------------
 // >> Boundary 0 (left boundary)
 const unsigned num_vertices_b0 = 3;

 Vector<Vector <double> > vertices(num_vertices_b0);

 for (unsigned i = 0; i < num_vertices_b0; i++)
   {
     vertices[i].resize(2);
   }

 vertices[0][0] = -width;
 vertices[0][1] = -length;

 vertices[1][0] = -width;
 vertices[1][1] = 0.0;

 vertices[2][0] = -width;
 vertices[2][1] = length;

 unsigned boundary_id = 0;
 Boundary0_pt =
   new TriangleMeshPolyLine(vertices, boundary_id);

 // ---------------------------------------------------------------------
 // >> Boundary 1
 const unsigned num_vertices_b1 = 3;
 vertices.resize(num_vertices_b1);
 for (unsigned i = 0; i < num_vertices_b1; i++)
   {
     vertices[i].resize(2);
   }

 vertices[0][0] = -width;
 vertices[0][1] = length;

 vertices[1][0] = 0.0;
 vertices[1][1] = length;

 vertices[2][0] = width;
 vertices[2][1] = length;

 boundary_id = 1;
 Boundary1_pt =
   new TriangleMeshPolyLine(vertices, boundary_id);

 // ---------------------------------------------------------------------
 // >> Boundary 2
 const unsigned num_vertices_b2 = 3;
 vertices.resize(num_vertices_b2);
 for (unsigned i = 0; i < num_vertices_b2; i++)
   {
     vertices[i].resize(2);
   }

 vertices[0][0] = width;
 vertices[0][1] = length;

 vertices[1][0] = width;
 vertices[1][1] = 0.0;

 vertices[2][0] = width;
 vertices[2][1] = -length;

 boundary_id = 2;
 Boundary2_pt =
   new TriangleMeshPolyLine(vertices, boundary_id);

 // ---------------------------------------------------------------------
 // >> Boundary 3
 const unsigned num_vertices_b3 = 3;
 vertices.resize(num_vertices_b3);
 for (unsigned i = 0; i < num_vertices_b3; i++)
   {
     vertices[i].resize(2);
   }

 vertices[0][0] = width;
 vertices[0][1] = -length;

 vertices[1][0] = 0.0;
 vertices[1][1] = -length;

 vertices[2][0] = -width;
 vertices[2][1] = -length;

 boundary_id = 3;
 Boundary3_pt =
   new TriangleMeshPolyLine(vertices, boundary_id);

 // ---------------------------------------------------------------------
 // >> Building the OUTER BOUNDARY
 // >> Setting up the domain with PolyLines
 Vector<TriangleMeshCurveSection*> outer_boundary_polyline_pt(4);

 outer_boundary_polyline_pt[0] = Boundary0_pt;
 outer_boundary_polyline_pt[1] = Boundary1_pt;
 outer_boundary_polyline_pt[2] = Boundary2_pt;
 outer_boundary_polyline_pt[3] = Boundary3_pt;

 // The outer polygon
 Outer_boundary_pt =
  new TriangleMeshClosedCurve(outer_boundary_polyline_pt);

 // --------------------------------------------------------------------
 //Create the mesh
 //---------------
 //Create mesh parameters object
 Triangle_mesh_parameters_pt = new TriangleMeshParameters(Outer_boundary_pt);

 // Set the element area into the mesh arguments
 Triangle_mesh_parameters_pt->element_area() = element_area;
}

template<class ELEMENT>
UnstructuredFvKProblem<ELEMENT>::UnstructuredFvKProblem(double element_area)
  : Element_area(element_area), Triangle_mesh_parameters_pt(0), Outer_boundary_pt(0),
    Boundary0_pt(0), Boundary1_pt(0), Boundary2_pt(0), Boundary3_pt(0)
{
 // Set the maximum residuals for Newton iterations
 Problem::Max_residuals = 1000;
 // Set up the mesh
 set_up_rectangular_mesh(element_area);
 Bulk_mesh_pt = new TriangleMesh<ELEMENT>(*Triangle_mesh_parameters_pt);

 // Store number of Poisson bulk elements (= number of elements so far).
 Bulk_mesh_pt->nelement();

 // Create "surface mesh" that will contain only the bc
 // elements. The constructor just creates the mesh without
 // giving it any elements, nodes, etc.
 Surface_mesh_pt = new Mesh;

 // Add the two sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_pt);

  // Build the Problem's global mesh from its various sub-meshes
  build_global_mesh();

  oomph_info<<"Number Elements "<<Bulk_mesh_pt->nelement()<<std::endl;
  oomph_info<<"Number Nodes "<<Bulk_mesh_pt->nnode()<<std::endl;

 // Open the trace and start recording
 char filename[100];
 sprintf(filename, "RESLT/trace.dat");
 Trace_file.open(filename);

 // Complete the problem setup
 complete_problem_setup();

 // Assign the equation numbers
 assign_equation_numbers();
}

//==start_of_complete======================================================
 /// Set boundary condition exactly, and complete the build of
 /// all elements
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::complete_problem_setup()
{
 // -------------------------------------------------------------------
 // Complete the build of all elements so they are fully functional
 const unsigned n_element = Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

   //Set the pressure function pointers and the physical constants
   el_pt->nu_pt() = &TestSoln::nu;
   el_pt->pressure_fct_pt() = &TestSoln::get_pressure;
  }
 // Apply boundary conditions
 apply_boundary_conditions();
}

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
   const unsigned num_nod=Bulk_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Get node
     Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
     // Initialise node_position and exact_w
     Vector<double> node_position(2,0.0), exact_w(6,0.0);
     node_position[0] =  nod_pt->x(0);
     node_position[1] =  nod_pt->x(1);

     // If the boundary is to be pinned
     if(TestSoln::is_boundary_pinned(ibound))
      {
       // Pin the value of w at nodes
       nod_pt->pin(0);   // w
       nod_pt->set_value(0,exact_w[0]); // w

       // Pin the tangent derivatives at nodes
       // On even boundaries y axis is the tangent direction
       // On odd boundaries the x axis is the tangent direction
       const unsigned idwdt = (ibound % 2 ? 2 : 1);
       const unsigned id2wdt2 = (ibound % 2 ? 5 : 3);
       // Pin tangent derivative dwdt
       nod_pt->pin(idwdt);
       nod_pt->set_value(idwdt,exact_w[idwdt]);
       // Pin second tangent derivative d2wdy2
       nod_pt->pin(id2wdt2);
       nod_pt->set_value(id2wdt2,exact_w[id2wdt2]);
      } // for (inod<num_nod)

     // If the angle is to be set
     if(TestSoln::is_boundary_held_at_zero_slope(ibound))
      {
       // Get node
       Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
       // Pin the cross derivative
       // d2wdydx
       nod_pt->pin(4);
       nod_pt->set_value(4,exact_w[4]);

       // Pin the normal derivatives
       // On odd boundaries normal is y
       if(ibound % 2 == 1)
        {
         nod_pt->pin(2);   // dwdy
         nod_pt->set_value(2,exact_w[2]);
        }
       // On even boundaries normal is x
       if (ibound % 2 == 0) // even boundaries
        {
         nod_pt->pin(1);   // dwdx
         nod_pt->set_value(1,exact_w[1]);
        }
      } // end if
    } // end for

   // Get number of nodes on ibound
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Get node
     Node* nod_pt=this->Bulk_mesh_pt->boundary_node_pt(ibound,inod);

     // Extract nodal coordinates from node:
     Vector<double> x(2);
     x[0]=nod_pt->x(0);
     x[1]=nod_pt->x(1);
    } // for (inod<num_nod)
  } // for (ibound<Number_of_boundaries)
} // end set bc

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

  ofstream some_file;
  char filename[100];

  // Number of plot points
  unsigned npts = 2;

  sprintf(filename,"%s/coarse_soln%i.dat",
          doc_info.directory().c_str(),
          doc_info.number());
  some_file.open(filename);
  this->Bulk_mesh_pt->output(some_file,npts);
  some_file.close();

  // Output with very high resolution
  npts = TestSoln::High_resolution ? 25 : 5;
  sprintf(filename,"%s/soln%i-%f.dat",
          doc_info.directory().c_str(),
          doc_info.number(),
          Element_area);
  some_file.open(filename);
  this->Bulk_mesh_pt->output(some_file,npts);
  some_file.close();

  //  Output exact solution
  sprintf(filename, "%s/exact_sol_%i.dat",
          doc_info.directory().c_str(),
          doc_info.number());
  ofstream exact_file(filename,ios::app);
  this->Bulk_mesh_pt->output_fct(exact_file,npts,TestSoln::get_exact_w);
  exact_file.close();

  sprintf(filename,"%s/bc%i.dat",
          doc_info.directory().c_str(),
          doc_info.number());
  some_file.open(filename);
  this->Surface_mesh_pt->output(some_file,npts);
  some_file.close();

 // HERE need to sever the link or fix the BellElement class - as this class
 // BREAKS locate zeta by intoducing unitialized value through n_position_type
 // which is  =  6 for use in interpolation
 // Find the solution at x = y = 0
 MeshAsGeomObject* Mesh_as_geom_obj_pt=
  new MeshAsGeomObject(Bulk_mesh_pt);
 Vector<double> s(2);
 GeomObject* geom_obj_pt=0;
 Vector<double> r(2,0.0);
 Mesh_as_geom_obj_pt->locate_zeta(r,geom_obj_pt,s);

 // The member function does not exist in this element
 // it is instead called interpolated_u_biharmonic and returns a vector of length
 // 6 - this may need tidying up
 Vector<double> u_0(6,0.0);
 dynamic_cast<ELEMENT*>(geom_obj_pt)->interpolated_u_biharmonic(s,u_0);
 oomph_info << "w in the middle: " << u_0[0] << std::endl;

 // Get the precision
 const unsigned default_precision = Trace_file.precision();
 // Print to trace file
 Trace_file.precision(15);
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
 feenableexcept(FE_INVALID | FE_DIVBYZERO);

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Swith timings on
 Global_timings::Doc_comprehensive_timings = true;

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define possible command line arguments and parse the ones that
 // were actually specified
 string output_dir="RESLT";
 CommandLineArgs::specify_command_line_flag("--dir", &output_dir);
 double element_area=.1;
 CommandLineArgs::specify_command_line_flag("--element_area", &element_area);
 CommandLineArgs::specify_command_line_flag("--high_resolution");
 CommandLineArgs::specify_command_line_flag("--p", &TestSoln::Pressure);
 CommandLineArgs::specify_command_line_flag("--nu", &TestSoln::nu);
 CommandLineArgs::specify_command_line_flag("--length", &TestSoln::Length_of_strip);

 // Parse command line
 CommandLineArgs::parse_and_assign();

 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Output in high resolution
 TestSoln::High_resolution = CommandLineArgs::
   command_line_flag_has_been_set("--high_resolution");

 // Label for output
 DocInfo doc_info;

 // Output directory
 doc_info.set_directory(output_dir);

 // Problem instance
 // Use a 3rd order curved Bell element, we don't need to upgrade it.
 UnstructuredFvKProblem<KirchhoffPlateBendingC1CurvedBellElement> problem(element_area);
 problem.max_newton_iterations() = 1;

 // Get some timings
 double tt_start = 0.0;
 if (Global_timings::Doc_comprehensive_timings)
  {
     tt_start=TimingHelpers::timer();
  }

 // Newton Solve
 problem.newton_solve();

 //Output solution
 problem.doc_solution(doc_info);

 oomph_info << std::endl;
 oomph_info << "---------------------------------------------" << std::endl;
 oomph_info << "Pressure (" << TestSoln::Pressure << ")" << std::endl;
 oomph_info << "Poisson ratio (" << TestSoln::nu << ")" << std::endl;
 oomph_info << "Solution number (" << doc_info.number()-1 << ")" << std::endl;
 oomph_info << "---------------------------------------------" << std::endl;
 oomph_info << std::endl;

 // Document the total timing
 if (Global_timings::Doc_comprehensive_timings)
  {
   // Total time for problem
   oomph_info
    << "Total problem time: "
    << TimingHelpers::timer()-tt_start << std::endl;
  }
} //End of main
