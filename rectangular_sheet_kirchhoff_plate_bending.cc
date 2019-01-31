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
 // Enumerated type describing different boundary cases we may want to test
 enum Boundary_case {
   all_pinned = 0, 
   three_pinned = 1, 
   opposite_pinned_opposite_free = 2,
   corners_pinned = 3, 
   opposite_pinned_opposite_clamped = 4,
   opposite_clamped_opposite_sliding = 5, 
   opposite_pinned_opposite_sliding = 6,
   all_clamped= 7,
   opposite_clamped_opposite_free = 8 ,
   opposite_pinned_free_clamped = 9 ,
   opposite_free_opposite_pinned = 10
   };

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

 // The problem we are solving 
 Boundary_case boundary_case = opposite_pinned_opposite_free;

 /// Nondim length of strip (width is lengthscale and therefore 1)
 // length is in the y-direction, width=1 in y
 double Length_of_strip=2.0; 

 /// Poisson's ratio
 double nu = 0.3; /// hierher FvK parameter (Poisson's ratio)
 
 /// The imposed angle on angle controlled boundaries
 double imposed_angle =0.0; 

 /// Current overall x extension (total increment in length imposed
 /// by moving the two ends of the strip apart -- in non-dim units!)
 double Displacement_x = 0.0;
 
 double Kirchhoff_shear = 0.0;

 double exact_soln_max_terms_in_expansion=50;
 
  /// Perturbation
 double P_cos=1.0;
  
 double Applied_normal_moment = 0.0;
 
 /// Wavenumber of perturbation (10 from Fig. 1 in PRL paper)
 unsigned N_cos = 10;

 /// sheet normals
 void sheet_normal_function(const Vector<double>& x, Vector<double>& n1,
  Vector<double>& n2)
  {
   // Test case 1 - just the same as the axis-aligned case
   n1[0]=1.0 ; n1[1]= 1.0;
   n2[0]=1.0 ; n2[1]=-1.0;
  }

 /// Bool to determine which boundaries require traction elements
 bool is_boundary_moment_free(unsigned& ibound)
 {
  switch (TestSoln::boundary_case)
   {
    // Cases with all four boundaries that are moment controlled.
    case TestSoln::all_pinned   : 
    case TestSoln::opposite_pinned_opposite_free : 
    case TestSoln::three_pinned : case TestSoln::corners_pinned :
    case TestSoln::opposite_free_opposite_pinned :
      return true ;  break;

    // Cases with Opposite angle controlled opposite moment controlled.
    case TestSoln::opposite_pinned_opposite_sliding :
    case TestSoln::opposite_pinned_opposite_clamped : 
      return ibound==1 || ibound==3; break;
    case TestSoln::opposite_clamped_opposite_free : 
      return ibound==0 || ibound==2; break;

    // Cases with all Angle imposed boundary conditions.
    case TestSoln::opposite_clamped_opposite_sliding :
    case::TestSoln::all_clamped :
      return false; break;

    // Default to all pinned
    default : 
      oomph_info <<"Boundary case not recognised. Defaulting to all pinned."; 
      return true; 
   }
 }

 /// Bool to determine which boundaries require traction elements
 bool is_boundary_at_zero_slope(unsigned& ibound)
 {
  // Its an angle boundaries if it isn't traction!
  return ! is_boundary_moment_free(ibound);  
 }

 /// Bool to determine which boundaries to pin
 bool is_boundary_pinned(unsigned& ibound)
 {
  switch (TestSoln::boundary_case)
   {
    // All boundaries pinned
    case TestSoln::opposite_pinned_opposite_clamped : case TestSoln::all_pinned :
    case TestSoln::all_clamped :
      return true;  break; 
    // All except boundary 1 pinned, 1 free
    case TestSoln::three_pinned : 
    case TestSoln::opposite_pinned_free_clamped : 
      return  ibound != 1 ; break;
    //  Opposite pinned (1 and 3) opposite sliding/free
    case TestSoln::opposite_clamped_opposite_sliding : 
    case TestSoln::opposite_pinned_opposite_sliding : 
    case TestSoln::opposite_pinned_opposite_free : 
    case TestSoln::opposite_clamped_opposite_free : 
      return ibound == 1 || ibound == 3; break;
    case TestSoln::opposite_free_opposite_pinned :
      return ibound == 0 || ibound == 2; break;
    // All free/sliding, pin corners in apply_boundary_conditions
    case TestSoln::corners_pinned : 
      return false ;  
    // Default to all_pinned
    return true; break;
    // Default to all pinned
    default :  
      oomph_info <<"Boundary case not recognised. Defaulting to all pinned."; 
      return true;
   }
 }

 // The usage cases with no analytic solution
 bool has_no_analytic_solution()
 {
  return boundary_case == three_pinned || 
         boundary_case ==  corners_pinned ||
         all_clamped == boundary_case      ||
         opposite_clamped_opposite_free==boundary_case ||
         opposite_pinned_free_clamped==boundary_case ||
         opposite_free_opposite_pinned==boundary_case;
 }

  // Assigns the value of pressure depending on the position (x,y)
  void get_pressure(const Vector<double>& x, double& pressure)
  {
   // Return the (constant) pressure
   pressure = P_cos; 
  }
 
  // Get the exact solution
  void get_exact_w(const Vector<double>& x, Vector<double>& exact_w)
  {
    // We assume the dimensions of the strip are a, b with b=1
    double a=Length_of_strip;
    // Copy of the position Vector
    Vector <double> xp(x);
    // Swap the coordinates, in line with analytic solution
    xp[0]=x[1]; xp[1]=x[0];

    // Locally use this namespace
    using namespace KirchhoffPlateBendingAnalyticSolutions;
    // Switch for theboundary cases
    switch(TestSoln::boundary_case){
     case all_pinned:
       // Get the analytic solution centred on zero-zero
       rectangle_pppp(xp,Length_of_strip,exact_w);  break;

     case opposite_pinned_opposite_free:
       // Get the solution
       rectangle_pfpf(xp,Length_of_strip,nu,exact_w); break;

     case opposite_pinned_opposite_clamped:
       // Get the solution
       rectangle_pcpc(xp,a,exact_w);  break;

     case opposite_clamped_opposite_sliding:
       // Get the solution
       rectangle_cscs(xp,a,exact_w);  break;

     case opposite_pinned_opposite_sliding:
       // Get the solution
       rectangle_psps(xp,a,exact_w);  break;

     default: break;
    }
  // Analytic solutions are for p=1, but because it is linear we can just rescale
  for(unsigned ideriv=0;ideriv<exact_w.size();++ideriv)
   { exact_w[ideriv]*=P_cos;}
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
   Trace_file.close();  
   delete this->Surface_mesh_pt;
   delete Bulk_mesh_pt;
   delete this->eigen_solver_pt();
  };
 
  /// Update after solve (empty)
  void actions_after_newton_solve() { }
  
  /// Update the problem specs before solve: Re-apply boundary conditons
  void actions_before_newton_solve()
  {
    complete_problem_setup();
    apply_boundary_conditions();
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
 
 void output_boundary_nodal_values() 
  {  
   // Loop over outer boundary and pin the nodes on the boundaries
   const unsigned nbound = 4;
   
   // loop over boundaries
   for(unsigned ibound=0;ibound<nbound;ibound++)
    {
     oomph_info<<"### Boundary "<<ibound<<" ###\n";
     const unsigned num_nod=Bulk_mesh_pt->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       // Get node
       Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
       oomph_info<<"Node at: ("<<nod_pt->x(0)<<","<<nod_pt->x(1)<<")\n";
       for(unsigned i=0;i<6;i++)
        { oomph_info<<nod_pt->value(i)<<(i==5?"\n":","); }
      } // for (inod<num_nod)
      
    } // end loop over boundaries
  }
private:
 
  /// Trace file to document norm of solution
  ofstream Trace_file;

  /// Helper function to apply boundary conditions
  void apply_boundary_conditions();

  /// \short Helper function to (re-)set boundary condition
  /// and //complete the build of  all elements
  void complete_problem_setup();
  
  /// Pointers to specific mesh
  TriangleMesh<ELEMENT>* Bulk_mesh_pt;
      
  /// Pointer to the "surface" mesh
  Mesh* Surface_mesh_pt;
  
protected:
    
  // The initial (and maximum) element area
  double Element_area;
 
}; // end_of_problem_class


template<class ELEMENT>
UnstructuredFvKProblem<ELEMENT>::UnstructuredFvKProblem(double element_area)
  :
Element_area(element_area)
{
 // Set the maximum residuals for Newton iterations
 Problem::Max_residuals = 1000;
 
 // Setup mesh
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
 TriangleMeshPolyLine *boundary0_pt =
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
 TriangleMeshPolyLine *boundary1_pt =
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
 TriangleMeshPolyLine *boundary2_pt =
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
 TriangleMeshPolyLine *boundary3_pt =
   new TriangleMeshPolyLine(vertices, boundary_id);
 
 // ---------------------------------------------------------------------
 // >> Building the OUTER BOUNDARY
 
 // >> Setting up the domain with PolyLines
 Vector<TriangleMeshCurveSection*> outer_boundary_polyline_pt(4);
 
 outer_boundary_polyline_pt[0] = boundary0_pt;
 outer_boundary_polyline_pt[1] = boundary1_pt;
 outer_boundary_polyline_pt[2] = boundary2_pt;
 outer_boundary_polyline_pt[3] = boundary3_pt;
 
 // The outer polygon
 TriangleMeshClosedCurve* outer_boundary_pt = 
  new TriangleMeshClosedCurve(outer_boundary_polyline_pt);
 
 // --------------------------------------------------------------------
 //Create the mesh
 //---------------
 
 //Create mesh parameters object
 TriangleMeshParameters triangle_mesh_parameters(outer_boundary_pt);

 // Set the element area into the mesh arguments
 triangle_mesh_parameters.element_area() = element_area;
 
 Bulk_mesh_pt = new TriangleMesh<ELEMENT>(triangle_mesh_parameters);
  
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
 
 // Open the trace and start recording 
 char filename[100];
 sprintf(filename, "RESLT/trace.dat");
 Trace_file.open(filename);

 // Complete the problem setup
 complete_problem_setup();
 
 // Assign equations numbers
 oomph_info << "Number of equations: "
            << this->assign_eqn_numbers() << '\n';  
 
 this->eigen_solver_pt() = new LAPACK_QZ; 
 
 // Document the number of elements in the mesh
 oomph_info << "Number of elements: " << Bulk_mesh_pt->nelement() << std::endl;

 // Clean up
 delete outer_boundary_pt;
 outer_boundary_pt = 0;
 delete boundary0_pt;
 boundary0_pt = 0;
 delete boundary1_pt;
 boundary1_pt = 0;
 delete boundary2_pt;
 boundary2_pt = 0;
 delete boundary3_pt;
 boundary3_pt = 0;
}

//==start_of_complete======================================================
 /// Set boundary condition exactly, and complete the build of 
 /// all elements
//========================================================================
template<class ELEMENT>
void UnstructuredFvKProblem<ELEMENT>::complete_problem_setup()
{  
 // Apply boundary conditions
 apply_boundary_conditions();
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
 // Print helpful message
 oomph_info<<"Applying boundary conditions to the sheet."<<std::endl;
 // Loop over outer boundary and pin the nodes on the boundaries
 const unsigned nbound = 4;
 
 // loop over boundaries
 for(unsigned ibound=0;ibound<nbound;ibound++)
  {
   // If the boundary is to be pinned
   if(TestSoln::is_boundary_pinned(ibound))
    {
     const unsigned num_nod=Bulk_mesh_pt->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       // Get node
       Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
       
       // Pin unknown values
       nod_pt->pin(0);   // w
       nod_pt->set_value(0,0.0); // w
       
       // Pin the tangent derivatives
       // On even boundaries y is tangent direction
       if(ibound % 2 == 0)         
        {
         // Pin tangent derivative dwdy
         nod_pt->pin(2);   
         nod_pt->set_value(2,0.0); 
         // Pin second tangent derivative d2wdy2
         nod_pt->pin(5);
         nod_pt->set_value(5,0.0); 
        }
       // On odd boundaries x is tangent direction
       else /*if (ibound % 2 ==1)*/
        {
         // Pin tangent derivative dwdx
         nod_pt->pin(1);
         nod_pt->set_value(1,0.0); 
         // Pin second tangent derivative d2wdx2
         nod_pt->pin(3);
         nod_pt->set_value(3,0.0); 
        }            
      } // for (inod<num_nod)
   }
 
   // If the angle is to be set
   if(TestSoln::is_boundary_at_zero_slope(ibound))
    {
     const unsigned num_nod=Bulk_mesh_pt->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       // Get node
       Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
       // Pin the cross derivative 
       // d2wdydx
       nod_pt->pin(4);   
       nod_pt->set_value(4,0.0);
       
       // Pin the normal derivatives
       // On odd boundaries normal is y
       if(ibound % 2 == 1) 
        {
         nod_pt->pin(2);   // dwdy
         nod_pt->set_value(2,0.0);
        }
       // On even boundaries normal is x
       if (ibound % 2 == 0) // even boundaries
        {
         nod_pt->pin(1);   // dwdx
         nod_pt->set_value(1,0.0);
        }
      } // for (inod<num_nod)
    } // end if 
    
   // Get number of nodes on ibound 
   const unsigned num_nod=this->Bulk_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Get node
     Node* nod_pt=this->Bulk_mesh_pt->boundary_node_pt(ibound,inod);
            
     // Extract nodal coordinates from node:
     Vector<double> x(2);
     x[0]=nod_pt->x(0);
     x[1]=nod_pt->x(1);

     // Corners - if we are on two boundaries we are a corner
     // don't count backward as well otherwise we will double count
     if(nod_pt->is_on_boundary(ibound) && nod_pt->is_on_boundary((ibound+1) % 4))
     {
      oomph_info<<"Found corner point\n";
      if( TestSoln::boundary_case == TestSoln::corners_pinned) // pin w on corners
       {
      nod_pt->pin(0);
      // Clamped boundary conditions w = 0 at boundary
      nod_pt->set_value(0,0.0); // w
       }
     }
    } // for (inod<num_nod)
  } // for (ibound<nbound)
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


  npts = 5;
  sprintf(filename,"%s/soln%i-%i-%f.dat",
          doc_info.directory().c_str(),
          doc_info.number(),TestSoln::boundary_case,
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
 Vector<double> u_0(12,0.0);
 dynamic_cast<ELEMENT*>(geom_obj_pt)->interpolated_u_biharmonic(s,u_0);
 oomph_info << "w in the middle: " << std::setprecision(15) << u_0[0] << std::endl;

 // Print to trace file
 Trace_file << std::setprecision(15) << u_0[0] << std::endl;

 // Doc error and return of the square of the L2 error
 //---------------------------------------------------
 if(!TestSoln::has_no_analytic_solution())
 {
 double dummy_error,zero_norm;
 sprintf(filename,"%s/error%i-%f.dat",
         doc_info.directory().c_str(),
         doc_info.number(),
         Element_area);
 some_file.open(filename);

 Bulk_mesh_pt->compute_error(some_file,TestSoln::get_exact_w,
                           dummy_error,zero_norm);

 oomph_info<<"L2 Error w : "<<sqrt(dummy_error)<<"\n";
 some_file.close();

 // Doc error and return of the square of the L2 error
 //---------------------------------------------------
 sprintf(filename,"%s/L2-norm%i-%f.dat",
         doc_info.directory().c_str(),
         doc_info.number(),
         Element_area);
 some_file.open(filename);
 
 some_file<<"### L2 Norm\n";
 some_file<<"##  Format: err^2 norm^2 log(err/norm) \n";
 // Print error in prescribed format
  some_file<< dummy_error <<" "<< zero_norm <<" ";

 // Only divide by norm if its nonzero
 if(dummy_error!=0 && zero_norm!=0)
  some_file<< 0.5*(log10(dummy_error)-log10(zero_norm))<<"\n";
 some_file.close();
 }
 else
 {
  oomph_info <<"No analytic solution defined, skipping error calculation."<<std::endl;
 }

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
 
 // The `usage' flag
 CommandLineArgs::specify_command_line_flag("--usage");
 // Output directory
 string output_dir="RESLT";
 CommandLineArgs::specify_command_line_flag("--dir",
                                            &output_dir);
 // 1Pcos
 CommandLineArgs::specify_command_line_flag("--p", &TestSoln::P_cos);

 // -------------------------------------------------------------------------
 // Adaptation parameters
 // -------------------------------------------------------------------------
 
 // The element size
 double element_area=.1;
 CommandLineArgs::specify_command_line_flag("--element_area", &element_area);
 
 // -------------------------------------------------------------------------
 // Problem/Test parameters
 // -------------------------------------------------------------------------
 
 // Number of steps
 unsigned n_step=1;
 CommandLineArgs::specify_command_line_flag("--n_step", &n_step);
 
 // Increase in pressure dp
 double dux=0.0;
 CommandLineArgs::specify_command_line_flag("--dp", &dux);

 // Set the Poisson's ratio
 CommandLineArgs::specify_command_line_flag("--nu", &TestSoln::nu);

 // Boundary case as an unsigned
 int boundary_case_num=0;
 CommandLineArgs::specify_command_line_flag("--case", &boundary_case_num);

 // -------------------------------------------------------------------------
 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();
 // -------------------------------------------------------------------------
 // Flag for bad user input
 bool invalid_input=false;

 // Assign the boundary case
 if (boundary_case_num==int(TestSoln::all_pinned))
  TestSoln::boundary_case=TestSoln::all_pinned;

 else if (boundary_case_num==int(TestSoln::three_pinned))
  TestSoln::boundary_case=TestSoln::three_pinned;

 else if (boundary_case_num==int(TestSoln::opposite_pinned_opposite_free))
  TestSoln::boundary_case=TestSoln::opposite_pinned_opposite_free;

 else if (boundary_case_num==int(TestSoln::corners_pinned))
  TestSoln::boundary_case=TestSoln::corners_pinned;

 else if (boundary_case_num==int(TestSoln::opposite_pinned_opposite_clamped))
  TestSoln::boundary_case=TestSoln::opposite_pinned_opposite_clamped;

 else if (boundary_case_num==int(TestSoln::opposite_clamped_opposite_sliding))
  TestSoln::boundary_case=TestSoln::opposite_clamped_opposite_sliding;

 else if (boundary_case_num==int(TestSoln::opposite_pinned_opposite_sliding))
  TestSoln::boundary_case=TestSoln::opposite_pinned_opposite_sliding;

 else if (boundary_case_num==int(TestSoln::all_clamped))
  TestSoln::boundary_case=TestSoln::all_clamped;

 else if (boundary_case_num==int(TestSoln::opposite_clamped_opposite_free))
  TestSoln::boundary_case=TestSoln::opposite_clamped_opposite_free;

 else if (boundary_case_num==int(TestSoln::opposite_free_opposite_pinned))
  TestSoln::boundary_case=TestSoln::opposite_free_opposite_pinned;

 else if (boundary_case_num==int(TestSoln::opposite_pinned_free_clamped))
  TestSoln::boundary_case=TestSoln::opposite_pinned_free_clamped;

 else // Default to what is set in TestSoln
  { 
   oomph_info<<"Boundary case \""<<boundary_case_num<<"\" not recognised.\n";
   invalid_input=true;
  }

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
  oomph_info<<"                      0 - all_pinned\n";
  oomph_info<<"                      1 - three_pinned\n";
  oomph_info<<"                      2 - opposite_pinned_opposite_free\n";
  oomph_info<<"                      3 - corners_pinned\n";
  oomph_info<<"                      4 - opposite_pinned_opposite_clamped\n";
  oomph_info<<"                      5 - opposite_clamped_opposite_sliding\n";
  oomph_info<<"                      6 - opposite_pinned_opposite_sliding\n";
  oomph_info<<"                      7 - all_clamped\n";
  oomph_info<<"                      8 - opposite_clamped_opposite_free\n";
  oomph_info<<"                      9 - opposite_pinned_free_clamped\n";
  oomph_info<<"                     10 - opposite_free_opposite_pinned\n";
                 
  // Terminate here
  return(0);
 }

 // Label for output
 DocInfo doc_info;
 
 // Output directory
 doc_info.set_directory(output_dir);

 // Problem instance
 // Use a 3rd order curved Bell element, we don't need to upgrade it.
 UnstructuredFvKProblem<KirchhoffPlateBendingC1CurvedBellElement<2,2,3> > problem(element_area);

 problem.max_newton_iterations()=1;
 // Get some timings
 double tt_start = 0.0;
 if (Global_timings::Doc_comprehensive_timings)
  {
     tt_start=TimingHelpers::timer();
  }

 //Output solution
 problem.doc_solution(doc_info);
 
 // Loop over the number of steps
 for (unsigned s = 0; s < n_step; s++)
  {
   // Newton Solve
   problem.newton_solve();
 
   //Output solution
   problem.doc_solution(doc_info);
   
   oomph_info << std::endl;
   oomph_info << "---------------------------------------------" << std::endl;
   oomph_info << "Doced with displacement Ux (" << TestSoln::Displacement_x/2 << ")" << std::endl;
   oomph_info << " Pcos (" << TestSoln::P_cos << ")" << std::endl;
   oomph_info << "Current dux (" << dux << ")" << std::endl;
   oomph_info << "Poisson ratio (" << TestSoln::nu << ")" << std::endl;
   oomph_info << "Solution number (" << doc_info.number()-1 << ")" << std::endl;
   oomph_info << "---------------------------------------------" << std::endl;
   oomph_info << std::endl;
   
   // Increase the displacement in x and go for the next one
   TestSoln::P_cos+=dux;   
  } // for (s < n_step)
 
 // Document the total timing
 if (Global_timings::Doc_comprehensive_timings)
  {
   // Total time for problem
   oomph_info
    << "Total problem time: " 
    << TimingHelpers::timer()-tt_start << std::endl;
  }
 
} //End of main
