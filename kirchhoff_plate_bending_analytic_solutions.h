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
// A file containing the analytic solutions to a number of bending problems

//Generic Routines
#include "generic.h"

// Solutions are centred on zero, with the typical size, b, set to unit length.
// Solutions are enumerated by the geometry type  and the following convention 
// for boundary conditions:
// As the operator is linear the solutions are for a unit pressure - which is
// equivalent to using a pressure and bending moment scaled deflection.
//
// For N--gons: ngon_bc1...bcn...bcN
// where bcn is one of following:
// p - pinned, zero moment (resting) 
// s - sliding, angle 0, zero shear
// c - clamped, angle 0, deflection 0
// f - free, zero moment, zero shear
//
// e.g eqtriangle_pfc is an equilateral triangle with a pinned edge a clamped
// edge and a free edge.
namespace oomph
{
namespace KirchhoffPlateBendingAnalyticSolutions  
{
 // Import Namespace MathematicalConstants (will not class in current namespace)
 using MathematicalConstants::Pi;

 // Analytic solution for opposite edges clamped and sliding repsectively
 inline void rectangle_cscs(const Vector<double>&x, const double& length_of_strip, 
  Vector<double>& w)
 {
  // Reset values just in case
  w=Vector<double>(6,0.0);
  // Variables
  double a=length_of_strip; // width b=1
  // Exact solution
  w[0]=(x[0]*x[0]-a*a/4.)*(x[0]*x[0]-a*a/4.)/24.;
  w[1]=x[0]*(x[0]*x[0]-a*a/4.)/6.;  
  w[2]=0;
  w[3]=(x[0]*x[0]/2.-a*a/24.); 
  w[4]=0;
  w[5]=0;
 } //the sliding/clamped rectangle solution

 // Analytic solution for opposite edges clamped and sliding repsectively
 inline void rectangle_psps(const Vector<double>&x, const double& length_of_strip, 
  Vector<double>& w)
 {
  // Reset values just in case
  w=Vector<double>(6,0.0);
  // Variables
  double a=length_of_strip; //Width b=1
  // Exact solution
  w[0]=(x[0]*x[0]-a*a/4)*(x[0]*x[0]-5*a*a/4)/24;
  w[1]=x[0]*(4*x[0]*x[0]-3*a*a)/24.;  
  w[2]=0;
  w[3]=(4*x[0]*x[0]-a*a)/8.; 
  w[4]=0;
  w[5]=0;
 } //the sliding/clamped rectangle solution

 // Analytic solution for a rectangle, fully pinned 
 inline void rectangle_pppp(const Vector<double>&x, const double& length_of_strip, 
  Vector<double>& w, const unsigned& terms=10)
 {
  // Reset values just in case
  w=Vector<double>(6,0.0);
  // Variables length and width
  double a=length_of_strip, b=1;
 
 for (unsigned j=0; j<terms;++j) //10 terms should be sufficient
  {
   unsigned n=2*j+1;
   for (unsigned k=0; k<terms;++k)
    {
     unsigned m=2*k+1;
     // w   - the deflection
     w[0] += 16/pow(Pi,6)*pow(-1,j+k) * cos(m*Pi*x[0]/a) * cos(n*Pi*x[1]/b)
                     / (m*n * pow((m*m/(a*a) + n*n/(b*b)),2));
     // g_x - the first x derivative
     w[1] -= 16/pow(Pi,6) * m*Pi/a * pow(-1,j+k)*sin(m*Pi*x[0]/a)
                    * cos(n*Pi*x[1]/b) / (m*n * pow((m*m/(a*a) + n*n/(b*b)),2));
     // g_y - the first y derivative
     w[2] -= 16/pow(Pi,6)* m*Pi/b * pow(-1,j+k) * cos(m*Pi*x[0]/a)
                    * sin(n*Pi*x[1]/b) / (m*n * pow((m*m/(a*a) + n*n/(b*b)),2));
  
     // g_xx - the 2nd x derivative
     w[3] -= 16/pow(Pi,6) * pow(m*Pi/a,2) * pow(-1,j+k)*cos(m*Pi*x[0]/a)
                    * cos(n*Pi*x[1]/b) / (m*n * pow((m*m/(a*a) + n*n/(b*b)),2));

     // g_xx - the 2nd y derivative
     w[5] -= 16/pow(Pi,6) * pow(m*Pi/b,2) * pow(-1,j+k)*cos(m*Pi*x[0]/a)
                    * cos(n*Pi*x[1]/b) / (m*n * pow((m*m/(a*a) + n*n/(b*b)),2));

     // g_xy - the cross derivative
     w[4] += 16/pow(Pi,6) * pow(m*Pi,2)/(b*a) * pow(-1,j+k)*sin(m*Pi*x[0]/a)
                    * sin(n*Pi*x[1]/b) / (m*n * pow((m*m/(a*a) + n*n/(b*b)),2));
/*     // q_x - the x component of the Querkraft
     w[3] -= 16/pow(Pi,6)*(pow(m/a,3) + m/a*pow(n/b,2))*pow(Pi/a,3)*pow(-1,j+k)
                    * sin(m*Pi*x[0]/a) * cos(n*Pi*x[1]/b)
                    / (m*n * pow((m*m/(a*a) + n*n/(b*b)),2));
     // q_y - the y component of the Querkraft
     w[4] -= 16/pow(Pi,6)*(pow(n/b,3) + pow(m/a,2)*n/b)*pow(Pi/a,3)*pow(-1,j+k)
                    * cos(m*Pi*x[0]/a) * sin(n*Pi*x[1]/b)
                    / (m*n * pow((m*m/(a*a) + n*n/(b*b)),2));*/
    }
  }
 } //the simply supported rectangle solution

 // Analytic solution for a rectangle, pinned on two opposing sides
 inline void rectangle_pfpf(const Vector<double>&x, const double& length_of_strip, 
  const double& nu, Vector<double>& w, const unsigned& terms=10)
 {
  //Reset values
  w=Vector<double>(6,0.0);
 
  // Variables
  double a=length_of_strip,b=1.;
  double alpha= Pi*b/(2*a); // parameters from Timoshenko
  double A_m=0, B_m=0; // initialise the horrible coefficients
 
  for (unsigned j=0; j<terms;++j)
  { 
   unsigned m=2*j+1;
   // Compute the coefficients
   A_m=4/pow(m*Pi,5)* (nu*(1+nu)*sinh(alpha*m) - nu*(1-nu)*alpha*m*cosh(alpha*m))
         / ((3+nu)*(1-nu)*sinh(alpha*m)*cosh(alpha*m) - (1-nu)*(1-nu)*alpha*m);
   B_m=4/pow(m*Pi,5)* (nu*(1-nu)*sinh(alpha*m))
         / ((3+nu)*(1-nu)*sinh(alpha*m)*cosh(alpha*m) - (1-nu)*(1-nu)*alpha*m);
  
   // w   - the deflection
   w[0]  += pow(a,4)*(4/pow(Pi*m,5) + A_m*cosh(m*Pi*x[1]/a)
              +  B_m*m*Pi*x[1]/a*sinh(m*Pi*x[1]/a))*cos(m*Pi*x[0]/a)*pow(-1,j);
   // gx  - the deflection
   w[1] -= pow(a,4)*m*Pi/a*(4/pow(Pi*m,5) + A_m*cosh(m*Pi*x[1]/a)
              +  B_m*m*Pi*x[1]/a*sinh(m*Pi*x[1]/a))*sin(m*Pi*x[0]/a)*pow(-1,j);
   // gy  - the deflection
   w[2] += pow(a,4)*m*Pi/a*((A_m + B_m)*sinh(m*Pi*x[1]/a) 
              +  B_m*m*Pi*x[1]/a*cosh(m*Pi*x[1]/a))*cos(m*Pi*x[0]/a)*pow(-1,j);
/*   // qx  - the deflection
   w[3] -= pow(a,4)*pow(m*Pi/a,3) * (-4/pow(Pi*m,5) + B_m*cosh(m*Pi*x[1]/a))
              *sin(m*Pi*x[0]/a)*pow(-1,j);
   // qy  - the deflection
   w[4] -= pow(a,4)*pow(m*Pi/a,3) * B_m*sinh(m*Pi*x[1]/a)
              *cos(m*Pi*x[0]/a)*pow(-1,j);
*/
  }
 }

 // Analytic solution for a rectangle, pinned and clamped on opposing sides 
 inline void rectangle_pcpc(const Vector<double>&x, const double& length_of_strip, 
  Vector<double>& w, const unsigned& terms=10)
 {
  // Reset values just in case
  w=Vector<double>(6,0.0);
  // Variables
  for (unsigned j=0; j<terms;++j)
   {
    double a=length_of_strip,b=1., A_m=0, B_m=0,Y_m=0;
    unsigned m = 2*j+1;
    // Local copies of coefficients nu, beta_m, and translated y to keep it short
    double am = Pi*b*m/(2*a), y = x[1];
    // Compute horrible coefficents
    A_m = - 2 *(am*cosh(am)+sinh(am))/(2*am+sinh(2*am));
    B_m = + 2 *sinh(am)/(2*am+sinh(2*am));
    // The Correction to the resting, strip solution
    Y_m = A_m*cosh(m*Pi*y/a) + B_m*m*Pi*y/a*sinh(m*Pi*y/a);
     
    // w   - the deflection
    w[0] += 4/pow(Pi,5)*pow(a,4)*(1+Y_m)*pow(m,-5)*cos(m*Pi*x[0]/a)
                  *pow(-1,j);
    // Fill in derivatives
  }
 }

// A circular sheet symmetrically pinned at k points subject to uniform pressure
// This was communicated in Bassali  Proc. Cambridge Phil. Soc. vol 53(2) p 525 
// - though there is an apparent transcription error in the final result for 
// uniformly forced circular plate.
// Include the radial derivatives
 void circle_pinned_symmetrically(const Vector<double>& polar_vector, const
unsigned& k, const double& nu, Vector<double>& w,  const unsigned& nterms=10)
 {
  // Throw for k<2
  if(k<2)
   {
   throw OomphLibError(
    "Pinned at k-points solution is valid only for k>1. Please choose a k>1.",
    OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
   }
   
  // The Polar Coordinates
  const double rho = polar_vector[0];
  const double theta = polar_vector[1];
  // A constant defined here for convenience
  const double beta = (1.-nu)/(1.+nu);
  // The fully pinned solution - as k -> infinity we will regain this solution
  w[0] = (1.-rho*rho)*((5.+nu)/(1.+nu)-rho*rho)/64.;
  w[1] = (-2*rho)*((5.+nu)/(1.+nu)-rho*rho)/64.+(1.-rho*rho)*(-2*rho)/64.;
  w[2] = 0.0;
  w[3] = (-2)*((5.+nu)/(1.+nu)-3*rho*rho)/64.+(1.-3*rho*rho)*(-2)/64.;
  w[4] = 0.0;
  w[5] = 0.0;;
  Vector<double> wcorr(6,0.0);
  
  // Now loop over the number of terms
  for(unsigned i=1;i<=nterms;++i)
   {
   // Now define the actual index (multiples of k)
   const unsigned n = k*i;
   // We need to add the terms on 
   wcorr[0] += ((n+n*beta+1.)/(n*n-1.) + beta/2.*((n*rho*rho)/(n+1.) - (n+2./beta)
        /(n-1.))*std::pow(rho,n)*std::cos(n*theta))/(n*n);

   // rho derivative 
   // Lowest k=2, then d rho^n / d rho = d rho ^2 / d rho = rho
   wcorr[1] += (beta/2.*(((n+2)*rho*rho)/(n+1.) - (n+2./beta)
        /(n-1.))*std::pow(rho,n-1)*std::cos(n*theta))/(n);

   // theta derivative 
   wcorr[2] -= (beta/2.*((n*rho*rho)/(n+1.) - (n+2./beta)
        /(n-1.))*std::pow(rho,n)*std::sin(n*theta))/(n*n);
 
   // second rho derivative - will cause arithmetic exception at rho=0 n=2
   if(n>2)
    {
    wcorr[3] += (beta/2.*(((n+2)*rho*rho) - (n-1)*(n+2./beta)
         /(n-1.))*std::pow(rho,n-2)*std::cos(n*theta))/n;
    }
   // For n=2 the second derivative should give a constant
   else
    {
    // wcorr[3] += (beta/2.*(((n+2))/(n-1.))*std::pow(rho,n)*std::cos(n*theta))/n;
    wcorr[3] += (beta/2.*(((n+2)*rho*rho) - (n-1)*(n+2./beta)
         /(n-1.))*std::cos(n*theta))/n;
    }

   // Rho theta derivative
   wcorr[4] -= (beta/2.*(((n+2)*rho*rho)/(n+1.) - (n+2./beta)
        /(n-1.))*std::pow(rho,n-1)*std::sin(n*theta));

   // theta theta derivative
   wcorr[5] -= (beta/2.*((n*rho*rho)/(n+1.) - (n+2./beta)
        /(n-1.))*std::pow(rho,n)*std::cos(n*theta));
   }

   // Construct the final solution
   w[0] += wcorr[0] / (beta *(3. + nu));
   w[1] += wcorr[1] / (beta *(3. + nu));
   w[2] += wcorr[2] / (beta *(3. + nu));
   w[3] += wcorr[3] / (beta *(3. + nu));
   w[4] += wcorr[4] / (beta *(3. + nu));
   w[5] += wcorr[5] / (beta *(3. + nu));
 }
} // end of namespace KpbAnalyticSolns
} //end of namespace oomph
