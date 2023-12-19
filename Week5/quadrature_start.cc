#include <cstdlib>
#include <iostream>
#include <cmath>
#include <cstdio>

using namespace std;

/*************************************************************/
/* function to be integrated, for input value x returns f(x) */
/*                                                           */
/*************************************************************/
double myfunction(const double x){

  double y(0.);
  return y;
}


/**********************************************************************/
/*                                                                    */
/* calculate next contribution for Trapez-rule from the new points    */
/*                                                                    */
/*   a,b: boundaries for integration                                  */
/*   m  : step                                                        */
/*                                                                    */
/*   in step m: add 2^(m-1) new points in [a,b]                       */
/**********************************************************************/
double new_trapez(const double a, const double b, const int m){

  double sum(0.);
  return sum;

}


/***************************************************************************/
/*                                                                         */
/* Quadrature using Trapez rule                                            */
/*                                                                         */
/* a,b: boundaries for integration                                         */
/* eps: precision                                                          */
/*                                                                         */
/* calculation done by subsequent calls to new_trapez                      */
/*                                                                         */
/***************************************************************************/
double quadrature_trapez(const double a, const double b, const double eps)
{
  // maximum number of steps
  const int MAXSTEP(20);

  // current and previous result
  double sum(0.), oldsum(0.);

  // start values
  //oldsum = sum = ????;

  // do steps
  for(int ii(1); ii<MAXSTEP; ++ii){
    
    // sum = ????; 

    // after minimum number of steps, check if precision reached
    if( ii > 5 )
      if( ( sum == 0.0 && oldsum == 0.0 ) || fabs(sum-oldsum) < eps*fabs(oldsum) ){
        printf("quadrature_trapez_stop stopped after %i steps\n",ii);
        return sum;
      }
    oldsum = sum;

  }

  printf("Maximum number of steps reached in quadrature_trapez_stop!\n");
  return sum;
  
}


/**************************************************************************/
/*                                                                        */
/* quadrature using Simpson's rule                                        */
/*                                                                        */
/*   a,b: boundaries for integration                                      */
/*   eps: precision                                                       */
/*                                                                        */
/*   using new_trapez to calculate new contribution                       */
/**************************************************************************/

double quadrature_simpson(const double a, const double b, const double eps)
{
  double sum(0.);
  return sum;

}



/**************************************************************************/
/*                                                                        */
/* quadrature using Romberg-Integration                                   */
/*                                                                        */
/*   uses  new_trapez for calculating new points of trapez-rule           */
/*    and  additional routine for extrapolation to zero h                 */
/*   e.g. an implementation of Neville's algorithm                        */
/*   a,b: boundaries for integration                                      */
/*   eps: precision                                                       */
/**************************************************************************/

double quadrature_romberg(const double a, const double b, const double eps){

  double sum(0.);
  return sum;
}

/*************************************************************************/

int main(){


  // boundaries for integration
  const double a (0.0), b(2.0);
  // precision for quadrature
  const double eps(1e-10);


  printf("quadrature (%e) with Trapez-rule: %e\n",
	 eps,quadrature_trapez(a,b,eps));

  printf("quadrature (%e) with Simpson's-rule: %e\n",
	 eps,quadrature_simpson(a,b,eps));

  printf("quadrature (%e) with Romberg-Integration:  %e\n",
	 eps,quadrature_romberg(a,b,eps));
  
}
