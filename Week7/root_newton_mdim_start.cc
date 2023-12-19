#include <cstdlib>
#include <iostream>
#include <cmath>
#include <cstdio>

using namespace std;

// dimension of the problem
const int dim(3);

/*******************************************************************/
/*   myfunction:                                                   */
/*     input                                                       */
/*       x: current estimate of the root (vector of length dim)    */
/*     output                                                      */
/*         F: vector of functions evaluated at x                   */
/*       Jac: Jacobi matrix of first derivatives evaluated at x    */
/*******************************************************************/
void myfunction(const double x[dim], double F[dim], double Jac[dim][dim]){

}
/*******************************************************************/
/*   solve:                                                        */
/*     input                                                       */
/*         F:  vector of functions evaluated at x                  */                         
/*       Jac: Jacobi matrix of first derivatives evaluated at x    */
/*     output                                                      */
/*       deltax: solution of J . deltax = -F                       */
/*******************************************************************/
void solve(const double F[dim], const double Jac[dim][dim], double deltax[dim]){

  //Use your solutions to previous exercises to solve
  //for deltax.

}

int main(int argc, char* argv[]){

  printf("usage %s [ start-value (%i floats) ]\n\n",argv[0],dim);
  
  // tolerance for deviation of F from zero
  const double tolf(1e-6);
  // precision for solution vector
  const double tolx(1e-6);
  // maximum number of steps
  const int MAXSTEPS(200);

  // start value for root
  double root[dim]={1,1,1};

  // alternative start values for the root from the command-line
  if(argc == dim+1){
    for(int dd(0); dd<dim; ++dd)
      root[dd] = atof(argv[dd+1]);
  }

  // value of Function, Jacobi-Matrix
  double Fvec[dim], Jac[dim][dim];

  // count steps
  int steps;

  printf("\nstart value for root:\n");
  for(int jj(0); jj<dim;++jj)
    printf("x[%i] = %e\n",jj,root[jj]);
  printf("\n");

  for( steps=0; steps<MAXSTEPS; ++steps){

    // get value of functions and Jacobi matrix at current root
    myfunction(root, Fvec, Jac);

    // ADD: calculate deviation of functions (|F|) from zero
    double errf(0.);

    printf("step %i: errf = %e", steps, errf);

    // ADD: exit loop if F approx zero

    double delta[dim];

    //solve for delta 
    solve(Fvec,Jac,delta);

    // ADD: calculate magnitude of deltax (errx) and update solution
    double errx(0.);

    printf(", errx = %e\n",errx);

    if(errx < tolx) break; // exit loop if small step
  }

  printf("\nfound root after %i steps:\n",steps);
  for(int jj(0); jj<dim;++jj)
    printf("x[%i] = %e\n",jj,root[jj]);
    
	

}
