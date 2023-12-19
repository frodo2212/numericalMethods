#include <cstdlib>
#include <iostream>
#include <cmath>
#include <cstdio>


using namespace std;


/******************************************************************/
/* function to be solved (f(x)=0), for input value x returns f(x) */
/*                                                                */
double myfunction(const double x){

  double y(0.);

  //calculate y=f(x)
  //e.g. y= x*x;

  y = 1.+cos(x)-sqrt(x);

  return y;

}

/*********************************************************************/
/* root with bisection                                               */

double bisection(const double x1, const double x2, const double xacc){

  double root(0);
  return root;

}


/*********************************************************************/
/* root with secant-method                                           */

double secant(const double x1, const double x2, const double xacc){

  const int MAXIT(40);

  double xl, root, fl, f;
  
  fl=myfunction(x1);
  f=myfunction(x2);

  if(fl*f>0) { printf("root not bracketed!\n"); exit(-1);}

  if(fl<0){ 
    xl=x1;
    root=x2;
  }
  else{
    xl=x2;
    root=x1;
    double tmp(fl);
    fl=f;
    f=tmp;
  }
  

  for(int ii(0); ii<MAXIT; ++ii){

    double dx = (root-xl)*f/(fl-f);

    xl = root;
    fl = f;

    root += dx;
    f=myfunction(root);

    if(fabs(dx) < xacc || f == 0.0) {
      printf("exit secant method in step %i\n",ii);
      return root;
    }

  }

  printf("maximum number of iterations reached in secant!\n");
  return root;

}

/*********************************************************************/
/* root with regula falsi                                            */

double regulafalsi(const double x1, const double x2, const double xacc){

  double root(0.);

  return root;
}



/***********************************************************************/

int main(){

  // start-interval (you must choose)
  double xl(-1.e-10), xh(1.e+10);
  
  // accuracy for root
  double xacc(1e-10);

  // find root with bisection
  double root_b = bisection(xl, xh, xacc);

  // find root with secant-method
  double root_s = secant(xl, xh, xacc);

  // find root with Regula Falsi
  double root_rf = regulafalsi(xl, xh, xacc);
  
  printf("root with bisection    : %e\n",root_b);
  printf("root with secant-method: %e\n",root_s);
  printf("root with regula falsi : %e\n",root_rf);
  

}


