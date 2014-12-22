/**
  Program: maxpower.c
  Author: Michael F. Hutt
  http://www.mikehutt.com
  Mar. 2, 2011
  gcc -Wall -o maxpower maxpower.c nmsimplex.c -DMACOSX
  
  Example to test out nelder-mead on a one dimensional optimization problem
 *
  The circuit used for this example is a voltage divider consisting of a 1k ohm
  and 470 ohm resistor with a 9V source. The output is take across the 470 ohm resistor.
  spice file circuit description
  Vs 1 0 9
  R1 1 2 1000
  R2 2 0 470
  RL 2 0 320
 *
  The Thevenin equivalent resistance is 319.73 ohms, which is the value we should find
  from the optimization.
 *
 *
  RL:320.000000 P:0.006474
  29 Function Evaluations
  13 Iterations through program

  Ported to Dlang by Laeeth Isharc in 2014.  Use at your own risk!
 * 
 */
 
import nmsimplex;
import std.stdio;
import std.math;
 //double simplex(double function(double[] x) objfunc, double[] start,int n, double EPSILON, double scale, void function(double[] x,int n) constrain);

 double objfun(double[] x)
 {
   /* find the power for this value of RL */
   double Req;
   double VL;
   double PL;
   double Vs = 9.0;
   double RL = x[0];
   
   Req = 470*RL/(470+RL);
   VL = Vs*Req/(1000+Req);
   PL = VL*VL/RL;
   writefln("RL:%f P:%f",RL,PL);
   return -PL; /* nm finds the minimum, so to find the max we return the negative value */
 }
 
 void my_constraints(double[] x, int n)
 {
   // resistance must be positive
   int i;
   for (i=0; i<n; i++) {
     if (x[i] < 0) {
       x[i] = fabs(x[i]);
     }
   }
 }

int main()
{
	double[] start = [100];
	double min;
	int i;
  int dim = 1;
  double eps = 1.0e-8;
  double scale = 1.0;

	min=simplex(&objfun,start,dim,eps,scale,&my_constraints);

	for (i=0; i<dim; i++) {
		writefln("%f",start[i]);
	}
	return 0;
}