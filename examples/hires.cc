# include <math.h>
/*	This is a sample file for systems of ODE's,
 *	written in C++. The meaning of the functions is 
 *	as follows:
 *		1. getx0():  Return the left boundary of the equation.
 *		2. getx1():  Return the right boundary of the equation.
 *		3. getnode():Return the amount of the equations.
 *		4. getnpoints(): Return the number points in which the system
 *				will try to solve the ODE.
 *		5. systemfun(node,x,y,yy): Return the system of equations to 
 *				be solved. In parameter y we put the value of y's
 *				and in the table yy the value of the derivatives of y's.
 *		6. systemf0(node, f0): It returns in table f0 the left boundaries of 
 *				the equations.
 * */
extern "C"
{

double	getx0()
{
	return 0.0;
}

double	getx1()
{
	return 321.8122;
}

int	getnode()
{
	return 8;
}

int	getnpoints()
{
	return 50;
}

double	systemfun(int node, double x, double *y, double *yy)
{
	double p1,p2,p3,p4,p5,p6,p7,p8;
	p1=yy[0]+1.71*y[0]-0.43*y[1]-8.32*y[2]-0.0007; //yy[0]-cos(x)-y[0]*y[0]-y[1]+(x*x+sin(x)*sin(x));
	p2=yy[1]-1.71*y[0]+8.75*y[1];
	p3=yy[2]+10.03*y[2]-0.43*y[3]-0.035*y[5];
	p4=yy[3]-8.32*y[1]-1.71*y[2]+1.12*y[3];
	p5=yy[4]+1.745 *y[4]-0.43*y[5]-0.43*y[6];
	p6=yy[5]+280.0*y[5]*y[7]-0.69*y[3]-1.71*y[4]+0.43*y[5]-0.69*y[6];
	p7=yy[6]-280*y[5]*y[7]+1.81*y[6];
	p8=yy[7]+280*y[5]*y[7]-1.81*y[7];
	return pow(p1,2.0)+pow(p2,2.0)+pow(p3,2.0)+pow(p4,2.0)+
		 pow(p5,2.0)+pow(p6,2.0)+pow(p7,2.0)+pow(p8,2.0);
}

void	systemf0(int node,double *f0)
{
	f0[0]=1.0;
	f0[1]=0.0;
	f0[2]=0.0;
	f0[3]=0.0;
	f0[4]=0.0;
	f0[5]=0.0;
	f0[6]=0.0;
	f0[7]=0.0057;
}

}
