# include <math.h>
/*	This is a sample file for ODE,
 *	written in C++. The meaning of the functions is 
 *	as follows:
 *		1. getx0():  Return the left boundary of the equation.
 *		2. getx1():  Return the right boundary of the equation.
 *		3. getkind():Return the kind of the equation.
 *		4. getnpoints(): Return the number points in which the system
 *				will try to solve the ODE.
 *		5. getf0():  Return the left boundary condition.
 *		6. getf1():  Return the right boundary condition.
 *		7. getff0(): Return the derivative of the left boundary condition.
 *		8. ode1ff(): Return the equation for first order equations.
 *		9. ode2ff(): Return the equation for second order equations.
 * */
/*	KINDS:
 *	=======================================================================
 * 	ODE1:	ODE of first order, with the boundary condition y(x0)=f0.
 *	ODE2:	ODE of second order, with the boundary conditions y(x0)=f0
 *		and y'(x0)=ff0.
 *	ODE3:	ODE of second order, with the boundary conditions y(x0)=f0
 *		and y(x1)=f1.
 * */
extern "C"
{

double	getx0()
{
	return 0.1;
}

double	getx1()
{
	return 1.0;
}

int	getkind()
{
	return 1;
}

int	getnpoints()
{
	return 10;
}

double	getf0()
{
	return 2.1/sin(0.1);
}

double	getff0()
{
	return 20.1;
}

double	getf1()
{
	return sin(10.0);
}

double	ode1ff(double x,double y,double yy)
{
	return yy-(1-y*cos(x))/sin(x);
}

double	ode2ff(double x,double y,double yy,double yyy)
{
	return yyy+100.0 * y;
}

}
