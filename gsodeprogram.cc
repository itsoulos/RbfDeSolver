# include <gsodeprogram.h>		
# include <math.h>
# include <QDebug>

typedef double(*DOUBLE_FUNCTION)();
typedef int(*INTEGER_FUNCTION)();
# define REDO_MAX       4
# define GINF           -1e+100
# define GLAMBDA        1000.0

GSodeProgram::GSodeProgram(Rbf *m,double X0,double X1,int Node,int Npoints)
    :GProgram(m)
{
	/*	The first constructor of the class. It sets the value of
	 *	x0, x1, node and npoints. The rest parameters bust be set
	 *	by the user, with the help of other methods.
	 * */
	x0 = X0;
	x1 = X1;
	node = Node;
	if(node<0) node=1;
	f0=new double[node];
	npoints=Npoints;
	if(npoints<0) npoints=10;
	systemfun=NULL;
	fsystemfun=NULL;
	systemf0=NULL;
	fsystemf0=NULL;

}

GSodeProgram::GSodeProgram(Rbf *m,QString filename)
    :GProgram(m)
{
	/*	The second constructor of the class. It tries to open
	 *	the dll filename. If the dll can not be opened, then
	 *	the constructor sets default values for the parameters.
	 *	If the opening of dll was sucessfull, it loads all
	 *	the parameters of the SODE from the dll.
	 * */
	int i;
    ptr=new QLibrary(filename);
    if(ptr==NULL)
	{
		x0=0.0;
		x1=1.0;
		npoints=10;
		node=1;
		f0=new double[node];
		systemfun=NULL;
		fsystemfun=NULL;
		systemf0=NULL;
		fsystemf0=NULL;

	}
	else
	{
		DOUBLE_FUNCTION X0, X1;
		INTEGER_FUNCTION NODE, NPOINTS;

        X0=(DOUBLE_FUNCTION)ptr->resolve("getx0");
        if(X0==NULL) X0=(DOUBLE_FUNCTION)ptr->resolve("_getx0");
        if(X0==NULL) X0=(DOUBLE_FUNCTION)ptr->resolve("getx0_");
        if(X0==NULL) X0=(DOUBLE_FUNCTION)ptr->resolve("_getx0_");
		if(X0==NULL) x0=0.0; else x0=X0();

        X1=(DOUBLE_FUNCTION)ptr->resolve("getx1");
        if(X1==NULL) X1=(DOUBLE_FUNCTION)ptr->resolve("_getx1");
        if(X1==NULL) X1=(DOUBLE_FUNCTION)ptr->resolve("getx1_");
        if(X1==NULL) X1=(DOUBLE_FUNCTION)ptr->resolve("_getx1_");
		if(X1==NULL) x1=1.0; else x1=X1();

        NODE=(INTEGER_FUNCTION)ptr->resolve("getnode");
        if(NODE==NULL) NODE=(INTEGER_FUNCTION)ptr->resolve("_getnode");
        if(NODE==NULL) NODE=(INTEGER_FUNCTION)ptr->resolve("getnode_");
        if(NODE==NULL) NODE=(INTEGER_FUNCTION)ptr->resolve("_getnode_");
		if(NODE==NULL) node=1;else node=NODE();

        NPOINTS=(INTEGER_FUNCTION)ptr->resolve("getnpoints");
		if(NPOINTS==NULL) 
            NPOINTS=(INTEGER_FUNCTION)ptr->resolve("_getnpoints");
		if(NPOINTS==NULL) 
            NPOINTS=(INTEGER_FUNCTION)ptr->resolve("getnpoints_");
		if(NPOINTS==NULL) 
            NPOINTS=(INTEGER_FUNCTION)ptr->resolve("_getnpoints_");
		if(NPOINTS==NULL) npoints=10; else npoints=NPOINTS();

		fsystemfun=NULL;
		fsystemf0=NULL;
        systemfun=(GSYSTEMFUN)ptr->resolve("systemfun");
		if(systemfun==NULL) 
            systemfun=(GSYSTEMFUN)ptr->resolve("_systemfun");
		if(systemfun==NULL) 
            fsystemfun=(GFSYSTEMFUN)ptr->resolve("systemfun_");
		if(fsystemfun==NULL) 
            fsystemfun=(GFSYSTEMFUN)ptr->resolve("systemfun_");

	setparam=(SETPARAM)ptr->resolve("setparam");
	if(setparam==NULL)
	setparam=(SETPARAM)ptr->resolve("_setparam");
	if(setparam==NULL)
	setparam=(SETPARAM)ptr->resolve("setparam_");
		
        systemf0=(GSYSTEMF0)ptr->resolve("systemf0");
		if(systemf0==NULL) 
            systemf0=(GSYSTEMF0)ptr->resolve("_systemf0");
		if(systemf0==NULL) 
		{
            fsystemf0=(GFSYSTEMF0)ptr->resolve("systemf0_");
			if(fsystemf0==NULL)
                fsystemf0=(GFSYSTEMF0)ptr->resolve("_systemf0_");
			f0=new double[node];
			fsystemf0(&node,f0);
		}
		else
		{
			f0=new double[node];
			systemf0(node,f0);
		}

	}
}

void	GSodeProgram::setGsystemFun(GSYSTEMFUN f)
{
	/*	Set the systemfun to f, if systemfun is NULL.
	 * */
	if(systemfun==NULL) systemfun=f;
}

void	GSodeProgram::setGsystemF0(GSYSTEMF0 f)
{
	/*	Set the systemf0 to f, if systemf0 is NULL.
	 * */
	if(systemf0==NULL) 
	{
		systemf0=f;
		systemf0(node,f0);
	}
}

double	GSodeProgram::getX0() const
{
	/*	Return the left boundary of the equations.
	 * */
	return x0;
}

double	GSodeProgram::getX1() const
{
	/*	Return the right boundary of the equations.
	 * */
	return x1;
}

int	GSodeProgram::getNode() const
{
	/*	Return the amount of ODE's in the system.
	 * */
	return node;
}

int	GSodeProgram::getNpoints() const
{
	/*	Return the amount of training points.
	 * */
	return npoints;
}


double GSodeProgram::fitness(vector<double> &genome)
{
	/*	Evaluation of the fitness value associated with the 
	 *	integer array genome. Each genome codes node functions,
	 *	in its elements.
	 * */
	if(systemfun==NULL && fsystemfun==NULL) return GINF;
    int i,j;

    double y_array[node];
    double yy_array[node];
    Data x;
    x.resize(1);
	double pt = 0.0;
    Data subgenome;
    subgenome.resize(genome.size()/node);
	double value=0.0;
	/*	Fitness evaluation at the npoints points.
	 * */

	double lastValue=1e+100;
	int sameValue=0;
	for(i=0;i<npoints;i++)
	{
	    x[0]=x0+i*1.0*(x1-x0)/(npoints-1.0);
	    for(int j=0;j<node;j++)
	    {
            for(int k=0;k<subgenome.size();k++)
            {
                subgenome[k]=genome[subgenome.size()*j+k];
            }
            rbf->setVariables(subgenome);
            yy_array[j]=rbf->getDerivative(x,0);
            y_array[j]=rbf->getValue(x);
	     }


	double pf;
	    if(fsystemfun==NULL)
	    		pf=systemfun(node,x[0],y_array,yy_array);
	    else
		    	pf=fsystemfun(&node,&x[0],y_array,yy_array);
	if(i && fabs(pf-lastValue)<1e-5) sameValue++;
	lastValue = pf;
	value = value + pf;

	}
	/*	Penalty evaluation.
	 * */
	double penalty=0.0;
	for(j=0;j<node;j++)
	{
            x[0]=x0;
            for(int k=0;k<subgenome.size();k++)
            {
                subgenome[k]=genome[subgenome.size()*j+k];
            }
            rbf->setVariables(subgenome);
            y_array[j]= rbf->getValue(x);
	    	penalty+=GLAMBDA * pow(y_array[j]-f0[j],2.0);
	}
                if(std::isnan(value) || std::isinf(value)) return GINF;
                if(std::isnan(penalty) || std::isinf(penalty)) return GINF;
	//printf("value: %lf penalty; %lf \n",value,penalty);
    value = value +(penalty);
	return -value*(1.0+pt);//sameValue*1.0/npoints);
}

GSodeProgram::~GSodeProgram()
{
	/*	The destructor of the class. It deallocates the memory of the
	 *	class.
	 * */
    if(ptr!=NULL) delete ptr;
	delete[] f0;
}
