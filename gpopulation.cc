# include <gpopulation.h>	
# include <math.h>
# include <stdio.h>
# include <iostream>
# include <omp.h>
# define MAX_RULE	256

GPopulation::GPopulation(int gcount,int gsize,vector<GProgram*> p)
{

    /*	Default values for the parameters
     *	of the population.
     * */
    selection_rate  = 0.1;
    mutation_rate   = 0.1;
    genome_count    = gcount;
    genome_size     = gsize;
    generation      = 0;
    ompProgram         = p;
    tournament_size =10;

    /*	Memory allocation.
     * */
    genome.resize(genome_count);
    children.resize(genome_count);

    /*	Random initialization of the population.
     * */
    for(int i=0;i<genome_count;i++)
    {
        genome[i].resize(genome_size);
        children[i].resize(genome_size);
        for(int j=0;j<genome_size;j++)
        {
            genome[i][j]=10* drand48();//0.1* rand()*1.0/RAND_MAX;
        }
    }
    fitness_array.resize(genome_count);
    calcFitnessArray();

}
GPopulation::GPopulation(int gcount,int gsize,GProgram *p)
{
	/*	Default values for the parameters
	 *	of the population.
	 * */
	selection_rate  = 0.1;
	mutation_rate   = 0.1;
	genome_count    = gcount;
	genome_size     = gsize;
	generation      = 0;
	program         = p;
	tournament_size =10;

    /*	Memory allocation.
     * */
    genome.resize(genome_count);
    children.resize(genome_count);
	
	/*	Random initialization of the population.
	 * */
	for(int i=0;i<genome_count;i++)
	{
        genome[i].resize(genome_size);
        children[i].resize(genome_size);
		for(int j=0;j<genome_size;j++)
		{
            genome[i][j]=0.01* drand48();//0.1* rand()*1.0/RAND_MAX;
	    
        }
	}
    fitness_array.resize(genome_count);
    calcFitnessArray();

}

class RbfProblem :public Problem
{
    public:
        GProgram *program;
        virtual double  funmin(Matrix);
};

double RbfProblem::funmin(Matrix x)
{
    return -program->fitness(x);
}


void	GPopulation::localSearch(int pos)
{
    vector<double> g;
	g.resize(genome_size);
	
    RbfProblem *pt=new RbfProblem();
    Data xl,xr;
    xl.resize(g.size());
    xr.resize(g.size());
    for(int i=0;i<g.size();i++)
    {
        xl[i]=0;
        xr[i]= 100;
    }
    pt->setLeftMargin(xl);
    pt->setRightMargin(xr);
    extern int threads;
    if(threads<=1)
    pt->program=program;
    else pt->program=ompProgram[omp_get_thread_num()];
    pt->setDimension(genome[0].size());
    MinInfo Info;
    Info.p= pt;
    Info.iters=2001;
    double ff=tolmin(genome[pos],Info);
    delete pt;
    fitness_array[pos]=-ff;
    return;
    /*for(int i=0;i<genome_size;i++)
    {
        int ipos = rand() % genome_size;
        int new_value;
        for(int k=0;k<20;k++)
        {
        int old_value = genome[pos][ipos];
        new_value =  old_value * (1.0+0.05 *(2.0*drand48()-1.0));
        genome[pos][ipos]=new_value;
        for(int j=0;j<genome_size;j++) g[j]=genome[pos][j];
        double trial_fitness=fitness(g);
        if(fabs(trial_fitness)<fabs(fitness_array[pos]))
        {
            fitness_array[pos]=trial_fitness;
    //		return;
        }
        else	genome[pos][ipos]=old_value;
        }
    }*/

    for(int iters=1;iters<=100;iters++)
    {
        int gpos=rand() % genome_count;
        int cutpoint=rand() % genome_size;
        for(int j=0;j<cutpoint;j++) g[j]=genome[pos][j];
        for(int j=cutpoint;j<genome_size;j++) g[j]=genome[gpos][j];
        double f=fitness(g);
        if(fabs(f)<fabs(fitness_array[pos]))
        {
            for(int j=0;j<genome_size;j++) genome[pos][j]=g[j];
            fitness_array[pos]=f;
//			return;
        }
        else
        {
            for(int j=0;j<cutpoint;j++) g[j]=genome[gpos][j];
            for(int j=cutpoint;j<genome_size;j++) g[j]=genome[pos][j];
            double f=fitness(g);
            if(fabs(f)<fabs(fitness_array[pos]))
            {
                for(int j=0;j<genome_size;j++) genome[pos][j]=g[j];
                fitness_array[pos]=f;
//			return;
            }
        }
    }

}

double 	GPopulation::fitness(vector<double> &g)
{
	/*	Return the fitness for the genome g. 
	 *	The value is taken from the program pointed
	 *	by program.
	 * */

    for(int i=0;i<g.size();i++)
        if(g[i]>100 || g[i]<-100) return -1e+100;
    extern int threads;
    if(threads<=1)
	return program->fitness(g);
    else
    {
        int tpos= omp_get_thread_num();
        return ompProgram[tpos]->fitness(g);
    }
}

void	GPopulation::select()
{
	/*	Sort the genome and the fitness_array table according
	 *	to the value of fitness_array. The sorting is performed
	 *	in a way that the best genomes go to the beginning of
	 *	the table genome. The sorting is performed with a 
	 *	simple bubble sort.
	 * */
    Data itemp;
    itemp.resize(genome_size);
	for(int i=0;i<genome_count;i++)
	{
		for(int j=0;j<genome_count-1;j++)
		{
			if(fitness_array[j+1]>fitness_array[j])
			{
				double dtemp;
				dtemp=fitness_array[j];
				fitness_array[j]=fitness_array[j+1];
				fitness_array[j+1]=dtemp;
				
                itemp=genome[j];
                genome[j]=genome[j+1];
                genome[j+1]=itemp;
			}
		}
	}
}

void	GPopulation::crossover()
{
	/*	parent:	A table that holds the indexes of the parents
	 *		in the genome array for every crossover.
	 *	nchildren: The amount of the individuals that the 
	 *		crossover procedure must create.
	 *	count_childrem: An auxiliary variable which counts 
	 *		the produced children.
	 * */
        int parent[2];
        int nchildren=(int)((1.0 - selection_rate) * genome_count);
	if(!(nchildren%2==0)) nchildren++;
        int count_children=0;
        while(1)
        {
		/*	Select two parents for crossover. The selection is
		 *	performed with the method of tournament selection.
		 * */
                for(int i=0;i<2;i++)
                {
			
                        double max_fitness=-1e+10;
                        int    max_index=-1;
			int r;
                        for(int j=0;j<tournament_size;j++)
                        {
				r=rand() % (genome_count);
                                if(j==0 || fitness_array[r]>max_fitness)
                                {
                                        max_index=r;
                                        max_fitness=fitness_array[r];
                                }
                        }
                        parent[i]=max_index;
			
                }




        for(int k=0;k<genome_size;k++)
        {
        double x1=genome[parent[0]][k];
        double x2=genome[parent[1]][k];
               double alfa=-0.5+2.0*drand48();
               double g1=alfa*x1+(1.0-alfa)*x2;
               double g2=alfa*x2+(1.0-alfa)*x1;
	       if(g1<0) g1=-g1;
	       if(g2<0) g2=-g2;
               children[count_children][k]=g1;
               children[count_children+1][k]=g2;
        }


		/*	Advance the children counter
		 * */
		count_children+=2;

		if(count_children>=nchildren) break;
	}
	
	
	/*	Replace the nchildren wors individuals with the children
	 *	produced from the crossover procedure.
	 * */
	for(int i=0;i<nchildren;i++)
	{
        genome[genome_count-i-1]=children[i];
	}
}

void	GPopulation::mutate()
{
	/*	pos:	Is the starting position in the population
	 *		for the application of the mutation procedure.
	 *		It depends on the value of selection_rate.
	 * */
	int pos=(int)(selection_rate*genome_count);
	/*	For every genome and for every element in the genome
	 *	we select a random number in the range[0,1]. If that
	 *	number is less 
	 * */
	for(int i=pos;i<genome_count;i++)
	{
		for(int j=0;j<genome_size;j++)
		{
			double r=rand()*1.0/RAND_MAX;
			if(r<mutation_rate)
			{
                genome[i][j]=fabs(genome[i][j]*(1.0+ 0.05*(2.0*drand48()-1.0)));
		
			}
		}
	}
}

void	GPopulation::calcFitnessArray()
{
	/*	g: auxiliary array, for the calculation of fitness values.
	 * */
    vector<double> g;
	g.resize(genome_size);

	/*	For every genome calculate its fitness value and store
	 *	it to the corresponding position in the fitness array.
     * */
extern int threads;
if(threads<=1)
{
    for(int i=0;i<genome_count;i++)
    {

        fitness_array[i]=fitness(genome[i]);
    }
}
else
{
#pragma omp parallel for num_threads(threads)
	for(int i=0;i<genome_count;i++)
    {
        fitness_array[i]=fitness(genome[i]);
	}
}
}

int	GPopulation::getGeneration() const
{
	/*	Return the number of generations.
	 * */
	return generation;
}

int	GPopulation::getCount() const
{
	/*	Return the amount of individuals in the population.
	 * */
	return genome_count;
}

int	GPopulation::getSize() const
{
	/*	Return the size of each genome in the population.
	 * */
	return genome_size;
}

void	GPopulation::nextGeneration()
{
	/*	Create the next generation,
	 *	by the application of 
	 *	genetic operators.
	 * */
	if(generation) mutate();
	calcFitnessArray();
	select();
	crossover();

   extern int  localSearchGenerations;
   extern int localSearchChromosomes;
    if((generation+1)%localSearchGenerations==0)
	{
		extern int threads ;
		if(threads<=1)
		{
        for(int i=0;i<localSearchChromosomes;i++)
                localSearch(rand() % genome_count);
		}
		else
		{
#pragma omp parallel for num_threads(threads)
        for(int i=0;i<localSearchChromosomes;i++)
                localSearch(rand() % genome_count);
		}
        select();
    }
	++generation;
}

void	GPopulation::setTournamentSize(int n)
{
	/*	If n is positive and less than genome_count,
	 *	then set tournament_size to n.
	 * */
	if(n<=0 || n>=genome_count) return;
	tournament_size=n;
}

int	GPopulation::getTournamentSize() const
{
	/*	Return the value of tournament_size.
	 * */
	return tournament_size;
}

void	GPopulation::setMutationRate(double r)
{
	/*	If r is between 0 and 1, set the value of
	 *	mutation_rate to r.
	 * */
	if(r>=0 && r<=1) mutation_rate = r;
}

void	GPopulation::setSelectionRate(double r)
{
	/*	If r is between 0 and 1, set the value of
	 *	selection_rate to r.
	 * */
	if(r>=0 && r<=1) selection_rate = r;
}

double	GPopulation::getSelectionRate() const
{
	/*	Return the value of selection rate.
	 * */
	return selection_rate;
}

double	GPopulation::getMutationRate() const
{
	/*	Return the value of mutation rate.
	 * */
	return mutation_rate;
}

double	GPopulation::getBestFitness() const
{
	/*	Return the best fitness in the population, 
	 *	without evaluation. The best fitness is 
	 *	located in the first position, because we
	 *	use elitism.
	 * */
	return  fitness_array[0];
}

vector<double> GPopulation::getBestGenome() const
{
	/*	Return the best genome in the population.
	 *	The best genome is located in the first 
	 *	position of the population, because we 
	 *	use elitism.
	 * */
    return genome[0];
}

double	GPopulation::evaluateBestFitness() 
{
	/*	It calculates the value of the best genome in the population.
     * */
    return fitness(genome[0]);
}

GPopulation::~GPopulation()
{

}
