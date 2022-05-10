#include <QCoreApplication>
# include <get_options.h>
# include <godeprogram.h>
# include <gpopulation.h>
# include <gsodeprogram.h>
# include <gpdeprogram.h>
Rbf *rbf = NULL;
GOdeProgram *ode_program  = NULL;
GPopulation *ode_population = NULL;
GPopulation *sode_population = NULL;
GSodeProgram *sode_program  =NULL;
GPdeProgram  *pde_program = NULL;
GPopulation *pde_population = NULL;
vector<GProgram*> ompODeProgram;
vector<Rbf*> ompRbf;
int chromosome_size = 0;
int dimension = 1;
void init_ode()
{
    rbf = new Rbf(dimension,weights);
    ode_program =new GOdeProgram(rbf,dll_name);
    chromosome_size  = dimension * weights + weights + weights;
    ompODeProgram.resize(threads);
    ompRbf.resize(threads);
    for(int i=0;i<threads;i++)
    {
        ompRbf[i] = new Rbf(dimension,weights);
        ompODeProgram[i]=new GOdeProgram(ompRbf[i],dll_name);
    }
    if(threads<=1)
        ode_population=new GPopulation(chromosome_count,chromosome_size,ode_program);
    else
    ode_population=new GPopulation(chromosome_count,chromosome_size,ompODeProgram);
    ode_population->setSelectionRate(selection_rate);
    ode_population->setMutationRate(mutation_rate);
}

void run_ode()
{
    int iters=0;
    vector<double> genome;
    genome.resize(chromosome_size);
    double fitness;
    for(iters=1;iters<=maxgenerations;iters++)
    {
        ode_population->nextGeneration();
        ode_population->evaluateBestFitness();
        genome=ode_population->getBestGenome();
        fitness=ode_population->getBestFitness();
        rbf->setVariables(genome);
        printf("GENERATION: %4d\t FITNESS: %.10lg RBF: %s \t \n",
                iters,fabs(fitness),rbf->toString().toStdString().c_str());
        if(fabs(fitness)<eps) break;
      }
}

void done_ode()
{
    for(int i=0;i<threads;i++)
    {
        delete ompODeProgram[i];
        delete ompRbf[i];
    }

    delete ode_program;
    delete ode_population;
    delete rbf;
}


void init_sode()
{
    rbf = new Rbf(dimension,weights);
    sode_program =new GSodeProgram(rbf,dll_name);
    ompODeProgram.resize(threads);

    chromosome_size  =sode_program->getNode()*( dimension * weights + weights + weights);
    ompODeProgram.resize(threads);
    ompRbf.resize(threads);
    for(int i=0;i<threads;i++)
    {
        ompRbf[i] = new Rbf(dimension,weights);
        ompODeProgram[i]=new GOdeProgram(ompRbf[i],dll_name);
    }

    if(threads<=1)
        sode_population=new GPopulation(chromosome_count,chromosome_size,sode_program);
    else
    sode_population=new GPopulation(chromosome_count,chromosome_size,ompODeProgram);
    sode_population->setSelectionRate(selection_rate);
    sode_population->setMutationRate(mutation_rate);

}

void run_sode()
{

    int iters=0,i,j;
    vector<double> genome;
    genome.resize(chromosome_size);


    double fitness;
    Data subgenome;
    subgenome.resize(genome.size()/sode_program->getNode());
    for(iters=1;iters<=maxgenerations;iters++)
    {
            sode_population->nextGeneration();
            genome=sode_population->getBestGenome();
            fitness=sode_population->getBestFitness();
            printf("GENERATION: %4d\t FITNESS: %.10lg \n",iters,fabs(fitness));
            for(int j=0;j<sode_program->getNode();j++)
            {
                for(int k=0;k<subgenome.size();k++)
                {
                    subgenome[k]=genome[subgenome.size()*j+k];
                }
                rbf->setVariables(subgenome);
                printf("Rbf[%d]=%s\n",j,rbf->toString().toStdString().c_str());
            }
            if(fabs(fitness)<eps) break;
    }
}

void done_sode()
{
    for(int i=0;i<threads;i++)
    {
        delete ompODeProgram[i];
        delete ompRbf[i];
    }
    delete sode_population;
    delete sode_program;
}


void init_pde()
{
    dimension = 2;
    rbf = new Rbf(dimension,weights);
    pde_program =new GPdeProgram(rbf,dll_name);

    chromosome_size  =( dimension * weights + weights + weights);
    ompODeProgram.resize(threads);
    ompRbf.resize(threads);
    for(int i=0;i<threads;i++)
    {
        ompRbf[i] = new Rbf(dimension,weights);
        ompODeProgram[i]=new GOdeProgram(ompRbf[i],dll_name);
    }

    if(threads<=1)
        pde_population=new GPopulation(chromosome_count,chromosome_size,pde_program);
    else
    pde_population=new GPopulation(chromosome_count,chromosome_size,ompODeProgram);
    pde_population->setSelectionRate(selection_rate);
    pde_population->setMutationRate(mutation_rate);
}

void run_pde()
{
    int iters=0;
    vector<double> genome;
    genome.resize(chromosome_size);
    string solution;
    double fitness;
    for(iters=1;iters<=maxgenerations;iters++)
    {
        pde_population->nextGeneration();
        genome=pde_population->getBestGenome();
        fitness=pde_population->getBestFitness();
        rbf->setVariables(genome);
        printf("GENERATION: %4d\t FITNESS: %.10lg RBF: %s \t \n",
                iters,fabs(fitness),rbf->toString().toStdString().c_str());
        if(fabs(fitness)<eps) break;
     }
}

void done_pde()
{
    for(int i=0;i<threads;i++)
    {
        delete ompODeProgram[i];
        delete ompRbf[i];
    }
    delete pde_population;
    delete pde_program;
}

int main(int argc, char *argv[])
{
    QCoreApplication app(argc,argv);
    setlocale(LC_ALL,"C");
    parseCmdLine(app.arguments());
    srand(genome_rand);
    if(dll_name=="")
    {
        print_usage();
        return 0;
    }
    if(ode_kind=="ode")
    {
        init_ode();
        run_ode();
        done_ode();
    }
    else
    if(ode_kind=="pde")
    {
      init_pde();
      run_pde();
      done_pde();
    }
    else
    if(ode_kind=="sode")
    {
      init_sode();
        run_sode();
        done_sode();
    }
    return 0;

}
