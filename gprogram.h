#ifndef GPROGRAM_H
#define GPROGRAM_H
# include <rbf.h>
# include <vector>
using namespace std;
typedef vector<double> Data;

class GProgram
{
protected:
    Rbf *rbf;
public:
    GProgram(Rbf *m);
    virtual double fitness(Data &x)=0;
};

#endif // GPROGRAM_H
