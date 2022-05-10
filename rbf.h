#ifndef RBF_H
#define RBF_H
# include <math.h>
# include <stdlib.h>
# include <string.h>
# include <stdio.h>
# include <QString>
# include <vector>
using namespace std;
typedef vector<double> Data;

class Rbf
{
private:
    int dimension;
    int weights;
    vector<Data>    center;
    Data            variance;
    Data            weight;
    double      getDistance(Data &x,Data &y);
    double      gaussian(Data &x,Data &m,double v);
    double      gaussianDerivative(Data &x,Data &m,double v,int pos);
    double      gaussianSecondDerivative(Data &x,Data &m,double v,int pos);

public:
    Rbf(int d,int w);
    int     getNodes() const;
    int     getDimension() const;
    void    setVariables(Data &w);
    double  getValue(Data &x);
    double  getDerivative(Data &x,int pos);
    double  getSecondDerivative(Data &x,int pos);
    QString toString();
    ~Rbf();
};

#endif // RBF_H
