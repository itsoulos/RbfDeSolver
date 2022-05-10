#include "rbf.h"
Rbf::Rbf(int d,int w)
{
    dimension = d;
    weights = w;
}

double      Rbf::getDistance(Data &x,Data &y)
{
    double sum = 0.0;
    for(unsigned int i=0;i<x.size();i++)
        sum+=(x[i]-y[i])*(x[i]-y[i]);
    return sqrt(sum);
}

double      Rbf::gaussianDerivative(Data &x,Data &m,double v,int pos)
{
      if(fabs(v)<1e-7) return 1e+100;
    double hx = gaussian(x,m,v);
    return hx * (-2.0/v)*(x[pos]-m[pos]);
}

double      Rbf::gaussianSecondDerivative(Data &x,Data &m,double v,int pos)
{
      if(fabs(v)<1e-7) return 1e+100;
    double hx = gaussian(x,m,v);
    double phixx = (-2.0/v)*(x[pos]-m[pos]);
    return hx *(phixx  * phixx +(-2.0/v));
}

double      Rbf::gaussian(Data &x,Data &m,double v)
{
    if(fabs(v)<1e-7) return 1e+100;
    double dist =getDistance(x,m);
    return exp(-dist * dist/v);
}

int     Rbf::getNodes() const
{
    return weights;
}

int     Rbf::getDimension() const
{
    return dimension;
}

void    Rbf::setVariables(Data &w)
{
    center.resize(weights);
    for(int i=0;i<weights;i++)
        center[i].resize(dimension);
    variance.resize(weights);
    weight.resize(weights);
    for(int i=0;i<weights;i++)
    {
        for(int j=0;j<dimension;j++)
        {
            center[i][j]=w[i*dimension+j];
        }
        variance[i]=w[weights * dimension+i];
        weight[i]=w[weights *dimension+weights+i];
    }
}

double  Rbf::getValue(Data &x)
{
    double sum = 0.0;
    for(int i=0;i<weights;i++)
    {
        double px = gaussian(x,center[i],variance[i]);
        sum+=weight[i] * px;
    }
    return sum;
}

double  Rbf::getDerivative(Data &x,int pos)
{
    double sum = 0.0;
    for(int i=0;i<weights;i++)
    {
        double px = gaussianDerivative(x,center[i],variance[i],pos);
        sum = sum + weight[i]*px;
    }
    return sum;
}

double  Rbf::getSecondDerivative(Data &x,int pos)
{
    double sum = 0.0;
    for(int i=0;i<weights;i++)
    {
        double px = gaussianSecondDerivative(x,center[i],variance[i],pos);
        sum = sum + weight[i]*px;
    }
    return sum;
}


QString doubletostring(double v)
{
    if(v>=0)
        return QString::number(v,'f',6);
    else
        return "("+QString::number(v,'f',6)+")";
}
QString Rbf::toString()
{
    QString s="";
    if(dimension == 1)
    {
    for(int i=0;i<weights;i++)
    {
        s = s + doubletostring(weight[i])+"*exp(-(x-"+
                doubletostring(center[i][0])+")*(x-"+
                doubletostring(center[i][0])+")/"+
                doubletostring(variance[i])+")";
        if(i!=weights-1)
            s=s+"+";
    }
    }
    else
    if(dimension == 2)
    {
        for(int i=0;i<weights;i++)
        {
            s = s + doubletostring(weight[i])+"*exp(-((x-"+
                    doubletostring(center[i][0])+")*(x-"+
                    doubletostring(center[i][0])+")+(y-"+
                    doubletostring(center[i][1])+")*(y-"+
                    doubletostring(center[i][1])+
                    "))/"+
                    doubletostring(variance[i])+")";
            if(i!=weights-1)
                s=s+"+";
        }
    }
    return s;
}

Rbf::~Rbf()
{

}
