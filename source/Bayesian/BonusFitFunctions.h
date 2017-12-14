#ifndef BONUSFITFUNCTIONS_H
#define BONUSFITFUNCTIONS_H

/*!
 * \file
 * \brief The class for user-defined functions build on MjjFitFunctionBase
 * \author Katherine Pachal
 * \date 2013
 * 
 * This header allows the user to specify their own function class
 * for use in the package. All functions must consist of an MjjTF1Wrapper
 * with the mathematical form of the function and an MjjFitFunction 
 * which uses it to create a TF1. Many standard dijet functions are
 * already implemented here, so a user can simply imitate these 
 * examples to add any new functions they wish.
 *
 * List of functions included (here and in MjjFitFunction.h):
 
 * ThreeParam2015FitFunction :
    f = par[0] * (1-x)^par[1] * x^par[2]
    
 * FourParamFitFunction :
    f = par[0] * (1-x)^par[1] * x^( - par[2] - par[3] * log(x))
    
 * FiveParamSqrtsFitFunction :
    f = par[0] * (1-x)^par[1] * x^( - par[2] - par[3] * log(x))  but x includes a free param

 * FiveParamLog2FitFunction :
    f = par[0] * (1-x)^par[1] * x^( - par[2] - par[3] * log(x) - par[4] * (log(x))^2)

 * UA2FitFunction :
    f = (par[0]/(x^par[1])) * e^(- par[2] * x - par[3] * x^2)

 * CDFFitFunction :
 * Yes, have cross checked the location of ECM scaling against http://arxiv.org/pdf/1110.5302v1.pdf
    f = (par[0]/(x^par[1])) * (1-x/ECM)^par[2]
    
 * CDF1997FitFunction :
 * Yes, have cross checked the location of ECM scaling against http://arxiv.org/pdf/1110.5302v1.pdf
    f = (par[0]/(x^par[1])) * (1 - x/ECM - par[3] * (x/ECM)^2)^par[2]

 * D0FitFunction :
    f = (par[0]/(x^par[1])) * log(par[2]/x) * log(par[3]/x^2) where x is just mass

 * TeV Gravity Analysis:
 
 * ThreeParamFitFunction :
    f = (1-x)^par[0] * x^( - par[1] - par[2]*log(x))

 * Our own experimental :

 * SixParamFitFunction :
    f = par[0] * (1-x)^par[1] * x^( - par[2] - par[3]*logx - par[4]*(log(x))^2)

 * SixParam2FitFunction :
    f = par[0] * (1-x)^par[1] * (1-x^2)^par[4] * x^( - par[2] - par[3]*log(x))

 * From multijet analysis:

 * MultijetFitFunction1 :
    f = par[0] * (1-x)^par[1] * x^( - par[2])
    NOT IMPLEMENTED: IDENTICAL TO ThreeParam2015FitFunction

 * MultijetFitFunction2 :
    f = par[0] * (1-x)^par[1] * e^(par[2]*x^2)
    
 * MultijetFitFunction3 :
    f = par[0] * (1-x)^par[1] * x^(par[2]* x)

 * MultijetFitFunction4 :
    f = par[0] * (1-x)^par[1] * x^(par[2]*log(x))
    
 * MultijetFitFunction5 :
    f = par[0] * (1-x)^par[1] * (1+x)^(par[2]*x)

 * MultijetFitFunction6 :
    f = par[0] * (1-x)^par[1] * (1+x)^(par[2]*log(x))
    
 * MultijetFitFunction7 :
    f = (par[0]/x) * (1-x)^(par[1] - par[2]*log(x))
 
 * MultijetFitFunction8 :
    f = (par[0]/x^2) * (1-x)^(par[1] - par[2]*log(x))
    
 * MultijetFitFunction9 :
    f = par[0] * (1 - x^(1/3))^par[1] * x^( - par[2])
    
 * MultijetFitFunction10 :
    f = par[0] * (1 - x^(1/3))^par[1] * x^(par[2] * log(x))

 * From gamma gamma high mass analysis:

 * TwoParGammaGamma :
    f = (1 - x)^par[0] * x^( - par[1])

 * ThreeParGammaGamma :
    f = (1 - x)^par[0] * x^( - par[1] - par[2] * log(x))
    NOT IMPLEMENTED: IDENTICAL TO ThreeParamFitFunction

 * FourParGammaGamma :
    f = (1 - x)^par[0] * x^( - par[1] - par[2] * (log(x))^2)

 * TwoParGammaGammaWithThird :
    f = (1 - x^(1/3))^par[0] * x^( - par[1])

 * ThreeParGammaGammaWithThird :
    f = (1 - x^(1/3))^par[0] * x^( - par[1] - par[2] * log(x))

 * FourParGammaGammaWithThird :
    f = (1 - x^(1/3))^par[0] * x^( - par[1] - par[2] * (log(x))^2)
 
 */

// ---------------------------------------------------------


#include <string>
#include <vector>
#include "assert.h"

#include "TH1.h"
#include "TF1.h"
#include "math.h"
#include "Bayesian/MjjHistogram.h"
#include "Bayesian/MjjFitFunctionBase.h"
#include "Bayesian/MjjFitFunction.h"

//////////////////////////////////////////////////////////
// Three parameter dijet (no norm)

// ---------------------------------------------------------
class ThreeParamDijetTF1 : public MjjTF1Wrapper {

  public:

  ThreeParamDijetTF1(double Ecm) {
    fCenterOfMassEnergy = Ecm;
    fUseWindowExclusion = false;
    fExcludeFromFitLow = 0;
    fExcludeFromFitHigh = 0;
  }
  ~ThreeParamDijetTF1() {}

  double operator()(double *m, double *par) const {

    double x = m[0] / fCenterOfMassEnergy;

    if (fUseWindowExclusion) {
      if (m[0] >= fExcludeFromFitLow && m[0] <= fExcludeFromFitHigh) {
        TF1::RejectPoint(1);
        return 0.;
      }
    }

    // Following the convention in https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2012-088/
    double f = pow(1-x , par[0]) * 1./pow(x, par[1] + par[2] * log(x));

    return f;
  }

};

// ---------------------------------------------------------
class ThreeParamFitFunction : public MjjFitFunction {
  
public:
  
  ThreeParamFitFunction(double & mjjLow, double & mjjHigh, double Ecm=8000.) {
    fNParameters = 3;
    fFitParameters.clear();
    for (int i=0; i<fNParameters; i++) {
      MjjFitParam * thisFitParam = new MjjFitParam();
      thisFitParam->SetParamDefault(1);
      fFitParameters.push_back(thisFitParam);
    }
    //default range
    fMjjLow = mjjLow;
    fMjjHigh = mjjHigh;
    ThreeParamDijetTF1 * fThreeParam = new ThreeParamDijetTF1(Ecm);
    fPersonalisedFunction = fThreeParam;
    fFitFunction = new TF1("TheThreeParamFitFunction", fThreeParam, fMjjLow, fMjjHigh, fNParameters, "ThreeParamFitFunction");
  }

  ~ThreeParamFitFunction() {
  }
   
};  



//////////////////////////////////////////////////////////
// Six parameter dijet (ln^2 x term)

// ---------------------------------------------------------
class SixParamTF1 : public MjjTF1Wrapper {
    
public:
    
    SixParamTF1(double Ecm) {
        fCenterOfMassEnergy = Ecm;
        fUseWindowExclusion = false;
        fExcludeFromFitLow = 0;
        fExcludeFromFitHigh = 0;
    }
    ~SixParamTF1() {}
    
    double operator()(double *m, double *par) const {
        
        double x = m[0] / (par[5] > 0. ? par[5] : 1e-15);
        
        if (fUseWindowExclusion) {
            if (m[0] >= fExcludeFromFitLow && m[0] <= fExcludeFromFitHigh) {
                TF1::RejectPoint(1);
                return 0.;
            }
        }
        
        // Following the convention in https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2012-088/
        double logx = log(x);
        double f = par[0] * pow(1-x , par[1]) * 1./pow(x, par[2] + par[3] * logx + (par[4] * logx * logx));
        
        return f;
    }
    
};

// ---------------------------------------------------------
class SixParamFitFunction : public MjjFitFunction {
    
public:
    
    SixParamFitFunction(double & mjjLow, double & mjjHigh, double Ecm=8000.) {
        fNParameters = 6;
        fFitParameters.clear();
        for (int i=0; i<fNParameters; i++) {
            MjjFitParam * thisFitParam = new MjjFitParam();
            thisFitParam->SetParamDefault(1);
            fFitParameters.push_back(thisFitParam);
        }
        //default range
        fMjjLow = mjjLow;
        fMjjHigh = mjjHigh;
        SixParamTF1 * fSixParam = new SixParamTF1(Ecm);
        fPersonalisedFunction = fSixParam;
        fFitFunction = new TF1("TheSixParamFitFunction", fSixParam, fMjjLow, fMjjHigh, fNParameters, "SixParamFitFunction");
    }
    
    ~SixParamFitFunction() {
    }
    
};

//////////////////////////////////////////////////////////
// Six parameter dijet, my way

// ---------------------------------------------------------
class SixParam2TF1 : public MjjTF1Wrapper {

public:

    SixParam2TF1(double Ecm) {
        fCenterOfMassEnergy = Ecm;
        fUseWindowExclusion = false;
        fExcludeFromFitLow = 0;
        fExcludeFromFitHigh = 0;
    }
    ~SixParam2TF1() {}

    double operator()(double *m, double *par) const {

        double x = m[0] / (par[5] > 0. ? par[5] : 1e-15);
//        double x = m[0] / fCenterOfMassEnergy;


        if (fUseWindowExclusion) {
            if (m[0] >= fExcludeFromFitLow && m[0] <= fExcludeFromFitHigh) {
                TF1::RejectPoint(1);
                return 0.;
            }
        }

        // Following the convention in https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2012-088/
        double logx = log(x);
        double f = par[0] * pow(1-x , par[1]) * pow(1-x*x, par[4]) * 1./pow(x, par[2] + par[3] * logx);

        return f;
    }

};

// ---------------------------------------------------------
class SixParam2FitFunction : public MjjFitFunction {
public:

    SixParam2FitFunction(double & mjjLow, double & mjjHigh, double Ecm=8000.) {
        fNParameters = 6;
        fFitParameters.clear();
        for (int i=0; i<fNParameters; i++) {
            MjjFitParam * thisFitParam = new MjjFitParam();
            thisFitParam->SetParamDefault(1);
            fFitParameters.push_back(thisFitParam);
        }
        //default range
        fMjjLow = mjjLow;
        fMjjHigh = mjjHigh;
        SixParam2TF1 * fSixParam2 = new SixParam2TF1(Ecm);
        fPersonalisedFunction = fSixParam2;
        fFitFunction = new TF1("TheSixParam2FitFunction", fSixParam2, fMjjLow, fMjjHigh, fNParameters, "SixParam2FitFunction");
    }

    ~SixParam2FitFunction() {
    }

};


//////////////////////////////////////////////////////////
// Multijet Function 2

// ---------------------------------------------------------
class MultijetFunction2TF1 : public MjjTF1Wrapper {

public:

    MultijetFunction2TF1(double Ecm) {
        fCenterOfMassEnergy = Ecm;
        fUseWindowExclusion = false;
        fExcludeFromFitLow = 0;
        fExcludeFromFitHigh = 0;
    }
    ~MultijetFunction2TF1() {}

    double operator()(double *m, double *par) const {

        double x = m[0] / fCenterOfMassEnergy;

        if (fUseWindowExclusion) {
            if (m[0] >= fExcludeFromFitLow && m[0] <= fExcludeFromFitHigh) {
                TF1::RejectPoint(1);
                return 0.;
            }
        }

        double f = par[0] * pow(1-x , par[1]) * exp(par[2] * x * x);

        return f;
    }

};

// ---------------------------------------------------------
class MultijetFitFunction2 : public MjjFitFunction {
public:

    MultijetFitFunction2(double & mjjLow, double & mjjHigh, double Ecm=8000.) {
        fNParameters = 3;
        fFitParameters.clear();
        for (int i=0; i<fNParameters; i++) {
            MjjFitParam * thisFitParam = new MjjFitParam();
            thisFitParam->SetParamDefault(1);
            fFitParameters.push_back(thisFitParam);
        }
        //default range
        fMjjLow = mjjLow;
        fMjjHigh = mjjHigh;
        MultijetFunction2TF1 * fMultijet2 = new MultijetFunction2TF1(Ecm);
        fPersonalisedFunction = fMultijet2;
        fFitFunction = new TF1("MultijetFitFunction2", fMultijet2, fMjjLow, fMjjHigh, fNParameters, "MultijetFitFunction2");
    }

    ~MultijetFitFunction2() {
    }

};

//////////////////////////////////////////////////////////
// Multijet Function 3

// ---------------------------------------------------------
class MultijetFunction3TF1 : public MjjTF1Wrapper {

public:

    MultijetFunction3TF1(double Ecm) {
        fCenterOfMassEnergy = Ecm;
        fUseWindowExclusion = false;
        fExcludeFromFitLow = 0;
        fExcludeFromFitHigh = 0;
    }
    ~MultijetFunction3TF1() {}

    double operator()(double *m, double *par) const {

        double x = m[0] / fCenterOfMassEnergy;

        if (fUseWindowExclusion) {
            if (m[0] >= fExcludeFromFitLow && m[0] <= fExcludeFromFitHigh) {
                TF1::RejectPoint(1);
                return 0.;
            }
        }

        double f = par[0] * pow(1-x , par[1]) * pow(x, par[2] * x);

        return f;
    }

};

// ---------------------------------------------------------
class MultijetFitFunction3 : public MjjFitFunction {
public:

    MultijetFitFunction3(double & mjjLow, double & mjjHigh, double Ecm=8000.) {
        fNParameters = 3;
        fFitParameters.clear();
        for (int i=0; i<fNParameters; i++) {
            MjjFitParam * thisFitParam = new MjjFitParam();
            thisFitParam->SetParamDefault(1);
            fFitParameters.push_back(thisFitParam);
        }
        //default range
        fMjjLow = mjjLow;
        fMjjHigh = mjjHigh;
        MultijetFunction3TF1 * fMultijet3 = new MultijetFunction3TF1(Ecm);
        fPersonalisedFunction = fMultijet3;
        fFitFunction = new TF1("MultijetFitFunction3", fMultijet3, fMjjLow, fMjjHigh, fNParameters, "MultijetFitFunction3");
    }

    ~MultijetFitFunction3() {
    }

};

//////////////////////////////////////////////////////////
// Multijet Function 4

// ---------------------------------------------------------
class MultijetFunction4TF1 : public MjjTF1Wrapper {

public:

    MultijetFunction4TF1(double Ecm) {
        fCenterOfMassEnergy = Ecm;
        fUseWindowExclusion = false;
        fExcludeFromFitLow = 0;
        fExcludeFromFitHigh = 0;
    }
    ~MultijetFunction4TF1() {}

    double operator()(double *m, double *par) const {

        double x = m[0] / fCenterOfMassEnergy;

        if (fUseWindowExclusion) {
            if (m[0] >= fExcludeFromFitLow && m[0] <= fExcludeFromFitHigh) {
                TF1::RejectPoint(1);
                return 0.;
            }
        }

        double logx = log(x);
        double f = par[0] * pow(1-x , par[1]) * pow(x, par[2] * logx);

        return f;
    }

};

// ---------------------------------------------------------
class MultijetFitFunction4 : public MjjFitFunction {
public:

    MultijetFitFunction4(double & mjjLow, double & mjjHigh, double Ecm=8000.) {
        fNParameters = 3;
        fFitParameters.clear();
        for (int i=0; i<fNParameters; i++) {
            MjjFitParam * thisFitParam = new MjjFitParam();
            thisFitParam->SetParamDefault(1);
            fFitParameters.push_back(thisFitParam);
        }
        //default range
        fMjjLow = mjjLow;
        fMjjHigh = mjjHigh;
        MultijetFunction4TF1 * fMultijet4 = new MultijetFunction4TF1(Ecm);
        fPersonalisedFunction = fMultijet4;
        fFitFunction = new TF1("MultijetFitFunction4", fMultijet4, fMjjLow, fMjjHigh, fNParameters, "MultijetFitFunction4");
    }

    ~MultijetFitFunction4() {
    }

};

//////////////////////////////////////////////////////////
// Multijet Function 5

// ---------------------------------------------------------
class MultijetFunction5TF1 : public MjjTF1Wrapper {

public:

    MultijetFunction5TF1(double Ecm) {
        fCenterOfMassEnergy = Ecm;
        fUseWindowExclusion = false;
        fExcludeFromFitLow = 0;
        fExcludeFromFitHigh = 0;
    }
    ~MultijetFunction5TF1() {}

    double operator()(double *m, double *par) const {

        double x = m[0] / fCenterOfMassEnergy;

        if (fUseWindowExclusion) {
            if (m[0] >= fExcludeFromFitLow && m[0] <= fExcludeFromFitHigh) {
                TF1::RejectPoint(1);
                return 0.;
            }
        }

        double f = par[0] * pow(1-x , par[1]) * pow(1 + x, par[2] * x);

        return f;
    }

};

// ---------------------------------------------------------
class MultijetFitFunction5 : public MjjFitFunction {
public:

    MultijetFitFunction5(double & mjjLow, double & mjjHigh, double Ecm=8000.) {
        fNParameters = 3;
        fFitParameters.clear();
        for (int i=0; i<fNParameters; i++) {
            MjjFitParam * thisFitParam = new MjjFitParam();
            thisFitParam->SetParamDefault(1);
            fFitParameters.push_back(thisFitParam);
        }
        //default range
        fMjjLow = mjjLow;
        fMjjHigh = mjjHigh;
        MultijetFunction5TF1 * fMultijet5 = new MultijetFunction5TF1(Ecm);
        fPersonalisedFunction = fMultijet5;
        fFitFunction = new TF1("MultijetFitFunction5", fMultijet5, fMjjLow, fMjjHigh, fNParameters, "MultijetFitFunction5");
    }

    ~MultijetFitFunction5() {
    }

};

//////////////////////////////////////////////////////////
// Multijet Function 6

// ---------------------------------------------------------
class MultijetFunction6TF1 : public MjjTF1Wrapper {

public:

    MultijetFunction6TF1(double Ecm) {
        fCenterOfMassEnergy = Ecm;
        fUseWindowExclusion = false;
        fExcludeFromFitLow = 0;
        fExcludeFromFitHigh = 0;
    }
    ~MultijetFunction6TF1() {}

    double operator()(double *m, double *par) const {

        double x = m[0] / fCenterOfMassEnergy;

        if (fUseWindowExclusion) {
            if (m[0] >= fExcludeFromFitLow && m[0] <= fExcludeFromFitHigh) {
                TF1::RejectPoint(1);
                return 0.;
            }
        }

        double logx = log(x);
        double f = par[0] * pow(1-x , par[1]) * pow(1 + x, par[2] * logx);

        return f;
    }

};

// ---------------------------------------------------------
class MultijetFitFunction6 : public MjjFitFunction {
public:

    MultijetFitFunction6(double & mjjLow, double & mjjHigh, double Ecm=8000.) {
        fNParameters = 3;
        fFitParameters.clear();
        for (int i=0; i<fNParameters; i++) {
            MjjFitParam * thisFitParam = new MjjFitParam();
            thisFitParam->SetParamDefault(1);
            fFitParameters.push_back(thisFitParam);
        }
        //default range
        fMjjLow = mjjLow;
        fMjjHigh = mjjHigh;
        MultijetFunction6TF1 * fMultijet6 = new MultijetFunction6TF1(Ecm);
        fPersonalisedFunction = fMultijet6;
        fFitFunction = new TF1("MultijetFitFunction6", fMultijet6, fMjjLow, fMjjHigh, fNParameters, "MultijetFitFunction6");
    }

    ~MultijetFitFunction6() {
    }

};

//////////////////////////////////////////////////////////
// Multijet Function 7

// ---------------------------------------------------------
class MultijetFunction7TF1 : public MjjTF1Wrapper {

public:

    MultijetFunction7TF1(double Ecm) {
        fCenterOfMassEnergy = Ecm;
        fUseWindowExclusion = false;
        fExcludeFromFitLow = 0;
        fExcludeFromFitHigh = 0;
    }
    ~MultijetFunction7TF1() {}

    double operator()(double *m, double *par) const {

        double x = m[0] / fCenterOfMassEnergy;

        if (fUseWindowExclusion) {
            if (m[0] >= fExcludeFromFitLow && m[0] <= fExcludeFromFitHigh) {
                TF1::RejectPoint(1);
                return 0.;
            }
        }

        double logx = log(x);
        double f = (par[0]/x) * pow(1-x , par[1] - par[2] * logx);

        return f;
    }

};

// ---------------------------------------------------------
class MultijetFitFunction7 : public MjjFitFunction {
public:

    MultijetFitFunction7(double & mjjLow, double & mjjHigh, double Ecm=8000.) {
        fNParameters = 3;
        fFitParameters.clear();
        for (int i=0; i<fNParameters; i++) {
            MjjFitParam * thisFitParam = new MjjFitParam();
            thisFitParam->SetParamDefault(1);
            fFitParameters.push_back(thisFitParam);
        }
        //default range
        fMjjLow = mjjLow;
        fMjjHigh = mjjHigh;
        MultijetFunction7TF1 * fMultijet7 = new MultijetFunction7TF1(Ecm);
        fPersonalisedFunction = fMultijet7;
        fFitFunction = new TF1("MultijetFitFunction7", fMultijet7, fMjjLow, fMjjHigh, fNParameters, "MultijetFitFunction7");
    }

    ~MultijetFitFunction7() {
    }

};

//////////////////////////////////////////////////////////
// Multijet Function 8

// ---------------------------------------------------------
class MultijetFunction8TF1 : public MjjTF1Wrapper {

public:

    MultijetFunction8TF1(double Ecm) {
        fCenterOfMassEnergy = Ecm;
        fUseWindowExclusion = false;
        fExcludeFromFitLow = 0;
        fExcludeFromFitHigh = 0;
    }
    ~MultijetFunction8TF1() {}

    double operator()(double *m, double *par) const {

        double x = m[0] / fCenterOfMassEnergy;

        if (fUseWindowExclusion) {
            if (m[0] >= fExcludeFromFitLow && m[0] <= fExcludeFromFitHigh) {
                TF1::RejectPoint(1);
                return 0.;
            }
        }

        double logx = log(x);
        double f = (par[0]/(x*x)) * pow(1-x , par[1] - par[2] * logx);

        return f;
    }

};

// ---------------------------------------------------------
class MultijetFitFunction8 : public MjjFitFunction {
public:

    MultijetFitFunction8(double & mjjLow, double & mjjHigh, double Ecm=8000.) {
        fNParameters = 3;
        fFitParameters.clear();
        for (int i=0; i<fNParameters; i++) {
            MjjFitParam * thisFitParam = new MjjFitParam();
            thisFitParam->SetParamDefault(1);
            fFitParameters.push_back(thisFitParam);
        }
        //default range
        fMjjLow = mjjLow;
        fMjjHigh = mjjHigh;
        MultijetFunction8TF1 * fMultijet8 = new MultijetFunction8TF1(Ecm);
        fPersonalisedFunction = fMultijet8;
        fFitFunction = new TF1("MultijetFitFunction8", fMultijet8, fMjjLow, fMjjHigh, fNParameters, "MultijetFitFunction8");
    }

    ~MultijetFitFunction8() {
    }

};

//////////////////////////////////////////////////////////
// Multijet Function 9

// ---------------------------------------------------------
class MultijetFunction9TF1 : public MjjTF1Wrapper {

public:

    MultijetFunction9TF1(double Ecm) {
        fCenterOfMassEnergy = Ecm;
        fUseWindowExclusion = false;
        fExcludeFromFitLow = 0;
        fExcludeFromFitHigh = 0;
    }
    ~MultijetFunction9TF1() {}

    double operator()(double *m, double *par) const {

        double x = m[0] / fCenterOfMassEnergy;

        if (fUseWindowExclusion) {
            if (m[0] >= fExcludeFromFitLow && m[0] <= fExcludeFromFitHigh) {
                TF1::RejectPoint(1);
                return 0.;
            }
        }

        double val = pow(x,1.0/3.0);
        double f = par[0] * pow((1 - val), par[1]) * pow(x, -1 * par[2]);

        return f;
    }

};

// ---------------------------------------------------------
class MultijetFitFunction9 : public MjjFitFunction {
public:

    MultijetFitFunction9(double & mjjLow, double & mjjHigh, double Ecm=8000.) {
        fNParameters = 3;
        fFitParameters.clear();
        for (int i=0; i<fNParameters; i++) {
            MjjFitParam * thisFitParam = new MjjFitParam();
            thisFitParam->SetParamDefault(1);
            fFitParameters.push_back(thisFitParam);
        }
        //default range
        fMjjLow = mjjLow;
        fMjjHigh = mjjHigh;
        MultijetFunction9TF1 * fMultijet9 = new MultijetFunction9TF1(Ecm);
        fPersonalisedFunction = fMultijet9;
        fFitFunction = new TF1("MultijetFitFunction9", fMultijet9, fMjjLow, fMjjHigh, fNParameters, "MultijetFitFunction9");
    }

    ~MultijetFitFunction9() {
    }

};

//////////////////////////////////////////////////////////
// Multijet Function 10

// ---------------------------------------------------------
class MultijetFunction10TF1 : public MjjTF1Wrapper {

public:

    MultijetFunction10TF1(double Ecm) {
        fCenterOfMassEnergy = Ecm;
        fUseWindowExclusion = false;
        fExcludeFromFitLow = 0;
        fExcludeFromFitHigh = 0;
    }
    ~MultijetFunction10TF1() {}

    double operator()(double *m, double *par) const {

        double x = m[0] / fCenterOfMassEnergy;

        if (fUseWindowExclusion) {
            if (m[0] >= fExcludeFromFitLow && m[0] <= fExcludeFromFitHigh) {
                TF1::RejectPoint(1);
                return 0.;
            }
        }
      
        double logx = log(x);
        double f = par[0] * pow(1 - pow(x , 1.0/3.0), par[1]) * pow(x, par[2] * logx);

        return f;
    }

};

// ---------------------------------------------------------
class MultijetFitFunction10 : public MjjFitFunction {
public:

    MultijetFitFunction10(double & mjjLow, double & mjjHigh, double Ecm=8000.) {
        fNParameters = 3;
        fFitParameters.clear();
        for (int i=0; i<fNParameters; i++) {
            MjjFitParam * thisFitParam = new MjjFitParam();
            thisFitParam->SetParamDefault(1);
            fFitParameters.push_back(thisFitParam);
        }
        //default range
        fMjjLow = mjjLow;
        fMjjHigh = mjjHigh;
        MultijetFunction10TF1 * fMultijet10 = new MultijetFunction10TF1(Ecm);
        fPersonalisedFunction = fMultijet10;
        fFitFunction = new TF1("MultijetFitFunction10", fMultijet10, fMjjLow, fMjjHigh, fNParameters, "MultijetFitFunction10");
    }

    ~MultijetFitFunction10() {
    }

};

//////////////////////////////////////////////////////////
// TwoParGammaGamma

// ---------------------------------------------------------
class TwoParGammaGammaTF1 : public MjjTF1Wrapper {

  public:

  TwoParGammaGammaTF1(double Ecm) {
    fCenterOfMassEnergy = Ecm;
    fUseWindowExclusion = false;
    fExcludeFromFitLow = 0;
    fExcludeFromFitHigh = 0;
  }
  ~TwoParGammaGammaTF1() {}

  double operator()(double *m, double *par) const {

    double x = m[0] / fCenterOfMassEnergy;

    if (fUseWindowExclusion) {
      if (m[0] >= fExcludeFromFitLow && m[0] <= fExcludeFromFitHigh) {
        TF1::RejectPoint(1);
        return 0.;
      }
    }

    // Following the convention in https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2012-088/
    double f = pow(1 - x, par[0]) * 1./pow(x, par[1]);

    return f;
  }

};

// ---------------------------------------------------------
class TwoParGammaGamma : public MjjFitFunction {
  
public:
  
    TwoParGammaGamma(double & mjjLow, double & mjjHigh, double Ecm=13000.) {
    fNParameters = 2;
    fFitParameters.clear();
    for (int i=0; i<fNParameters; i++) {
      MjjFitParam * thisFitParam = new MjjFitParam();
      thisFitParam->SetParamDefault(1);
      fFitParameters.push_back(thisFitParam);
    }
    //default range
    fMjjLow = mjjLow;
    fMjjHigh = mjjHigh;
    TwoParGammaGammaTF1 * fFourParam = new TwoParGammaGammaTF1(Ecm);
    fPersonalisedFunction = fFourParam;
    fFitFunction = new TF1("TheFourParamFitFunction", fFourParam, fMjjLow, fMjjHigh, fNParameters, "FourParamFitFunction");
  }

  TwoParGammaGamma() {
  }
   
};


//////////////////////////////////////////////////////////
// FourParGammaGamma

// ---------------------------------------------------------
class FourParGammaGammaTF1 : public MjjTF1Wrapper {

  public:

  FourParGammaGammaTF1(double Ecm) {
    fCenterOfMassEnergy = Ecm;
    fUseWindowExclusion = false;
    fExcludeFromFitLow = 0;
    fExcludeFromFitHigh = 0;
  }
  ~FourParGammaGammaTF1() {}

  double operator()(double *m, double *par) const {

    double x = m[0] / fCenterOfMassEnergy;

    if (fUseWindowExclusion) {
      if (m[0] >= fExcludeFromFitLow && m[0] <= fExcludeFromFitHigh) {
        TF1::RejectPoint(1);
        return 0.;
      }
    }

    // Following the convention in https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2012-088/
    double f = pow(1 - x, par[0]) * 1./pow(x, par[1] + par[2] * log(x) * log(x));

    return f;
  }

};

// ---------------------------------------------------------
class FourParGammaGamma : public MjjFitFunction {
  
public:
  
    FourParGammaGamma(double & mjjLow, double & mjjHigh, double Ecm=13000.) {
    fNParameters = 4;
    fFitParameters.clear();
    for (int i=0; i<fNParameters; i++) {
      MjjFitParam * thisFitParam = new MjjFitParam();
      thisFitParam->SetParamDefault(1);
      fFitParameters.push_back(thisFitParam);
    }
    //default range
    fMjjLow = mjjLow;
    fMjjHigh = mjjHigh;
    FourParGammaGammaTF1 * fFourParam = new FourParGammaGammaTF1(Ecm);
    fPersonalisedFunction = fFourParam;
    fFitFunction = new TF1("TheFourParamFitFunction", fFourParam, fMjjLow, fMjjHigh, fNParameters, "FourParamFitFunction");
  }

  FourParGammaGamma() {
  }
   
};

//////////////////////////////////////////////////////////
// TwoParGammaGammaWithThird

// ---------------------------------------------------------
class TwoParGammaGammaWithThirdTF1 : public MjjTF1Wrapper {

  public:

  TwoParGammaGammaWithThirdTF1(double Ecm) {
    fCenterOfMassEnergy = Ecm;
    fUseWindowExclusion = false;
    fExcludeFromFitLow = 0;
    fExcludeFromFitHigh = 0;
  }
  ~TwoParGammaGammaWithThirdTF1() {}

  double operator()(double *m, double *par) const {

    double x = m[0] / fCenterOfMassEnergy;

    if (fUseWindowExclusion) {
      if (m[0] >= fExcludeFromFitLow && m[0] <= fExcludeFromFitHigh) {
        TF1::RejectPoint(1);
        return 0.;
      }
    }

    // Following the convention in https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2012-088/
    double f = pow(1 - pow(x,1.0/3.0), par[0]) * 1./pow(x, par[1]);

    return f;
  }

};

// ---------------------------------------------------------
class TwoParGammaGammaWithThird : public MjjFitFunction {
  
public:
  
    TwoParGammaGammaWithThird(double & mjjLow, double & mjjHigh, double Ecm=13000.) {
    fNParameters = 2;
    fFitParameters.clear();
    for (int i=0; i<fNParameters; i++) {
      MjjFitParam * thisFitParam = new MjjFitParam();
      thisFitParam->SetParamDefault(1);
      fFitParameters.push_back(thisFitParam);
    }
    //default range
    fMjjLow = mjjLow;
    fMjjHigh = mjjHigh;
    TwoParGammaGammaWithThirdTF1 * fFourParam = new TwoParGammaGammaWithThirdTF1(Ecm);
    fPersonalisedFunction = fFourParam;
    fFitFunction = new TF1("TheFourParamFitFunction", fFourParam, fMjjLow, fMjjHigh, fNParameters, "FourParamFitFunction");
  }

  TwoParGammaGammaWithThird() {
  }
   
};

//////////////////////////////////////////////////////////
// ThreeParGammaGammaWithThird

// ---------------------------------------------------------
class ThreeParGammaGammaWithThirdTF1 : public MjjTF1Wrapper {

  public:

  ThreeParGammaGammaWithThirdTF1(double Ecm) {
    fCenterOfMassEnergy = Ecm;
    fUseWindowExclusion = false;
    fExcludeFromFitLow = 0;
    fExcludeFromFitHigh = 0;
  }
  ~ThreeParGammaGammaWithThirdTF1() {}

  double operator()(double *m, double *par) const {

    double x = m[0] / fCenterOfMassEnergy;

    if (fUseWindowExclusion) {
      if (m[0] >= fExcludeFromFitLow && m[0] <= fExcludeFromFitHigh) {
        TF1::RejectPoint(1);
        return 0.;
      }
    }

    // Following the convention in https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2012-088/
    double f = pow(1 - pow(x,1.0/3.0), par[0]) * 1./pow(x, par[1] + par[2] * log(x));

    return f;
  }

};

// ---------------------------------------------------------
class ThreeParGammaGammaWithThird : public MjjFitFunction {
  
public:
  
    ThreeParGammaGammaWithThird(double & mjjLow, double & mjjHigh, double Ecm=13000.) {
    fNParameters = 3;
    fFitParameters.clear();
    for (int i=0; i<fNParameters; i++) {
      MjjFitParam * thisFitParam = new MjjFitParam();
      thisFitParam->SetParamDefault(1);
      fFitParameters.push_back(thisFitParam);
    }
    //default range
    fMjjLow = mjjLow;
    fMjjHigh = mjjHigh;
    ThreeParGammaGammaWithThirdTF1 * fFourParam = new ThreeParGammaGammaWithThirdTF1(Ecm);
    fPersonalisedFunction = fFourParam;
    fFitFunction = new TF1("TheFourParamFitFunction", fFourParam, fMjjLow, fMjjHigh, fNParameters, "FourParamFitFunction");
  }

  ThreeParGammaGammaWithThird() {
  }
   
};  

//////////////////////////////////////////////////////////
// FourParGammaGammaWithThird

// ---------------------------------------------------------
class FourParGammaGammaWithThirdTF1 : public MjjTF1Wrapper {

  public:

  FourParGammaGammaWithThirdTF1(double Ecm) {
    fCenterOfMassEnergy = Ecm;
    fUseWindowExclusion = false;
    fExcludeFromFitLow = 0;
    fExcludeFromFitHigh = 0;
  }
  ~FourParGammaGammaWithThirdTF1() {}

  double operator()(double *m, double *par) const {

    double x = m[0] / fCenterOfMassEnergy;

    if (fUseWindowExclusion) {
      if (m[0] >= fExcludeFromFitLow && m[0] <= fExcludeFromFitHigh) {
        TF1::RejectPoint(1);
        return 0.;
      }
    }

    // Following the convention in https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2012-088/
    double f = pow(1 - pow(x,1.0/3.0), par[0]) * 1./pow(x, par[1] + par[2] * log(x) * log(x));

    return f;
  }

};

// ---------------------------------------------------------
class FourParGammaGammaWithThird : public MjjFitFunction {
  
public:
  
    FourParGammaGammaWithThird(double & mjjLow, double & mjjHigh, double Ecm=13000.) {
    fNParameters = 4;
    fFitParameters.clear();
    for (int i=0; i<fNParameters; i++) {
      MjjFitParam * thisFitParam = new MjjFitParam();
      thisFitParam->SetParamDefault(1);
      fFitParameters.push_back(thisFitParam);
    }
    //default range
    fMjjLow = mjjLow;
    fMjjHigh = mjjHigh;
    FourParGammaGammaWithThirdTF1 * fFourParam = new FourParGammaGammaWithThirdTF1(Ecm);
    fPersonalisedFunction = fFourParam;
    fFitFunction = new TF1("TheFourParamFitFunction", fFourParam, fMjjLow, fMjjHigh, fNParameters, "FourParamFitFunction");
  }

  FourParGammaGammaWithThird() {
  }
   
};  


// ---------------------------------------------------------

#endif
