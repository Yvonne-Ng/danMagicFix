#ifndef MJJFITFUNCTION_H
#define MJJFITFUNCTION_H

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


//////////////////////////////////////////////////////////
// 3 parameter dijet Function (2015 version)

// ---------------------------------------------------------
class ThreeParam2015TF1 : public MjjTF1Wrapper {

  public:

  ThreeParam2015TF1(double Ecm) {
    fCenterOfMassEnergy = Ecm;
    fUseWindowExclusion = false;
    fExcludeFromFitLow = 0;
    fExcludeFromFitHigh = 0;
  }
  ~ThreeParam2015TF1() {}

  double operator()(double *m, double *par) const {

    double x = m[0];

    if (fUseWindowExclusion) {
      if (m[0] >= fExcludeFromFitLow && m[0] <= fExcludeFromFitHigh) {
        TF1::RejectPoint(1);
        return 0.;
      }
    }

    double f = par[0]*pow((1-x/fCenterOfMassEnergy),par[1])*pow((x/fCenterOfMassEnergy),par[2]);
    return f;
  }

};

// ---------------------------------------------------------
class ThreeParam2015FitFunction : public MjjFitFunction {

public:

  ThreeParam2015FitFunction(double & mjjLow, double & mjjHigh, double Ecm=13000.) {
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
    ThreeParam2015TF1 * fThreeparam2015 = new ThreeParam2015TF1(Ecm);
    fPersonalisedFunction = fThreeparam2015;
    fFitFunction = new TF1("ThreeParam2015FitFunction", fThreeparam2015, fMjjLow, fMjjHigh, fNParameters, "ThreeParam2015FitFunction");
  }

  ~ThreeParam2015FitFunction() {
  }

};


//////////////////////////////////////////////////////////
// Four parameter dijet

// ---------------------------------------------------------
class FourParamTF1 : public MjjTF1Wrapper {

  public:

  FourParamTF1(double Ecm) {
    fCenterOfMassEnergy = Ecm;
    fUseWindowExclusion = false;
    fExcludeFromFitLow = 0;
    fExcludeFromFitHigh = 0;
  }
  ~FourParamTF1() {}

  double operator()(double *m, double *par) const {

    double x = m[0] / fCenterOfMassEnergy;

    if (fUseWindowExclusion) {
      if (m[0] >= fExcludeFromFitLow && m[0] <= fExcludeFromFitHigh) {
        TF1::RejectPoint(1);
        return 0.;
      }
    }

    // Following the convention in https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2012-088/
    double f = par[0] * pow(1-x , par[1]) * 1./pow(x, par[2] + par[3] * log(x));

    return f;
  }

};

// ---------------------------------------------------------
class FourParamFitFunction : public MjjFitFunction {
  
public:
  
    FourParamFitFunction(double & mjjLow, double & mjjHigh, double Ecm=13000.) {
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
    FourParamTF1 * fFourParam = new FourParamTF1(Ecm);
    fPersonalisedFunction = fFourParam;
    fFitFunction = new TF1("TheFourParamFitFunction", fFourParam, fMjjLow, fMjjHigh, fNParameters, "FourParamFitFunction");
  }

  ~FourParamFitFunction() {
  }
   
};  


//////////////////////////////////////////////////////////
// Five parameter dijet (sqrt(s)->free parameter)

// ---------------------------------------------------------
class FiveParamSqrtsTF1 : public MjjTF1Wrapper {
    
public:
    
    FiveParamSqrtsTF1() {
        fCenterOfMassEnergy = 13000; // Caution this is hardcoded to 13000 GeV
        fUseWindowExclusion = false;
        fExcludeFromFitLow = 0;
        fExcludeFromFitHigh = 0;
    }
    ~FiveParamSqrtsTF1() {}
    
    double operator()(double *m, double *par) const {
        
        double x = m[0] / (par[4] > 0. ? par[4] : 1e-15);
        
        if (fUseWindowExclusion) {
            if (m[0] >= fExcludeFromFitLow && m[0] <= fExcludeFromFitHigh) {
                TF1::RejectPoint(1);
                return 0.;
            }
        }
        
        // Following the convention in https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2012-088/
        double f = par[0] * pow(1-x , par[1]) * 1./pow(x, par[2] + par[3] * log(x));
        
        return f;
    }
    
};

// ---------------------------------------------------------
class FiveParamSqrtsFitFunction : public MjjFitFunction {
    
public:

    FiveParamSqrtsFitFunction(double & mjjLow, double & mjjHigh, double Ecm=13000.) {
        fNParameters = 5;
        fFitParameters.clear();
        for (int i=0; i<fNParameters; i++) {
            MjjFitParam * thisFitParam = new MjjFitParam();
            thisFitParam->SetParamDefault(1);
            fFitParameters.push_back(thisFitParam);
        }
        //default range
        fMjjLow = mjjLow;
        fMjjHigh = mjjHigh;
        FiveParamSqrtsTF1 * fFiveParam = new FiveParamSqrtsTF1();
        fPersonalisedFunction = fFiveParam;
        fFitFunction = new TF1("TheFiveParamSqrtsFitFunction", fFiveParam, fMjjLow, fMjjHigh, fNParameters, "FiveParamSqrtsFitFunction");
    }
    
    ~FiveParamSqrtsFitFunction() {
    }
    
};  

//////////////////////////////////////////////////////////
// Five parameter dijet (ln^2 x term)

// ---------------------------------------------------------
class FiveParamLog2TF1 : public MjjTF1Wrapper {
    
public:
    
    FiveParamLog2TF1(double Ecm) {
        fCenterOfMassEnergy = Ecm;
        fUseWindowExclusion = false;
        fExcludeFromFitLow = 0;
        fExcludeFromFitHigh = 0;
    }
    ~FiveParamLog2TF1() {}
    
    double operator()(double *m, double *par) const {
        
        double x = m[0] / fCenterOfMassEnergy;
        
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
class FiveParamLog2FitFunction : public MjjFitFunction {
    
public:
    FiveParamLog2FitFunction(double & mjjLow, double & mjjHigh, double Ecm=13000.) {
        fNParameters = 5;
        fFitParameters.clear();
        for (int i=0; i<fNParameters; i++) {
            MjjFitParam * thisFitParam = new MjjFitParam();
            thisFitParam->SetParamDefault(1);
            fFitParameters.push_back(thisFitParam);
        }
        //default range
        fMjjLow = mjjLow;
        fMjjHigh = mjjHigh;
        FiveParamLog2TF1 * fFiveParam = new FiveParamLog2TF1(Ecm);
        fPersonalisedFunction = fFiveParam;
        fFitFunction = new TF1("TheFiveParamLog2FitFunction", fFiveParam, fMjjLow, fMjjHigh, fNParameters, "FiveParamLog2FitFunction");
    }
    
    ~FiveParamLog2FitFunction() {
    }
    
};

//////////////////////////////////////////////////////////
// UA2 Function

// ---------------------------------------------------------
class UA2TF1 : public MjjTF1Wrapper {

  public:

  UA2TF1(double Ecm) {
    fCenterOfMassEnergy = Ecm;
    fUseWindowExclusion = false;
    fExcludeFromFitLow = 0;
    fExcludeFromFitHigh = 0;
  }

  ~UA2TF1() {}

  double operator()(double *m, double *par) const {

    double x = m[0]/fCenterOfMassEnergy ;

    if (fUseWindowExclusion) {
      if (m[0] >= fExcludeFromFitLow && m[0] <= fExcludeFromFitHigh) {
        TF1::RejectPoint(1);
        return 0.;
      }
    }

    double f = par[0] * pow(x,-1*par[1]) * exp(-1*(par[2]*x + par[3]*pow(x,2)));
    return f;
  }

};

// ---------------------------------------------------------
class UA2FitFunction : public MjjFitFunction {

public:

  UA2FitFunction(double & mjjLow, double & mjjHigh, double Ecm=13000.) {
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
    UA2TF1 * fUa2 = new UA2TF1(Ecm);
    fPersonalisedFunction = fUa2;
    fFitFunction = new TF1("TheFourParamFitFunction", fUa2, fMjjLow, fMjjHigh, fNParameters, "FourParamFitFunction");
  }

  ~UA2FitFunction() {
  }

};

//////////////////////////////////////////////////////////
// CDF Function (1995 version)

// ---------------------------------------------------------
class CDFTF1 : public MjjTF1Wrapper {

  public:

  CDFTF1(double Ecm) {
    fCenterOfMassEnergy = Ecm;
    fUseWindowExclusion = false;
    fExcludeFromFitLow = 0;
    fExcludeFromFitHigh = 0;
  }
  ~CDFTF1() {}

  double operator()(double *m, double *par) const {

    double x = m[0];

    if (fUseWindowExclusion) {
      if (m[0] >= fExcludeFromFitLow && m[0] <= fExcludeFromFitHigh) {
        TF1::RejectPoint(1);
        return 0.;
      }
    }

    // Following the convention in https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2012-088/
    double f = par[0]*pow(x,-1*par[1])*pow((1-x/fCenterOfMassEnergy),par[2]);
    return f;
  }

};

// ---------------------------------------------------------
class CDFFitFunction : public MjjFitFunction {

public:

  CDFFitFunction(double & mjjLow, double & mjjHigh, double Ecm=13000.) {
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
    CDFTF1 * fCdfparam = new CDFTF1(Ecm);
    fPersonalisedFunction = fCdfparam;
    fFitFunction = new TF1("CDFFitFunction", fCdfparam, fMjjLow, fMjjHigh, fNParameters, "CDFFitFunction");
  }

  ~CDFFitFunction() {
  }

};


//////////////////////////////////////////////////////////
// CDF Function (1997 version)

// ---------------------------------------------------------
class CDF1997TF1 : public MjjTF1Wrapper {

  public:

  CDF1997TF1(double Ecm) {
    fCenterOfMassEnergy = Ecm;
    fUseWindowExclusion = false;
    fExcludeFromFitLow = 0;
    fExcludeFromFitHigh = 0;
  }
  ~CDF1997TF1() {}

  double operator()(double *m, double *par) const {

    double x = m[0];

    if (fUseWindowExclusion) {
      if (m[0] >= fExcludeFromFitLow && m[0] <= fExcludeFromFitHigh) {
        TF1::RejectPoint(1);
        return 0.;
      }
    }

    double f = par[0]*pow(x,-1*par[1])*pow((1-x/fCenterOfMassEnergy-par[3]*pow((x/fCenterOfMassEnergy),2)),par[2]);
    return f;
  }

};

// ---------------------------------------------------------
class CDF1997FitFunction : public MjjFitFunction {

public:

  CDF1997FitFunction(double & mjjLow, double & mjjHigh, double Ecm=13000.) {
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
    CDF1997TF1 * fCdf1997param = new CDF1997TF1(Ecm);
    fPersonalisedFunction = fCdf1997param;
    fFitFunction = new TF1("CDF1997FitFunction", fCdf1997param, fMjjLow, fMjjHigh, fNParameters, "CDF1997FitFunction");
  }

  ~CDF1997FitFunction() {
  }

};


//////////////////////////////////////////////////////////
// D0 Fitting Function

// ---------------------------------------------------------
class D0TF1 : public MjjTF1Wrapper {

  public:

  D0TF1(double Ecm) {
    fCenterOfMassEnergy = Ecm;
    fUseWindowExclusion = false;
    fExcludeFromFitLow = 0;
    fExcludeFromFitHigh = 0;
  }
  ~D0TF1() {}

  double operator()(double *m, double *par) const {

    double x = m[0];

    if (fUseWindowExclusion) {
      if (m[0] >= fExcludeFromFitLow && m[0] <= fExcludeFromFitHigh) {
        TF1::RejectPoint(1);
        return 0.;
      }
    }

    double f = par[0]*pow((x),-1*par[1])*log(par[2]/x)*log(par[3]/pow(x,2));
    return f;
  }

};

// ---------------------------------------------------------
class D0FitFunction : public MjjFitFunction {

public:

  D0FitFunction(double & mjjLow, double & mjjHigh, double Ecm=13000.) {
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
    D0TF1 * fD0func = new D0TF1(Ecm);
    fPersonalisedFunction = fD0func;
    fFitFunction = new TF1("TheD0FitFunction", fD0func, fMjjLow, fMjjHigh, fNParameters, "D0FitFunction");
  }

  ~D0FitFunction() {
  }

};

//////////////////////////////////////////////////////////
// David's Decorrelated Fitting Madness!

// ---------------------------------------------------------
class DecorrelatedTF1 : public MjjTF1Wrapper {

  public:

  typedef enum { THREE_PARAMETER=3, FOUR_PARAMETER=4, SIX_PARAMETER=6 } Variant;
    
  const Variant fVariant;
    
  DecorrelatedTF1(double Ecm,
                  Variant variant=SIX_PARAMETER)
    : fVariant(variant)
  {
    fCenterOfMassEnergy = Ecm;
    fUseWindowExclusion = false;
    fExcludeFromFitLow = 0;
    fExcludeFromFitHigh = 0;
  }
  ~DecorrelatedTF1() {}

  double operator()(double *m, double *par) const {

    double x = m[0];

    if (fUseWindowExclusion) {
      if (m[0] >= fExcludeFromFitLow && m[0] <= fExcludeFromFitHigh) {
        TF1::RejectPoint(1);
        return 0.;
      }
    }

    if (fVariant == SIX_PARAMETER ) {
      /* Latest function from Caterina. Convention wrt p19 of the note
       * par[0] --> p1
       * par[1] --> p2
       * par[2] --> p3
       * par[3] --> p5
       * par[4] --> p6
       * par[5] --> p9
       * */
        double f = (par[0]*pow(1 - par[4]/fCenterOfMassEnergy,-par[1] - par[2]*par[5])*pow(1 - (x/fCenterOfMassEnergy),par[1] + par[2]*par[5])*pow(x/fCenterOfMassEnergy,par[2]))/pow(par[4]/fCenterOfMassEnergy,par[2]);
      return f;
    } else if (fVariant == THREE_PARAMETER) {
        const double par3 = fCenterOfMassEnergy;
        const double par4 = 1220.;
        const double par5 = 4.3;
        double f = par[0]/(pow((1-(par4/par3)),par[1]+par5*par[2])*pow(par4/par3,par[2]))*pow(par4/par3,par[2])*pow(1-(x/par3),par[1]+par5*par[2]);
        return f;
    } else {
      /* Kate's four-parameter variant */
      //double f = par[0] * pow((1-x),par[1]) * pow(x,par[2]) * (pow((1-par[3]),-1*par[1]) * pow(par[3],-1*par[2]));
      double f = par[0] * pow(((1-x/par[3])/(1-1000./par[3])),par[1]) * pow((x/1000.),par[2]);
      return f;
    }
      
  }

};


// ---------------------------------------------------------
class DavidsDecorrelatedFitFunction : public MjjFitFunction {

public:

  DavidsDecorrelatedFitFunction(double & mjjLow, double & mjjHigh, double Ecm=13000., DecorrelatedTF1::Variant variant=DecorrelatedTF1::THREE_PARAMETER) {
    switch (variant) {
    case DecorrelatedTF1::SIX_PARAMETER:
      fNParameters = 6;
      break;
    case DecorrelatedTF1::THREE_PARAMETER:
      fNParameters = 3;
      break;
    case DecorrelatedTF1::FOUR_PARAMETER:
      fNParameters = 4;
      break;
    }
    fFitParameters.clear();
    for (int i=0; i<fNParameters; i++) {
      MjjFitParam * thisFitParam = new MjjFitParam();
      thisFitParam->SetParamDefault(1);
      fFitParameters.push_back(thisFitParam);
    }
    //default range
    fMjjLow = mjjLow;
    fMjjHigh = mjjHigh;
    DecorrelatedTF1 * fDecorrfunc = new DecorrelatedTF1(Ecm,variant);
    fPersonalisedFunction = fDecorrfunc;
    const char* uniqueFuncName = "";
    switch (variant) {
    case DecorrelatedTF1::SIX_PARAMETER:
      uniqueFuncName = "DavidsDecorrelatedFitFunction6";
      break;
    case DecorrelatedTF1::THREE_PARAMETER:
      uniqueFuncName = "DavidsDecorrelatedFitFunction3";
      break;
    case DecorrelatedTF1::FOUR_PARAMETER:
      uniqueFuncName = "DavidsDecorrelatedFitFunction4";
      break;
    }
    fFitFunction = new TF1(uniqueFuncName, fDecorrfunc, fMjjLow, fMjjHigh, fNParameters, uniqueFuncName);
  }

  ~DavidsDecorrelatedFitFunction() {
  }

};

// ---------------------------------------------------------
class FourParamTimesXTF1 : public MjjTF1Wrapper {

  public:

  FourParamTimesXTF1(double Ecm) {
    fCenterOfMassEnergy = Ecm;
    fUseWindowExclusion = false;
    fExcludeFromFitLow = 0;
    fExcludeFromFitHigh = 0;
  }
  ~FourParamTimesXTF1() {}

  double operator()(double *m, double *par) const {

    double x = m[0] / fCenterOfMassEnergy;

    if (fUseWindowExclusion) {
      if (m[0] >= fExcludeFromFitLow && m[0] <= fExcludeFromFitHigh) {
        TF1::RejectPoint(1);
        return 0.;
      }
    }

    // Following the convention in https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2012-088/
    double f = x * (par[0] * pow(1-x , par[1]) * 1./pow(x, par[2] + par[3] * log(x)));

    return f;
  }

};

// ---------------------------------------------------------
class FourParamFitFunctionTimesX : public MjjFitFunction {

public:

  FourParamFitFunctionTimesX(double & mjjLow, double & mjjHigh, double Ecm=13000.) {
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
    FourParamTimesXTF1 * fFourParamTimesX = new FourParamTimesXTF1(Ecm);
    fPersonalisedFunction = fFourParamTimesX;
    fFitFunction = new TF1("TheFourParamFitFunctionTimesX", fFourParamTimesX, fMjjLow, fMjjHigh, fNParameters, "FourParamFitFunctionTimesX");
  }

  ~FourParamFitFunctionTimesX() {
  }

};

// ---------------------------------------------------------

#endif
