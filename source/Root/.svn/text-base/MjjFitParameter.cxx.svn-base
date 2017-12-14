// ---------------------------------------------------------

#include "Bayesian/MjjFitParameter.h"

// ---------------------------------------------------------
MjjFitParam::MjjFitParam() 
{
   fDefault = 0;
   fIsFixed = false;
   fLowerLimit = -999;
   fUpperLimit = -999;
   fHasLowerLimit = false;
   fHasUpperLimit = false;
}

// ---------------------------------------------------------
MjjFitParam::~MjjFitParam() 
{

}

// ---------------------------------------------------------
void MjjFitParam::SetParamDefault(double value) 
{
   fDefault = value; 
   return;
}

// ---------------------------------------------------------
void MjjFitParam::SetParamLimitLow(double low) 
{
   fLowerLimit = low;
   fHasLowerLimit = true;
   return;
}

// ---------------------------------------------------------
void MjjFitParam::SetParamLimitHigh(double high) 
{
   fUpperLimit = high;
   fHasUpperLimit = true;
}

// ---------------------------------------------------------
void MjjFitParam::ActivateParamLimits(bool low, bool high) 
{
   fHasLowerLimit = low;
   fHasUpperLimit = high;
   return;
}

// ---------------------------------------------------------
void MjjFitParam::SetFixParameter(bool doFix)
{
   fIsFixed = doFix;
   return;
}

// ---------------------------------------------------------
std::pair<bool,bool> MjjFitParam::GetHasParamLimits() 
{ 
   return std::make_pair(fHasLowerLimit,fHasUpperLimit);
}
