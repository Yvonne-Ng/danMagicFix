#ifndef MJJFITFUNCTIONBASE_H
#define MJJFITFUNCTIONBASE_H

// ---------------------------------------------------------

#include <iostream>
#include <string>
#include <vector>
#include "assert.h"

#include "TH1.h"
#include "TF1.h"
#include "math.h"
#include "Bayesian/MjjHistogram.h"
#include "Bayesian/MjjFitParameter.h"

using std::string;

/*!
 * \class MjjTF1Wrapper
 * \brief A class providing behavioural controls for a user-defined TF1.
 * \author Katherine Pachal
 * \date 2013
 * 
 * One of the two classes which make the base upon which all 
 * FitFunctions used by the package are built. MjjTF1Wrapper
 * is used to define the function itself, with built-in
 * methods for exluding particular windows from the fit
 * and parameters for the window exclusion and center of mass
 * energy. 
 */

// ---------------------------------------------------------
class MjjTF1Wrapper
{

   public :

      /** \name Constructors and destructors */
      /** @{ */

      /**
       * The default constructor. */
      MjjTF1Wrapper() {
        fExcludeFromFitLow = 0;
        fExcludeFromFitHigh = 0;
      }

      /**
       * The default destructor. */
      ~MjjTF1Wrapper() {}


      /** @} */
      /** \name Member functions (get) */
      /** @{ */


      /**
       * Determine whether the window exclusion
       * is activated.
       * @return True if window exclusion activated,
       * false if not.  */
      bool GetDoWindowExclusion()
        { return fUseWindowExclusion; }

      /**
       * @return Low and high edges of window exclusion as
       * pair of doubles.  */
      std::pair <double,double> GetWindowBoundaries()
          { return std::make_pair(fExcludeFromFitLow,fExcludeFromFitHigh); }

      /** @} */

      /** \name Member functions (set) */
      /** @{ */

      /* Set window within which to exclude points from
       * fit. Does NOT automatically activate window exclusion.
       * Simply sets range.
       * @param excLow X-value of low edge of region to exclude.
       * @param excHigh X-value of high edge of region to exclude. */
      void SetExclusionWindowFromRange(double excLow, double excHigh) {
        fExcludeFromFitLow = excLow;
        fExcludeFromFitHigh = excHigh;
        return;
      }

      /* Specify whether or not to exclude window
       * from fit. 
       * @param yesOrNo True if window should be 
       * excluded, else false. */
      void SetDoWindowExclusion(bool yesOrNo) {
        fUseWindowExclusion = yesOrNo;
        return;
      }

   protected:

      /**
       * The center of mass energy of the data to fit. */
      double fCenterOfMassEnergy;

      /**
       * Specifies whether or not to exclude window from fit. */
      bool fUseWindowExclusion;

      /**
       * x-value of the lower edge of window to exclude. */
      double fExcludeFromFitLow;

      /**
       * x-value of the upper edge of window to exclude. */
      double fExcludeFromFitHigh;

};

/*!
 * \class MjjFitFunction
 * \brief A class for manipulating and accessing fit functions.
 * \author Katherine Pachal
 * \date 2013
 * 
 * One of the two classes which make the base upon which all 
 * FitFunctions used by the package are built.
 * Once the function is constructed, MjjFitFunction 
 * defines how it can be used within in the larger code,
 * with methods for getting and setting parameters, controlling
 * the range, translating the function into a histogram, etc.
 * Please note that parameter ranges can be specified here, but
 * are not set using methods of the TF1. All parameters are left
 * unbounded in the TF1 and their behavour in the fit is
 * controlled through Minuit in the class MjjFitter.
 */

// ---------------------------------------------------------
class MjjFitFunction 
{

   public:

      /** \name Constructors and destructors */
      /** @{ */

      /**
       * The default constructor. */
      MjjFitFunction() {}

      /**
       * The default destructor. */
      ~MjjFitFunction() {}


      /** @} */
      /** \name Member functions (get) */
      /** @{ */

      /**
       * @return The TF1 * which underlies this MjjFitFunction. */
      TF1 * GetFitFunction() 
          { return fFitFunction; };

      /**
       * Determine whether the window exclusion
       * is activated.
       * @return True if window exclusion activated,
       * false if not.  */
      bool GetDoWindowExclusion() 
          { return fPersonalisedFunction->GetDoWindowExclusion(); };

      /**
       * @return Low and high edges of window exclusion, if active, as
       * pair of doubles.  */
      std::pair <double,double> GetWindowBoundaries()
          { return fPersonalisedFunction->GetWindowBoundaries(); };

      /**
       * @return Minimum x-value at which FitFunction is defined. */
      double GetMinXVal() 
          { return fMjjLow; };

      /**
       * @return Maximum x-value at which FitFunction is defined. */
      double GetMaxXVal() 
          { return fMjjHigh; };

      /**
       * Returns number of parameters in the TF1 which defines this
       * function. If a signal parameter is being used in the fit,
       * it is NOT counted here. This is parameters in the function alone.
       * @return Number of parameters. */
      int GetNParams() 
          { return fFitFunction->GetNpar(); };

      /**
       * @param index The index of the parameter to retrieve
       * @return A pointer to the MjjFitParam at the specified index. */
      MjjFitParam * GetParameter(int index) 
          { return fFitParameters.at(index); };

      /**
       * Retrieves vector of the current values of all function
       * parameters in order defined.
       * @return The parameters. */
      vector<double> GetCurrentParameterValues();

      /** @} */
      /** \name Member functions (set) */
      /** @{ */

      /**
       * Replaces stored default values for function parameters with
       * values from given vector. Fails with error if the passed
       * vector is not the same length as the number of parameters
       * in the function.
       * @param newParameterDefaults Desired defaults for fit parameters,
       * in the same order as the parameters are in the function. */
      void SetParameterDefaults(vector<double> newParameterDefaults);

      /**
       * Sets current values of function parameters to the values 
       * passed in the vector. Fails with error if the passed
       * vector is not the same length as the number of parameters
       * in the function.
       * @param desiredParameterValues Desired values for fit parameters,
       * in the same order as the parameters are in the function. */
      void SetCurrentParameterValues(vector<double> desiredParameterValues);

      /**
       * Sets the specified lower boundary on the range of the 
       * specified function parameter. Please note the input
       * limit must be a number. To remove limits, use
       * GetParameter()->ActivateParamLimits(false,false)
       * @param paramIndex The index of the parameter to limit.
       * @param limit The new lower bound for the parameter. */
      void SetParameterLimitLow(int paramIndex, double limit);

      /**
       * Sets the specified upper boundary on the range of the 
       * specified function parameter. Please note the input
       * limit must be a number. To remove limits, use
       * GetParameter()->ActivateParamLimits(false,false)
       * @param paramIndex The index of the parameter to limit.
       * @param limit The new upper bound for the parameter. */
      void SetParameterLimitHigh(int paramIndex, double limit);

      /**
       * Specify whether to exclude window from fit
       * @param yesOrNo True if intend to exclude window,
       * false otherwise. */
      void SetDoWindowExclusion(bool yesOrNo);

      /**
       * Specify range in x to exclude from fit, if 
       * SetDoWindowExclusion is activated. Does not 
       * activate automatically when this is called!
       * @param xlow Lower boundary of region to exclude.
       * @param xhigh Upper boundary of region to exclude. */
      void SetExclusionWindowFromRange(double xlow, double xhigh);

      /**
       * Specify window to exclude from fit as most central
       * region of the given histogram containing the given
       * percentage of its contents. Window exclusion does not 
       * activate automatically when this is called!
       * @param signalTemplate Histogram from which to 
       * calculate exclusion window.
       * @param percentage Percentage of histogram used 
       * to calculate exclusion window. */
      void SetExclusionWindowFromHisto(MjjHistogram & signalTemplate, double percentage=0.68);

      /** @} */
      /** \name Member functions (miscellaneous methods) */
      /** @{ */

      /**
       * Reset all parameters to their default values. */
      void RestoreParameterDefaults();

      /**
       * Fill the passed histogram such that the y-value in each bin corresponds
       * to the height of the function in the center of that bin. 
       * Fills histogram within specified range. If no range given, uses range
       * of function.
       * @param blankHist The TH1D* to fill.
       * @param xmin The lower end of the region to fill.
       * @param xmax The upper end of the region to fill.  */
      void MakeHistFromFunction(TH1D* blankHist, double xmin=-1, double xmax=-1, bool overwrite = true);

      /** @} */

   protected:

      /**
       * The TF1 underlying this MjjFitFunction. */
      TF1 * fFitFunction;

      /**
       * The MjjTF1Wrapper specifying behaviour of this MjjFitFunction. */
      MjjTF1Wrapper * fPersonalisedFunction;

      /**
       * The lower boundary of definition of the function. */
      double fMjjLow;

      /**
       * The upper boundary of definition of the function. */
      double fMjjHigh;

      /**
       * The number of parameters in the TF1. */
      double fNParameters;

      /**
       * The vector of MjjFitParam objects which store the default
       * values and ranges of the fit parameters. */
      vector<MjjFitParam * > fFitParameters;

};

#endif
