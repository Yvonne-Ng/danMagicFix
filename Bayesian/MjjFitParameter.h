#ifndef MJJFITPARAMETER_H
#define MJJFITPARAMETER_H

/*!
 * \class MjjFitParam
 * \brief A class to store parameter limits and defaults
 * for input into fits
 * \author Katherine Pachal
 * \date 2013
 * 
 * This class is used by MjjFitFunction to store
 * default values, limits, and their state of activation
 * for fit parameters. Note that these constraints are
 * almost entirely implemented via the MjjFitter.
 * Parameters in the TF1 stored in MjjFitFunction objects
 * are never limited or fixed; such behaviour is
 * controlled only in the fitter.
 */

// ---------------------------------------------------------

#include <iostream>
#include <vector>

// ---------------------------------------------------------
class MjjFitParam 
{

   public:

      /** \name Constructors and destructors */
      /** @{ */

      /**
       * The default constructor. */
      MjjFitParam();

      /**
       * The default destructor. */
      ~MjjFitParam(); 

      /** @} */
      /** \name Member functions (get) */
      /** @{ */

      /**
       * @return The default value of this parameter. */
      double GetParamDefault()
          { return fDefault; }

      /**
       * @return A bool indicating whether the parameter
       * should be fixed to its default value (true) or
       * left to float (false) */
      bool GetIsFixed()
          { return fIsFixed; };

      /**
       * Tells the user whether the parameter has 
       * activated lower and upper limits.
       * @return A pair of bools indicating state
       * of lower limit, upper limit respectively.
       * True if limit active, false otherwise. */
      std::pair<bool,bool> GetHasParamLimits();

      /**
       * @return The lower limit of this parameter. 
       * Returns -999 if no lower limit set. */
      double GetParamLimitLow()
          { return fLowerLimit; }

      /**
       * @return The upper limit of this parameter. 
       * Returns -999 if no upper limit set. */
      double GetParamLimitHigh()
          { return fUpperLimit; }

      /** @} */
      /** \name Member functions (set) */
      /** @{ */

      /**
       * Set the default value for this parameter.
       * @param value Desired parameter default. */
      void SetParamDefault(double value);

      /**
       * Set the lower limit for this parameter,
       * and activates lower limit.
       * @param low Desired lower parameter limit. */
      void SetParamLimitLow(double low);

      /**
       * Set the upper limit for this parameter,
       * and activates upper limit.
       * @param high Desired upper parameter limit. */
      void SetParamLimitHigh(double high); 

      /** @} */
      /** \name Member functions (miscellaneous methods) */
      /** @{ */

      /**
       * Activate or deactivate the limits on a parameter.
       * Use false to deactivate and true to activate.
       * This method permits one-sided limits: just leave
       * the appropriate edge of the parameter range free.
       * @param low Do or do not activate lower bound on parameter.
       * @param high Do or do not activate upper bound on parameter. */
      void ActivateParamLimits(bool low, bool high);

      /**
       * @param doFix Specifies whether to fix or unfix
       * this parameter in the fit. */
      void SetFixParameter(bool doFix);

      /** @} */

   private:

     /**
      * The default parameter value. */
     double fDefault;

     /**
      * True if parameter should be fixed to 
      * default value, else false. */
     bool fIsFixed;

     /**
      * True if parameter has lower limit;
      * false if unbounded. */
     bool fHasLowerLimit;

     /**
      * True if parameter has upper limit;
      * false if unbounded. */
     bool fHasUpperLimit;

     /**
      * Value of parameter lower limit. */
     double fLowerLimit;

     /**
      * Value of parameter upper limit. */
     double fUpperLimit;

};

// ---------------------------------------------------------

#endif

