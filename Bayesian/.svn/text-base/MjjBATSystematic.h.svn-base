#ifndef MJJBATSYSTEMATIC_H
#define MJJBATSYSTEMATIC_H

/*!
 * \class MjjBATSystematic
 * \brief A class desribing a systematic uncertainty.
 * Based on BCMTFSystematicVariation by Daniel Kollar and Kevin Kroeninger.
 * \author Katherine Pachal
 * \date 2013
 * 
 * This class describes a systematic uncertainty.
 * It does not specify the effects of the systematic, as it may
 * have different effects on different processes. It simply
 * defines the common systematic object whose effects can be
 * specified using an MjjBATSystematicVariation class.
 */

// ---------------------------------------------------------

#include <string>

// ---------------------------------------------------------
class MjjBATSystematic 
{

    public:

      /** \name Constructors and destructors */
      /** @{ */

      /**
       * The default constructor. 
       * @param name The name of the source of systematic uncertainty. */
      MjjBATSystematic(std::string name);

      /** 
       * The default destructor. */
      ~MjjBATSystematic();

      /** @} */
      /** \name Member functions (get) */
      /** @{ */

      /**
       * @return The name of the systematic uncertainty. */
      std::string GetName()
         { return fSystematicName; };

      /**
       * @return Index of parameter corresponding to systematic.*/
      int GetParamIndex()
         { return fParamIndex; };

      /** 
       * @return A flag defining if this uncertainty is active or not. */
      bool GetFlagSystematicActive()
         { return fFlagSystematicActive; };

      /** @} */
      /** \name Member functions (set) */
      /** @{ */

      /**
       * Store index of parameter that was created for this process.
       * @param index The parameter index. */
      void SetSystematicParameterIndex(int index)
         { fParamIndex = index; };


      /**
       * Set a flag defining if this uncertainty is active or not.
       * @param flag The flag. */ 
      void SetFlagSystematicActive(bool flag)
         { fFlagSystematicActive = flag; };

      /** @}*/

   private:

      /**
       * The name of the source of the systematic uncertainty. */
      std::string fSystematicName;

      /**
       * Index of parameter in model which corresponds to this systematic. */
      int fParamIndex;

      /**
       * A flag defining if this uncertainty is active or not. */
      bool fFlagSystematicActive;

};
// ---------------------------------------------------------

#endif

