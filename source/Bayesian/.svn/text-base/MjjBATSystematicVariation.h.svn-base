#ifndef MJJBATSYSTEMATICVARIATION_H
#define MJJBATSYSTEMATICVARIATION_H

/*!
 * \class MjjBATSystematicVariation
 * \brief A class describing a systematic variation.
 * Based on BCMTFSystematicVariation by Daniel Kollar and Kevin Kroeninger.
 * \author Katherine Pachal
 * \date 2013
 * 
 * This class describes the impact of a systematic uncertainty on a 
 * specific process. It is a base class on which several different sorts
 * of systematic variations with different kinds of impacts are defined.
 */

// ---------------------------------------------------------

#include <string>
#include <vector>
#include <TAxis.h>

class TH1D;

// ---------------------------------------------------------
class MjjBATSystematicVariation
{

   public:

      /** \name Constructors and destructors */
      /** @{ */

      /**
       * The default constructor. 
       * Empty, as derived classes need their own constructors. */
      MjjBATSystematicVariation();

      /** 
       * The default destructor. */
      virtual ~MjjBATSystematicVariation() = 0;

      /** @} */
      /** \name Member functions (get) */
      /** @{ */

      /**
       * Get the name of the systematic which corresponds to this variation. */
      std::string GetParentSystematic()
          { return fSystematicName; };

      /**
       * Get bin structure of template or matrix that defines this systematic variation. */
      TArrayD * GetBinStructure()
          { return &fBinStructure; };

      /** @} */
      /** \name Member functions (set) */
      /** @{ */

      /** 
       * Record which systematic corresponds to this variation. */
      void SetParentSystematic(std::string name)
          { fSystematicName = name; };

      /** @} */

   private:

      /**
       * The name of the corresponding source of systematic
       * uncertainty. */
      std::string fSystematicName;

   protected:

      /**
       * The bin structure of all histograms for this systematic. */
      TArrayD fBinStructure;

};
// ---------------------------------------------------------

#endif

