#ifndef MJJBATTEMPLATESYST_H
#define MJJBATTEMPLATESYST_H

// ---------------------------------------------------------

#include <iostream>
#include <algorithm>
#include <TH1D.h>

#include <Bayesian/MjjBATSystematicVariation.h>

using std::vector;

/*!
 * \struct sort_pairs_TH1D
 * \brief A struct for sorting pairs of doubles and TH1Ds used by 
 * MjjBATTemplateSyst.
 * \author Katherine Pachal
 * \date 2013
 * 
 * This struct sorts a vector of pairs of doubles and TH1Ds according
 * to the size of the double, from smallest to largest. 
 */

// ---------------------------------------------------------

struct sort_pairs_TH1D {

   /**
    * Operator actually used in sorting vector of pairs.
    * \return true if first pair has smaller double value than second. */
   bool operator()(const std::pair<double,TH1D*> &left, const std::pair<double,TH1D*> &right) {
      return left.first < right.first;
   }

};

/*!
 * \class MjjBATTemplateSyst
 * \brief A class for defining systematic variations that
 * change the shape of a template using another template.
 * Based on BCMTFSystematicVariation by Daniel Kollar and Kevin Kroeninger.
 * \author Katherine Pachal
 * \date 2013
 * 
 * This class defines an MjjBATSystematicVariation
 * for the case that the nominal template is varied by some scaled
 * version of another template. This cannot be mixed with the 
 * transfer matrix MjjBATShapeChangingSyst class, and only 
 * one MjjBATTemplateSyst can be defined per process.
 */

// ---------------------------------------------------------
class MjjBATTemplateSyst : public MjjBATSystematicVariation {

   public :

      /** \name Constructors and destructors */
      /** @{ */

      /**
       * The default constructor. */
      MjjBATTemplateSyst(bool isScale);

      /** 
       * The default destructor. */
      ~MjjBATTemplateSyst();

      /** @} */
      /** \name Member enums and structs */
      /** @{ */

      /** @} */
      /** \name Member functions (get) */
      /** @{ */

      /**
       * @param param Parameter value for calculation.
       * @return A histogram constituting the adjustment for the given parameter value. */
      vector<double> GetAdjustment(double param);

      /**
       * @param bin The bin in which to find adjustment.
       * @param param Parameter value for calculation.
       * @return The adjustment for the specified bin and the given parameter value. */
      double GetBinAdjustment(int bin, double param);

      bool GetIsScale()
          { return fIsScale; };

      /**
       * @return A vector of the shifted spectra used to define the 
       * effects of this systematic. */
      vector<TH1D> GetSpectra()
          { return fDiscreteSpectra; };

      /** @} */
      /** \name Member functions (set) */
      /** @{ */

      /**
       * Set collection of spectra from which to extrapolate for adjustments
       * from systematic.
       * @param sigmasAndSpectra Vector of pairs with first being signed sigma 
       * of fluctuation and second being spectrum at that fluctuation. */
      void SetSpectra(vector<std::pair<double, TH1D*> > sigmasAndSpectra);

   private :

      /**
       * Vector of spectra for different values of systematic. */
      vector<TH1D> fDiscreteSpectra;

      /**
       * Vector of fluctuation sizes in sigmas corresponding to spectra. */
       vector<double> fCorrespondingSigmas;

      /**
       * Index of the spectrum corresponding to no variation (plain data) */
      int fOriginalIndex;

      /**
       * Number of bins in spectra. */
      int fNBinsInSpectra;

      bool fIsScale;
};
// ---------------------------------------------------------

#endif
