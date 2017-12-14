#ifndef MJJBATSCALECHANGINGSYST_H
#define MJJBATSCALECHANGINGSYST_H

/*!
 * \class MjjBATScaleChangingSyst
 * \brief A class for defining systematic variations that
 * change the scale of a template.
 * Based on BCMTFSystematicVariation by Daniel Kollar and Kevin Kroeninger.
 * \author Katherine Pachal
 * \date 2013
 * 
 * This class defines an MjjBATSystematicVariation
 * for the case that the scale of the signal parameter is changed by
 * the systematic. It uses either bin contents or bin errors multiplied
 * by a specified scale.
 */

// ---------------------------------------------------------

#include <iostream>
#include <cmath>
#include <TVectorD.h>
#include <TH1D.h>
#include <Bayesian/MjjBATSystematicVariation.h>

using std::vector;

// ---------------------------------------------------------
class MjjBATScaleChangingSyst : public MjjBATSystematicVariation {

  public :

   /** \name Constructors and destructors */
   /** @{ */

   /**
    * The constructor. 
    * @param scale is the scale by which the  
    * systematic should change the contents of the template. */
    MjjBATScaleChangingSyst(double scale);

   /** 
    * The default destructor. */
   ~MjjBATScaleChangingSyst();

   /** @} */
   /** \name Member functions (set) */
   /** @{ */

   /**
    * Ensures that the content of each bin in the final shifted
    * template will be multiplied by the scale in the loglikelihood. 
    * This is the default setting. */
   void SetScaleFromBinContent();

   /**
    * Ensures that the error of each bin in the final shifted
    * template will be multiplied by the scale in the loglikelihood. 
    * This is not the default setting and must be specified if desired. */
   void SetScaleFromBinError();

   /** @} */
   /** \name Member functions (get) */
   /** @{ */

   /**
    * Caution: being non-generic for the sake of speed. If at any point background 
    * prediction is coming from weighted MC or anything where the error is not
    * the root of the bin content then this will not hold, but for now it's faster
    * than re-making histograms. 
    * @param bins The bin contents before adjustment.
    * @param param The value of the parameter to use.
    * @return the bin values after the adjustment for the given parameter value. */
   TVectorD GetAdjusted(TVectorD bins, double param);

   /**
    * Caution: being non-generic for the sake of speed. If at any point background 
    * prediction is coming from weighted MC or anything where the error is not
    * the root of the bin content then this will not hold, but for now it's faster
    * than re-making histograms. 
    * @param bincontent The content of the specific bin before adjustment.
    * @param param The value of the parameter to use.
    * @return The bin value after the adjustment for the given parameter value. */
   double GetAdjustedBin(double bincontent, double param);

  private :

   /**
    * To hold the scaled value in one bin. */
   double fBinObserved;

   /**
    * The scale by which to alter bin contents. */
   double fScale;

   /**
    * If true, use fScale with bin content to make adjustments.  */
   bool fUseBinContent;

   /**
    * If true, use fScale with bin errors to make adjustments.  */
   bool fUseBinError;

};
// ---------------------------------------------------------

#endif
