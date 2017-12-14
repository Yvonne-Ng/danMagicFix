#ifndef MJJBATSHAPECHANGINGSYST_H
#define MJJBATSHAPECHANGINGSYST_H

// ---------------------------------------------------------

#include <iostream>
#include <algorithm>
#include <TH2D.h>
#include <TMatrixDSparse.h>
#include <TArrayD.h>

#include <Bayesian/MjjBATSystematicVariation.h>

using std::vector;

/*!
 * \struct sort_pairs_TH2D
 * \brief A struct for sorting pairs of doubles and TH2Ds used by 
 * MjjBATShapeChangingSyst.
 * \author Katherine Pachal
 * \date 2013
 * 
 * This struct sorts a vector of pairs of doubles and TH2Ds according
 * to the size of the double, from smallest to largest. 
 */

// ---------------------------------------------------------

struct sort_pairs_TH2D {

   /**
    * Operator actually used in sorting vector of pairs.
    * \return true if first pair has smaller double value than second. */
   bool operator()(const std::pair<double,TH2D*> &left, const std::pair<double,TH2D*> &right) {
      return left.first < right.first;
   }

}; 

/*!
 * \class MjjBATShapeChangingSyst
 * \brief A class for defining systematic variations that
 * change the shape of a template.
 * Based on BCMTFSystematicVariation by Daniel Kollar and Kevin Kroeninger.
 * \author Katherine Pachal
 * \date 2013
 * 
 * This class defines an MjjBATSystematicVariation
 * for the case that the shape of the signal parameter is changed by
 * the systematic. It is based on transfer matrices.
 */

// ---------------------------------------------------------
class MjjBATShapeChangingSyst : public MjjBATSystematicVariation {

   public :

      /** \name Constructors and destructors */
      /** @{ */

      /**
       * The default constructor. */ 
      MjjBATShapeChangingSyst();

      /** 
       * The default destructor. */
      ~MjjBATShapeChangingSyst();

      /** @} */
      /** \name Member functions (get) */
      /** @{ */

      /**
       * @param param The value of the parameter for which to find a matrix.
       * @return The adjustment for the specified bin and the given parameter value. */
      TMatrixDSparse GetMatrix(double param);

      /**
       * @return A vector of the matrices used to define the systematic. */
      vector<TMatrixDSparse> GetMatrices()
          { return fTransferMatrices; };

      /** @} */
      /** \name Member functions (set) */
      /** @{ */

      /**
       * Set collection of spectra from which to extrapolate for adjustments
       * from systematic.
       * @param sigmasAndSpectra Vector of pairs with first being signed sigma 
       * of fluctuation and second being spectrum at that fluctuation. */
      void SetSpectra(vector<std::pair<double, TH2D*> > sigmasAndSpectra);

   private :

      /**
       * Vector of spectra for different values of systematic. */
      vector<TMatrixDSparse> fTransferMatrices;

      /**
       * Vector of arrays to hold rows in calculations. */
      vector<TArrayI> fStoreRows;

      /**
       * Vector of arrays to hold columns in calculations. */
       vector<TArrayI> fStoreColumns;

      /**
       * Vector of arrays to hold data in calculations. */
      vector<TArrayD> fStoreData;

      /**
       * Vector of degree of fluctuation in sigmas corresponding to spectra. */
      vector<double> fCorrespondingSigmas;

      /**
       * Vector of differences between neighbouring sigma intervals, calculated
       * in advance to speed up processing time. */
      vector<double> fSigmaDifferences;

      /**
       * Vector of differences between neighbouring matrices, calculated
       * in advance to speed up processing time. */
      vector<TMatrixDSparse> fMatrixDifferences;

};
// ---------------------------------------------------------

#endif
