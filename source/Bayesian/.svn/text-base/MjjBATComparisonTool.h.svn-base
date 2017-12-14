#ifndef MJJBATCOMPARISONTOOL_H
#define MJJBATCOMPARISONTOOL_H

/*!
 * \class MjjBATComparisonTool
 * \brief A helper class for MjjBATAnalysisFacility storing information.
 * Unmodified from BCMTFComparisonTool.
 * \author Daniel Kollar
 * \author Kevin Kr&ouml;ninger
 * 
 * This is a helper class for MjjBATAnalysisFacility storing information.
 * Identical to the class BCMTFComparisonTool, but generated for the MjjBATModel. 
 */

// ---------------------------------------------------------
#include <string>
#include <vector>
#include <iostream>

#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TLatex.h>


// ---------------------------------------------------------
class MjjBATComparisonTool
{

   public:

      /** \name Constructors and destructors */
      /** @{ */

      /**
       * The default constructor. 
       * @param name The name of the class. */
      MjjBATComparisonTool(std::string name);

      /** 
       * The defaul destructor. */
      ~MjjBATComparisonTool();

      /** @} */
      /** \name Member functions (get) */
      /** @{ */

      /**
       * @return The name of the class. */
      std::string GetName()
         { return fName; };

      /**
       * @return The number of contributions. */
      int GetNContributions()
         { return (int) fHistogramContainer.size(); };

      /** @} */
      /** \name Member functions (miscellaneous methods) */
      /** @{ */

      /**
       * Add a constribution. 
       * @param name The name of the contribution.
       * @param hist The histogram. */
      void AddContribution(std::string name, TH1D hist);

      /**
       * Add a constribution. 
       * @param name The name of the contribution.
       * @param centralvalue The central value.
       * @param uncertainty The uncertainty. */
      void AddContribution(std::string name, double centralvalue, double uncertainty);

      /** 
       * Draw an overview. */ 
      void DrawOverview();

      /** @} */
      /** \name Member functions (output methods) */
      /** @{ */

      /**
       * Print all histograms to a file. 
       * @param filename The name of the file. */
      void PrintHistograms(std::string filename);

      /**
       * Print an overview to a file.
       * @param filename The name of the file. */
      void PrintOverview(std::string filename);

      /** @} */

   private:

      /**
       * The name of the class. */
      std::string fName;

      /**
       * The names of the contributions. */
      std::vector<std::string> fNameContainer;

      /** 
       * A container of TH1D histograms. */
      std::vector<TH1D *> fHistogramContainer;

      /**
       * A container of central values. */
      std::vector<double> fCentralValueContainer;

      /**
       * A container of uncertainties. */
      std::vector<double> fUncertaintyContainer;

};
// ---------------------------------------------------------

#endif

