#ifndef MJJBATMODEL_H
#define MJJBATMODEL_H

/*!
 * \class MjjBATModel
 * \brief A class for fitting several templates to a data set.
 * Based on BCMTF by Daniel Kollar and Kevin Kroeninger.
 * \author Katherine Pachal
 * \date 2013
 * 
 * This class can be used for fitting several template
 * histograms to a data histogram. The templates are assumed to have
 * no statistical uncertainty whereas the data are assumed to have
 * Poissonian fluctuations in each bin. This method takes MjjHistogram
 * formatted data and can therefore handle weights. There are two
 * available methods for expressing systematic uncertainties: templates
 * or matrices. In the case that more than one systematic which 
 * changes the shape of the signal template is used, all such systematics
 * must be expressed as matrices. Templates are only permitted when
 * there is a single shape-changing systematic.
 */

// ---------------------------------------------------------

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <iomanip>

#include <TCanvas.h>
#include <THStack.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TAxis.h>
#include <TF1.h>
#include <TMath.h>
#include <TGraphAsymmErrors.h>

#include "BAT/BCMath.h"
#include "BAT/BCLog.h"
#include "BAT/BCModel.h"
#include "BAT/BCParameter.h"

#include "Bayesian/MjjHistogram.h"
#include "Bayesian/MjjBATProcess.h"
#include "Bayesian/MjjBATSystematic.h"
#include "Bayesian/MjjBATSystematicVariation.h"
#include "Bayesian/MathFunctions.h"

enum LikeType { DEFAULT=1, BUMPHUNTER=2 };


// ---------------------------------------------------------
class MjjBATModel : public BCModel
{

   public:

      /** \name Constructors and destructors */
      /** @{ */

      /**
       * The default constructor. */
      MjjBATModel();

      /**
       * A constructor.
       * @param name The name of the model */
      MjjBATModel(std::string name);

      /**
       * The default destructor. */
      ~MjjBATModel();

      /** @} */
      /** \name Member functions (get) */
      /** @{ */

      /**
       * @return The number of processes. */
      int GetNProcesses()
         { return fNProcesses; };

      /**
       * @return The number of systematics. */
      int GetNSystematics()
         { return fNSystematics; };

      /**
       * @param name The name of the process.
       * @return The process index. */
      int GetProcessIndex(std::string name);

      /**
       * @param name The name of the systematic.
       * @return The systematic uncertainty index. */
      int GetSystematicIndex(std::string name);

      /**
       * @param index The parameter index (mtf counting) .
       * @return The parameter number corresponding to the parameter index (BAT counting). */
      vector<int> GetParIndicesProcess(int index)
         { return GetProcess(index)->GetParamIndices(); };

      /**
       * @param index The systematic uncertainty index (mtf counting).
       * @return The parameter number corresponding to the systematic uncertainty (BAT counting). */
      int GetParIndexSystematic(int index)
         { return GetSystematic(index)->GetParamIndex(); };

      /**
       * @param index The process index.
       * @return The process object. */			
      MjjBATProcess * GetProcess(int index)
         { return fProcessContainer.at(index); };

      /**
       * @param index The systematic uncertainty index.
       * @return The systematic uncertainty object. */			
      MjjBATSystematic * GetSystematic(int index)
         { return fSystematicContainer.at(index); };

      /**
       * @return The MjjHistogram defining the data. */
      MjjHistogram GetMjjHistogram() { return fMjjHist; };

      /**
       * @return The data. */
      TH1D * GetData() { return fData; };

      /**
       * @return The effective data. */
      TH1D * GetEffectiveData() { return fEffectiveData; };

      /**
       * @return The histogram of weights. */
      TH1D * GetWeights() { return fWeights; };

      /**
       * @return The total normalization of each process (in order)
       * for the specified parameter values. */
      vector<double> GetProcessNorms(const std::vector<double> params);

      /** @} */
      /** \name Member functions (set) */
      /** @{ */

      /**
       * Set the data histogram
       * @param hist The Mjj histogram.
       * @param firstBinToUse The first bin to include in comparisons. 
       * If left blank, sets automatically to first bin with data.
       * @param lastBinToUse The last bin to include in comparisons.
       * If left blank, sets automatically to last bin with data.
       * @param minimum The minimum number of expected events (used for calculation of uncertainty bands).
       * If left blank, sets automatically to zero.
       * @param maximum The maximum number of expected events (used for calculation of uncertainty bands).
       * If left blank, sets automatically to the maximum value of the histogram.
       * @return An error code. */			
      int SetData(MjjHistogram hist, int firstBinToUse=-1, int lastBinToUse=-1, double minimum=-1, double maximum=-1);

      /**
       * Set the template for a specific process using a central histogram.
       * This method can only be used when no normalisation nuisance parameter
       * is to be considered.
       * @param processname The name of the process.
       * @param centralhist The central (nominal) hist for this process.
       * @return An error code. */
      int SetTemplate(std::string processname, TH1D centralhist);

      /**
       * Set the template for a specific process using a central and a variation histogram.
       * @param processname The name of the process.
       * @param centralhist The central (nominal) hist for this process.
       * @param variationhist The histogram of the difference between nominal and nominal+1sigma
       * @return An error code. */
      int SetTemplate(std::string processname, TH1D centralhist, TH1D variationhist);

      /**
       * Set the template for a specific process using a template histogram and a function.
       * This is in development -- I don't trust it and nor should you.
       * @param processname The name of the process.
       * @param function An MjjFitFunction which will be used to define the template shape. 
       * NParams for process = NParams for function.
       * @param templatehist A histogram with the bin structure desired for the process.
       * @param central The central or nominal value for each parameter in the function
       * @param error The error on each parameter in the function, which will set their ranges
       * @param luminosity The luminosity to use.
       * @return An error code. */
      int SetTemplate(std::string processname, MjjFitFunction & function, TH1D templatehist, vector<double> central, vector<double> error, double luminosity);

      /**
       * Set the impact of a source of systematic uncertainty for a
       * particular systematic and process. The input
       * MjjBATSystematicVariation specified how it will be processed.
       * @param processname The name of the process.
       * @param systematicname The name of the source of systematic uncertainty.
       * @param systVariation An MjjBATSystematicVariation object that defines the
       * effect of the systematic on this process 
       * @return An error code. */			
      int SetSystematicVariation(std::string processname, std::string systematicname, MjjBATSystematicVariation * systVariation);
  
      /**
       * Specify a range of bins smaller than the overall range to consider in the likelihood
       * calculation. This allows tools such as the BumpHunter to examine only sub-ranges at
       * any time.
       * @param firstBin The lowest bin that will be included.
       * @param lastBin The highest bin that will be included.
       * @return An error code. */
      int SetInterestingBins(int firstBin, int lastBin);
  
      void SetLikelihoodChoice(LikeType type)
          { fLikelihood = type; };
  
      void SetDebug(bool doDebug)
          { fDebug = doDebug; };


      /** @} */
      /** \name Member functions (miscellaneous methods) */
      /** @{ */

      /**
       * Create a process and add the associated BAT parameter (s). Use this
       * one when defining a process from a function rather than a template.
       * @param process The process name
       * @param min Vector of minimum sigmas for each parameter associated with process
       * @param max Vector of maximum sigmas for each parameter associated with process
       * @return An error code. */			
     int AddProcess(std::string process, vector<double> min, vector<double> max);

      /**
       * Create a process and add the associated BAT parameter(s). Use this
       * one when the process is defined by a histogram acting as a nominal
       * template. If the parameter range is included, it will be associated
       * to a nusiance parameter for the normalisation of the histogram. This
       * uncertainty in the normalisation must then be specified by a histogram
       * to be added in a later step.
       * @param name The process name
       * @param npar The number of parameters associated with this process. 
       * Use 1 if defining process from template.
       * @param nmin Minimum number of sigmas in variation template or each function parameter uncertainty
       * @param nmax Maximum number of sigmas in variation template or each function parameter uncertainty
       * @return An error code. */
     int AddProcess(std::string name, bool doUnc=true, double nmin = -3., double nmax = 3.);

      /**
       * Add a source of systematic uncertainty and the associated BAT (nuisance) parameter.
       * @param name The systematic uncertainty name.
       * @param min The lower limit on the BAT parameter values, typically -5 sigma if Gaussian constraint is used.
       * @param max The upper limit on the BAT parameter values, typically +5 sigma if Gaussian constraint is used.
       * @return An error code. */			
      int AddSystematic(std::string name, double min = -3., double max = 3.);

      /**
       * Return the expected number of events across all bins using the given parameters.
       * @param parameters A reference to the parameters used to calculate the expectation.
       * @return The expectation value as a vector of bin contents. Starts at bin 0 (underflow) */			
      vector<double> Expectation(const std::vector<double> & parameters);

      /**
       * Return the expected number of events across all bins for the specified process
       * using the given parameters.
       * @param processindex The index of the process in question.
       * @param parameters A reference to the parameters used to calculate the expectation.
       * @return The expectation value as a vector of bin contents. Starts at bin 0 (underflow) */
      std::vector<std::pair<double,double> > ProcessExpectation(int processindex, const std::vector<double> & parameters);

      /** 
       * Compare binning of two histograms. 
       * @param firstBins The bins of the first histogram.
       * @param secondBins The bins of the second histogram.
       * @return True if the structure is identical, else false. */
      bool DoBinsMatch(const TArrayD * firstBins, const TArrayD * secondBins);

      /** 
       * Compare binning of histogram to established binning. 
       * @param firstBins The bins of the histogram to test.
       * @return True if the structure is identical, else false. */
      bool DoBinsMatchStored(const TArrayD * firstBins);

      /**
       * Calculate a chi2 given a set of parameters.
       * @param parameters A reference to the parameters used to calculate the chi2.
       * @return A chi2 value. */
      double CalculateChi2(const std::vector<double> & parameters);

      /**
       * Calculate the Cash statistic
       * @param parameters A reference to the parameters used to calculate the Cash statistic.
       * @return The Cash statistic. */
      double CalculateCash(const std::vector<double> & parameters);

      /** @} */
      /** \name Member functions (output methods) */
      /** @{ */

      /**
       * Print a summary of the fit into an ASCII file.
       * @param filename The name of the file.
       * @return An error code */			
      int PrintSummary(std::string filename = "summary.txt");

      /**
       * Print the stack of templates together with the data.
       * Several plot options are available:\n
       * "logx" : plot the x-axis on a log scale \n
       * "logy" : plot the y-axis on a log scale \n
       * "bw"   : plot in black and white \n
       * "sum"  : draw a line corresponding to the sum of all templates \n
       * "stack" : draw the templates as a stack \n
       * "e0" : do not draw error bars \n
       * "e1" : draw error bars corresponding to sqrt(n) \n
       * "b0" : draw an error band on the expectation corresponding to the central 68% probability \n
       * "b1" : draw bands showing the probability to observe a certain number of events given the expectation. The green (yellow, red) bands correspond to the central 68% (95%, 99.8%) probability \n
       * @param parameters A reference to the parameters used to scale the templates.
       * @param filename The output filename. By default, "stack.eps."
       * @param options The plotting options. By default, "e1b0stack."
       * @return An error code. */			
      int PrintStack(const std::vector<double> & parameters, std::string filename = "stack.eps", std::string options = "e1b0stack");

      /**
       * Print the templates in this channel.
       * @param filename The name of the file. */
      void PrintTemplates(std::string filename);

      /** 
       * Print histogram for uncertainty band calculation.
       * @param filename The name of the file. */
      void PrintHistUncertaintyBandExpectation(std::string filename);

      /** 
       * Print histogram for uncertainty band calculation.
       * @param filename The name of the file. */
      void PrintHistUncertaintyBandPoisson(std::string filename);

      /**
       * Print cumulative histogram for uncertainty band calculation.
       * @param filename The name of the file. */
      void PrintHistCumulativeUncertaintyBandPoisson(std::string filename);

      /**
       * Set the y-ranges for printing. 
       * @param min The minimum range. 
       * @param max The maximum range. */
      void SetRangeY(double min, double max) {
         fRangeYMin = min; fRangeYMax = max; };

      /** @} */
      /** \name Member functions (overloaded from BCModel) */
      /** @{ */

      /**
       * Calculates natural logarithm of the likelihood.
       * This model defined it using weighted data and a
       * set range of bins.
       * @param parameters A set of parameter values
       * @return Natural logarithm of the likelihood */
      double LogLikelihood(const std::vector<double> & parameters);

      /**
       * Calculates natural logarithm of the likelihood.
       * This model is prepared to calculate the Poisson p-value 
       * used in the BumpHunter for one background process and no signal.
       * @param parameters A set of parameter values
       * @return Natural logarithm of the likelihood */
      double BHLogLikelihood(const std::vector<double> & parameters);

      /**
       * Calculates log(Poisson) in fewer steps to retain
       * sensitivity to very large numbers.
       * For use in LikeMin. */
      double SimpleLogPoisson(double x, double par);

      /**
       * Method executed for every iteration of the MCMC. User's code should be
       * provided via overloading in the derived class*/
      void MCMCUserIterationInterface();

      /** @} */

 private:

      /** 
       * Calculate histogram for uncertainty band calculation. */
      void CalculateHistUncertaintyBandPoisson();

      /**
       * Calculate histogram for uncertainty band calculation and
       * return a TH1D.
       * @param minimum The minimum value on the expectation.
       * @param maximum The maximum value on the expectation.
       * @param color The color scheme.
       * @return A TH1D histogram. */
      TH1D* CalculateUncertaintyBandPoisson(double minimum, double maximumm, int color);

      /**
       * A container of processes. */
      std::vector<MjjBATProcess *> fProcessContainer;

      /** 
       * A container of sources of systematic uncertainy. */
      std::vector<MjjBATSystematic *> fSystematicContainer;
  
      /**
       * Specifies the likelihood to use, when there are several. */
      LikeType fLikelihood;


      /**
       * The number of processes. */
      int fNProcesses;

      /**
       * The number of systematics uncertainties. */
      int fNSystematics;

      /**
       * The data histogram. */
      TH1D * fData;

      /**
       * The effective data histogram. */
      TH1D * fEffectiveData;

      /**
       * The histogram of weights relating data to effective data. */
      TH1D * fWeights;

      /**
       * The MjjHistogram containing all of the above. */
      MjjHistogram fMjjHist;

      /**
       * The first bin used in calculations. */
      int fFirstBin;

      /**
       * The last bin used in calculations. */
      int fLastBin;
  
      /**
       * The first interesting bin to be considered in upcoming likelihood
       * calculations. If only a subrange is to be considered this can be
       * smaller than the overall range of the spectrum.*/
      int fFirstInterestingBin;
  
      /**
       * The last interesting bin to be considered in upcoming likelihood
       * calculations. If only a subrange is to be considered this can be
       * smaller than the overall range of the spectrum.*/
      int fLastInterestingBin;

      /**
       * The bin structure which all signals must have in common with the data.*/
      TArrayD fBinStructure;

      /** 
       * The minimal y-range for printing. */
      double fRangeYMin;

      /**
       * The maximal y-range for printing. */
      double fRangeYMax;

      /**
       * A histogram for the calculation of uncertainty bands. */
      TH2D* fHistUncertaintyBandExpectation;

      /** 
       * A histogram for the calculation of uncertainty bands. */
      TH2D* fHistUncertaintyBandPoisson;

      /** 
       * A TVectorD for holding the expectation calculated in 
       * individual processes. */
      TVectorD fProcessExpectation;

      /**
       * Control print statements. */
      bool fDebug;


};
// ---------------------------------------------------------

#endif

