#ifndef MJJBATPROCESS_H
#define MJJBATPROCESS_H

/*!
 * \class MjjBATProcess
 * \brief A class describing a process.
 * Based on BCMTFProcess by Daniel Kollar and Kevin Kroeninger.
 * \author Katherine Pachal
 * \date 2013
 * 
 * This class describes a process: signal or one 
 * background, etc, which will be added to other processes to
 * define the precition for data in the MjjBATModel.
 */

// ---------------------------------------------------------

#include <string>
#include <assert.h>

#include <TVectorD.h>

#include "Bayesian/MjjBATSystematicVariation.h"
#include "Bayesian/MjjBATShapeChangingSyst.h"
#include "Bayesian/MjjBATScaleChangingSyst.h"
#include "Bayesian/MjjBATTemplateSyst.h"
#include "Bayesian/MjjFitFunction.h"

// ---------------------------------------------------------
class MjjBATProcess
{
   public:

      /** \name Constructors and destructors */
      /** @{ */

      /**
       * The default constructor. 
       * name The name of the process. 
       * doNormalisationUnc Sets whether or not a nuisance
       * parameter should be associated with the normalisation of the
       * process. This should be set by including a histogram with bin
       * contents set to the relevant 1 sigma normalisation uncertainty.
       * doStatUnc Sets whether or not statistical uncertainties should
       * be considered in the process. If so, they will be calculated from
       * the uncertainty on the input nominal histogram, so make sure the
       * errors are appropriate! */
      MjjBATProcess(std::string name, bool doUnc = true); // , bool doStatUnc = false); // left over from old tests

      /**
       * The default destructor. */
      ~MjjBATProcess();

      /** @} */
      /** \name Member functions (get) */
      /** @{ */

      /**
       * @return The name of the process. */
      std::string GetName()
         { return fName; };

      /**
       * @return Indices of parameters corresponding to process.*/
      vector<int> GetParamIndices()
         { return fParamIndices; };

      /**
       * Returns index of the parameter corresponding to the normalisation of
       * the sample. If there is no such parameter this returns -1.
       * The normalisation parameter is also included in the list returned by
       * GetParamIndices.
       * @return Normalisation parameter index. */
      int GetNormalisationParamIndex()
         { return fNormalisationParameter; };
       

      /**
       * @return Names of parameters corresponding to process.*/
      vector<std::string> GetParamNames()
         { return fParamNames; };

      /**
       * @return The nominal TH1D histogram. */
      TH1D * GetNominalHistogram()
         { return fNominal; };

      /**
       * @return The +1sigma TH1D histogram. */
      TH1D * GetVariationHistogram()
         { return fVariation; };

      /**
       * @return All the systematic variations which change the shape of this process. */
      vector<MjjBATShapeChangingSyst*> GetShapeChangingSysts() 
         { return shapeVariations; };

      /**
       * @return All the systematic variations which change the scale of this process. */
      vector<MjjBATScaleChangingSyst*> GetScaleChangingSysts() 
         { return scaleVariations; };

      /**
       * @return All the systematic variations which change the scale of this process. */
      vector<MjjBATTemplateSyst*> GetTempScaleChangingSysts()
         { return templateScaleVariations; };

      /**
       * @return A systematic variation which changes the process using template-based 
       * shifts, if one exists. */
      MjjBATTemplateSyst * GetTemplateSyst() 
         { return templateVariation; };

      /**
       * @param bin The bin in question.
       * @param params The parameter values to use.
       * @return Content of given bin for process at specified parameter value(s). */
      double GetBinValue(int bin, vector<double> params);

      /**
       * @param params The parameter values to use.
       * @return Array of all bin values for process at specified parameter value(s). */
      TVectorD GetBins(vector<double> params);

      /**
       * @param params The parameter values to use.
       * @return Overall normalisation of process given specified parameter value(s). */
      double GetNEvents(vector<double> params);

      /**
       * @return The function (if it exists). */
      MjjFitFunction * GetFunction()
         { return fFunction; };

      /**
       * @return true if the process expectation should be trimmed to central specified amount. */
      bool GetTrimProcess()
         { return fTrimProcess; };

      /**
       * @return fraction of event to be kept in trimmed process expectation. */
      double GetTrimPercentage()
         { return fTrimPercentage; };
  
      /**
       * @return whether or not statistical uncertainties are considered. These are
       * decorrelated uncertainties bin-to-bin. */
/**      bool GetDoStatUnc()
	     { return fUseStatUnc; }; */
  
      /**
       * @return whether or not a normalisation uncertainty is considered. */
      bool GetDoNormalisationUnc()
         { return fUseNormUnc; };

      /** @} */
      /** \name Member functions (set) */
      /** @{ */

      /** 
       * Set the name of the process.
       * @param name The name of the process. */
      void SetName(std::string name)
         { fName = name; };

      /**
       * Store index of parameter that was created for this process. 
       * @param index The index to store. */
      void AddProcessParameterIndex(int index)
         { fParamIndices.push_back(index); };
  
      /**
       * Add index of normalisation parameter. */
      void AddNormalisationParIndex(int index)
         { fNormalisationParameter = index; };

      /**
       * Store name of parameter that was created for this process. 
       * @param name The name of the parameter. */
      void AddProcessParameterName(std::string name)
         { fParamNames.push_back(name); };

      /**
       * Set the nominal histogram and the +1-sigma histogram.
       * @param centralTemplate The nominal TH1D histogram. 
       * @param variationTemplate The TH1D histogram corresponding
       * to a +1 sigma variation. */
      void SetTemplateFromHistograms(TH1D * centralTemplate, TH1D * variationTemplate = 0);

      /** 
       * Set a function to define shape and scale of template.
       * The function will define a histogram whose first and last bins 
       * will be the first and last of blankHist which contain data.
       * @param functionTemplate The function.
       * @param blankHist A histogram with the bin structure desired. 
       * @param nominalFuncParams Nominal value for each param in fitted function.
       * @param errorFuncParams One-sigma error on each param in fitted function. 
       * @param luminosity The luminosity to normalise the process to. */
      void SetTemplateFromFunction(MjjFitFunction * functionTemplate, TH1D * blankHist, vector<double> nominalFuncParams, vector<double> errorFuncParams, double luminosity);

      /**
       * Set a shape-changing systematic variation for the given systematic and this process.
       * @param systematic Name of the corresponding systematic
       * @param systvar The shape-changing systematic variation */
      int SetSystematicShapeVariation(std::string systematic, MjjBATShapeChangingSyst * systvar);

      /**
       * Set a scale-changing systematic variation for the given systematic and this process.
       * @param systematic Name of the corresponding systematic
       * @param systvar The scale-changing systematic variation */
      void SetSystematicScaleVariation(std::string systematic, MjjBATScaleChangingSyst * systvar);

      /**
       * Set a template-based shape-changing systematic variation 
       * for the given systematic and this process.
       * @param systematic Name of the corresponding systematic
       * @param systvar The template-based systematic variation */
      int SetSystematicTemplateVariation(std::string systematic, MjjBATTemplateSyst * systvar);

      /**
       * If activated without specifying amount via SetTrimPercentage, uses default 0.95
       * @param doTrim Set true if the process expectation should be trimmed to central specified amount. */
      void SetTrimProcess(bool doTrim)
         { fTrimProcess = doTrim; };

      /**
       * Sets fraction of event to be kept in trimmed process expectation. 
       * @param percentage Desired fraction of events to keep. */
      void SetTrimPercentage(double percentage)
         {  fTrimProcess = true;
            fTrimPercentage = percentage; };

      /** @} */

 private:

      /**
        * The name of the process. */
      std::string fName;

      /**
        * Luminosity of the dataset, for function normalization. */
      double fLuminosity;

      /** 
       * The nominal TH1D histogram. */
      TH1D * fNominal;

      /** 
       * The +1sigma TH1D histogram. */
      TH1D * fVariation;

      /** 
       * The number of bins in the process definition. */
      int fNBinsInSpectrum;

      /** 
       * The bin content of the nominal histogram as a TVectorD. */
      TVectorD fNominalBins;

      /** 
       * The bin content of the variation histogram as a TVectorD. */
      TVectorD fVariationBins;

      /**
       * A TF1 alternative for defining template
       * shape and scale. */
      MjjFitFunction * fFunction;

      /**
       * Nominal values of function parameters */
      vector<double> fFuncParamsNominal;

      /**
       * 1-sigma uncertainty on parameters
       * where order is the same as in fFuncParamsNominal. */
      vector<double> fFuncParamsError;

      /**
       * A vector of names of systematics which affect this process. */
      vector<std::string> systematics;

      /**
       * A vector of shape-changing systematic variations on this process. */
      vector<MjjBATShapeChangingSyst*> shapeVariations;

      /**
       * A vector of scale-changing systematic variations on this process. */
      vector<MjjBATScaleChangingSyst*> scaleVariations;

      /**
       * Scale-changing variations coming from template systematics on this process. */
      vector<MjjBATTemplateSyst*> templateScaleVariations;

      /**
       * A single systematic variation based on templates. */
      MjjBATTemplateSyst * templateVariation;

      /**
       * Indices of parameters in model which correspond to this process. */
      vector<int> fParamIndices;
  
      /**
      * Index of the parameter corresponding to normalisation. Returns -1 if none. */
      int fNormalisationParameter;

      /**
       * Names of parameters in model which correspond to this process. */
      vector<std::string> fParamNames;

      /**
       * The number of parameters in BCModel which define this process. */
      int fNParams;

      /**
       * Do we trim the signal expectation to a central range? */
      bool fTrimProcess;

      /**
       * Set percentage to which to trim expectation */
      double fTrimPercentage;
  
      /**
       * Determines whether to perform a summation over an uncorrelated nuisance parameter
       * representing the statistical uncertainty in each bin. */
//      bool fUseStatUnc;
  
      /**
       * Determines whether to introduce one overall parameter for the normalisation of the 
       * sample. The size will be set from the variation histogram, so this can be used
       * for either correlated errors (as from a function) or for the normalisation of an MC
       * sample, as long as it is set correctly. */
      bool fUseNormUnc;
  
};
// ---------------------------------------------------------

#endif

