#ifndef ADAPTIVE_LIKELIHOODS_H
#define ADAPTIVE_LIKELIHOODS_H


#include "waveform_util.h"
#include <vector>


using cpl = std::complex<double>;
using vectordbl = std::vector<double>;
using vectorint = std::vector<int>;
using vectorcpl = std::vector<cpl>;

const vectordbl gammasDefault = {-5./3., -2./3., 1., 5./3., 7./3.};

// Struct for an interferometer's strain and PSD at specific frequencies
struct ifo_data_struct
{
	// Interferometer strain
	vectorcpl strain;
	// Frequency array
	vectordbl freqs;
	// PSD
	vectordbl psd;
	// Signal duration
	double duration;
};

// Struct for each interferometer.
// Holds their binned strain and summary data.
struct fiducial_data_struct
{
	vectorcpl strain;
	// Summary data
	vectorcpl A0, A1, B0, B1;
};


/** 
 * Base class for adaptive likelihoods
*/
class AdaptiveLikelihood
{
public:
    AdaptiveLikelihood() = default;
    virtual ~AdaptiveLikelihood() = default;

	// Compute the log likelihood -1/2[(h|h) - 2(d|h)] over all detectors
	virtual double log_likelihood(
		std::string *detectors, int num_detectors,
		gen_params_base<double> *params, std::string generation_method,
		bool reuse_WF
	) = 0;
};


// RELATIVE BINNING

class RelativeBinningLikelihood: public AdaptiveLikelihood
{
public:
	RelativeBinningLikelihood(
		double chi, double epsilon,
		const std::vector<ifo_data_struct> &ifos_data,
		const cpl *const *fiducial_data, int num_detectors,
		const double *frequencies, int data_length,
		const vectordbl gammas = gammasDefault
	);

	double log_likelihood(
		std::string *detectors, int num_detectors,
		gen_params_base<double> *params, std::string generation_method,
		bool reuse_WF
	);

	vectordbl get_bin_freqs() { return bin_freqs; }

private:
	double chi, epsilon;
	// Signal duration for tc purposes
	double duration;
	int number_of_bins;
	// Check that bins have been setup
	bool bins_are_setup = false;
	// Interferometers
	std::vector<fiducial_data_struct> ifos_fiducial_data;

	// Maximum frequency 
	double max_frequency;
	// Array of frequency powers for binning criterion
	vectordbl gammas;
	// Frequency array bin edge indices
	vectorint bin_inds;
	// Frequencies at bin edges
	vectordbl bin_freqs;
	// Distances between bin indices
	vectordbl bin_sizes;
	// Distances between edge frequencies
	vectordbl bin_widths;
	// Frequency at bin centers
	vectordbl bin_centers;

	double find_max_frequency(const cpl *const *fiducial_data_in,
		const double *frequencies, const int data_length,
		const int num_detectors);
	void setup_bins(const double *frequencies, const int data_length);
	void setup_fiducial_data(
		const cpl *const *fiducial_data_in,
		const int num_detectors
	);
	void compute_summary_data(
		const cpl *const *fiducial_data,
		const std::vector<ifo_data_struct> &ifos_data
	);


	void compute_waveform_ratios(
		vectorcpl &r0,
		vectorcpl &r1,
		const cpl *h,
		const fiducial_data_struct *fiducial
	);
	double log_likelihood_per_detector(
		const cpl *h,
		const fiducial_data_struct *fiducial
	);
};

#endif // ADAPTIVE_LIKELIHOODS_H