#include "adaptivelikelihoods.h"
#include <algorithm>
#include <cmath>
#include <assert.h>
#include <numeric>


using std::conj;


/* Print statement */
void RelativeBinningPrinter(std::string message)
{
    std::cout << "RELATIVE BINNING:\t";
    std::cout << message << "\n";
}

/**
 * Initialize relative binning
*/
RelativeBinningLikelihood::RelativeBinningLikelihood(
    double chi, /**< size of PN perturbation */
    double epsilon, /**< size of dephasing */
    const std::vector<ifo_data_struct> &ifos_data, /**< Vector of signal data */
    const cpl *const *fiducial_data, /**< Strain of fiducial data for each interferometer */
    int num_detectors, /**< Number of detectors */
    const double *frequencies, /**< Frequency array shared among all interferometers */
    int data_length, /**< length of each fiducial array */
    const vectordbl gammas /**< Vector of gammas for binning criterion */
)
{
    std::cout << "RELATIVE BINNING INITIALIZING\n";
    this->chi = chi;
    this->epsilon = epsilon;
    this->gammas = gammas;

    std::cout << "\tUsing gamma = { ";
    for (auto &gamma : gammas)
    {
        std::cout << gamma << " ";
    }
    std::cout << "}\n";

    duration = ifos_data[0].duration;
    max_frequency = find_max_frequency(
        fiducial_data, frequencies, data_length, num_detectors
    );
    RelativeBinningPrinter("Max frequency: " + std::to_string(max_frequency));

    setup_bins(frequencies, data_length);
    RelativeBinningPrinter(std::to_string(number_of_bins) + " bins setup");

    setup_fiducial_data(fiducial_data, num_detectors);
    RelativeBinningPrinter("Fiducial data setup");

    compute_summary_data(fiducial_data, ifos_data);
    RelativeBinningPrinter("Summary data setup");

    // Calculate likelihood
    double logL = 0.;
    cpl *data;
    for (int i = 0; i < num_detectors; i++)
    {
        data = ifos_fiducial_data.at(i).strain.data();
        logL += log_likelihood_per_detector(
            data, &(ifos_fiducial_data.at(i))
        );
    }
    RelativeBinningPrinter("Fiducial likelihood: " + std::to_string(logL));
}

/**
 * Find maximum frequency at which we can evaluate the likelihood,
 * which will determine the binning.
 * Defined as the last frequency where the fiducial data is non-zero
*/
double RelativeBinningLikelihood::find_max_frequency(
    const cpl *const *fiducial_data_in,
    const double *frequencies, const int data_length,
    const int num_detectors
)
{
    vectorcpl fiducial_data;
    vectorcpl::reverse_iterator ref;
    int idx;
    double maxFreq = frequencies[data_length-1];
    
    for (int d = 0; d < num_detectors; d++)
    {
        fiducial_data = vectorcpl(
            fiducial_data_in[d], fiducial_data_in[d] + data_length
        );

        ref = std::find_if(fiducial_data.rbegin(), fiducial_data.rend(),
            [](cpl &dat) {return dat != cpl(0.);});
        if (ref == fiducial_data.rbegin())
        {
            // If the last data point is non-zero
            continue;
        }
        
        idx = std::distance(ref, fiducial_data.rend())-1;
        maxFreq = fmin(frequencies[idx], maxFreq);
    }

    return maxFreq;
}

/**
 * Setup frequency bins based on PN ansatz.
 * See Binning Criterion in arxiv:2312.06009
*/
void RelativeBinningLikelihood::setup_bins(
    const double *frequencies, /**< Frequency array */
    const int data_length /**< Size of frequency array */
)
{
    double min_frequency;
    double perturb_size = GWAT_TWOPI * chi;

    // Find the frequency range over all detectors
    // Assumes array is in sequential order
    min_frequency = frequencies[0];
    
    // Eqs. 15-16 of 2312.06009
    // To optimize computing the sum in Eq. 16 we do not solve for \Delta_\alpha
    // but rather find the f_{k*} for each k
    vectordbl freqs_gamma;
    vectorint signs_gamma;
    double freq_gamma;
    bool negative_gamma;
    for (const double &gamma : gammas)
    {
        negative_gamma = gamma < 0;
        freq_gamma = negative_gamma ? min_frequency : max_frequency;
        freqs_gamma.push_back(freq_gamma);
        signs_gamma.push_back(negative_gamma ? -1 : 1);
    }

    // Eq. 16, second line instead of first
    vectordbl d_phis;
    double d_phi_f, f;
    size_t k;
    for (int i = 0; i < data_length; i++)
    {
        f = frequencies[i];
        if (f > max_frequency) { break; }

        d_phi_f = 0.;

        for (k = 0; k < gammas.size(); k++)
        {
            d_phi_f += signs_gamma[k] * pow(
                f/freqs_gamma[k], gammas[k]
            );
        }

        d_phis.push_back(perturb_size*d_phi_f);
    }

    // \Delta\Psi(f) - \Delta\Psi(f_{\min})
    vectordbl d_phi_from_start;
    for (const double &d_phi : d_phis)
    {
        d_phi_from_start.push_back(d_phi-d_phis[0]);
    }

    // Number of bins
    int num_bins = static_cast<int> (
        std::floor(d_phi_from_start.back()/epsilon)
    );

    // Find bin edges
    vectordbl::iterator bin_itr;
    // last_* variables to avoid starting the find_if searches from the beginning
    const double *last_base_ptr = frequencies;
    vectordbl::iterator last_itr = d_phi_from_start.begin();
    int bin_ind, last_ind = -1;
    double d_phi_lower, bin_freq;
    double d_phi_lower_den = d_phi_from_start.back() / (double)num_bins;
	// Find in integer increments of epsilon where a bin will have reached a 
	// dephasing on the order of epsilon
    for (int i = 0; i < num_bins+1; i++)
    {
        d_phi_lower = i * d_phi_lower_den;
        
        bin_itr = std::find_if(last_itr, d_phi_from_start.end(),
            [d_phi_lower](double &d_phi) {return d_phi >= d_phi_lower;});
        bin_ind = std::distance(d_phi_from_start.begin(), bin_itr);
        if (bin_ind == last_ind)
        {
            continue;
        }

        last_itr = bin_itr;
        last_ind = bin_ind;
        bin_freq = frequencies[bin_ind];
        last_base_ptr = std::find_if(last_base_ptr,
            frequencies + data_length,
            [bin_freq](double f) {return f >= bin_freq;});
        bin_ind = std::distance(frequencies, last_base_ptr);

        bin_inds.push_back(bin_ind);
        bin_freqs.push_back(bin_freq);
    }

    // Set bin info
    number_of_bins = bin_inds.size()-1;
    for (int i = 1; i < bin_inds.size(); i++)
    {
        bin_sizes.push_back(bin_inds[i]-bin_inds[i-1]);
        bin_widths.push_back(bin_freqs[i]-bin_freqs[i-1]);
        bin_centers.push_back((bin_freqs[i]+bin_freqs[i-1])/2);
    }

    bins_are_setup = true;
}

/**
 * Store fiducial waveforms at bin edges
 * Saved within ifos_fiducial_data
*/
void RelativeBinningLikelihood::setup_fiducial_data(
    const cpl *const *fiducial_data_in, /**< Fiducial data of each interferometer */
    const int num_detectors /** Number of detectors */
)
{
    assert(bins_are_setup);

    const cpl *det_wf;
    int d;

    for (d = 0; d < num_detectors; d++)
    {
        fiducial_data_struct ifo;
        det_wf = fiducial_data_in[d];

        for (int &ind : bin_inds)
        {
            ifo.strain.push_back(det_wf[ind]);
        }
        
        ifos_fiducial_data.push_back(ifo);
    }
}

/**
 * Compute the A0, A1, B0, and B2 summary data for each bin in each interferometer.
 * See Eq. 5 of 2312.06009
*/
void RelativeBinningLikelihood::compute_summary_data(
    const cpl *const *fiducial_data,
    const std::vector<ifo_data_struct> &ifos_data
)
{
    const ifo_data_struct *ifo;
    const cpl *fiducial;
    fiducial_data_struct* ifo_fiducial;
    vectorint ifo_bin_ends;
    vectorcpl data, h0;
    vectordbl psd, freqs;
    int start_ind, end_ind;
    cpl A_fac, B_fac, A0_b, A1_b, B0_b, B1_b;
    double delta_f, inner_product_weight;

    // Compute summary data for each interferometer
    for (int i = 0; i < ifos_data.size(); i++)
    {
        // The ifo holding the data
        ifo = &(ifos_data.at(i));
        // The fiducial data
        fiducial = fiducial_data[i];
        // The ifo to store the summary data
        ifo_fiducial = &ifos_fiducial_data.at(i);

        inner_product_weight = 4./ifo->duration;
        
        ifo_bin_ends = vectorint(bin_inds);
        // Cover the last frequency in the bin array
        ifo_bin_ends.end() += 1;

        // Compute summary data per bin
        for (int b = 0; b < number_of_bins; b++)
        {
            // Get the starting and ending indices of the bin
            start_ind = ifo_bin_ends[b];
            end_ind = ifo_bin_ends[b+1];

            // Grab the data for each bin
            data = vectorcpl(ifo->strain.begin() + start_ind,
                ifo->strain.begin() + end_ind);
            psd = vectordbl(ifo->psd.begin() + start_ind,
                ifo->psd.begin() + end_ind);
            h0 = vectorcpl(fiducial + start_ind,
                fiducial + end_ind);
            freqs = vectordbl(ifo->freqs.begin() + start_ind,
                ifo->freqs.begin() + end_ind);

            // Compute the Riemann sums
            A0_b = A1_b = B0_b = B1_b = 0.;
            for (int f = 0; f < freqs.size(); f++)
            {
                A_fac = conj(h0[f]) / psd[f];
                B_fac = h0[f] * A_fac;
                A_fac *= data[f];
                delta_f = freqs[f] - bin_centers[b];

                A0_b += A_fac;
                A1_b += A_fac*delta_f;
                B0_b += B_fac;
                B1_b += B_fac*delta_f;
            }
            // Store the data
            ifo_fiducial->A0.push_back(inner_product_weight*A0_b);
            ifo_fiducial->A1.push_back(inner_product_weight*A1_b);
            ifo_fiducial->B0.push_back(inner_product_weight*B0_b);
            ifo_fiducial->B1.push_back(inner_product_weight*B1_b);
        }
    }
}

void RelativeBinningLikelihood::compute_waveform_ratios(
    vectorcpl &r0, /**< [out] Array of r0 coefficients in each bin */
    vectorcpl &r1, /**< [out] Array of r1 coefficients in each bin */
    const cpl *h, /**< Template waveform evaluated at bin edges */
    const fiducial_data_struct *fiducial /**< Interferometer data */
)
{
    // Ratios at left edge
    cpl ratio_left = h[0]/fiducial->strain.front();
    // Ratio at right edge
    cpl ratio_right;

    for (int i = 0; i < number_of_bins; i++)
    {
        ratio_right = h[i+1]/fiducial->strain[i+1];

        r0.push_back( 0.5*(ratio_right + ratio_left) );
        r1.push_back( (ratio_right - ratio_left)/bin_widths[i] );

        // For next bin
        ratio_left = ratio_right;
    }
}

double RelativeBinningLikelihood::log_likelihood_per_detector(
    const cpl *h, /**< Template waveform. Must be evaluated at the bin edges only */
    const fiducial_data_struct *fiducial /**< Interferometer data */
)
{
    double d_h = 0.;
    double h_h = 0.;

    // Obtain the waveform ratios
    vectorcpl r0, r1;
    compute_waveform_ratios(r0, r1, h, fiducial);

    for (int b = 0; b < number_of_bins; b++)
    {
        // Inner product of the data with the template
        d_h += real(
            fiducial->A0[b] * conj(r0[b]) + fiducial->A1[b] * conj(r1[b])
        );
        // Inner product of the template
        h_h += real(
            fiducial->B0[b] * r0[b]*conj(r0[b]) + 
            2.*fiducial->B1[b] * real(r0[b]*conj(r1[b]))
        );
    }

    return -0.5*h_h + d_h;
}

double RelativeBinningLikelihood::log_likelihood(
    std::string *detectors, /**< Detector names */
    int num_detectors, /**< Number of detectors */
    gen_params_base<double> *params, /**< Template parameters */
    std::string generation_method, /**< Template model name */
    bool reuse_WF /**< Option to obtain detector responses using the same waveform */
)
{
    // Set up template arrays a là MCMC_likelihood_extrinsic (mcmc_gw.cpp)
    int i;
    int *data_lengths = new int[num_detectors];
    cpl **responses = new cpl*[num_detectors];
    double **frequencies = new double*[num_detectors];

    for (int i = 0; i < num_detectors; i++)
    {
        data_lengths[i] = bin_freqs.size();
        responses[i] = new cpl[data_lengths[i]];
        frequencies[i] = bin_freqs.data();
    }

    // Time shift
    // TODO: How important is shifting tc as in MCMC_likelihood_extrinsic? If not, delete this block.
    //double T = duration;
    //double tc_ref = T - params->tc;
    //params->tc = tc_ref;

    // Obtain the data
    if (num_detectors == 1)
    {
        create_single_GW_detection(
            responses[0],
            detectors[0],
            frequencies[0],
            data_lengths[0],
            params,
            generation_method
        );
    }
    else
    {
        create_coherent_GW_detection(
            detectors, num_detectors, frequencies, data_lengths, reuse_WF,
            params, generation_method, responses
        );
    }

    // Calculate likelihood
    double logL = 0.;
    for (i = 0; i < num_detectors; i++)
    {
        logL += log_likelihood_per_detector(
            responses[i], &(ifos_fiducial_data[i])
        );
    }

    // Clean up
    for (i = 0; i < num_detectors; i++)
    {
        delete [] responses[i];
    }
    delete [] frequencies;
    delete [] responses;
    delete [] data_lengths;

    return logL;
}
