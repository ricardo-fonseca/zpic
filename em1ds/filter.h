#ifndef __FILTER__
#define __FILTER__

/**
 * @brief Spectral filter types
 * 
 */
enum filter_type {
	FILTER_NONE,	///< No filtering
	FILTER_GAUSS,	///< Gaussian filtering
	FILTER_SHARP	///< Sharp low-pass filtering
};

/**
 * @brief Spectral filter
 * 
 */
typedef struct {
	/// Type of filtering
	enum filter_type type;

	/// Number of points
	int nk;

	/// Filter response datqa
	float *Sk;
} t_filter;

/**
 * @brief Initializes filter object
 * 
 * @param filter Filter object
 * @param nk     Number of k-space points (assumes k at nk-1 is Nyquist wavenumber)
 */
void filter_new( t_filter *filter, int nk );

/**
 * @brief Cleanup filter object
 * 
 * @param filter Filter object
 */
void filter_delete( t_filter *filter );

/**
 * @brief Updates filter values for the specified parameters
 * 
 * @param filter Filter object
 * @param type   Type of filter to use
 * @param ck     Cutoff wavenumber for the filter (filter response = 1/2)
 * @return       0 on success, -1 on error
 */
int filter_set( t_filter* const filter, enum filter_type const type, float const ck );

#endif
