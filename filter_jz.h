/***** filter.h *****/
/**
 * Includes all the necessary code to implement a second order IIR filter
 */

#include <cmath>

#define FILTER_TYPE_LOW_PASS 0
#define FILTER_TYPE_BAND_PASS 1
#define FILTER_TYPE_HIGH_PASS 2

/**
 * Second order filter of variable type.
 */
class Filter {
	public:
		//Constructors
		Filter();
		Filter(int filterType_, float cutoffFrequency_, float T);
		~Filter();
		
		//Public functions
		float processNew(float x_);
		
	private:
		int filterType = FILTER_TYPE_LOW_PASS;
		float cutoffFrequency = 440 * 2 * M_PI;
		
		//Filter coefficients
		float a_0, a_1, a_2, b_0, b_1, b_2;
		
		//Previous filter state;
		float x0=0;
		float x1=0;
		float x2=0;
		float y0=0;
		float y1=0;
		float y2=0;
		
		void calculate_coefficients(float T);
};
