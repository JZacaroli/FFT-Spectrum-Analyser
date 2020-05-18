/***** filter.h *****/
/**
 * Includes all the necessary code to implement a second order IIR filter
 */

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

Filter::Filter() {
	//Do nothing.
}

Filter::Filter(int filterType_, float cutoffFrequency_, float T) {
	filterType = filterType_;
	cutoffFrequency = cutoffFrequency_;
	calculate_coefficients(T);
	
}

// Calculate filter coefficients given a T. (Usually 1/sampleRate)
void Filter::calculate_coefficients(float T)
{
	//Do the pre-warping here.
	float w = (2/T)*atan(M_PI*cutoffFrequency*T);
	
	// Could make this variable in the future?
	float q = 1/sqrt(2);
	
	switch(filterType) {
		case FILTER_TYPE_LOW_PASS:
			//LPF coefficients:
			a_0 = 1 + (2/(w*q*T)) + (4/(w*w*T*T));
			a_1 = 2 - (8/(w*w*T*T));
			a_2 = 1 - (2/(w*q*T)) + (4/(w*w*T*T));
			b_0 = 1;
			b_1 = 2;
			b_2 = 1;
			break;
		case FILTER_TYPE_BAND_PASS:
			//BPF coefficients:
			a_0 = (4/(T*T)) + ((2*w)/(q*T)) + (w*w);
			a_1 = (2*w*w) - (8/(T*T));
			a_2 = (4/(T*T)) - ((2*w)/(q*T)) + (w*w);
			b_0 = (2*w)/(q*T);
			b_1 = 0;
			b_2 = -(2*w)/(q*T);
			break;
		case FILTER_TYPE_HIGH_PASS:
			//HPF coefficients:
			a_0 = (4/(T*T)) + ((2*w)/(q*T)) + (w*w);
			a_1 = (2*w*w) - (8/(T*T));
			a_2 = (4/(T*T)) - ((2*w)/(q*T)) + (w*w);
			b_0 = 4/(T*T);
			b_1 = -8/(T*T);
			b_2 = 4/(T*T);
			break;
	}
}

float Filter::processNew(float x_) {
	// Process new input
	x0 = x_;
	y0 = ((x0*b_0 + x1*b_1 + x2*b_2) - (y1*a_1 + y2*a_2))/a_0;
	
	// Update the filter state for the next time we come through processNew
	x2 = x1;
	x1 = x0;
	y2 = y1;
	y1 = y0;
	
	return y0;
}

Filter::~Filter() {
	//Do nothing here.
}