/***** filter.cpp *****/
#include <filter_jz.h>

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