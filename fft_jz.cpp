/***** fft_jz.cpp *****/
#include <fft_jz.h>

MyFFT::MyFFT(int N_) {
	mN = N_;
	mlog2N = log2(N_);
	
	// Calculate the array of twiddle factors used in the FFT.
	calculateTwiddleArray();
}

MyFFT::~MyFFT() {
	for (int i=0; i<mlog2N; i++) {
		delete[] mTwiddleArray[i];
	}
    delete[] mTwiddleArray;
}

/**
 * At each stage of the FFT, the twiddles rotate by a certain factor
 * Stage 1: pi/2, Stage 2: pi/4, Stage 3: pi/8
 * Angle = -pi/(2^stage)
 */
Complex MyFFT::calculateTwiddleRotationFactor(int stage) {
	Complex rotationFactor;
	float angle = -M_PI/(1<<stage);
	rotationFactor.re = cos(angle);
	rotationFactor.im = sin(angle);
	return rotationFactor;
}

/**
 * This is the smallest building block of the DFT
 */
void MyFFT::butterfly(Complex *inputArray, int inputStartPointer, Complex *outputArray, Complex twiddle_1) {
	Complex input_1 = inputArray[inputStartPointer++];
	Complex input_2 = inputArray[inputStartPointer--];
	
	input_2 = complexMultiply(input_2, twiddle_1);
	
	outputArray[inputStartPointer++] = complexAdd(input_1, input_2);
	outputArray[inputStartPointer] = complexSubtract(input_1, input_2);
}

/**
 * De-Interleave two arrays of Complex numbers
 */
void MyFFT::de_interleave(Complex *interleaved, Complex *de_interleaved, int N_) {
	for (int i=0; i<N_/2; i++) {
		de_interleaved[i] = interleaved[2*i];
	}
	for (int i=N_/2; i<N_; i++) {
		de_interleaved[i] = interleaved[(2*i)-N_+1];
	}
}


/**
 * My FFT implementation (just returning the magnitude as that's all I'm interested in for spectral visualisation)
 * input is an array of floats
 * output is an array of floats that is written to
 * N is the size of the FFT (a power of 2)
 */
void MyFFT::processFft(float *input, float *output) {
	
	// Copy input into a new (Complex) array in bit reversed order.
	Complex *fftTempBuffer = new Complex[mN];
	Complex *fftTempBuffer_interleaved = new Complex[mN];
	for (unsigned int i=0; i<mN; i++) {
		Complex c;
		c.re = input[bitReversal(i, mlog2N)];
		fftTempBuffer[i] = c;
	}
	
	// For each stage... (log2(N) number of stages), we perform N/2 butterflies down the array.
	for (int stage=0; stage<mlog2N; stage++) {
		switch (mfft_type) {
			case Cpp:
			{
				//---------------------------------------------------------------------//
				// Perform N/2 butterflies down the array with pre-calculated Twiddles //
				//---------------------------------------------------------------------//
				int butterfly_input_pointer = 0;
				for (int n=0; n<mN/2; n++) {
					Complex twiddles = mTwiddleArray[stage][n];
					butterfly(fftTempBuffer, butterfly_input_pointer, fftTempBuffer_interleaved, twiddles);
					butterfly_input_pointer++;
					butterfly_input_pointer++;
				}
				break;
			}
			case Asm:
			{
				do_butterflies_asm_precalculated_twiddles(fftTempBuffer, fftTempBuffer_interleaved, mTwiddleArray[stage]);
				break;
			}
			case AsmNeon:
			{
				do_butterflies_neon_precalculated_twiddles(fftTempBuffer, fftTempBuffer_interleaved, mTwiddleArray[stage]);
				break;
			}
		}
		
		//De-interleave the array from fftTempBuffer_interleaved into fftTempBuffer
		de_interleave(fftTempBuffer_interleaved, fftTempBuffer, mN);
	
	}
	
	
	//Put the magnitude into the output.
	for (int i=0; i<mN; i++) {
		float bin = sqrt((fftTempBuffer[i].re*fftTempBuffer[i].re) + (fftTempBuffer[i].im*fftTempBuffer[i].im))/(mN/4.0); //Dividing by N/4 to normalise amplitude to 1.
		output[i] = bin;
	}
	
	// Free up these buffers.
	delete[] fftTempBuffer;
	delete[] fftTempBuffer_interleaved;
}

void MyFFT::calculateTwiddleArray() {
	mTwiddleArray = new Complex*[mlog2N];
	
	for (int stage=0; stage<mlog2N; stage++) {
		
		int twiddleIndexForStage=0;
		Complex twiddleFactor;
		Complex twiddleRotationFactor = calculateTwiddleRotationFactor(stage);
		mTwiddleArray[stage] = new Complex[mN/2];
		
		for (int butterflyType=0; butterflyType<1<<stage; butterflyType++) {
			
			for (int butterflyNum=0; butterflyNum < mN/(2<<stage); butterflyNum++) {
				mTwiddleArray[stage][twiddleIndexForStage] = twiddleFactor;
				twiddleIndexForStage++;
			}
			
			twiddleFactor = complexMultiply(twiddleFactor, twiddleRotationFactor);
		}
	}
}
