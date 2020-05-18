/***** fft_jz.h *****/
#include <complex_jz.h>

/* Assembly-language filtering function */
extern "C" {
	// Fast bit reversal using rotating registers.
	// Only using 32 bit registers so second argument's maximum value is 32 if using unsigned ints, 31 if signed.
	int bitReversal(int intToReverse, int numOfBitsToReverse);
	
	// Assembly implementation of all stages of butterflies in the FFT.
	void do_butterflies_asm(Complex* fftbuffer, Complex* fftbuffer_interleaved, Complex* twiddleRotationFactor, int stage);
	
	// Assembly implementations of the first log2n-2 stages of butterflies in the FFT.
	void do_butterflies_neon(Complex* fftbuffer, Complex* fftbuffer_interleaved, Complex* twiddleRotationFactor, int stage);
	void do_butterflies_neon_pipeline(Complex* fftbuffer, Complex* fftbuffer_interleaved, Complex* twiddleRotationFactor, int stage);
}

/**
 * At each stage of the FFT, the twiddles rotate by a certain factor
 * Stage 1: pi/2, Stage 2: pi/4, Stage 3: pi/8
 * Angle = -pi/(2^stage)
 */
Complex calculateTwiddleRotationFactor(int stage) {
	Complex rotationFactor;
	float angle = -M_PI/(1<<stage);
	rotationFactor.re = cos(angle);
	rotationFactor.im = sin(angle);
	return rotationFactor;
}

/**
 * This is the smallest building block of the DFT
 */
void butterfly(Complex *inputArray, int inputStartPointer, Complex *outputArray, Complex twiddle_1) {
	Complex input_1 = inputArray[inputStartPointer++];
	Complex input_2 = inputArray[inputStartPointer--];
	
	input_2 = complexMultiply(input_2, twiddle_1);
	
	outputArray[inputStartPointer++] = complexAdd(input_1, input_2);
	outputArray[inputStartPointer] = complexSubtract(input_1, input_2);
}

/**
 * De-Interleave two arrays of Complex numbers
 */
void de_interleave(Complex *interleaved, Complex *de_interleaved, int N) {
	for (int i=0; i<N/2; i++) {
		de_interleaved[i] = interleaved[2*i];
	}
	for (int i=N/2; i<N; i++) {
		de_interleaved[i] = interleaved[(2*i)-N+1];
	}
}


/**
 * My FFT implementation (just returning the magnitude as that's all I'm interested in for spectral visualisation)
 * input is an array of floats
 * output is an array of floats that is written to
 * N is the size of the FFT (a power of 2)
 * log2N is log2(N)
 */
void my_fft(float *input, int N, int log2N, float *output) {
	// Copy input into a new (Complex) array in bit reversed order.
	Complex *fftTempBuffer = new Complex[N];
	Complex *fftTempBuffer_interleaved = new Complex[N];
	for (unsigned int i=0; i<N; i++) {
		Complex c;
		c.re = input[bitReversal(i, log2N)];
		//c.im = (float)i;
		fftTempBuffer[i] = c;
	}
	
	// For each stage... (log2(N) number of stages), we perform N/2 butterflies down the array.
	for (int stage=0; stage<log2N; stage++) {
		
		// Initialises to a complex value of (1 + j0), which is the initial twiddle factor. Handy dandy.
		Complex twiddleFactor;
		
		// Each twiddle factor rotates by a certain amount for each stage. (Zeroth: No rotation, First stage is pi/2, second stage is pi/4, pi/8 etc...)
		Complex twiddleRotationFactor;
		if (stage!=0) {
			twiddleRotationFactor = calculateTwiddleRotationFactor(stage);
		}
		
		//----------------------------------------//
		// Perform N/2 butterflies down the array //
		//----------------------------------------//
		// int butterfly_input_pointer = 0;
		// //There's 2^stage different types of butterfly in each stage
		// //rt_printf("Types of butterfly: %d\n",(1<<stage));
		// //rt_printf("Number of each butterfly: %d\n", N/(2<<stage));
		// for (int butterflyType=0; butterflyType < (1<<stage); butterflyType++) {
			
		// 	//There's N/(2^(stage+1)) of each butterfly type.
		// 	for (int butterflyNum=0; butterflyNum < N/(2<<stage); butterflyNum++) {
		// 		butterfly(fftTempBuffer, butterfly_input_pointer, fftTempBuffer_interleaved, twiddleFactor);
		// 		//rt_printf("FftTempInterleaved: {%.3f+j%.3f}, {%.3f+j%.3f}, {%.3f+j%.3f}, {%.3f+j%.3f}\n", fftTempBuffer_interleaved[0].re, fftTempBuffer_interleaved[0].im, fftTempBuffer_interleaved[1].re, fftTempBuffer_interleaved[1].im,
		// 		//fftTempBuffer_interleaved[2].re, fftTempBuffer_interleaved[2].im, fftTempBuffer_interleaved[3].re, fftTempBuffer_interleaved[0].im);
		// 		butterfly_input_pointer++;
		// 		butterfly_input_pointer++;
		// 		//rt_printf("%d\n", butterfly_input_pointer);
		// 	}
			
		// 	twiddleFactor = complexMultiply(twiddleFactor, twiddleRotationFactor);
		// }
		
		//----------------------------------------//
		// Perform N/2 butterflies using assembly //
		//----------------------------------------//
		// do_butterflies_asm(fftTempBuffer, fftTempBuffer_interleaved, &twiddleRotationFactor, stage);
		
		//-------------------------//
		// Using the NEON unit now //
		//-------------------------//
		if (stage < log2N-2) {
			do_butterflies_neon(fftTempBuffer, fftTempBuffer_interleaved, &twiddleRotationFactor, stage);
			// The 'pipelining' I tried makes no difference in the end.
			// do_butterflies_neon_pipeline(fftTempBuffer, fftTempBuffer_interleaved, &twiddleRotationFactor, stage);
		} else {
			do_butterflies_asm(fftTempBuffer, fftTempBuffer_interleaved, &twiddleRotationFactor, stage);
		}
		
		//De-interleave the array from fftTempBuffer_interleaved into fftTempBuffer
		de_interleave(fftTempBuffer_interleaved, fftTempBuffer, N);
		
	}
	
	
	//Put the magnitude into the output.
	for (int i=0; i<N; i++) {
		float bin = sqrt((fftTempBuffer[i].re*fftTempBuffer[i].re) + (fftTempBuffer[i].im*fftTempBuffer[i].im))/(N/4.0); //Dividing by N/4 to normalise amplitude to 1.
		output[i] = bin;
	}
	
	// Free up these buffers.
	delete[] fftTempBuffer;
	delete[] fftTempBuffer_interleaved;
}
