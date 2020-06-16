/***** fft_jz.h *****/
#include <complex_jz.h>
#include <cmath>
#include <Bela.h>

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

	void do_butterflies_asm_precalculated_twiddles(Complex* fftbuffer, Complex* fftbuffer_interleaved, Complex* twiddleArrayForStage);
	void do_butterflies_neon_precalculated_twiddles(Complex* fftbuffer, Complex* fftbuffer_interleaved, Complex* twiddleArrayForStage);
}

enum FFT_TYPE { Cpp, Asm, AsmNeon };

class MyFFT {
		private:
		//Variables
		int mlog2N;
		Complex **mTwiddleArray;
		FFT_TYPE mfft_type = Cpp;
		
		//Methods
		Complex calculateTwiddleRotationFactor(int s);
		
		void calculateTwiddleArray();
		void de_interleave(Complex *interleaved, Complex *de_interleaved, int N_);
		void butterfly(Complex *inputArray, int inputStartPointer, Complex *outputArray, Complex twiddle_1);
	
	public:
		//Constructors
		MyFFT(int N_);
		~MyFFT();
		
		//Variables
		int mN;
		
		//Methods
		void processFft(float *input, float *output);
		
		FFT_TYPE getmFftType() { return mfft_type; };
		void setmFftType(FFT_TYPE type) {mfft_type = type;};
};
