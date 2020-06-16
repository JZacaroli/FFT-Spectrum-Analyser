// ECS7012P Music and Audio Programming
// School of Electronic Engineering and Computer Science
// Queen Mary University of London
// Spring 2020
// Final Project - Joe Zacaroli

/**
 *			************************README***************************************
 *			*You can play around with gFFTSize and gHopSize if you fancy.       *
 *			*Or swap in different wav files to analyze/listen to with gFileName.*
 *			*********************************************************************
 */

#include <Bela.h>
#include <SampleLoader.h>
#include <libraries/Gui/Gui.h>
#include <stdlib.h>
#include <cmath>

#include <chrono>
#include <thread>

//NE10's implementation STUFF
#include <libraries/ne10/NE10.h>
ne10_fft_cpx_float32_t* timeDomainIn;
ne10_fft_cpx_float32_t* frequencyDomain;
ne10_fft_cfg_float32_t cfg;
bool TESTING_NE10 = false;

//My include files.
#include <fft_jz.h> //fft_jz includes complex_jz.h so we have access to complex numbers in this file too.
#include <filter_jz.h>
#include <complex_jz.h>

// GUI object declaration
Gui gui;

// Determines the size of the FFT processed.
int gFFTSize = 32768;
// Determines how long we should wait between processing FFTs.
int gHopSize = gFFTSize/10;
//How often we send to the GUI. (Every SEND_WAIT_TIME FFTs)
int SEND_WAIT_TIME=1;
int sendCounter=0;


// Input circular buffer. Make sure this is bigger than FFT size.
#define BUFFER_SIZE 35000
float gFftInputBuffer[BUFFER_SIZE] = {0};
int gFftInputBufferPointer = 0;

MyFFT *fftObject;

// Output frequency linear Buffer
float *gFftOutputBuffer;

// Hop counter use to determine when we should start the next FFT.
int gHopCounter = 0;

// Buffer for holding a precalculated window function
// Allocated dynamically in setup()
float *gWindowBuffer;

// Used to rt_printf() when we're not keeping up with real-time.
int gLastNumberOfFftsCurrentlyProcessing=0;
int gNumberOfFftsCurrentlyProcessing=0;


// Sample info
// string gFilename = "Chopin.wav"; // Name of the sound file (in project folder)
string gFilename = "dubstep.wav";
// string gFilename = "2000hz.wav";
// string gFilename = "DrumBeat.wav";
float *gSampleBuffer;			 // Buffer that holds the sound file
int gSampleBufferLength;		 // The length of the buffer in frames
int gReadPointer = 0;			 // Position of the last frame we played 

// Thread for FFT processing
AuxiliaryTask gFFTTask;
int gCachedFftInputBufferPointer = 0; //Cache the pointer into the input circular buffer so the background thread can start at the right place.
void process_fft_background(void *);

//Audio Processing Variables
Filter LowsLowPass_0, LowsLowPass_1, LowsLowPass_2, LowsLowPass_3, MidsHighPass_0, MidsHighPass_1, MidsHighPass_2, MidsHighPass_3;
Filter MidsLowPass_0, MidsLowPass_1, MidsLowPass_2, MidsLowPass_3, HighsHighPass_0, HighsHighPass_1, HighsHighPass_2, HighsHighPass_3;
bool filtersOn = false;
float gains[3] = {1};

bool setup(BelaContext *context, void *userData)
{
	//--------------------//
	// Setup Filter state //
	//--------------------//
	// Edges of the low / mid / high bands
	float Fc1 = 200;
	float Fc2 = 2000;
	float T = 1/(context->audioSampleRate);
	LowsLowPass_0 = Filter(FILTER_TYPE_LOW_PASS, Fc1, T);
	LowsLowPass_1 = Filter(FILTER_TYPE_LOW_PASS, Fc1, T);
	LowsLowPass_2 = Filter(FILTER_TYPE_LOW_PASS, Fc1, T);
	LowsLowPass_3 = Filter(FILTER_TYPE_LOW_PASS, Fc1, T);
	
	//Mids
	MidsHighPass_0 = Filter(FILTER_TYPE_HIGH_PASS, Fc1, T);
	MidsHighPass_1 = Filter(FILTER_TYPE_HIGH_PASS, Fc1, T);
	MidsHighPass_2 = Filter(FILTER_TYPE_HIGH_PASS, Fc1, T);
	MidsHighPass_3 = Filter(FILTER_TYPE_HIGH_PASS, Fc1, T);
	MidsLowPass_0 = Filter(FILTER_TYPE_LOW_PASS, Fc2, T);
	MidsLowPass_1 = Filter(FILTER_TYPE_LOW_PASS, Fc2, T);
	MidsLowPass_2 = Filter(FILTER_TYPE_LOW_PASS, Fc2, T);
	MidsLowPass_3 = Filter(FILTER_TYPE_LOW_PASS, Fc2, T);
	
	//Highs
	HighsHighPass_0 = Filter(FILTER_TYPE_HIGH_PASS, Fc2, T);
	HighsHighPass_1 = Filter(FILTER_TYPE_HIGH_PASS, Fc2, T);
	HighsHighPass_2 = Filter(FILTER_TYPE_HIGH_PASS, Fc2, T);
	HighsHighPass_3 = Filter(FILTER_TYPE_HIGH_PASS, Fc2, T);
	
	//----------------------//
	// Load audio samples   //
	//----------------------//
	// Check the length of the audio file and allocate memory
    gSampleBufferLength = getNumFrames(gFilename);
    
    if(gSampleBufferLength <= 0) {
    	rt_printf("Error loading audio file '%s'\n", gFilename.c_str());
    	return false;
    }
    gSampleBuffer = new float[gSampleBufferLength];
    
    // Make sure the memory allocated properly
    if(gSampleBuffer == 0) {
    	rt_printf("Error allocating memory for the audio buffer.\n");
    	return false;
    }
    
    // Load the sound into the file (note: this example assumes a mono audio file)
    getSamples(gFilename, gSampleBuffer, 0, 0, gSampleBufferLength);

    rt_printf("Loaded the audio file '%s' with %d frames (%.1f seconds)\n", 
    			gFilename.c_str(), gSampleBufferLength,
    			gSampleBufferLength / context->audioSampleRate);
    
    //------------------------------------//
    // Sort out window buffer for the FFT //
    //------------------------------------//
	// Allocate the window buffer based on the FFT size
	gWindowBuffer = new float[gFFTSize];
	if(gWindowBuffer == 0)
		return false;
	// Calculate a Hann window
	for(int n = 0; n < gFFTSize; n++) {
		gWindowBuffer[n] = 0.5f * (1.0f - cosf(2.0 * M_PI * n / (float)(gFFTSize - 1)));
	}
	
	// Create the FFT object.
	fftObject = new MyFFT(gFFTSize);
	fftObject->setmFftType(Cpp); //Options are Cpp, Asm, AsmNeon
	
	// Set up the thread for the FFT
	gFFTTask = Bela_createAuxiliaryTask(process_fft_background, 75, "bela-process-fft");
	
	// Allocate the FftOutputBuffer. (Linear buffer)
	gFftOutputBuffer = new float[gFFTSize];
	
	//------------------//
	// Sort out the GUI //
	//------------------//
	// Setup GUI. By default, the Bela GUI runs on port 5555 and address 'gui'
	gui.setup(context->projectName);
	// Setup buffer of ints (holding a maximum of 1 int)
	gui.setBuffer('d', 1); // Buffer Index #0
	// Setup buffer of floats (holding a maximum of 3 floats)
	gui.setBuffer('f', 3); // Buffer Index #1
	
	// NE10 FFT-implementation setup
	timeDomainIn = (ne10_fft_cpx_float32_t*) NE10_MALLOC (gFFTSize * sizeof (ne10_fft_cpx_float32_t));
	frequencyDomain = (ne10_fft_cpx_float32_t*) NE10_MALLOC (gFFTSize * sizeof (ne10_fft_cpx_float32_t));
	cfg = ne10_fft_alloc_c2c_float32_neon (gFFTSize);
	rt_printf("Looking good for setup\n");
    	
	return true;
}

// This function handles the FFT processing once the buffer has
// been assembled.
void process_fft(float *inBuffer, int inWritePointer)
{
	//I could probably have done the bit reversal and complex Struct creation here. Oh well.
	float fftInput[gFFTSize];
	
	// Copy buffer into FFT input, unwrapping the circular buffer
	int pointer = (inWritePointer - gFFTSize + BUFFER_SIZE) % BUFFER_SIZE;
	for(int n = 0; n < gFFTSize; n++) {
		
		fftInput[n] = inBuffer[pointer] * gWindowBuffer[n];
		
		pointer++;
		if(pointer >= BUFFER_SIZE)
			pointer = 0;
	}

	// Run the FFT
	// rt_printf("Calculating an FFT!\n");
	fftObject->processFft(fftInput, gFftOutputBuffer);
	
}

/**
 * NEON FFT Implementation. Taken from the fft_overlap_add_threads lecture example.
 */
void process_fft_ne10(float *inBuffer, int inWritePointer) {
	// Copy buffer into FFT input, unwrapping the circular buffer
	int pointer = (inWritePointer - gFFTSize + BUFFER_SIZE) % BUFFER_SIZE;
	for(int n = 0; n < gFFTSize; n++) {
		timeDomainIn[n].r = (ne10_float32_t) inBuffer[pointer] * gWindowBuffer[n];
		timeDomainIn[n].i = 0;

		pointer++;
		if(pointer >= BUFFER_SIZE)
			pointer = 0;
	}
	
	// Run the FFT
	ne10_fft_c2c_1d_float32_neon (frequencyDomain, timeDomainIn, cfg, 0);
	
	// Put it back into the global output buffer.
	for (int i=0; i<gFFTSize/2; i++) {
		gFftOutputBuffer[i] = sqrt(frequencyDomain[i].r*frequencyDomain[i].r + frequencyDomain[i].i*frequencyDomain[i].i) / (gFFTSize/4.0);
	}
}

// This function runs in an auxiliary task on Bela, calling process_fft
void process_fft_background(void *)
{
	if (TESTING_NE10) {
		process_fft_ne10(gFftInputBuffer, gCachedFftInputBufferPointer);
	} else {
		process_fft(gFftInputBuffer, gCachedFftInputBufferPointer);
	}
	
	//This SEND_WAIT_TIME was useful in changing how quickly we send spectral data to the GUI.
	if (sendCounter++ == SEND_WAIT_TIME);
	{
		//For some reason I don't think I'm able to send a dynamically allocated array to the buffer.
		std::vector<float> arrayToSend;
		for (int i=0; i<gFFTSize/2; i++) {
			arrayToSend.push_back(gFftOutputBuffer[i]);
		}
		
		sendCounter=0;
		//Send the array on Buffer #3.
		gui.sendBuffer(2, arrayToSend);
	}
	
	gNumberOfFftsCurrentlyProcessing--;
}


void render(BelaContext *context, void *userData)
{
	//----------------------------------------//
	// Read in any Buffer values from the GUI //
	//----------------------------------------//
	// Get buffer 0 (filtersOn)
	DataBuffer buffer0 = gui.getDataBuffer(0);
	// Retrieve contents of the buffer as ints
	int* data0 = buffer0.getAsInt();
	if ((*data0) == 1) {
		filtersOn = true;
		
	} else if ((*data0) == 0) {
		
		filtersOn = false;
	}
	// Get buffer 1 (gains)
	DataBuffer buffer1 = gui.getDataBuffer(1);
	float* data1 = buffer1.getAsFloat();
	for (int i=0; i<3; i++) {
		gains[i] = data1[i];
	}
	
	
	for(unsigned int n = 0; n < context->audioFrames; n++) {
        // Read the next sample from the stored sound
        float in = gSampleBuffer[gReadPointer];
        if(++gReadPointer >= gSampleBufferLength)
        	gReadPointer = 0;
        
        
		//------------------------------//
		// Do any audio processing here //
		//------------------------------//
		float out = in;
		
		if (filtersOn) {
			//Filter into different bands
			float outLows = LowsLowPass_0.processNew(in);
			outLows = LowsLowPass_1.processNew(outLows);
			outLows = LowsLowPass_2.processNew(outLows);
			outLows = LowsLowPass_3.processNew(outLows);
			
			float outMids = MidsLowPass_0.processNew(in);
			outMids = MidsLowPass_1.processNew(outMids);
			outMids = MidsLowPass_2.processNew(outMids);
			outMids = MidsLowPass_3.processNew(outMids);
			outMids = MidsHighPass_0.processNew(outMids);
			outMids = MidsHighPass_1.processNew(outMids);
			outMids = MidsHighPass_2.processNew(outMids);
			outMids = MidsHighPass_3.processNew(outMids);
			
			float outHighs = HighsHighPass_0.processNew(in);
			outHighs = HighsHighPass_1.processNew(outHighs);
			outHighs = HighsHighPass_2.processNew(outHighs);
			outHighs = HighsHighPass_3.processNew(outHighs);
			
			//Then apply gains to each band separately. Using a log scaling so it's perceptually nice.
			out =  pow(10.0,(60.0*(gains[0]-1))/20.0)*outLows + pow(10.0,(60.0*(gains[1]-1))/20.0)*outMids + pow(10.0,(60.0*(gains[2]-1))/20.0)*outHighs;
		}
		//Noise generator. Just an idea for an extra GUI component.
		//out += 0.03 * ((float)rand())/((float)RAND_MAX);
		
		//----------------------------//
		// Stop audio processing here //
		//----------------------------//
        
        
        //---------------------------------//
        // FFT processing code starts here //
		//---------------------------------//
		// Store the input in the input circular buffer
		gFftInputBuffer[gFftInputBufferPointer] = out;
		gFftInputBufferPointer++;
		if(gFftInputBufferPointer >= BUFFER_SIZE)
			gFftInputBufferPointer = 0;
		
		// Increment the hop counter, and start a new FFT when we reach the hop size
		gHopCounter++;
		if(gHopCounter >= gHopSize) {
			gHopCounter = 0;
			
			gNumberOfFftsCurrentlyProcessing++;
			if(gNumberOfFftsCurrentlyProcessing != gLastNumberOfFftsCurrentlyProcessing) {
				rt_printf("%d FFT(s) processing..\n", gNumberOfFftsCurrentlyProcessing);
				gLastNumberOfFftsCurrentlyProcessing = gNumberOfFftsCurrentlyProcessing;
			}
			
			//Schedule the auxiliary task to process the FFT.
			gCachedFftInputBufferPointer = gFftInputBufferPointer;
			Bela_scheduleAuxiliaryTask(gFFTTask);
		}
		//-------------------------------//
        // FFT processing code ends here //
		//-------------------------------//
		
		
		// Play the output sound (mono).
		for(unsigned int channel = 0; channel < context->audioOutChannels; channel++) {
			audioWrite(context, n, channel, out);
		}
	}
	
}

void cleanup(BelaContext *context, void *userData)
{
	// Free up any buffers allocated
    delete[] gSampleBuffer;
    delete[] gWindowBuffer;
    delete[] gFftOutputBuffer;
    delete[] fftObject;
}
