
#include <stdlib.h>
#include <string.h>
#include "WAVheader.h"
#include "expander.h"

#define BLOCK_SIZE 16
#define MAX_NUM_CHANNEL 8
#define NUM_OF_CHANNELS 5

double sampleBuffer[MAX_NUM_CHANNEL][BLOCK_SIZE];


void processing(double pInbuff[][BLOCK_SIZE], double pOutbuff[][BLOCK_SIZE])
{
	int mode = 1;
	double minus6db_gain = 0.25;
	double minus16db_gain = 0.025;
	double minus3db_gain = 0.5;
	double minus5db_gain = 0.32;

	AudioExpander_t * expander;


	for (int i = 0; i < BLOCK_SIZE; i++)
	{
		pInbuff[0][i] = pInbuff[0][i] * minus6db_gain; // L -> input gain = -6db
	}

	for (int i = 0; i < BLOCK_SIZE; i++)
	{
		pInbuff[1][i] = pInbuff[1][i] * minus6db_gain; //R -> input gain = -6db
	}

	for (int i = 0; i < BLOCK_SIZE; i++)
	{
		pInbuff[2][i] = pInbuff[0][i] * minus16db_gain; //L[0]-> input gain = -16db
	}

	for (int i = 0; i < BLOCK_SIZE; i++)
	{
		pInbuff[3][i] = pInbuff[0][i] * minus6db_gain;  //L[1]-> input gain = -6db
	}

	for (int i = 0; i < BLOCK_SIZE; i++)
	{
		pInbuff[4][i] = pInbuff[0][i] * minus5db_gain;  //L[2]-> input gain = -5db
	}

	for (int i = 0; i < BLOCK_SIZE; i++)
	{
		pInbuff[5][i] = pInbuff[0][i] * minus3db_gain; //L[3]-> input gain = -3db
	}

	for (int i = 0; i < BLOCK_SIZE; i++)
	{
		pInbuff[6][i] = pInbuff[1][i] * (-1); //R[0]-> (-1)
	}

	gst_audio_dynamic_transform_expander_double(expander, pInbuff[0], BLOCK_SIZE);
	gst_audio_dynamic_transform_expander_double(expander, pInbuff[1], BLOCK_SIZE);

	for (int i = 0; i < BLOCK_SIZE; i++)
	{
		if (mode == 1)
		{
			pOutbuff[0][i] = pInbuff[2][i];
		}
		else
		{
			pOutbuff[0][i] = pInbuff[3][i];
		}

		pOutbuff[1][i] = pInbuff[0][i];

		if (mode == 1)
		{
			pOutbuff[2][i] = pInbuff[5][i];
		}
		else
		{
			pOutbuff[2][i] = pInbuff[4][i];
		}

		pOutbuff[3][i] = pInbuff[6][i];

		pOutbuff[4][i] = pInbuff[1][i];
	}
}

int main(int argc, char* argv[])
{
	FILE *wav_in=NULL;
	FILE *wav_out=NULL;
	char WavInputName[256];
	char WavOutputName[256];
	WAV_HEADER inputWAVhdr,outputWAVhdr;	

	AudioExpander_t expander;

	// Init channel buffers
	for(int i=0; i<MAX_NUM_CHANNEL; i++)
		memset(&sampleBuffer[i],0,BLOCK_SIZE);

	// Open input and output wav files
	//-------------------------------------------------
	strcpy(WavInputName,argv[1]);
	wav_in = OpenWavFileForRead (WavInputName,"rb");
	strcpy(WavOutputName,argv[2]);
	wav_out = OpenWavFileForRead (WavOutputName,"wb");
	//-------------------------------------------------

	// Read input wav header
	//-------------------------------------------------
	ReadWavHeader(wav_in,inputWAVhdr);
	//-------------------------------------------------
	
	// Set up output WAV header
	//-------------------------------------------------	
	outputWAVhdr = inputWAVhdr;
	outputWAVhdr.fmt.NumChannels = inputWAVhdr.fmt.NumChannels; // change number of channels

	int oneChannelSubChunk2Size = inputWAVhdr.data.SubChunk2Size/inputWAVhdr.fmt.NumChannels;
	int oneChannelByteRate = inputWAVhdr.fmt.ByteRate/inputWAVhdr.fmt.NumChannels;
	int oneChannelBlockAlign = inputWAVhdr.fmt.BlockAlign/inputWAVhdr.fmt.NumChannels;
	
	outputWAVhdr.data.SubChunk2Size = oneChannelSubChunk2Size*outputWAVhdr.fmt.NumChannels;
	outputWAVhdr.fmt.ByteRate = oneChannelByteRate*outputWAVhdr.fmt.NumChannels;
	outputWAVhdr.fmt.BlockAlign = oneChannelBlockAlign*outputWAVhdr.fmt.NumChannels;


	// Write output WAV header to file
	//-------------------------------------------------
	WriteWavHeader(wav_out,outputWAVhdr);

	//Init audio expander
	audio_expander_init(&expander);


	// Processing loop
	//-------------------------------------------------	
	{
		int sample;
		int BytesPerSample = inputWAVhdr.fmt.BitsPerSample/8;
		const double SAMPLE_SCALE = -(double)(1 << 31);		//2^31
		int iNumSamples = inputWAVhdr.data.SubChunk2Size/(inputWAVhdr.fmt.NumChannels*inputWAVhdr.fmt.BitsPerSample/8);
		
		// exact file length should be handled correctly...
		for(int i=0; i<iNumSamples/BLOCK_SIZE; i++)
		{	
			for(int j=0; j<BLOCK_SIZE; j++)
			{
				for(int k=0; k<inputWAVhdr.fmt.NumChannels; k++)
				{	
					sample = 0; //debug
					fread(&sample, BytesPerSample, 1, wav_in);
					sample = sample << (32 - inputWAVhdr.fmt.BitsPerSample); // force signextend
					sampleBuffer[k][j] = sample / SAMPLE_SCALE;				// scale sample to 1.0/-1.0 range		
				}
			}

			processing(sampleBuffer, sampleBuffer);

			for(int j=0; j<BLOCK_SIZE; j++)
			{
				for(int k=0; k<outputWAVhdr.fmt.NumChannels; k++)
				{	
					sample = sampleBuffer[k][j] * SAMPLE_SCALE ;	// crude, non-rounding 			
					sample = sample >> (32 - inputWAVhdr.fmt.BitsPerSample);
					fwrite(&sample, outputWAVhdr.fmt.BitsPerSample/8, 1, wav_out);		
				}
			}		
		}
	}
	
	// Close files
	//-------------------------------------------------	
	fclose(wav_in);
	fclose(wav_out);
	//-------------------------------------------------	

	return 0;
}