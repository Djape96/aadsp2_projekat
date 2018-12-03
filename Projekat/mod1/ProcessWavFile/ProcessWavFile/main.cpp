
#include <stdlib.h>
#include <string.h>
#include "WAVheader.h"
#include "common.h"
#include <cmath>
#include "expander.h"

#define BLOCK_SIZE 16
#define MAX_NUM_CHANNEL 8
#define NUM_OF_CHANNELS 5

//enum output_mode o_mode = ( mode1,mode2 );

double sampleBuffer[MAX_NUM_CHANNEL][BLOCK_SIZE];
AudioExpander_t expander;

void processing()
{
	int i;
	int mode = 1;
	const double tap_gain_minus6db = 0;
	const double tap_gain_5db = 0;
	const double tap_gain_minus16db = 0;
	const double tap_gain_3db = 0;
	

	double* pChannel0 = sampleBuffer[0];
	double* pChannel1 = sampleBuffer[1];
	double* pChannel2 = sampleBuffer[2];
	double* pChannel3 = sampleBuffer[3];
	double* pChannel4 = sampleBuffer[4];

	double* tmp_channel0 = pChannel0;
	double* tmp_channel1 = pChannel1;

	for (i = 0; i < BLOCK_SIZE; i++)
	{ 
		*(tmp_channel0) *= tap_gain_minus6db;
		*(tmp_channel1) *= tap_gain_minus6db;
	}
	for (i = 0; i < BLOCK_SIZE; i++) {
		*pChannel3 = *(tmp_channel1) * (-1);
	}
	if (mode == 1)
	{
		for (i = 0; i < BLOCK_SIZE; i++)
		{
			*pChannel0 = *(tmp_channel0)* tap_gain_minus6db;
			*pChannel2 = *(tmp_channel0)* tap_gain_3db;
		}
	}
	else {
		for (i = 0; i < BLOCK_SIZE; i++)
		{
			*pChannel0 = *(tmp_channel0)* tap_gain_minus16db;
			*pChannel2 = *(tmp_channel0)* tap_gain_5db;
		}
	}
	
	gst_audio_dynamic_transform_expander_double(&expander, tmp_channel0);
	gst_audio_dynamic_transform_expander_double(&expander, tmp_channel1);

}

int main(int argc, char* argv[])
{
	FILE *wav_in=NULL;
	FILE *wav_out=NULL;
	char WavInputName[256];
	char WavOutputName[256];
	WAV_HEADER inputWAVhdr,outputWAVhdr;	

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

			//processing();

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