#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <signal.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <curl/curl.h>
#include <complex.h>
#include "rtl-sdr.h"

using namespace std;

#define DEFAULT_FREQUENCY			433900000
#define DEFAULT_SAMPLE_RATE			1003008
#define DEFAULT_LEVEL_LIMIT			30
#define DEFAULT_MAX_SAMPLES			80000
#define DEFAULT_MIN_PROCESS_SAMPLES	77000
#define DEFAULT_BIT_SAMPLE_LENGTH	256
#define DEFAULT_SIGNAL_LENGTH		350
#define DEFAULT_BUF_LENGTH			(16 * 16384)
#define DEFAULT_FFT_LENGTH			131072		/* It has to be a power of 2, 65536 samples of i + j  */
#define PI	M_PI	/* pi to machine precision, defined in math.h */
#define TWOPI	(2.0*PI)

/*#define HIFREQ					157640
#define LOFREQ						82500*/

/*#define HIFREQ						154500
#define LOFREQ						85500*/
#define HIFREQ						157750
#define LOFREQ						82500

static int do_exit = 0;
static int do_exit_async = 0, frequencies = 0, events = 0;
static rtlsdr_dev_t *dev = NULL;
uint32_t samp_rate = DEFAULT_SAMPLE_RATE;
uint32_t freq = DEFAULT_FREQUENCY;

typedef double complex cplx;

struct sample
{
	int8_t Real;
	int8_t Imag;
};

struct dm_state {
	sample *samples;
	float *magnitudes;
	int sampleCount;

	uint8_t *signal;
	int bitCount;

	double *fft;
};

static void send_to_emon(uint16_t watts)
{
	CURL *curl;
	CURLcode res;

	char url[256];

	sprintf(url, "%s?json={power:{%u}}&apikey=%s&time=%u", "http://emoncms.org/input/post.json", watts, "e944f669563901f53cfadf1869804832", time(NULL));

	curl = curl_easy_init();
	if (curl) {
		curl_easy_setopt(curl, CURLOPT_URL, url);

		curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
		curl_easy_setopt(curl, CURLOPT_VERBOSE, 0L);
		curl_easy_setopt(curl, CURLOPT_TIMEOUT, 5L);

		res = curl_easy_perform(curl);

		if (res != CURLE_OK)
			fprintf(stderr, "curl_easy_perform() failed: %s\n",
			curl_easy_strerror(res));

		curl_easy_cleanup(curl);
	}
}

static float goertzel_mag(int numSamples, int TARGET_FREQUENCY, int SAMPLING_RATE, float* data)
{
	int     k, i;
	float   floatnumSamples;
	float   omega, sine, cosine, coeff, q0, q1, q2, magnitude, real, imag;

	float   scalingFactor = numSamples / 2.0;

	floatnumSamples = (float)numSamples;
	k = (int)(0.5 + ((floatnumSamples * TARGET_FREQUENCY) / SAMPLING_RATE));
	omega = (2.0 * M_PI * k) / floatnumSamples;
	sine = sin(omega);
	cosine = cos(omega);
	coeff = 2.0 * cosine;
	q0 = 0;
	q1 = 0;
	q2 = 0;

	for (i = 0; i<numSamples; i++)
	{
		q0 = coeff * q1 - q2 + data[i];
		q2 = q1;
		q1 = q0;
	}

	// calculate the real and imaginary results
	// scaling appropriately
	real = (q1 - q2 * cosine) / scalingFactor;
	imag = (q2 * sine) / scalingFactor;

	magnitude = sqrtf(real*real + imag*imag);
	return magnitude;
}

static bool demanchesterise(uint8_t *signal, int length)
{
	/*int bufflen = length / 2;
	if (length % 2 == 1) bufflen++;

	uint8_t *buffer = (uint8_t *)malloc(bufflen);*/

	int finalBitPos = 0;

	for (int i = 0; i < length; i++)
	{
		uint8_t b = signal[i];

		for (int bitPos = 0; bitPos < 7; bitPos += 2)
		{
			uint8_t n = (b >> (6 - bitPos)) & 3;
			if (n == 0 || n == 3)
				return false;
			else if (n == 1)
			{
				signal[finalBitPos / 8] &= 0xff ^ (1 << (7 - (finalBitPos % 8)));
				finalBitPos++;
			}
			else if (n == 2)
			{
				signal[finalBitPos / 8] |= 1 << (7 - (finalBitPos % 8));
				finalBitPos++;
			}
		}
	}

	return true;
}

static void log_to_file(uint16_t watts)
{
	time_t rawtime;
	struct tm *info;
	char buffer[80];

	time(&rawtime);

	info = localtime(&rawtime);

	FILE *f = fopen("watts.log", "a");

	strftime(buffer, 80, "%d/%m/%Y %H:%M:%S", info);

	fprintf(f, "%s - %u watts\n", buffer, watts);

	fclose(f);
}

void _fft(cplx buf[], cplx out[], int n, int step)
{
	if (step < n) {
		_fft(out, buf, n, step * 2);
		_fft(out + step, buf + step, n, step * 2);

		for (int i = 0; i < n; i += 2 * step) {
			cplx t = cexp(-I * M_PI * i / n) * out[i + step];
			buf[i / 2] = out[i] + t;
			buf[(i + n) / 2] = out[i] - t;
		}
	}
}

void fft(cplx buf[], int n)
{
	cplx out[n];
	for (int i = 0; i < n; i++) out[i] = buf[i];

	_fft(buf, out, n, 1);
}

static void process_samples(void *ctx)
{
	struct dm_state *demod = (dm_state *)ctx;
	int offset = 25088;

	/* Debug */
	/*FILE *f = fopen("packet.csv", "w");
	for (int i = 0; i < demod->sampleCount; i++)
	{
		fprintf(f, "%f\n", demod->magnitudes[i]);
	}
	fclose(f);*/

	/* Work out FFT to find 2 most common frequencies */
	//memset(demod->fft, 0, DEFAULT_FFT_LENGTH);
	//
	//int fftpos = 0;
	//for (int sPos = offset; sPos < demod->sampleCount && sPos < offset + DEFAULT_FFT_LENGTH; sPos++)
	//{
	//	demod->fft[fftpos] = (double)demod->samples[sPos];		// This is the i component
	//	
	//	//demod->fft[fftpos + 1] = 0.0;			// This is what this does but is not required as memset does this above (this is the j component)
	//	fftpos += 2;
	//}

	//// We don't have to worry about padding to the nearest power of 2, because memset does this above (which fills the 131072)
	//// Calculate FFT
	/*four1((int8_t*)demod->samples, DEFAULT_FFT_LENGTH / 2, 1);
	FILE *f = fopen("fft.csv", "w");
	for (int i = 0; i < DEFAULT_FFT_LENGTH / 2; i++)
	{
		float mag = sqrt(demod->samples[i].Imag * demod->samples[i].Imag + demod->samples[i].Real * demod->samples[i].Real);
		fprintf(f, "%f\n", mag);
	}
	fclose(f);*/

	///* Find the first mode frequency */
	//int hipos = 0;
	//double himag = 0;
	//for (int sPos = 0; sPos < DEFAULT_FFT_LENGTH / 2; sPos += 2)
	//{
	//	double mag = sqrt((demod->fft[sPos] * demod->fft[sPos]) + (demod->fft[sPos + 1] * demod->fft[sPos + 1]));
	//	if (mag > himag)
	//	{
	//		hipos = sPos / 2;
	//		himag = mag;
	//	}
	//}

	///* Find the second mode frequency */
	//int lopos = 0;
	//double lomag = 0;
	//for (int sPos = 0; sPos < DEFAULT_FFT_LENGTH / 2; sPos += 2)
	//{
	//	double mag = sqrt((demod->fft[sPos] * demod->fft[sPos]) + (demod->fft[sPos + 1] * demod->fft[sPos + 1]));
	//	if (mag > himag && abs((sPos / 2)-hipos) > 256)
	//	{
	//		lopos = sPos / 2;
	//		lomag = mag;
	//	}
	//}

	///* The should give us the frequencies */
	//hipos = (hipos / DEFAULT_FFT_LENGTH / 2) * DEFAULT_SAMPLE_RATE;
	//lopos = (lopos / DEFAULT_FFT_LENGTH / 2) * DEFAULT_SAMPLE_RATE;

	/* Reset bit count */
	demod->bitCount = 0;
	memset(demod->signal, 0, DEFAULT_SIGNAL_LENGTH);

	for (int sPos = offset; sPos < demod->sampleCount - DEFAULT_BIT_SAMPLE_LENGTH; sPos += DEFAULT_BIT_SAMPLE_LENGTH)
	{
		float lofreq = goertzel_mag(DEFAULT_BIT_SAMPLE_LENGTH, HIFREQ, DEFAULT_SAMPLE_RATE, &demod->magnitudes[sPos]);
		float hifreq = goertzel_mag(DEFAULT_BIT_SAMPLE_LENGTH, LOFREQ, DEFAULT_SAMPLE_RATE, &demod->magnitudes[sPos]);

		/* Not sure why the bits are shifted at this position, this corrects it */
		if (demod->bitCount == 72)
		{
			demod->signal[demod->bitCount / 8] |= (uint8_t)(1 << (7 - (demod->bitCount % 8)));
			demod->bitCount++;
		}

		if (hifreq > lofreq)
		{
			if (demod->bitCount == 0)
			{
				demod->bitCount = 2;       // Add two 0 bits at the beginning
			}

			/* Add 1 bit */
			demod->signal[demod->bitCount / 8] |= (uint8_t)(1 << (7 - (demod->bitCount % 8)));
			demod->bitCount++;
		}
		else if (demod->bitCount > 0)  // Make sure we have a 1 bit first
		{
			demod->bitCount++;         // Low bit, add 0 bit (no computation required)
		}
	}

	/* We do not have a valid signal*/
	if (demod->signal[4] != 0x2d || demod->signal[5] != 0xd4)
		return;

	if (demanchesterise(&demod->signal[10], 4))
	{
		if ((demod->signal[10] & 0x80) == 0x80)
		{
			// We have a valid wattage!! :)
			uint16_t wattage = demod->signal[11] | ((demod->signal[10] & 0x7f) << 8);
			//fprintf(stderr, "Wattage: %u W\n", wattage);
			log_to_file(wattage);
			send_to_emon(wattage);
		}
	}
}

static void rtlsdr_callback(unsigned char *buf, uint32_t len, void *ctx) {
	struct dm_state *demod = (dm_state *)ctx;
	//signed char *buffer = (signed char *)buf;
	//uint16_t* sbuf = (uint16_t*)buf;

	if (do_exit || do_exit_async)
		return;

	for (int bPos = 0; bPos < len - 1; bPos += 2)
	{
		int8_t im = buf[bPos] - 127;
		int8_t re = buf[bPos + 1] - 127;
		float mag = sqrtf(re * re + im * im);

		/*demod->fft[bPos] = (double)im;
		demod->fft[bPos+1] = (double)re;*/

		if (mag > DEFAULT_LEVEL_LIMIT)
		{
			demod->samples[demod->sampleCount].Real = re;
			demod->samples[demod->sampleCount].Imag = im;
			demod->magnitudes[demod->sampleCount] = mag;

			demod->sampleCount++;
		}
		else
		{
			if (demod->sampleCount >= DEFAULT_MIN_PROCESS_SAMPLES)
			{
				process_samples(ctx);
			}

			demod->sampleCount = 0;
		}
	}
}

static void sighandler(int signum) {
	if (signum == SIGPIPE) {
		signal(SIGPIPE, SIG_IGN);
	}
	else {
		fprintf(stderr, "Signal caught, exiting!\n");
	}
	do_exit = 1;
	rtlsdr_cancel_async(dev);
}

int main(int argc, char *argv[])
{
	int device_count, i, r;
	char vendor[256], product[256], serial[256];
	uint32_t dev_index = 0;
	struct sigaction sigact;
	uint32_t out_block_size = DEFAULT_BUF_LENGTH;
	uint8_t *buffer;
	struct dm_state* demod;
	
	demod = (dm_state*)malloc(sizeof(struct dm_state));
	memset(demod, 0, sizeof(struct dm_state));

	demod->samples = (sample*)malloc(sizeof(sample) * DEFAULT_MAX_SAMPLES);
	demod->magnitudes = (float*)malloc(sizeof(float) * DEFAULT_MAX_SAMPLES);
	demod->signal = (uint8_t*)malloc(sizeof(uint8_t) * DEFAULT_SIGNAL_LENGTH);
	demod->fft = (double*)malloc(sizeof(double) * DEFAULT_FFT_LENGTH);

	cplx buf[] = { 1, 1, 1, 1, 0, 0, 0, 0 };

	/*char sz[] = "Hello, World!";	//Hover mouse over "sz" while debugging to see its contents
	cout << sz << endl;	//<================= Put a breakpoint here*/

	/*char sz[128];
	gethostname(sz, sizeof(sz));
	printf("The hostname is %s\n", sz);*/

	buffer = (uint8_t *)malloc(out_block_size * sizeof(uint8_t));

	device_count = rtlsdr_get_device_count();
	if (!device_count) {
		fprintf(stderr, "No supported devices found.\n");
		exit(1);
	}
	
	fprintf(stderr, "Found %d device(s):\n", device_count);
	for (i = 0; i < device_count; i++) {
		rtlsdr_get_device_usb_strings(i, vendor, product, serial);
		fprintf(stderr, "  %d:  %s, %s, SN: %s\n", i, vendor, product, serial);
	}
	fprintf(stderr, "\n");

	fprintf(stderr, "Using device %d: %s\n",
		dev_index, rtlsdr_get_device_name(dev_index));

	r = rtlsdr_open(&dev, dev_index);
	if (r < 0) {
		fprintf(stderr, "Failed to open rtlsdr device #%d.\n", dev_index);

		exit(1);
	}

	sigact.sa_handler = sighandler;
	sigemptyset(&sigact.sa_mask);
	sigact.sa_flags = 0;
	sigaction(SIGINT, &sigact, NULL);
	sigaction(SIGTERM, &sigact, NULL);
	sigaction(SIGQUIT, &sigact, NULL);
	sigaction(SIGPIPE, &sigact, NULL);

	/* Set the sample rate */
	r = rtlsdr_set_sample_rate(dev, samp_rate);
	if (r < 0)
		fprintf(stderr, "WARNING: Failed to set sample rate.\n");
	else
		fprintf(stderr, "Sample rate set to %d.\n", rtlsdr_get_sample_rate(dev)); // Unfortunately, doesn't return real rate

	/* Set the tuner gain */
	r = rtlsdr_set_tuner_gain_mode(dev, 0);
	if (r < 0)
		fprintf(stderr, "WARNING: Failed to enable automatic gain.\n");
	else
		fprintf(stderr, "Tuner gain set to Auto.\n");

	r = rtlsdr_set_freq_correction(dev, 0);

	/* Reset endpoint before we start reading from it (mandatory) */
	r = rtlsdr_reset_buffer(dev);
	if (r < 0)
		fprintf(stderr, "WARNING: Failed to reset buffers.\n");

	while (!do_exit) {
		/* Set the frequency */
		r = rtlsdr_set_center_freq(dev, freq);
		if (r < 0)
			fprintf(stderr, "WARNING: Failed to set center freq.\n");
		else
			fprintf(stderr, "Tuned to %u Hz.\n", rtlsdr_get_center_freq(dev));
		r = rtlsdr_read_async(dev, rtlsdr_callback, (void *)demod,
			0, out_block_size);
		do_exit_async = 0;
	}

	if (do_exit)
		fprintf(stderr, "\nUser cancel, exiting...\n");
	else
		fprintf(stderr, "\nLibrary error %d, exiting...\n", r);

	rtlsdr_close(dev);


	return 0;
}

