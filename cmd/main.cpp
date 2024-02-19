#include <iostream>
#define LOG_EXPR(expr) std::cout << #expr << " = " << (expr) << "\n";

#include <ctime>

#include "../signalsmith-stretch.h"
#include "util/simple-args.h"
#include "util/wav.h"

int main(int argc, char* argv[]) {
	SimpleArgs args(argc, argv);
	
	std::string inputWav = args.arg<std::string>("input.wav", "16-bit WAV file");
	std::string outputWav = args.arg<std::string>("output.wav", "output WAV file");
	
	double semitones = args.flag<double>("semitones", "pitch-shift amount", 0);
	double tonality = args.flag<double>("tonality", "tonality limit (Hz)", 8000);
	double time = args.flag<double>("time", "time-stretch factor", 1);
	bool exactLength = args.hasFlag("exact", "trims the start/end so the output has the correct length");
	args.errorExit();

	Wav inWav;
	if (!inWav.read(inputWav).warn()) args.errorExit("failed to read WAV");
	size_t inputLength = inWav.samples.size()/inWav.channels;

	Wav outWav;
	outWav.channels = inWav.channels;
	outWav.sampleRate = inWav.sampleRate;
	int outputLength = std::round(inputLength*time);

	signalsmith::stretch::SignalsmithStretch<float> stretch;
	stretch.presetDefault(inWav.channels, inWav.sampleRate);
	stretch.setTransposeSemitones(semitones, tonality/inWav.sampleRate);

	// pad the input at the end, since we'll be reading slightly ahead
	size_t paddedInputLength = inputLength + stretch.inputLatency();
	inWav.samples.resize(paddedInputLength*inWav.channels);
	// pad the output at the end, since we have output latency as well
	int tailSamples = exactLength ? stretch.outputLatency() : (stretch.outputLatency() + stretch.inputLatency()); // if we don't need exact length, add a bit more output to catch any wobbles past the end
	int paddedOutputLength = outputLength + tailSamples;
	outWav.samples.resize(paddedOutputLength*outWav.channels);

	// The simplest way to deal with input latency is to always be slightly ahead in the input
	stretch.seek(inWav, stretch.inputLatency(), 1/time);
	
	// Process it all in one call, although it works just the same if we split into smaller blocks
	inWav.offset += stretch.inputLatency();

	// These lengths in the right ratio to get the time-stretch
	stretch.process(inWav, inputLength, outWav, outputLength);

	// Read the last bit of output without giving it any more input
	outWav.offset += outputLength;
	stretch.flush(outWav, tailSamples);
	outWav.offset -= outputLength;
	
	if (exactLength) {
		// The start has some extra output - we could just trim it, but we might as well fold it back into the output
		for (int c = 0; c < outWav.channels; ++c) {
			for (int i = 0; i < stretch.outputLatency(); ++i) {
				double trimmed = outWav[stretch.outputLatency() - 1 - i][c];
				outWav[stretch.outputLatency() + i][c] -= trimmed; // reversed in time and negated
			}
		}
		// Skips the output
		outWav.offset += stretch.outputLatency();

		// the `.flush()` call already handled foldback stuff at the end (since we asked for a shorter `tailSamples`)
	}

	if (!outWav.write(outputWav).warn()) args.errorExit("failed to write WAV");
}
