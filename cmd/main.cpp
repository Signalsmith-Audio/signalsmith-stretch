#include "util/stopwatch.h"
#include "util/memory-tracker.hxx"

#include <iostream>
#define LOG_EXPR(expr) std::cout << #expr << " = " << (expr) << "\n";

#include "../signalsmith-stretch.h"
#include "util/simple-args.h"
#include "util/wav.h"

int main(int argc, char* argv[]) {
	signalsmith::stretch::SignalsmithStretch<float/*, std::mt19937*/> stretch; // optional cheaper RNG for performance comparison

	SimpleArgs args(argc, argv);
	
	if (args.hasFlag("v", "prints the version")) {
		std::cout << stretch.version[0] << "." << stretch.version[1] << "." << stretch.version[2] << "\n";
		return 0;
	}
	
	std::string inputWav = args.arg<std::string>("input.wav", "16-bit WAV file");
	std::string outputWav = args.arg<std::string>("output.wav", "output WAV file");
	
	double semitones = args.flag<double>("semitones", "pitch-shift amount", 0);
	double tonality = args.flag<double>("tonality", "tonality limit (Hz)", 8000);
	double time = args.flag<double>("time", "time-stretch factor", 1);
	bool exactLength = args.hasFlag("exact", "trims the start/end so the output has the correct length");
	args.errorExit();
	
	std::cout << Console::Bright << inputWav << Console::Reset;
	std::cout << " -> ";
	std::cout << Console::Bright << outputWav << Console::Reset << "\n";
	std::cout << "\tsemitones: " << semitones << "\n\t     time: " << time << "x" << (exactLength ? " (exact)" : "") << "\n\t tonality: " << tonality << "Hz\n";

	Wav inWav;
	if (!inWav.read(inputWav).warn()) args.errorExit("failed to read WAV");
	size_t inputLength = inWav.samples.size()/inWav.channels;
	
	Wav prevWav; // Used during development, to compare against known-good previous render
	bool compareReference = (time <= 1.6);
	if (compareReference && !prevWav.read(outputWav + "-reference.wav")) {
		if (prevWav.read(outputWav)) {
			prevWav.write(outputWav + "-reference.wav");
		}
	}

	Wav outWav;
	outWav.channels = inWav.channels;
	outWav.sampleRate = inWav.sampleRate;
	int outputLength = std::round(inputLength*time);

	signalsmith::MemoryTracker initMemory;
	signalsmith::Stopwatch stopwatch;

	stopwatch.start();
	stretch.presetDefault(inWav.channels, inWav.sampleRate);
	stretch.setTransposeSemitones(semitones, tonality/inWav.sampleRate);
	double initSeconds = stopwatch.seconds(stopwatch.lap());

	initMemory = initMemory.diff();
	std::cout << "Setup:\n\t" << initSeconds << "s\n";
	if (initMemory.implemented) {
		std::cout << "\tallocated " << (initMemory.allocBytes/1000) << "kB, freed " << (initMemory.freeBytes/1000) << "kB\n";
	}

	// pad the input at the end, since we'll be reading slightly ahead
	size_t paddedInputLength = inputLength + stretch.inputLatency();
	inWav.samples.resize(paddedInputLength*inWav.channels);
	// pad the output at the end, since we have output latency as well
	int tailSamples = exactLength ? stretch.outputLatency() : (stretch.outputLatency() + stretch.inputLatency()); // if we don't need exact length, add a bit more output to catch any wobbles past the end
	int paddedOutputLength = outputLength + tailSamples;
	outWav.samples.resize(paddedOutputLength*outWav.channels);

	signalsmith::MemoryTracker processMemory;

	stopwatch.start();
	// The simplest way to deal with input latency (when have access to the audio buffer) is to always be slightly ahead in the input
	stretch.seek(inWav, stretch.inputLatency(), 1/time);
	inWav.offset += stretch.inputLatency();
	// Process it all in one call, although it works just the same if we split into smaller blocks
	stretch.process(inWav, inputLength, outWav, outputLength);
	// Read the last bit of output without giving it any more input
	outWav.offset += outputLength;
	stretch.flush(outWav, tailSamples);
	outWav.offset -= outputLength;

	double processSeconds = stopwatch.seconds(stopwatch.lap());
	double processRate = (inWav.length()/inWav.sampleRate)/processSeconds;
	double processPercent = 100/processRate;
	processMemory = processMemory.diff();
	std::cout << "Process:\n\t" << processSeconds << "s, " << processRate << "x realtime, " << processPercent << "% CPU\n";
	if (processMemory.implemented) {
		std::cout << "\tallocated " << (processMemory.allocBytes/1000) << "kB, freed " << (processMemory.freeBytes/1000) << "kB\n";
		if (processMemory) args.errorExit("allocated during process()");
	}
	
	if (exactLength) {
		// The start has some extra output - we could just trim it, but we might as well fold it back into the output
		for (size_t c = 0; c < outWav.channels; ++c) {
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
	
	if (compareReference && prevWav.result) {
		outWav.read(outputWav);
		if (prevWav.length() != outWav.length()) args.errorExit("lengths differ");
		double diff2 = 0;
		for (size_t i = 0; i < prevWav.samples.size(); ++i) {
			double diff = prevWav.samples[i] - outWav.samples[i];
			diff2 += diff*diff;
		}
		diff2 /= prevWav.samples.size();
		double diffDb = 10*std::log10(diff2);
		std::cout << "Reference:\n\tdifference: ";
		if (diff2 < 1e-6) std::cout << Console::Red;
		if (diff2 < 1e-8) std::cout << Console::Yellow;
		if (diff2 < 1e-10) std::cout << Console::Green;
		std::cout << Console::Bright << diffDb << Console::Reset << " dB\n";
		if (diffDb > -60) args.errorExit("too much difference\n");
	}
}
