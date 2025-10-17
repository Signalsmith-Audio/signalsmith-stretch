// helper for debugging
#include <iostream>
#define LOG_EXPR(expr) std::cout << #expr << " = " << (expr) << "\n";

#include "signalsmith-stretch/signalsmith-stretch.h"
using SignalsmithStretch = signalsmith::stretch::SignalsmithStretch<float>;

#include "./util/simple-args.h"
#include "./util/wav.h"
#include "./util/stopwatch.h"

#include "plot/plot.h"

int main(int argc, char* argv[]) {
	SimpleArgs args(argc, argv);
	
	if (args.hasFlag("v", "prints the version")) {
		auto &version = SignalsmithStretch::version;
		std::cout << version[0] << "." << version[1] << "." << version[2] << "\n";
		return 0;
	}

	std::string inputWav = args.arg<std::string>("input.wav", "16-bit WAV file");
	std::string outputWav = args.arg<std::string>("output.wav", "output WAV file");
	double time = args.flag<double>("time", "time-stretch factor", 1);
	double semitones = args.flag<double>("semitones", "pitch-shift amount", 0);
	double formants = args.flag<double>("formant", "formant-shift amount (semitones)", 0);
	bool formantComp = args.hasFlag("formant-comp", "formant compensation");
	double formantBase = args.flag<double>("formant-base", "formant base frequency (Hz, 0=auto)", 100);
	double tonality = args.flag<double>("tonality", "tonality limit (Hz)", 8000);
	double asymmetry = args.flag<double>("asymmetry", "asymmetrical STFT analysis (0-1)", 0);
	bool splitComputation = args.hasFlag("split-computation", "distributes the computation more evenly (but higher latency)");
	int processChunkSize = args.flag<int>("process-chunk", "process chunk size in samples", -1);
	args.errorExit(); // exits on error, or with `--help`

	std::cout << inputWav << " -> " << outputWav << "\n";

	Wav inWav;
	if (!inWav.read(inputWav).warn()) args.errorExit("failed to read WAV");
	
	size_t inputLength = inWav.length();
	size_t outputLength = std::round(inputLength*time);
	
	Wav outWav;
	outWav.channels = inWav.channels;
	outWav.sampleRate = inWav.sampleRate;
	outWav.resize(outputLength);

	SignalsmithStretch stretch;
	stretch.configure(int(inWav.channels), inWav.sampleRate*0.12, inWav.sampleRate*0.03, splitComputation, asymmetry);
	stretch.setTransposeSemitones(semitones, tonality/inWav.sampleRate);
	stretch.setFormantSemitones(formants, formantComp);
	stretch.setFormantBase(formantBase/inWav.sampleRate);

	signalsmith::plot::Figure figure;
	auto writeLater = figure.writeLater(outputWav + "-blocks.svg");
	size_t plotBlockCounter = 0;
	{
		auto &inputPlot = figure(0, 0).plot(800, 150);
		inputPlot.x.major(0);
		inputPlot.y.major(0);
		auto &outputPlot = figure(0, 1).plot(800, 150);
		outputPlot.x.major(0);
		outputPlot.y.major(0);

		stretch.debugAnalysis = [&](int inputOffset, const float *window, size_t windowOffset, bool isPrevious){
			if (plotBlockCounter > 10) return;
			int blockSamples = stretch.blockSamples();
			auto &line = inputPlot.line(plotBlockCounter);
			for (int i = 0; i < blockSamples; ++i) {
				line.add(inputOffset - blockSamples + i + int(inWav.offset), window[i]);
			}
			line.marker(inputOffset - blockSamples + int(windowOffset) + int(inWav.offset), window[windowOffset]);
		};
		stretch.debugSynthesis = [&](int outputOffset, const float *window, size_t windowOffset){
			if (plotBlockCounter > 10) return;
			int blockSamples = stretch.blockSamples();
			auto &line = outputPlot.line(plotBlockCounter);
			for (int i = 0; i < blockSamples; ++i) {
				line.add(outputOffset + i + int(outWav.offset), window[i]);
			}
			line.marker(outputOffset + windowOffset + int(outWav.offset), window[windowOffset]);

			++plotBlockCounter;
		};
	}

	/* Since the WAV helper allows sample access like `wav[c][index]`, we could just call:
	
		stretch.exact(inWav, int(inputLength), outWav, int(outputLength));
		
	However, we'll do it in separate stages to show more of the API. */
	
	// First, an "output seek", where we provide a chunk of input.
	// This is suitable for starting playback of a sample at a given playback rate.
	auto seekLength = stretch.outputSeekLength(1/time);
	signalsmith::Stopwatch stopwatch;
	stretch.outputSeek(inWav, seekLength);
	double seekTime = stopwatch.seconds(stopwatch.lap());
	// At this point, the next output samples we get will correspond to the beginning of the audio file.

	// We're going to process until *just* before the end of the audio file (so we can get a tidier end using `.flush()`.
	int outputIndex = outputLength - stretch.intervalSamples();

	// Stretch's internal output position is slightly ahead of the output samples we get
	int outputPos = outputIndex + stretch.outputLatency();
	// Time-map: where do we want the input position to be at that moment?
	int inputPos = std::round(outputPos/time);
	// And therefore which input samples do we need to be supplying?
	int inputIndex = inputPos + stretch.inputLatency();
	
	// In this particular case, our `inputPos` will be at the end of the file
	// and `inputIndex` will be beyond the end, so we pad with 0s to have enough input
	inWav.resize(inputIndex);

	// OK, go for it
	inWav.offset = seekLength;
	if (processChunkSize <= 0) {
		stretch.process(inWav, inputIndex - seekLength, outWav, outputIndex);
	} else {
		signalsmith::plot::Plot2D timePlot(500, 200);
		timePlot.x.major(0);
		timePlot.y.major(0);
		timePlot.y.minor(0.01*processChunkSize/inWav.sampleRate, "1%");
		timePlot.y.minor(0.02*processChunkSize/inWav.sampleRate, "2%");
		auto &timeLine = timePlot.line();
		auto &timeLineSeek = timePlot.line().fillToY(0);
		timeLine.add(outWav.offset, 0); // output seek
		timeLineSeek.add(0, 0);
		timeLineSeek.add(0, seekTime);
		timeLineSeek.add(inWav.offset, seekTime);
		timeLineSeek.add(inWav.offset, 0);
	
		float residue = 0.f;
		while (inWav.offset < size_t(inputIndex)) {
			int toProcess = std::min<int>(processChunkSize, inputIndex - inWav.offset);
			float outputPrecise = toProcess * time + residue;
			int outputSamples = std::round(outputPrecise);
			residue = outputPrecise - outputSamples;

			stopwatch.startLap();
			stretch.process(inWav, toProcess, outWav, outputSamples);
			double time = stopwatch.seconds(stopwatch.lap());
			timeLine.add(outWav.offset, time);
			timeLine.add(outWav.offset + toProcess, time);
			
			inWav.offset += toProcess;
			outWav.offset += outputSamples;
		}
		
		timeLine.add(outWav.offset, 0);
		timePlot.write(outputWav + ".svg");
	}
	
	// And as promised, get the last bits using `.flush()`, which does some extra stuff to avoid introducing clicks.
	outWav.offset = outputIndex;
	stretch.flush(outWav, outputLength - outputIndex);
	outWav.offset = 0;

	if (!outWav.write(outputWav).warn()) args.errorExit("failed to write WAV");
}
