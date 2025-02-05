#ifndef SIGNALSMITH_STRETCH_H
#define SIGNALSMITH_STRETCH_H

//#include "dsp/spectral.h"
//#include "dsp/delay.h"
#include "dsp/perf.h"
SIGNALSMITH_DSP_VERSION_CHECK(1, 6, 0); // Check version is compatible

#include "linear/stft.h"
#include <vector>
#include <algorithm>
#include <functional>
#include <random>

namespace signalsmith { namespace stretch {

template<typename Sample=float, class RandomEngine=std::default_random_engine>
struct SignalsmithStretch {
	static constexpr size_t version[3] = {1, 1, 1};

	SignalsmithStretch() : randomEngine(std::random_device{}()) {}
	SignalsmithStretch(long seed) : randomEngine(seed) {}

	int blockSamples() const {
		return stft.blockSamples();
	}
	int intervalSamples() const {
		return stft.defaultInterval();
	}
	int inputLatency() const {
		return stft.blockSamples() - stft.analysisOffset();
	}
	int outputLatency() const {
		return stft.synthesisOffset();
	}
	
	void reset() {
		stft.reset(0.1);
		prevInputOffset = -1;
		channelBands.assign(channelBands.size(), Band());
		silenceCounter = 0;
		didSeek = false;
	}

	// Configures using a default preset
	void presetDefault(int nChannels, Sample sampleRate) {
		configure(nChannels, sampleRate*0.12, sampleRate*0.03);
	}
	void presetCheaper(int nChannels, Sample sampleRate) {
		configure(nChannels, sampleRate*0.1, sampleRate*0.04);
	}

	// Manual setup
	void configure(int nChannels, int blockSamples, int intervalSamples) {
		channels = nChannels;
		stft.configure(channels, channels, blockSamples, intervalSamples + 1);
		stft.setInterval(intervalSamples, stft.kaiser);
		stft.reset(0.1);
		tmpBuffer.resize(blockSamples + intervalSamples);

		bands = stft.bands();
		channelBands.assign(bands*channels, Band());
		
		peaks.reserve(bands/2);
		energy.resize(bands);
		smoothedEnergy.resize(bands);
		outputMap.resize(bands);
		channelPredictions.resize(channels*bands);
	}

	/// Frequency multiplier, and optional tonality limit (as multiple of sample-rate)
	void setTransposeFactor(Sample multiplier, Sample tonalityLimit=0) {
		freqMultiplier = multiplier;
		if (tonalityLimit > 0) {
			freqTonalityLimit = tonalityLimit/std::sqrt(multiplier); // compromise between input and output limits
		} else {
			freqTonalityLimit = 1;
		}
		customFreqMap = nullptr;
	}
	void setTransposeSemitones(Sample semitones, Sample tonalityLimit=0) {
		setTransposeFactor(std::pow(2, semitones/12), tonalityLimit);
		customFreqMap = nullptr;
	}
	// Sets a custom frequency map - should be monotonically increasing
	void setFreqMap(std::function<Sample(Sample)> inputToOutput) {
		customFreqMap = inputToOutput;
	}

	// Provide previous input ("pre-roll"), without affecting the speed calculation.  You should ideally feed it one block-length + one interval
	template<class Inputs>
	void seek(Inputs &&inputs, int inputSamples, double playbackRate) {
		tmpBuffer.resize(0);
		tmpBuffer.resize(stft.blockSamples() + stft.defaultInterval());

		Sample totalEnergy = 0;
		for (int c = 0; c < channels; ++c) {
			auto &&inputChannel = inputs[c];
			int startIndex = std::max<int>(0, inputSamples - int(tmpBuffer.size()));
			for (int i = startIndex; i < inputSamples; ++i) {
				Sample s = inputChannel[i];
				totalEnergy += s*s;
				tmpBuffer[i - startIndex] = s;
			}
			
			stft.writeInput(c, 0, tmpBuffer.size(), tmpBuffer.data());
		}
		if (totalEnergy >= noiseFloor) {
			silenceCounter = 0;
			silenceFirst = true;
		}
		didSeek = true;
		seekTimeFactor = (playbackRate*stft.defaultInterval() > 1) ? 1/playbackRate : stft.defaultInterval();
	}

	template<class Inputs, class Outputs>
	void process(Inputs &&inputs, int inputSamples, Outputs &&outputs, int outputSamples) {
		int prevCopiedInput = 0;
		auto copyInput = [&](int toIndex){
			int length = std::min<int>(stft.blockSamples() + stft.defaultInterval(), toIndex - prevCopiedInput);
			tmpBuffer.resize(length);
			int offset = toIndex - length;
			for (int c = 0; c < channels; ++c) {
				auto &&inputBuffer = inputs[c];
				for (int i = 0; i < length; ++i) {
					tmpBuffer[i] = inputBuffer[i + offset];
				}
				stft.writeInput(c, length, tmpBuffer.data());
			}
			stft.moveInput(toIndex - prevCopiedInput);
			prevCopiedInput = toIndex;
		};

		Sample totalEnergy = 0;
		for (int c = 0; c < channels; ++c) {
			auto &&inputChannel = inputs[c];
			for (int i = 0; i < inputSamples; ++i) {
				Sample s = inputChannel[i];
				totalEnergy += s*s;
			}
		}
		if (totalEnergy < noiseFloor) {
			if (silenceCounter >= 2*stft.blockSamples()) {
				if (silenceFirst) { // first block of silence processing
					silenceFirst = false;
					//stft.reset();
					for (auto &b : channelBands) {
						b.input = b.prevInput = b.output = 0;
						b.inputEnergy = 0;
					}
				}
			
				if (inputSamples > 0) {
					// copy from the input, wrapping around if needed
					for (int outputIndex = 0; outputIndex < outputSamples; ++outputIndex) {
						int inputIndex = outputIndex%inputSamples;
						for (int c = 0; c < channels; ++c) {
							outputs[c][outputIndex] = inputs[c][inputIndex];
						}
					}
				} else {
					for (int c = 0; c < channels; ++c) {
						auto &&outputChannel = outputs[c];
						for (int outputIndex = 0; outputIndex < outputSamples; ++outputIndex) {
							outputChannel[outputIndex] = 0;
						}
					}
				}

				// Store input in history buffer
				copyInput(inputSamples);
				return;
			} else {
				silenceCounter += inputSamples;
			}
		} else {
			silenceCounter = 0;
			silenceFirst = true;
		}
		
		for (int outputIndex = 0; outputIndex < outputSamples; ++outputIndex) {
			if (stft.samplesSinceSynthesis() >= stft.defaultInterval()) {
				// Time to process a spectrum!  Where should it come from in the input?
				int inputOffset = std::round(outputIndex*Sample(inputSamples)/outputSamples);
				int inputInterval = inputOffset - prevInputOffset;
				prevInputOffset = inputOffset;
				
				copyInput(inputOffset);

				bool newSpectrum = didSeek || (inputInterval > 0);
				if (newSpectrum) {
					if (didSeek || inputInterval != int(stft.samplesSinceAnalysis())) { // make sure the previous input is the correct distance in the past
						stft.analyse(stft.defaultInterval());
						// Copy previous analysis to our band objects
						for (int c = 0; c < channels; ++c) {
							auto channelBands = bandsForChannel(c);
							auto *spectrumBands = stft.spectrum(c);
							for (int b = 0; b < bands; ++b) {
								channelBands[b].prevInput = spectrumBands[b];
							}
						}
					}

					stft.analyse();
					// Copy analysed spectrum into our band objects
					for (int c = 0; c < channels; ++c) {
						auto channelBands = bandsForChannel(c);
						auto *spectrumBands = stft.spectrum(c);
						for (int b = 0; b < bands; ++b) {
							channelBands[b].input = spectrumBands[b];
						}
					}
				}
				
				Sample timeFactor = didSeek ? seekTimeFactor : stft.defaultInterval()/std::max<Sample>(1, inputInterval);
				processSpectrum(newSpectrum, timeFactor);
				didSeek = false;

				for (int c = 0; c < channels; ++c) {
					auto channelBands = bandsForChannel(c);
					auto *spectrumBands = stft.spectrum(c);
					for (int b = 0; b < bands; ++b) {
						spectrumBands[b] = channelBands[b].output;
					}
				}
				stft.synthesise();
			};

			for (int c = 0; c < channels; ++c) {
				auto &&outputChannel = outputs[c];
				Sample v = 0;
				stft.readOutput(c, 1, &v);
				outputChannel[outputIndex] = v;
			}
			stft.moveOutput(1);
		}
		
		copyInput(inputSamples);
		prevInputOffset -= inputSamples;
	}

	// Read the remaining output, providing no further input.  `outputSamples` should ideally be at least `.outputLatency()`
	template<class Outputs>
	void flush(Outputs &&outputs, int outputSamples) {
		int plainOutput = std::min<int>(outputSamples, stft.blockSamples());
		int foldedBackOutput = std::min<int>(outputSamples, int(stft.blockSamples()) - plainOutput);
		for (int c = 0; c < channels; ++c) {
			tmpBuffer.resize(plainOutput);
			stft.readOutput(c, plainOutput, tmpBuffer.data());
			auto &&outputChannel = outputs[c];
			for (int i = 0; i < plainOutput; ++i) {
				// TODO: plain output should be gain-
				outputChannel[i] = tmpBuffer[i];
			}
			tmpBuffer.resize(foldedBackOutput);
			stft.readOutput(c, plainOutput, foldedBackOutput, tmpBuffer.data());
			for (int i = 0; i < foldedBackOutput; ++i) {
				outputChannel[outputSamples - 1 - i] -= tmpBuffer[i];
			}
		}
		stft.reset(0.1);

		// Reset the phase-vocoder stuff, so the next block gets a fresh start
		for (int c = 0; c < channels; ++c) {
			auto channelBands = bandsForChannel(c);
			for (int b = 0; b < bands; ++b) {
				channelBands[b].prevInput = channelBands[b].output = 0;
			}
		}
	}
private:
	using Complex = std::complex<Sample>;
	static constexpr Sample noiseFloor{1e-15};
	static constexpr Sample maxCleanStretch{2}; // time-stretch ratio before we start randomising phases
	size_t silenceCounter = 0;
	bool silenceFirst = true;

	Sample freqMultiplier = 1, freqTonalityLimit = 0.5;
	std::function<Sample(Sample)> customFreqMap = nullptr;

	signalsmith::linear::DynamicSTFT<Sample, false, true> stft;
	std::vector<Sample> tmpBuffer;

	int channels = 0, bands = 0;
	int prevInputOffset = -1;
	bool didSeek = false;
	Sample seekTimeFactor = 1;

	Sample bandToFreq(Sample b) const {
		return stft.binToFreq(b);
	}
	Sample freqToBand(Sample f) const {
		return stft.freqToBin(f);
	}
	
	struct Band {
		Complex input, prevInput{0};
		Complex output{0};
		Sample inputEnergy;
	};
	std::vector<Band> channelBands;
	Band * bandsForChannel(int channel) {
		return channelBands.data() + channel*bands;
	}
	template<Complex Band::*member>
	Complex getBand(int channel, int index) {
		if (index < 0 || index >= bands) return 0;
		return channelBands[index + channel*bands].*member;
	}
	template<Complex Band::*member>
	Complex getFractional(int channel, int lowIndex, Sample fractional) {
		Complex low = getBand<member>(channel, lowIndex);
		Complex high = getBand<member>(channel, lowIndex + 1);
		return low + (high - low)*fractional;
	}
	template<Complex Band::*member>
	Complex getFractional(int channel, Sample inputIndex) {
		int lowIndex = std::floor(inputIndex);
		Sample fracIndex = inputIndex - lowIndex;
		return getFractional<member>(channel, lowIndex, fracIndex);
	}
	template<Sample Band::*member>
	Sample getBand(int channel, int index) {
		if (index < 0 || index >= bands) return 0;
		return channelBands[index + channel*bands].*member;
	}
	template<Sample Band::*member>
	Sample getFractional(int channel, int lowIndex, Sample fractional) {
		Sample low = getBand<member>(channel, lowIndex);
		Sample high = getBand<member>(channel, lowIndex + 1);
		return low + (high - low)*fractional;
	}
	template<Sample Band::*member>
	Sample getFractional(int channel, Sample inputIndex) {
		int lowIndex = std::floor(inputIndex);
		Sample fracIndex = inputIndex - lowIndex;
		return getFractional<member>(channel, lowIndex, fracIndex);
	}

	struct Peak {
		Sample input, output;
	};
	std::vector<Peak> peaks;
	std::vector<Sample> energy, smoothedEnergy;
	struct PitchMapPoint {
		Sample inputBin, freqGrad;
	};
	std::vector<PitchMapPoint> outputMap;
	
	struct Prediction {
		Sample energy = 0;
		Complex input;

		Complex makeOutput(Complex phase) {
			Sample phaseNorm = std::norm(phase);
			if (phaseNorm <= noiseFloor) {
				phase = input; // prediction is too weak, fall back to the input
				phaseNorm = std::norm(input) + noiseFloor;
			}
			return phase*std::sqrt(energy/phaseNorm);
		}
	};
	std::vector<Prediction> channelPredictions;
	Prediction * predictionsForChannel(int c) {
		return channelPredictions.data() + c*bands;
	}

	RandomEngine randomEngine;

	void processSpectrum(bool newSpectrum, Sample timeFactor) {
		timeFactor = std::max<Sample>(timeFactor, 1/maxCleanStretch);
		bool randomTimeFactor = (timeFactor > maxCleanStretch);
		std::uniform_real_distribution<Sample> timeFactorDist(maxCleanStretch*2*randomTimeFactor - timeFactor, timeFactor);
		
		if (newSpectrum) {
			for (int c = 0; c < channels; ++c) {
				auto bins = bandsForChannel(c);

				Complex rot = std::polar(Sample(1), bandToFreq(0)*stft.defaultInterval()*Sample(2*M_PI));
				Sample freqStep = bandToFreq(1) - bandToFreq(0);
				Complex rotStep = std::polar(Sample(1), freqStep*stft.defaultInterval()*Sample(2*M_PI));
				
				for (int b = 0; b < bands; ++b) {
					auto &bin = bins[b];
					bin.output = signalsmith::perf::mul(bin.output, rot);
					bin.prevInput = signalsmith::perf::mul(bin.prevInput, rot);
					rot = signalsmith::perf::mul(rot, rotStep);
				}
			}
		}

		Sample smoothingBins = Sample(stft.fftSamples())/stft.defaultInterval();
		int longVerticalStep = std::round(smoothingBins);
		if (customFreqMap || freqMultiplier != 1) {
			findPeaks(smoothingBins);
			updateOutputMap();
		} else { // we're not pitch-shifting, so no need to find peaks etc.
			for (int c = 0; c < channels; ++c) {
				Band *bins = bandsForChannel(c);
				for (int b = 0; b < bands; ++b) {
					bins[b].inputEnergy = std::norm(bins[b].input);
				}
			}
			for (int b = 0; b < bands; ++b) {
				outputMap[b] = {Sample(b), 1};
			}
		}

		// Preliminary output prediction from phase-vocoder
		for (int c = 0; c < channels; ++c) {
			Band *bins = bandsForChannel(c);
			auto *predictions = predictionsForChannel(c);
			for (int b = 0; b < bands; ++b) {
				auto mapPoint = outputMap[b];
				int lowIndex = std::floor(mapPoint.inputBin);
				Sample fracIndex = mapPoint.inputBin - lowIndex;

				Prediction &prediction = predictions[b];
				Sample prevEnergy = prediction.energy;
				prediction.energy = getFractional<&Band::inputEnergy>(c, lowIndex, fracIndex);
				prediction.energy *= std::max<Sample>(0, mapPoint.freqGrad); // scale the energy according to local stretch factor
				prediction.input = getFractional<&Band::input>(c, lowIndex, fracIndex);

				auto &outputBin = bins[b];
				Complex prevInput = getFractional<&Band::prevInput>(c, lowIndex, fracIndex);
				Complex freqTwist = signalsmith::perf::mul<true>(prediction.input, prevInput);
				Complex phase = signalsmith::perf::mul(outputBin.output, freqTwist);
				outputBin.output = phase/(std::max(prevEnergy, prediction.energy) + noiseFloor);
			}
		}

		// Re-predict using phase differences between frequencies
		for (int b = 0; b < bands; ++b) {
			// Find maximum-energy channel and calculate that
			int maxChannel = 0;
			Sample maxEnergy = predictionsForChannel(0)[b].energy;
			for (int c = 1; c < channels; ++c) {
				Sample e = predictionsForChannel(c)[b].energy;
				if (e > maxEnergy) {
					maxChannel = c;
					maxEnergy = e;
				}
			}

			auto *predictions = predictionsForChannel(maxChannel);
			auto &prediction = predictions[b];
			auto *bins = bandsForChannel(maxChannel);
			auto &outputBin = bins[b];

			Complex phase = 0;
			auto mapPoint = outputMap[b];

			// Upwards vertical steps
			if (b > 0) {
				Sample binTimeFactor = randomTimeFactor ? timeFactorDist(randomEngine) : timeFactor;
				Complex downInput = getFractional<&Band::input>(maxChannel, mapPoint.inputBin - binTimeFactor);
				Complex shortVerticalTwist = signalsmith::perf::mul<true>(prediction.input, downInput);

				auto &downBin = bins[b - 1];
				phase += signalsmith::perf::mul(downBin.output, shortVerticalTwist);
				
				if (b >= longVerticalStep) {
					Complex longDownInput = getFractional<&Band::input>(maxChannel, mapPoint.inputBin - longVerticalStep*binTimeFactor);
					Complex longVerticalTwist = signalsmith::perf::mul<true>(prediction.input, longDownInput);

					auto &longDownBin = bins[b - longVerticalStep];
					phase += signalsmith::perf::mul(longDownBin.output, longVerticalTwist);
				}
			}
			// Downwards vertical steps
			if (b < bands - 1) {
				auto &upPrediction = predictions[b + 1];
				auto &upMapPoint = outputMap[b + 1];

				Sample binTimeFactor = randomTimeFactor ? timeFactorDist(randomEngine) : timeFactor;
				Complex downInput = getFractional<&Band::input>(maxChannel, upMapPoint.inputBin - binTimeFactor);
				Complex shortVerticalTwist = signalsmith::perf::mul<true>(upPrediction.input, downInput);

				auto &upBin = bins[b + 1];
				phase += signalsmith::perf::mul<true>(upBin.output, shortVerticalTwist);
				
				if (b < bands - longVerticalStep) {
					auto &longUpPrediction = predictions[b + longVerticalStep];
					auto &longUpMapPoint = outputMap[b + longVerticalStep];

					Complex longDownInput = getFractional<&Band::input>(maxChannel, longUpMapPoint.inputBin - longVerticalStep*binTimeFactor);
					Complex longVerticalTwist = signalsmith::perf::mul<true>(longUpPrediction.input, longDownInput);

					auto &longUpBin = bins[b + longVerticalStep];
					phase += signalsmith::perf::mul<true>(longUpBin.output, longVerticalTwist);
				}
			}

			outputBin.output = prediction.makeOutput(phase);
			
			// All other bins are locked in phase
			for (int c = 0; c < channels; ++c) {
				if (c != maxChannel) {
					auto &channelBin = bandsForChannel(c)[b];
					auto &channelPrediction = predictionsForChannel(c)[b];
					
					Complex channelTwist = signalsmith::perf::mul<true>(channelPrediction.input, prediction.input);
					Complex channelPhase = signalsmith::perf::mul(outputBin.output, channelTwist);
					channelBin.output = channelPrediction.makeOutput(channelPhase);
				}
			}
		}

		if (newSpectrum) {
			for (auto &bin : channelBands) {
				bin.prevInput = bin.input;
			}
		}
	}
	
	// Produces smoothed energy across all channels
	void smoothEnergy(Sample smoothingBins) {
		Sample smoothingSlew = 1/(1 + smoothingBins*Sample(0.5));
		for (auto &e : energy) e = 0;
		for (int c = 0; c < channels; ++c) {
			Band *bins = bandsForChannel(c);
			for (int b = 0; b < bands; ++b) {
				Sample e = std::norm(bins[b].input);
				bins[b].inputEnergy = e; // Used for interpolating prediction energy
				energy[b] += e;
			}
		}
		for (int b = 0; b < bands; ++b) {
			smoothedEnergy[b] = energy[b];
		}
		Sample e = 0;
		for (int repeat = 0; repeat < 2; ++repeat) {
			for (int b = bands - 1; b >= 0; --b) {
				e += (smoothedEnergy[b] - e)*smoothingSlew;
				smoothedEnergy[b] = e;
			}
			for (int b = 0; b < bands; ++b) {
				e += (smoothedEnergy[b] - e)*smoothingSlew;
				smoothedEnergy[b] = e;
			}
		}
	}
	
	Sample mapFreq(Sample freq) const {
		if (customFreqMap) return customFreqMap(freq);
		if (freq > freqTonalityLimit) {
			Sample diff = freq - freqTonalityLimit;
			return freqTonalityLimit*freqMultiplier + diff;
		}
		return freq*freqMultiplier;
	}
	
	// Identifies spectral peaks using energy across all channels
	void findPeaks(Sample smoothingBins) {
		smoothEnergy(smoothingBins);

		peaks.resize(0);
		
		int start = 0;
		while (start < bands) {
			if (energy[start] > smoothedEnergy[start]) {
				int end = start;
				Sample bandSum = 0, energySum = 0;
				while (end < bands && energy[end] > smoothedEnergy[end]) {
					bandSum += end*energy[end];
					energySum += energy[end];
					++end;
				}
				Sample avgBand = bandSum/energySum;
				Sample avgFreq = bandToFreq(avgBand);
				peaks.emplace_back(Peak{avgBand, freqToBand(mapFreq(avgFreq))});

				start = end;
			}
			++start;
		}
	}
	
	void updateOutputMap() {
		if (peaks.empty()) {
			for (int b = 0; b < bands; ++b) {
				outputMap[b] = {Sample(b), 1};
			}
			return;
		}
		Sample bottomOffset = peaks[0].input - peaks[0].output;
		for (int b = 0; b < std::min<int>(bands, std::ceil(peaks[0].output)); ++b) {
			outputMap[b] = {b + bottomOffset, 1};
		}
		// Interpolate between points
		for (size_t p = 1; p < peaks.size(); ++p) {
			const Peak &prev = peaks[p - 1], &next = peaks[p];
			Sample rangeScale = 1/(next.output - prev.output);
			Sample outOffset = prev.input - prev.output;
			Sample outScale = next.input - next.output - prev.input + prev.output;
			Sample gradScale = outScale*rangeScale;
			int startBin = std::max<int>(0, std::ceil(prev.output));
			int endBin = std::min<int>(bands, std::ceil(next.output));
			for (int b = startBin; b < endBin; ++b) {
				Sample r = (b - prev.output)*rangeScale;
				Sample h = r*r*(3 - 2*r);
				Sample outB = b + outOffset + h*outScale;
				
				Sample gradH = 6*r*(1 - r);
				Sample gradB = 1 + gradH*gradScale;
				
				outputMap[b] = {outB, gradB};
			}
		}
		Sample topOffset = peaks.back().input - peaks.back().output;
		for (int b = std::max<int>(0, peaks.back().output); b < bands; ++b) {
			outputMap[b] = {b + topOffset, 1};
		}
	}
};

}} // namespace
#endif // include guard
