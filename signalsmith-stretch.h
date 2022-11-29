#ifndef SIGNALSMITH_STRETCH_H
#define SIGNALSMITH_STRETCH_H

#include "dsp/spectral.h"
#include "dsp/delay.h"
#include "dsp/curves.h"
#include "dsp/perf.h"
SIGNALSMITH_DSP_VERSION_CHECK(1, 3, 3); // Check version is compatible
#include <vector>
#include <algorithm>
#include <functional>

namespace signalsmith { namespace stretch {

template<typename Sample=float>
struct SignalsmithStretch {

	int blockSamples() const {
		return stft.windowSize();
	}
	int intervalSamples() const {
		return stft.interval();
	}
	int inputLatency() const {
		return stft.windowSize()/2;
	}
	int outputLatency() const {
		return stft.windowSize() - inputLatency();
	}
	
	void reset() {
		stft.reset();
		inputBuffer.reset();
		prevInputOffset = -1;
		channelBands.assign(channelBands.size(), Band());
		silenceCounter = 2*stft.windowSize();
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
		stft.resize(channels, blockSamples, intervalSamples);
		inputBuffer.resize(channels, blockSamples);
		timeBuffer.assign(stft.fftSize(), 0);
		channelBands.assign(stft.bands()*channels, Band());
		
		// Various phase rotations
		rotCentreSpectrum.resize(stft.bands());
		rotPrevInput.assign(stft.bands(), 0);
		rotPrevInputShift = -1;
		rotPrevOutput.resize(stft.bands());
		timeShiftPhases(blockSamples*Sample(-0.5), rotCentreSpectrum);
		timeShiftPhases(-intervalSamples, rotPrevOutput);
		peaks.reserve(stft.bands());
		energy.resize(stft.bands());
		smoothedEnergy.resize(stft.bands());
		outputMap.resize(stft.bands());
		channelPredictions.resize(channels*stft.bands());
		maxEnergyChannel.resize(stft.bands());
	}
	
	template<class Inputs, class Outputs>
	void process(Inputs &&inputs, int inputSamples, Outputs &&outputs, int outputSamples) {
		Sample totalEnergy = 0;
		for (int c = 0; c < channels; ++c) {
			auto &&inputChannel = inputs[c];
			for (int i = 0; i < inputSamples; ++i) {
				Sample s = inputChannel[i];
				totalEnergy += s*s;
			}
		}
		if (totalEnergy < noiseFloor) {
			if (silenceCounter >= 2*stft.windowSize()) {
				if (silenceFirst) {
					silenceFirst = false;
					for (auto &b : channelBands) {
						b.input = b.prevInput = b.output = b.prevOutput = 0;
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
				for (int c = 0; c < channels; ++c) {
					auto &&inputChannel = inputs[c];
					auto &&bufferChannel = inputBuffer[c];
					int startIndex = std::max<int>(0, inputSamples - stft.windowSize());
					for (int i = startIndex; i < inputSamples; ++i) {
						bufferChannel[i] = inputChannel[i];
					}
				}
				inputBuffer += inputSamples;
				return;
			} else {
				silenceCounter += inputSamples;
			}
		} else {
			silenceCounter = 0;
			silenceFirst = true;
		}

		for (int outputIndex = 0; outputIndex < outputSamples; ++outputIndex) {
			stft.ensureValid(outputIndex, [&](int outputOffset) {
				// Time to process a spectrum!  Where should it come from in the input?
				int inputOffset = (outputOffset*inputSamples)/outputSamples - stft.windowSize();
				int inputInterval = inputOffset - prevInputOffset;
				prevInputOffset = inputOffset;

				if (inputInterval > 0) {
					if (inputInterval != rotPrevInputShift) { // Only recompute if needed
						timeShiftPhases(-inputInterval, rotPrevInput);
						rotPrevInputShift = inputInterval;
					}
					for (int c = 0; c < channels; ++c) {
						// Copy from the history buffer, if needed
						auto &&bufferChannel = inputBuffer[c];
						for (int i = 0; i < -inputOffset; ++i) {
							timeBuffer[i] = bufferChannel[i + inputOffset];
						}
						// Copy the rest from the input
						auto &&inputChannel = inputs[c];
						for (int i = std::max<int>(0, -inputOffset); i < stft.windowSize(); ++i) {
							timeBuffer[i] = inputChannel[i + inputOffset];
						}
						stft.analyse(c, timeBuffer);
					}
				}
				
				for (int c = 0; c < channels; ++c) {
					auto bands = bandsForChannel(c);
					auto &&spectrumBands = stft.spectrum[c];
					for (int b = 0; b < stft.bands(); ++b) {
						bands[b].input = signalsmith::perf::mul(spectrumBands[b], rotCentreSpectrum[b]);
					}
				}
				
				processSpectrum(inputInterval);

				for (int c = 0; c < channels; ++c) {
					auto bands = bandsForChannel(c);
					auto &&spectrumBands = stft.spectrum[c];
					for (int b = 0; b < stft.bands(); ++b) {
						spectrumBands[b] = signalsmith::perf::mul<true>(bands[b].output, rotCentreSpectrum[b]);
					}
				}
			});

			for (int c = 0; c < channels; ++c) {
				auto &&outputChannel = outputs[c];
				auto &&stftChannel = stft[c];
				outputChannel[outputIndex] = stftChannel[outputIndex];
			}
		}

		// Store input in history buffer
		for (int c = 0; c < channels; ++c) {
			auto &&inputChannel = inputs[c];
			auto &&bufferChannel = inputBuffer[c];
			int startIndex = std::max<int>(0, inputSamples - stft.windowSize());
			for (int i = startIndex; i < inputSamples; ++i) {
				bufferChannel[i] = inputChannel[i];
			}
		}
		inputBuffer += inputSamples;
		stft += outputSamples;
		prevInputOffset -= inputSamples;
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
	
private:
	static constexpr Sample noiseFloor{1e-15};
	int silenceCounter = 0;
	bool silenceFirst = true;

	using Complex = std::complex<Sample>;
	Sample freqMultiplier = 1, freqTonalityLimit = 0.5;
	std::function<Sample(Sample)> customFreqMap = nullptr;

	signalsmith::spectral::STFT<Sample> stft{0, 1, 1};
	signalsmith::delay::MultiBuffer<Sample> inputBuffer;
	int channels = 0;
	int prevInputOffset = -1;
	std::vector<Sample> timeBuffer;

	std::vector<Complex> rotCentreSpectrum, rotPrevOutput, rotPrevInput;
	int rotPrevInputShift = -1;
	Sample bandToFreq(int b) const {
		return (b + Sample(0.5))/stft.fftSize();
	}
	void timeShiftPhases(Sample shiftSamples, std::vector<Complex> &output) const {
		for (int b = 0; b < stft.bands(); ++b) {
			Sample phase = bandToFreq(b)*shiftSamples*Sample(-2*M_PI);
			output[b] = {std::cos(phase), std::sin(phase)};
		}
	}
	
	struct Band {
		Complex input, prevInput{0};
		Complex output, prevOutput{0};
		Sample inputEnergy;
	};
	std::vector<Band> channelBands;
	Band * bandsForChannel(int channel) {
		return channelBands.data() + channel*stft.bands();
	}
	template<Complex Band::*member>
	Complex getBand(int channel, int index) {
		if (index >= stft.bands()) return 0;
		if (index < 0) {
			return std::conj(getBand<member>(channel, -1 - index));
		}
		return channelBands[index + channel*stft.bands()].*member;
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
		Sample fracIndex = inputIndex - std::floor(inputIndex);
		return getFractional<member>(channel, lowIndex, fracIndex);
	}
	template<Sample Band::*member>
	Sample getBand(int channel, int index) {
		if (index < 0) index = -1 - index;
		if (index >= stft.bands()) return 0;
		return channelBands[index + channel*stft.bands()].*member;
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
		Sample fracIndex = inputIndex - std::floor(inputIndex);
		return getFractional<member>(channel, lowIndex, fracIndex);
	}

	struct Peak {
		Sample input, output;
		
		bool operator< (const Peak &other) const {
			return output < other.output;
		}
	};
	std::vector<Peak> peaks;
	std::vector<Sample> energy, smoothedEnergy;
	struct PitchMapPoint {
		Sample inputBin, freqGrad;
	};
	std::vector<PitchMapPoint> outputMap;
	
	struct Prediction {
		Sample energy;
		Complex input;
		Complex freqPrediction;
		Complex shortVerticalTwist, longVerticalTwist;
	};
	std::vector<Prediction> channelPredictions;
	Prediction * predictionsForChannel(int c) {
		return channelPredictions.data() + c*stft.bands();
	}
	std::vector<int> maxEnergyChannel;

	void processSpectrum(int inputInterval) {
		int outputInterval = stft.interval();
		int bands = stft.bands();
		
		Sample rate = outputInterval/std::max<Sample>(1, inputInterval);
		rate = std::min<Sample>(2, rate); // For now, limit the intra-block time stretching to 2x
		
		if (inputInterval > 0) {
			for (int c = 0; c < channels; ++c) {
				auto bins = bandsForChannel(c);
				for (int b = 0; b < stft.bands(); ++b) {
					auto &bin = bins[b];
					bin.prevOutput = signalsmith::perf::mul(bin.prevOutput, rotPrevOutput[b]);
					bin.prevInput = signalsmith::perf::mul(bin.prevInput, rotPrevInput[b]);
				}
			}
		}

		Sample smoothingBins = Sample(stft.fftSize())/stft.interval();
		findPeaks(smoothingBins);
		updateOutputMap(smoothingBins);

		int longVerticalStep = std::round(smoothingBins);
		auto *predictions0 = predictionsForChannel(0);
		for (auto &c : maxEnergyChannel) c = -1;
		for (int c = 0; c < channels; ++c) {
			Band *bins = bandsForChannel(c);
			auto *predictions = predictionsForChannel(c);
			for (int b = 0; b < stft.bands(); ++b) {
				auto &outputBin = bins[b];
				auto mapPoint = outputMap[b];
				int lowIndex = std::floor(mapPoint.inputBin);
				Sample fracIndex = mapPoint.inputBin - std::floor(mapPoint.inputBin);

				Prediction prediction;

				prediction.energy = getFractional<&Band::inputEnergy>(c, lowIndex, fracIndex);
				prediction.energy *= std::max<Sample>(0, mapPoint.freqGrad); // scale the energy according to local stretch factor
				prediction.input = getFractional<&Band::input>(c, lowIndex, fracIndex);
				
				Complex prevInput = getFractional<&Band::prevInput>(c, lowIndex, fracIndex);
				Complex freqTwist = signalsmith::perf::mul<true>(prediction.input, prevInput);
				prediction.freqPrediction = signalsmith::perf::mul(outputBin.prevOutput, freqTwist);

				if (b > 0) {
					Complex downInput = getFractional<&Band::input>(c, mapPoint.inputBin - rate);
					prediction.shortVerticalTwist = signalsmith::perf::mul<true>(prediction.input, downInput);
					if (b > longVerticalStep) {
						Complex longDownInput = getFractional<&Band::input>(c, mapPoint.inputBin - longVerticalStep*rate);
						prediction.longVerticalTwist = signalsmith::perf::mul<true>(prediction.input, longDownInput);
					} else {
						prediction.longVerticalTwist = 0;
					}
				} else {
					prediction.shortVerticalTwist = prediction.longVerticalTwist = 0;
				}

				predictions[b] = prediction;

				// Rough output prediction based on phase-vocoder, sensitive to previous input/output magnitude
				outputBin.output = prediction.freqPrediction/(prediction.energy + noiseFloor);
			}
		}
		for (int b = 0; b < stft.bands(); ++b) {
			// Find maximum-energy channel and calculate that
			int maxChannel = 0;
			Sample maxEnergy = predictions0[b].energy;
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
			auto mapPoint = outputMap[b];

			Complex phase = 0;

			// Short steps
			if (b > 0) {
				auto &otherBin = bins[b - 1];
				phase += signalsmith::perf::mul(otherBin.output, prediction.shortVerticalTwist);
			}
			if (b < stft.bands() - 1) {
				auto &otherBin = bins[b + 1];
				auto &otherPrediction = predictions[b + 1];
				phase += signalsmith::perf::mul<true>(otherBin.output, otherPrediction.shortVerticalTwist);
			}
			// longer verticals
			if (b > longVerticalStep) {
				auto &otherBin = bins[b - longVerticalStep];
				phase += signalsmith::perf::mul(otherBin.output, prediction.longVerticalTwist);
			}
			if (b < stft.bands() - longVerticalStep) {
				auto &otherBin = bins[b + longVerticalStep];
				auto &otherPrediction = predictions[b + longVerticalStep];
				phase += signalsmith::perf::mul<true>(otherBin.output, otherPrediction.longVerticalTwist);
			}

			Sample phaseNorm = std::norm(phase);
			if (phaseNorm > noiseFloor) {
				outputBin.output = phase*std::sqrt(prediction.energy/phaseNorm);
			} else {
				outputBin.output = prediction.input;
			}
			
			// All other bins are locked in phase
			for (int c = 0; c < channels; ++c) {
				if (c != maxChannel) {
					auto &channelBin = bandsForChannel(c)[b];
					auto &channelPrediction = predictionsForChannel(c)[b];
					
					Complex channelTwist = signalsmith::perf::mul<true>(channelPrediction.input, prediction.input);
					Complex channelPhase = signalsmith::perf::mul(outputBin.output, channelTwist);
					
					Sample channelPhaseNorm = std::norm(channelPhase);
					if (channelPhaseNorm > noiseFloor) {
						channelBin.output = channelPhase*std::sqrt(channelPrediction.energy/channelPhaseNorm);
					} else {
						channelBin.output = channelPrediction.input;
					}
				}
			}
		}

		if (inputInterval > 0) {
			for (auto &bin : channelBands) {
				bin.prevOutput = bin.output;
				bin.prevInput = bin.input;
			}
		} else {
			for (auto &bin : channelBands) bin.prevOutput = bin.output;
		}
	}
	
	// Produces smoothed energy across all channels
	void smoothEnergy(Sample smoothingBins) {
		Sample smoothingSlew = 1/(1 + smoothingBins*Sample(0.5));
		for (auto &e : energy) e = 0;
		for (int c = 0; c < channels; ++c) {
			Band *bins = bandsForChannel(c);
			for (int b = 0; b < stft.bands(); ++b) {
				Sample e = std::norm(bins[b].input);
				bins[b].inputEnergy = e; // Used for interpolating prediction energy
				energy[b] += e;
			}
		}
		for (int b = 0; b < stft.bands(); ++b) {
			smoothedEnergy[b] = energy[b];
		}
		Sample e = 0;
		for (int repeat = 0; repeat < 2; ++repeat) {
			for (int b = stft.bands() - 1; b >= 0; --b) {
				e += (smoothedEnergy[b] - e)*smoothingSlew;
				smoothedEnergy[b] = e;
			}
			for (int b = 0; b < stft.bands(); ++b) {
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
		while (start < stft.bands()) {
			if (energy[start] > smoothedEnergy[start]) {
				int end = start + 1;
				while (end < stft.bands() && energy[end] > smoothedEnergy[end]) {
					++end;
				}
				// Take the average frequency and energy across the peak range
				Sample freqSum = 0, energySum = 0;
				for (int b = start; b < end; ++b) {
					Sample e = energy[b];
					freqSum += (b + 0.5)*e;
					energySum += e;
				}
				Sample avgFreq = freqSum/(stft.fftSize()*energySum);
				Sample avgEnergy = energySum/(end - start);
				peaks.emplace_back(Peak{avgFreq*stft.fftSize(), mapFreq(avgFreq)*stft.fftSize()});

				start = end;
			}
			++start;
		}
	}
	
	void updateOutputMap(Sample peakWidthBins) {
		if (peaks.empty()) {
			for (int b = 0; b < stft.bands(); ++b) {
				outputMap[b] = {Sample(b), 1};
			}
			return;
		}
		Sample linearZoneBins = peakWidthBins*Sample(0.5);
		Sample bottomOffset = peaks[0].input - peaks[0].output;
		for (int b = 0; b < std::min<int>(stft.bands(), peaks[0].output); ++b) {
			outputMap[b] = {b + bottomOffset, 1};
		}
		for (size_t p = 1; p < peaks.size(); ++p) {
			const Peak &prev = peaks[p - 1], &next = peaks[p];
			Sample prevEnd = prev.output + linearZoneBins;
			Sample nextStart = next.output - linearZoneBins;
			if (nextStart < prevEnd) nextStart = prevEnd = (nextStart + prevEnd)*Sample(0.5);
			signalsmith::curves::Linear<Sample> segment(prevEnd, nextStart, prev.input + linearZoneBins, next.input - linearZoneBins);
			Sample segmentGrad = ((prev.input + linearZoneBins) - (next.input - linearZoneBins))/(prevEnd - nextStart + noiseFloor);

			prevEnd = std::max<Sample>(0, std::min<Sample>(stft.bands(), prevEnd));
			nextStart = std::max<Sample>(0, std::min<Sample>(stft.bands(), nextStart));

			for (int b = std::max<int>(0, std::ceil(prev.output)); b < prevEnd; ++b) {
				outputMap[b] = {b + prev.input - prev.output, 1};
			}
			for (int b = std::ceil(prevEnd); b < nextStart; ++b) {
				outputMap[b] = {segment(b), segmentGrad};
			}
			for (int b = std::ceil(nextStart); b < std::min<int>(stft.bands(), std::ceil(next.output)); ++b) {
				outputMap[b] = {b + next.input - next.output, 1};
			}
		}
		Sample topOffset = peaks.back().input - peaks.back().output;
		for (int b = std::max<int>(0, peaks.back().output); b < stft.bands(); ++b) {
			outputMap[b] = {b + topOffset, 1};
		}
	}
};

}} // namespace
#endif // include guard
