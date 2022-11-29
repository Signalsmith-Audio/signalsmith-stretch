#ifndef SIGNALSMITH_STRETCH_H
#define SIGNALSMITH_STRETCH_H

#include "dsp/spectral.h"
#include "dsp/delay.h"
#include "dsp/curves.h"
#include <vector>
#include <algorithm>

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
	}

	/// Configures using a default preset
	void presetDefault(int nChannels, Sample sampleRate) {
		configure(nChannels, sampleRate*0.12, sampleRate*0.03);
		freqWeight = 1;
		timeWeight = 2;
	}

	// manual parameters
	Sample freqWeight = 1, timeWeight = 2, maxWeight = 2;

	/// Manual setup
	void configure(int nChannels, int blockSamples, int intervalSamples) {
		channels = nChannels;
		stft.resize(channels, blockSamples, intervalSamples);
		inputBuffer.resize(channels, blockSamples);
		timeBuffer.assign(stft.fftSize(), 0);
		channelBands.assign(stft.bands()*channels, Band());
		
		// Various phase rotations
		rotCentreSpectrum.resize(stft.bands());
		rotPrevInput.assign(stft.bands(), 0);
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
		Sample timeScaling = Sample(inputSamples)/outputSamples;

		for (int outputIndex = 0; outputIndex < outputSamples; ++outputIndex) {
			stft.ensureValid(outputIndex, [&](int outputOffset) {
				// Time to process a spectrum!  Where should it come from in the input?
				int inputOffset = (outputOffset*inputSamples)/outputSamples - stft.windowSize();
				int inputInterval = inputOffset - prevInputOffset;
				prevInputOffset = inputOffset;

				if (inputInterval > 0) {
					timeShiftPhases(-inputInterval, rotPrevInput);
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
						bands[b].input = spectrumBands[b]*rotCentreSpectrum[b];
					}
				}
				
				processSpectrum(inputInterval);

				for (int c = 0; c < channels; ++c) {
					auto bands = bandsForChannel(c);
					auto &&spectrumBands = stft.spectrum[c];
					for (int b = 0; b < stft.bands(); ++b) {
						spectrumBands[b] = bands[b].output*std::conj(rotCentreSpectrum[b]);
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
	}
	void setTransposeSemitones(Sample semitones, Sample tonalityLimit=0) {
		setTransposeFactor(std::pow(2, semitones/12), tonalityLimit);
	}
	
private:
	using Complex = std::complex<Sample>;
	Sample freqMultiplier = 1, freqTonalityLimit = 0.5;

	signalsmith::spectral::STFT<Sample> stft{0, 1, 1};
	signalsmith::delay::MultiBuffer<Sample> inputBuffer;
	int channels = 0;
	int prevInputOffset = -1;
	std::vector<Sample> timeBuffer;

	std::vector<Complex> rotCentreSpectrum, rotPrevOutput, rotPrevInput;
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
		Sample input, output, energy;
		
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
		
		Sample rate = Sample(inputInterval)/outputInterval;
		
		if (inputInterval > 0) {
			for (int c = 0; c < channels; ++c) {
				auto bins = bandsForChannel(c);
				for (int b = 0; b < stft.bands(); ++b) {
					auto &bin = bins[b];
					bins[b].prevOutput *= rotPrevOutput[b];
					bins[b].prevInput *= rotPrevInput[b];
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
				Complex freqTwist = prediction.input*std::conj(prevInput);
				prediction.freqPrediction = outputBin.prevOutput*freqTwist;

				if (b > 0) {
					Complex downInput = getFractional<&Band::input>(c, mapPoint.inputBin - rate);
					prediction.shortVerticalTwist = prediction.input*std::conj(downInput);
					if (b > longVerticalStep) {
						Complex longDownInput = getFractional<&Band::input>(c, mapPoint.inputBin - longVerticalStep*rate);
						prediction.longVerticalTwist = prediction.input*std::conj(longDownInput);
					} else {
						prediction.longVerticalTwist = 0;
					}
				} else {
					prediction.shortVerticalTwist = prediction.longVerticalTwist = 0;
				}

				predictions[b] = prediction;

				// Rough output prediction based on phase-vocoder, sensitive to previous input/output magnitude
				outputBin.output = prediction.freqPrediction/(prediction.energy + Sample(1e-10));
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

			Complex phase = prediction.freqPrediction*freqWeight;

			// Track the strongest prediction
			Complex maxPrediction = prediction.freqPrediction;
			Sample maxPredictionNorm = std::norm(maxPrediction);

			// Short steps
			if (b > 0) {
				auto &otherBin = bins[b - 1];
				Complex newPrediction = otherBin.output*prediction.shortVerticalTwist;
				phase += newPrediction*timeWeight;
				if (std::norm(newPrediction) > maxPredictionNorm) {
					maxPredictionNorm = std::norm(newPrediction);
					maxPrediction = newPrediction;
				}
			}
			if (b < stft.bands() - 1) {
				auto &otherBin = bins[b + 1];
				auto &otherPrediction = predictions[b + 1];
				Complex newPrediction = otherBin.output*std::conj(otherPrediction.shortVerticalTwist);
				phase += newPrediction*timeWeight;
				if (std::norm(newPrediction) > maxPredictionNorm) {
					maxPredictionNorm = std::norm(newPrediction);
					maxPrediction = newPrediction;
				}
			}
			// longer verticals
			if (b > longVerticalStep) {
				auto &otherBin = bins[b - longVerticalStep];
				Complex newPrediction = otherBin.output*prediction.longVerticalTwist;
				phase += newPrediction*timeWeight;
				if (std::norm(newPrediction) > maxPredictionNorm) {
					maxPredictionNorm = std::norm(newPrediction);
					maxPrediction = newPrediction;
				}
			}
			if (b < stft.bands() - longVerticalStep) {
				auto &otherBin = bins[b + longVerticalStep];
				auto &otherPrediction = predictions[b + longVerticalStep];
				Complex newPrediction = otherBin.output*std::conj(otherPrediction.longVerticalTwist);
				phase += newPrediction*timeWeight;
				if (std::norm(newPrediction) > maxPredictionNorm) {
					maxPredictionNorm = std::norm(newPrediction);
					maxPrediction = newPrediction;
				}
			}

			phase += maxPrediction*maxWeight;
			
			Sample phaseNorm = std::norm(phase);
			if (phaseNorm > 1e-15) {
				outputBin.output = phase*std::sqrt(prediction.energy/phaseNorm);
			} else {
				outputBin.output = prediction.input;
			}
			
			// All other bins are locked in phase
			for (int c = 0; c < channels; ++c) {
				if (c != maxChannel) {
					auto &channelBin = bandsForChannel(c)[b];
					auto &channelPrediction = predictionsForChannel(c)[b];
					
					Complex channelTwist = prediction.input*std::conj(channelPrediction.input);
					Complex channelPhase = outputBin.output*channelTwist;
					
					Sample channelPhaseNorm = std::norm(channelPhase);
					if (channelPhaseNorm > 1e-15) {
						channelBin.output = channelPhase*std::sqrt(prediction.energy/channelPhaseNorm);
					} else {
						channelBin.output = channelPrediction.input;
					}
				}
			}
		}

		for (auto &bin : channelBands) {
			bin.prevOutput = bin.output;
			bin.prevInput = bin.input;
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
	
	Sample defaultFreqMap(Sample freq) const {
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
				peaks.emplace_back(Peak{avgFreq*stft.fftSize(), defaultFreqMap(avgFreq)*stft.fftSize(), avgEnergy});

				start = end;
			}
			++start;
		}
	}
	
	void updateOutputMap(Sample peakWidthBins) {
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
			Sample segmentGrad = ((prev.input + linearZoneBins) - (next.input - linearZoneBins))/(prevEnd - nextStart + Sample(1e-10));

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
