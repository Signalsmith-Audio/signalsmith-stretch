#ifndef SIGNALSMITH_STRETCH_H
#define SIGNALSMITH_STRETCH_H

#include "dsp/spectral.h"
#include "dsp/delay.h"
#include "dsp/curves.h"
#include <vector>
#include <functional>
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
		return stft.windowSize() - inputLatencySamples();
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
		channelWeight = 0.5;
	}

	/// Manual setup
	Sample freqWeight = 1, timeWeight = 2, channelWeight = 0.5;
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
		smoothedEnergy.resize(stft.bands());
		inputBinMap.resize(stft.bands());
		outputGainMap.resize(stft.bands());
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

	struct Peak {
		Sample input, output, energy;
		
		bool operator< (const Peak &other) const {
			return output < other.output;
		}
	};
	std::vector<Peak> peaks;
	/// This function is called once per channel, from inside `.process()`, so that you can alter the mapping in `.peaks`
	void setMap(std::function<void(int)> freqMap) {
		frequencyMapFn = freqMap;
	}
	
private:
	using Complex = std::complex<Sample>;
	Sample freqMultiplier = 1, freqTonalityLimit = 0.5;
	std::function<void(int)> frequencyMapFn;

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
		Complex timeChange{0};
		Sample energy, prevEnergy;
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

	Sample peakThreshold = 1;
	std::vector<Sample> smoothedEnergy, inputBinMap, outputGainMap;

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
		Band *bins0 = bandsForChannel(0);
		for (int c = 0; c < channels; ++c) {
			Band *bins = bandsForChannel(c);

			findPeaks(bins, smoothingBins);
			if (frequencyMapFn) frequencyMapFn(c);
	
			// Scale so they map bins, not frequency
			for (auto &p : peaks) {
				p.input *= stft.fftSize();
				p.output *= stft.fftSize();
			}
			// Create the input/output bin map
			updateBinMap(smoothingBins);
			
			for (int b = 0; b < stft.bands(); ++b) {
				Sample inputIndex = inputBinMap[b];
				int lowIndex = std::floor(inputIndex);
				Sample fracIndex = inputIndex - std::floor(inputIndex);
				Sample outputEnergy = getFractional<&Band::energy>(c, lowIndex, fracIndex);
				
				Band &outputBin = bins[b];
				Complex input = getFractional<&Band::input>(c, lowIndex, fracIndex);
				Complex prevInput = getFractional<&Band::prevInput>(c, lowIndex, fracIndex);
				Complex timeChange = input*std::conj(prevInput);

				Complex prediction = outputBin.prevOutput*timeChange*freqWeight;

				if (b > 0) {
					Sample downIndex = inputIndex - rate;
					int downLowIndex = std::floor(downIndex);
					Sample fracDownIndex = downIndex - std::floor(downIndex);
					Complex downInput = getFractional<&Band::input>(c, downLowIndex, fracDownIndex);
					Complex freqChange = input*std::conj(downInput);
					Complex outputDown = bins[b - 1].output;
					prediction += outputDown*freqChange*timeWeight;
				}
				int longStep = std::round(smoothingBins);
				if (b > longStep) {
					Sample downIndex = inputIndex - longStep*rate;
					int downLowIndex = std::floor(downIndex);
					Sample fracDownIndex = downIndex - std::floor(downIndex);
					Complex downInput = getFractional<&Band::input>(c, downLowIndex, fracDownIndex);
					Complex freqChange = input*std::conj(downInput);
					Complex outputDown = bins[b - longStep].output;
					prediction += outputDown*freqChange*timeWeight;
				}
				
				if (c > 0) {
					Complex ch0Input = getFractional<&Band::input>(0, lowIndex, fracIndex);
					Complex ch0Output = bins0[b].output;
					Complex channelRot = input*std::conj(ch0Input);
					prediction += ch0Output*channelRot*channelWeight;
				}
				
				Sample predictionNorm = std::norm(prediction);
				if (predictionNorm > 1e-15) {
					outputBin.output = prediction*std::sqrt(outputEnergy/predictionNorm);
				} else {
					outputBin.output = input;
				}
				outputBin.output *= outputGainMap[b];
			}
		}

		for (auto &bin : channelBands) {
			bin.prevOutput = bin.output;
			bin.prevInput = bin.input;
			bin.prevEnergy = bin.energy;
		}
	}
	
	void smoothEnergy(Band *bins, Sample smoothingBins) {
		Sample smoothingSlew = 1/(1 + smoothingBins*Sample(0.5));
		for (int b = 0; b < stft.bands(); ++b) {
			auto &bin = bins[b];
			Sample e = std::norm(bin.input);
			bin.energy = e;
			smoothedEnergy[b] = e*peakThreshold;
		}
		Sample e = 0;
		for (int repeat = 0; repeat < 2; ++repeat) {
			for (int b = stft.bands() - 1; b >= 0; --b) {
				auto &bin = bins[b];
				e += (smoothedEnergy[b] - e)*smoothingSlew;
				smoothedEnergy[b] = e;
			}
			for (int b = 0; b < stft.bands(); ++b) {
				auto &bin = bins[b];
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
	
	void findPeaks(Band *bins, Sample smoothingBins) {
		smoothEnergy(bins, smoothingBins);

		peaks.resize(0);
		// Artificial peak at 0
		peaks.emplace_back(Peak{0, 0, 0});
		
		int start = 0;
		while (start < stft.bands()) {
			if (bins[start].energy > smoothedEnergy[start]) {
				int end = start + 1;
				while (end < stft.bands() && bins[end].energy > smoothedEnergy[end]) {
					++end;
				}
				// Take the average frequency and energy across the peak range
				Sample freqSum = 0, energySum = 0;
				for (int b = start; b < end; ++b) {
					Sample e = bins[b].energy;
					freqSum += (b + 0.5)*e;
					energySum += e;
				}
				Sample avgFreq = freqSum/(stft.fftSize()*energySum);
				Sample avgEnergy = energySum/(end - start);
				peaks.emplace_back(Peak{avgFreq, defaultFreqMap(avgFreq), avgEnergy});

				start = end;
			}
			++start;
		}
		// Artificial peak at Nyquist
		peaks.emplace_back(Peak{0.5, defaultFreqMap(freqMultiplier), 0});
	}
	
	void updateBinMap(Sample peakWidthBins) {
		std::stable_sort(peaks.begin(), peaks.end());
		Sample linearZoneBins = peakWidthBins*Sample(0.5);
		for (auto &g : outputGainMap) g = 1; // reset gains
		for (int b = 0; b < std::min<int>(stft.bands(), peaks[0].output); ++b) {
			inputBinMap[b] = peaks[0].input;
			outputGainMap[b] = 0;
		}
		for (size_t p = 1; p < peaks.size(); ++p) {
			const Peak &prev = peaks[p - 1], &next = peaks[p];
			Sample prevEnd = prev.output + linearZoneBins;
			Sample nextStart = next.output - linearZoneBins;
			signalsmith::curves::Linear<Sample> segment(prevEnd, nextStart, prev.input + linearZoneBins, next.input - linearZoneBins);

			if (nextStart < prevEnd) nextStart = prevEnd = (nextStart + prevEnd)*Sample(0.5);
			prevEnd = std::max<Sample>(0, std::min<Sample>(stft.bands(), prevEnd));
			nextStart = std::max<Sample>(0, std::min<Sample>(stft.bands(), nextStart));

			for (int b = std::max<int>(0, std::ceil(prev.output)); b < prevEnd; ++b) {
				inputBinMap[b] = b + prev.input - prev.output;
			}
			for (int b = std::ceil(prevEnd); b < nextStart; ++b) {
				inputBinMap[b] = segment(b);
			}
			for (int b = std::ceil(nextStart); b < std::min<int>(stft.bands(), std::ceil(next.output)); ++b) {
				inputBinMap[b] = b + next.input - next.output;
			}
		}
		for (int b = std::max<int>(0, peaks.back().output); b < stft.bands(); ++b) {
			inputBinMap[b] = peaks.back().input;
			outputGainMap[b] = 0;
		}
	}
};

}} // namespace
#endif // include guard
