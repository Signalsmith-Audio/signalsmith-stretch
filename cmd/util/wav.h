#ifndef RIFF_WAVE_H_
#define RIFF_WAVE_H_

#include <vector>
#include <iostream>
#include <fstream>

// TODO: something better here that doesn't assume little-endian architecture
template<bool littleEndian=true>
struct BigEndian {
	static uint32_t read16(std::istream& in) {
		unsigned char a[2];
		in.read((char*)a, sizeof(a));
		return ((uint32_t)a[0]) + ((uint32_t)a[1])*256;
	}
	static uint32_t read32(std::istream& in) {
		unsigned char a[4];
		in.read((char*)a, sizeof(a));
		return ((uint32_t)a[0]&0xff) + ((uint32_t)a[1])*256 + ((uint32_t)a[2])*65536 + ((uint32_t)a[3])*256*65536;
	}
	
	static void write16(std::ostream& out, uint16_t value) {
		char a[2] = {(char)(value>>0), (char)(value>>8)};
		out.write(a, sizeof(a));
	}
	static void write32(std::ostream& out, uint32_t value) {
		char a[4] = {(char)(value>>0), (char)(value>>8), (char)(value>>16), (char)(value>>24)};
		out.write(a, sizeof(a));
	}
};

class Wav : BigEndian<true> {
	// Little-endian versions of text values
	uint32_t value_RIFF = 0x46464952;
	uint32_t value_WAVE = 0x45564157;
	uint32_t value_fmt = 0x20746d66;
	uint32_t value_data = 0x61746164;

	using BigEndian<true>::read16;
	using BigEndian<true>::read32;
	using BigEndian<true>::write16;
	using BigEndian<true>::write32;

public:
	struct Result {
		enum class Code {
			OK = 0,
			IO_ERROR,
			FORMAT_ERROR,
			UNSUPPORTED,
			WEIRD_CONFIG
		};
		Code code = Code::OK;
		std::string reason;
		
		Result(Code code, std::string reason="") : code(code), reason(reason) {};
		// Used to neatly test for success
		explicit operator bool () const {
			return code == Code::OK;
		};
		const Result & warn(std::ostream& output=std::cerr) const {
			if (!(bool)*this) {
				output << "WAV error: " << reason << std::endl;
			}
			return *this;
		}
	};
	
	size_t sampleRate = 48000;
	size_t channels = 1, offset = 0;
	std::vector<double> samples;
	size_t length() const {
		size_t perChannel = samples.size()/channels;
		return (perChannel >= offset) ? perChannel - offset : 0;
	}
	void resize(size_t length) {
		samples.resize((offset + length)*channels, 0);
	}
	template<bool isConst>
	class ChannelReader {
		using CSample = typename std::conditional<isConst, const double, double>::type;
		CSample *data;
		size_t stride;
	public:
		ChannelReader(CSample *samples, size_t channels) : data(samples), stride(channels) {}
		
		CSample & operator [](size_t i) {
			return data[i*stride];
		}
	};
	ChannelReader<false> operator [](size_t c) {
		return ChannelReader<false>(samples.data() + offset*channels + c, channels);
	}
	ChannelReader<true> operator [](size_t c) const {
		return ChannelReader<true>(samples.data() + offset*channels + c, channels);
	}
	
	Result result = Result(Result::Code::OK);

	Wav() {}
	Wav(double sampleRate, size_t channels) : sampleRate(sampleRate), channels(channels) {}
	Wav(double sampleRate, size_t channels, const std::vector<double> &samples) : sampleRate(sampleRate), channels(channels), samples(samples) {}
	Wav(std::string filename) {
		result = read(filename).warn();
	}
	
	enum class Format {
		PCM=1
	};
	bool formatIsValid(uint16_t format, uint16_t bits) const {
		if (format == (uint16_t)Format::PCM) {
			if (bits == 16) {
				return true;
			}
		}
		return false;
	}
	
	Result read(std::string filename) {
		std::ifstream file;
		file.open(filename, std::ios::binary);
		if (!file.is_open()) return result = Result(Result::Code::IO_ERROR, "Failed to open file: " + filename);

		// RIFF chunk
		if (read32(file) != value_RIFF) return result = Result(Result::Code::FORMAT_ERROR, "Input is not a RIFF file");
		read32(file); // File length - we don't check this
		if (read32(file) != value_WAVE) return result = Result(Result::Code::FORMAT_ERROR, "Input is not a plain WAVE file");
		
		auto blockStart = file.tellg(); // start of the blocks - we will seek back to here periodically
		bool hasFormat = false, hasData = false;
		
		Format format = Format::PCM; // Shouldn't matter, we should always read the `fmt ` chunk before `data`
		while (!file.eof()) {
			auto blockType = read32(file), blockLength = read32(file);
			if (file.eof()) break;
			if (!hasFormat && blockType == value_fmt) {
				auto formatInt = read16(file);
				format = (Format)formatInt;
				channels = read16(file);
				if (channels < 1) return result = Result(Result::Code::FORMAT_ERROR, "Cannot have zero channels");
				
				sampleRate = read32(file);
				if (sampleRate < 1) return result = Result(Result::Code::FORMAT_ERROR, "Cannot have zero sampleRate");

				size_t expectedBytesPerSecond = read32(file);
				size_t bytesPerFrame = read16(file);
				size_t bitsPerSample = read16(file);
				if (!formatIsValid(formatInt, bitsPerSample)) return result = Result(Result::Code::UNSUPPORTED, "Unsupported format:bits: " + std::to_string(formatInt) + ":" + std::to_string(bitsPerSample));
				// Since it's plain WAVE, we can do some extra checks for consistency
				if (bitsPerSample*channels != bytesPerFrame*8) return result = Result(Result::Code::FORMAT_ERROR, "Format sizes don't add up");
				if (expectedBytesPerSecond != sampleRate*bytesPerFrame) return result = Result(Result::Code::FORMAT_ERROR, "Format sizes don't add up");

				hasFormat = true;
				file.clear();
				file.seekg(blockStart);
			} else if (hasFormat && blockType == value_data) {
				std::vector<double> samples(0);
				switch (format) {
				case Format::PCM:
					samples.reserve(blockLength/2);
					for (size_t i = 0; i < blockLength/2; ++i) {
						uint16_t value = read16(file);
						if (file.eof()) break;
						if (value >= 32768) {
							samples.push_back(((double)value - 65536)/32768);
						} else {
							samples.push_back((double)value/32768);
						}
					}
				}
				while (samples.size()%channels != 0) {
					samples.push_back(0);
				}
				this->samples = samples;
				offset = 0;
				hasData = true;
			} else {
				// We either don't recognise
				file.ignore(blockLength);
			}
		}
		if (!hasFormat) return result = Result(Result::Code::FORMAT_ERROR, "missing `fmt ` block");
		if (!hasData) return result = Result(Result::Code::FORMAT_ERROR, "missing `data` block");
		return result = Result(Result::Code::OK);
	}
	
	Result write(std::string filename, Format format=Format::PCM) {
		if (channels == 0 || channels > 65535) return result = Result(Result::Code::WEIRD_CONFIG, "Invalid channel count");
		if (sampleRate <= 0 || sampleRate > 0xFFFFFFFFu) return result = Result(Result::Code::WEIRD_CONFIG, "Invalid sample rate");
		
		std::ofstream file;
		file.open(filename, std::ios::binary);
		if (!file.is_open()) return result = Result(Result::Code::IO_ERROR, "Failed to open file: " + filename);
		
		size_t bytesPerSample;
		switch (format) {
		case Format::PCM:
			bytesPerSample = 2;
			break;
		}
		
		// File size - 44 bytes is RIFF header, "fmt" block, and "data" block header
		size_t dataLength = (samples.size() - offset*channels)*bytesPerSample;
		size_t fileLength = 44 + dataLength;

		// RIFF chunk
		write32(file, value_RIFF);
		write32(file, uint32_t(fileLength - 8)); // File length, excluding the RIFF header
		write32(file, value_WAVE);
		// "fmt " block
		write32(file, value_fmt);
		write32(file, 16); // block length
		write16(file, uint16_t(format));
		write16(file, uint16_t(channels));
		write32(file, uint32_t(sampleRate));
		size_t expectedBytesPerSecond = sampleRate*channels*bytesPerSample;
		write32(file, uint32_t(expectedBytesPerSecond));
		write16(file, uint16_t(channels*bytesPerSample)); // Bytes per frame
		write16(file, uint16_t(bytesPerSample*8)); // bist per sample
		
		// "data" block
		write32(file, value_data);
		write32(file, uint32_t(dataLength));
		switch (format) {
		case Format::PCM:
			for (size_t i = offset*channels; i < samples.size(); i++) {
				double value = samples[i]*32768;
				if (value > 32767) value = 32767;
				if (value <= -32768) value = -32768;
				if (value < 0) value += 65536;
				write16(file, (uint16_t)value);
			}
			break;
		}
		return result = Result(Result::Code::OK);
	}
	
	void makeMono() {
		std::vector<double> newSamples(samples.size()/channels, 0);
		
		for (size_t channel = 0; channel < channels; ++channel) {
			for (size_t i = 0; i < newSamples.size(); ++i) {
				newSamples[i] += samples[i*channels + channel];
			}
		}
		for (size_t i = 0; i < newSamples.size(); ++i) {
			newSamples[i] /= channels;
		}
		
		channels = 1;
		samples = newSamples;
	}
};

#endif // RIFF_WAVE_H_
