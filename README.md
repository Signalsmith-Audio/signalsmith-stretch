# Signalsmith Stretch: C++ pitch/time library

This is a C++11 library for pitch and time stretching, using the final approach from the ADC22 presentation _Four Ways To Write A Pitch-Shifter_.

It can handle a wide-range of pitch-shifts (multiple octaves) but time-stretching sounds best for more modest changes (between 0.75x and 1.5x).  There are some audio examples on the [main project page](https://signalsmith-audio.co.uk/code/stretch/).

## How to use it

```cpp
#include "signalsmith-stretch.h"

signalsmith::stretch::SignalsmithStretch<float> stretch;
```

### Configuring

The easiest way to configure is a `.preset???()` method:

```cpp
stretch.presetDefault(channels, sampleRate);
stretch.presetCheaper(channels, sampleRate);
```

If you want to test out different block-sizes etc. then you can use `.configure()` manually.

### Processing (and resetting)

```cpp
// Clears internal buffers
stretch.reset();

float **inputBuffers, **outputBuffers;
int inputSamples, outputSamples;
stretch.process(inputBuffers, inputSamples, outputBuffers, outputSamples);

// Inspect latency
int totalLatency = stretch.inputLatency() + stretch.outputLatency();
```

The `.process()` method takes anything where `buffer[channel][index]` gives you a sample.  This could be a `float **` or a `double **` or some custom object.

To get a time-stretch, just hand it differently-sized input/output buffers.

### Pitch-shifting

```cpp
stretch.setTransposeFactor(2); // up one octave

stretch.setTransposeSemitones(12); // also one octave
```

You can set a "tonality limit", which uses a non-linear frequency map to preserve a bit more of the timbre:

```cpp
stretch.setTransposeSemitones(4, 8000/sampleRate);
```

Alternatively, you can set a custom frequency map, mapping input frequencies to output frequencies (both normalised against the sample-rate): 

```cpp
stretch.setFreqMap([](float inputFreq) {
	return inputFreq*2; // up one octave
});
```

## Compiling

Just include `signalsmith-stretch.h` in your build.

It's pretty slow if optimisation is disabled though, so you might want to enable optimisation just where it's used.

### DSP Library

This uses the Signalsmith DSP library for FFTs and other bits and bobs.

For convenience, a copy of the license is included (with its own `LICENSE.txt`) in `dsp/`, but if you're already using this elsewhere then you should remove this copy to avoid versioning issues.

## License

[MIT License](LICENSE.txt) for now - get in touch if you need anything else.
