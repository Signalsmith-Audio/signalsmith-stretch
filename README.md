# Signalsmith Stretch: pitch/time library

This is a C++11 library for pitch and time stretching, using the final approach from the ADC22 presentation _Four Ways To Write A Pitch-Shifter_.

## How to use it


```cpp
#include "signalsmith-stretch.h"

signalsmith::stretch::SignalsmithStretch<float> stretch;
```

### Configuring

The easiest way to configure is `.presetDefault()`:

```cpp
stretch.presetDefault(channels, sampleRate);
```

If you want to test out different block-sizes etc. then you can use `.configure()` manually, and even change `.freqWeight`/`.timeWeight`/`.channelWeight`.

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

### Custom pitch map

This stretcher does (fairly rough) peak-detection, and creates a non-linear frequency map based on that.

You can hook into this to define your own pitch-map, by providing a callback which is called once per channel, for every FFT block:

```cpp
stretch.setMap([&](int channel) {
	for (auto &peak : stretch.peaks) {
		peak.output = peak.input*2; // up one octave
	}
});
```

The input/output frequencies are relative to Nyquist.  It's not currently-tested what happens if your map is non-monotonic.

## Compiling

Just include `signalsmith-stretch.h` in your build.

It's pretty slow if optimisation is disabled though, so you might want to enable optimisation just where it's used.

### DSP Library

This uses the Signalsmith DSP library for FFTs and other bits and bobs.

For convenience, a copy of the license is included (with its own `LICENSE.txt`) in `dsp/`, but if you're already using this elsewhere then you should remove this copy to avoid versioning issues.

## License

[MIT License](LICENSE.txt) for now - get in touch if you need anything else.
