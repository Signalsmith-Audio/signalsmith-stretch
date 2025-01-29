#ifndef SIGNALSMITH_STOPWATCH_UTIL_H
#define SIGNALSMITH_STOPWATCH_UTIL_H

#include <limits>
#include <cmath>
#include <atomic>
#include <algorithm>

// We want CPU time, not wall-clock time, so we can't use `std::chrono::high_resolution_clock`
#ifdef WINDOWS
#	include <windows.h>
namespace signalsmith {
class Stopwatch {
	using Time = __int64;
	inline Time now() {
		LARGE_INTEGER result;
		QueryPerformanceCounter(&result);
		return result.QuadPart;
	}
	static double timeToSeconds(double t) {
		LARGE_INTEGER freq;
		QueryPerformanceFrequency(&freq);
		return t/double(freq);
	}
#else
#	include <ctime>
namespace signalsmith {
class Stopwatch {
	using Time = std::clock_t;
	inline Time now() {
		return std::clock();
	}
	static double timeToSeconds(double t) {
		return t/double(CLOCKS_PER_SEC);
	}
#endif

	std::atomic<Time> lapStart; // the atomic store/load should act as barriers for reordering operations
	Time lapBest, lapTotal, lapTotal2;
	double lapOverhead = 0;
	int lapCount = 0;
	
public:
	Stopwatch(bool compensate=true) {
		if (compensate) {
			start();
			const int repeats = 1000;
			for (int i = 0; i < repeats; ++i) {
				startLap();
				lap();
			}
			lapOverhead = (double)lapTotal/lapCount;
		}
		start();
	}

	static double seconds(double time) {
		return timeToSeconds(time);
	}

	void start() {
		lapCount = 0;
		lapTotal = lapTotal2 = 0;
		lapBest = std::numeric_limits<Time>::max();
		startLap();
	}
	void startLap() {
		lapStart.store(now());
	}
	double lap() {
		auto start = lapStart.load();
		auto diff = now() - start;

		if (diff < lapBest) lapBest = diff;
		lapCount++;
		lapTotal += diff;
		lapTotal2 += diff*diff;

		startLap();
		return diff;
	}
	double total() const {
		return std::max(0.0, lapTotal - lapCount*lapOverhead);
	}
	double mean() const {
		return total()/lapCount;
	}
	double var() const {
		double m = (double)lapTotal/lapCount, m2 = (double)lapTotal2/lapCount;
		return std::max(0.0, m2 - m*m);
	}
	double std() const {
		return sqrt(var());
	}
	double best() const {
		return std::max(0.0, lapBest - lapOverhead);
	}
	double optimistic(double deviations=1) const {
		return std::max(best(), mean() - std()*deviations);
	}
};

} // namespace
#endif // include guard
