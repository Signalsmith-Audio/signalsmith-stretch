/* Currently only working/tested on Mac.  You need to compile in `memory-tracker.cpp` as well, which does the actual stuff */
#ifndef SIGNALSMITH_UTIL_MEMORY_TRACKER_H
#define SIGNALSMITH_UTIL_MEMORY_TRACKER_H

#include <cstddef>

namespace signalsmith {

struct MemoryTracker {
	size_t allocBytes, freeBytes, currentBytes;
	MemoryTracker();
	
	MemoryTracker diff() const {
		MemoryTracker now;
		return {now.allocBytes - allocBytes, now.freeBytes - freeBytes};
	}
private:
	MemoryTracker(size_t allocBytes, size_t freeBytes) : allocBytes(allocBytes), freeBytes(freeBytes), currentBytes(allocBytes - freeBytes) {}
};

} // namespace
#endif // include guard
