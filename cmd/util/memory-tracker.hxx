#include "./memory-tracker.h"

static size_t memoryTrackerAllocCounter = 0;
static size_t memoryTrackerFreeCounter = 0;

signalsmith::MemoryTracker::MemoryTracker() : signalsmith::MemoryTracker::MemoryTracker(memoryTrackerAllocCounter, memoryTrackerFreeCounter) {}

#include <cstdlib>
#include <cstddef>
#include <dlfcn.h>
#include <cassert>
#include <utility>

static void * (*originalCalloc)(size_t, size_t) = nullptr;
static void * (*originalMalloc)(size_t) = nullptr;
static void * (*originalRealloc)(void*, size_t) = nullptr;
static void (*originalFree)(void*) = nullptr;

template<class Fn>
static void cacheOriginal(Fn& fn, const char *symbolName) {
	if (!fn) {
		fn = (Fn)dlsym(RTLD_NEXT, symbolName);
		if (!fn) exit(1);
	}
}

template<class Fn, typename ...Args>
auto callOriginal(Fn& fn, const char *symbolName, Args &&...args)
		-> decltype(fn(std::forward<Args>(args)...)) {
	cacheOriginal(fn, symbolName);
	return fn(std::forward<Args>(args)...);
}

static constexpr size_t EXTRA_SIZE = sizeof(std::max_align_t)*2;
void * storeAllocInfo(void *offsetPointer, void *originalPointer, size_t size) {
	if (!originalPointer) return nullptr;
	memoryTrackerAllocCounter += size;
	
	assert(!((size_t(offsetPointer))%sizeof(size_t))); // make sure it's aligned to size_t
	size_t *sizePtr = (size_t *)offsetPointer;
	sizePtr[-1] = size_t(originalPointer);
	sizePtr[-2] = size;
	return offsetPointer;
}
size_t getAllocSize(void *ptr) {
	assert(!(size_t(ptr)%sizeof(size_t)));
	size_t *sizePtr = (size_t *)ptr;
	return sizePtr[-2];
}
void * getAllocPointer(void *ptr) {
	assert(!(size_t(ptr)%sizeof(size_t)));
	size_t *sizePtr = (size_t *)ptr;
	return (void *)sizePtr[-1];
}

extern "C" {
	void * malloc(size_t size) {
		void *ptr = callOriginal(originalMalloc, "malloc", size + EXTRA_SIZE);
		return storeAllocInfo((unsigned char *)ptr + EXTRA_SIZE, ptr, size);
	}

	void * calloc(size_t size, size_t count) {
		size_t extraCount = (EXTRA_SIZE + size - 1)/size; // enough extra items to store what we need
		void *ptr = callOriginal(originalCalloc, "calloc", size, count + extraCount);
		return storeAllocInfo((unsigned char *)ptr + size*extraCount, ptr, size*count);
	}

	void * realloc(void *ptr, size_t size) {
		void *originalPtr = getAllocPointer(ptr);
		auto pointerOffset = (unsigned char *)ptr - (unsigned char *)originalPtr;
		size_t originalSize = getAllocSize(ptr);
		memoryTrackerFreeCounter += originalSize;
		
		ptr = callOriginal(originalRealloc, "realloc", originalPtr, size + pointerOffset);
		return storeAllocInfo((unsigned char *)ptr + pointerOffset, ptr, size);
	}

	void free(void *ptr) {
		void *originalPtr = getAllocPointer(ptr);
		size_t originalSize = getAllocSize(ptr);
		memoryTrackerFreeCounter += originalSize;

		callOriginal(originalFree, "free", originalPtr);
	}
}

#include <new>

void * operator new(size_t size) {
	return malloc(size);
}

void * operator new[](size_t size) {
	return malloc(size);
}

void operator delete(void *ptr) noexcept {
	free(ptr);
}

void operator delete[](void *ptr) noexcept {
	free(ptr);
}
