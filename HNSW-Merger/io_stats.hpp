// io_stats.hpp
#pragma once
#include <atomic>
#include <cstdint>
#include <cstdio>

struct IoStats {
  std::atomic<unsigned long long> rchar{0};      // bytes read via iostreams
  std::atomic<unsigned long long> wchar{0};      // bytes written via iostreams
  std::atomic<unsigned long long> read_calls{0}; // number of read() calls
  std::atomic<unsigned long long> write_calls{0};
};

inline IoStats& io_stats() {
  static IoStats s;
  return s;
}

// helpers
inline void io_add_read(std::size_t n) {
  if (n) {
    io_stats().rchar.fetch_add((unsigned long long)n, std::memory_order_relaxed);
  }
  io_stats().read_calls.fetch_add(1ULL, std::memory_order_relaxed);
}
inline void io_add_write(std::size_t n) {
  if (n) {
    io_stats().wchar.fetch_add((unsigned long long)n, std::memory_order_relaxed);
  }
  io_stats().write_calls.fetch_add(1ULL, std::memory_order_relaxed);
}

// optional pretty print
inline void io_print_stats(const char* title = "C++ iostream I/O") {
  auto R = io_stats().rchar.load(std::memory_order_relaxed);
  auto W = io_stats().wchar.load(std::memory_order_relaxed);
  auto RC = io_stats().read_calls.load(std::memory_order_relaxed);
  auto WC = io_stats().write_calls.load(std::memory_order_relaxed);

  auto to_gib = [](unsigned long long x)->double { return x / (1024.0*1024.0*1024.0); };
  std::printf("[%s]\n", title);
  std::printf("  Read : %llu bytes (%.3f GiB) in %llu calls\n", R, to_gib(R), RC);
  std::printf("  Write: %llu bytes (%.3f GiB) in %llu calls\n", W, to_gib(W), WC);
}
