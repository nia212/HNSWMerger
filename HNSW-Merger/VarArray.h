#pragma once
#include <vector>
#include <map>
#include <unordered_map>
#include <list>
#include <cstring>
#include <stdexcept>
#include <algorithm>
#if defined(__GLIBC__)
  #include <malloc.h>   // declares int malloc_trim(size_t);
#endif

namespace hnswlib {

class VarArray {
public:
    VarArray() : data_size_(0), size_(0), capacity_(0), max_items_(0) {}

    explicit VarArray(size_t data_size) : data_size_(data_size), size_(0), capacity_(0), max_items_(0) {}

    // Set/Get the maximum number of elements (0 means unlimited)
    void set_max_items(size_t max_items) {
        max_items_ = max_items;
        enforce_limit_after_setting();
    }

    size_t get_max_items() const { return max_items_; }

    size_t cnt = 0;
    // Push single element with external id = `start`
    void push_back(const char *src, unsigned int start) {
        // Ensure space by LRU first so target size respects max_items_
        ensure_space_by_lru(1);

        // After eviction, ensure capacity for the target size (bounded by max_items_)
        ensure_capacity_for_target_size(size_ + 1);

        const size_t pos = size_;
        std::memcpy(data_.data() + pos * data_size_, src, data_size_);

        id_data_map_[start] = static_cast<unsigned int>(pos);
        ensure_index_to_key_capacity();
        index_to_key_[pos] = start;

        ++size_;
        touch_lru(start);
        cnt++;
    }

    // Push a batch: external ids are [start, start+count)
    void push_back(const char *src, unsigned int start, size_t count) {
        if (count == 0)
            return;

        // Ensure space by LRU so that eventual size respects max_items_
        ensure_space_by_lru(count);

        // Compute the actual number we can accept after eviction
        size_t target_size = size_ + count;
        if (max_items_ > 0 && target_size > max_items_)
            target_size = max_items_;
        if (target_size <= size_)
            return;
        const size_t actual = target_size - size_;

        // Ensure capacity for the bounded target size (single reserve path)
        ensure_capacity_for_target_size(target_size);

        // Copy payload
        std::memcpy(data_.data() + size_ * data_size_, src, actual * data_size_);

        ensure_index_to_key_capacity();
        for (size_t i = 0; i < actual; ++i) {
            const unsigned int key = start + static_cast<unsigned int>(i);
            const size_t pos = size_ + i;
            id_data_map_[key] = static_cast<unsigned int>(pos);
            index_to_key_[pos] = key;
            touch_lru(key);
        }
        size_ += actual;
    }

    bool find(unsigned int idx) const { return id_data_map_.find(idx) != id_data_map_.end(); }

    // Non-const access by external id (updates LRU)
    char *operator[](unsigned int idx) {
        auto it = id_data_map_.find(idx);
        if (it == id_data_map_.end())
            throw std::out_of_range("VarArray: key not found");
        touch_lru(idx);
        return data_.data() + static_cast<size_t>(it->second) * data_size_;
    }

    // Const access by external id (updates LRU)
    const char *operator[](unsigned int idx) const {
        auto it = id_data_map_.find(idx);
        if (it == id_data_map_.end())
            throw std::out_of_range("VarArray: key not found");
        touch_lru(idx);
        return data_.data() + static_cast<size_t>(it->second) * data_size_;
    }

    // Non-const access by physical position (updates LRU)
    char *operator()(size_t idx) {
        if (idx >= size_)
            throw std::out_of_range("VarArray: index out of range");
        unsigned int key = index_to_key_[idx];
        touch_lru(key);
        return data_.data() + idx * data_size_;
    }

    // Const access by physical position (updates LRU)
    const char *operator()(size_t idx) const {
        if (idx >= size_)
            throw std::out_of_range("VarArray: index out of range");
        unsigned int key = index_to_key_[idx];
        touch_lru(key);
        return data_.data() + idx * data_size_;
    }

    size_t size() const { return size_; }

    void clear() {
        size_ = 0;
        id_data_map_.clear();
        lru_list_.clear();
        lru_pos_.clear();
        // keep capacity and buffers for reuse
    }

    void release(bool aggressive_trim = false) {
        size_ = 0;
        capacity_ = 0;
        cnt = 0;

        std::vector<char>().swap(data_);

        std::vector<unsigned int>().swap(index_to_key_);

        id_data_map_.clear();
        std::map<unsigned int, unsigned int>().swap(id_data_map_);

        lru_list_.clear();
        std::list<unsigned int>().swap(lru_list_);

        lru_pos_.clear();
        lru_pos_.rehash(0);
        decltype(lru_pos_)().swap(lru_pos_);

#if defined(__GLIBC__)
    if (aggressive_trim) { ::malloc_trim(0); }
#endif
    }

    ~VarArray() { release(true); }

    // Optional: shrink buffers to reclaim memory after setting a tight max_items_
    // By default, trims to max(size_, max_items_) to avoid immediate re-growth.
    void shrink_to_limit() {
        size_t target = (max_items_ > 0) ? std::max(size_, max_items_) : size_;
        if (capacity_ > target) {
            std::vector<char> new_data(target * data_size_);
            if (size_ > 0) {
                std::memcpy(new_data.data(), data_.data(), size_ * data_size_);
            }
            data_.swap(new_data);
            capacity_ = target;

            if (index_to_key_.size() > capacity_)
                index_to_key_.resize(capacity_);
        }
    }

    // private:
    void reserve(size_t new_capacity) {
        std::vector<char> new_data(new_capacity * data_size_);
        if (!data_.empty()) {
            std::memcpy(new_data.data(), data_.data(), size_ * data_size_);
        }
        data_.swap(new_data);
        capacity_ = new_capacity;

        if (index_to_key_.size() < capacity_) {
            index_to_key_.resize(capacity_);
        }
    }

    // Ensure capacity for a given target logical size (after eviction), bounded by max_items_
    void ensure_capacity_for_target_size(size_t target_size) {
        if (max_items_ > 0 && target_size > max_items_)
            target_size = max_items_;
        if (target_size > capacity_) {
            size_t new_cap = (capacity_ == 0) ? 1 : capacity_;
            while (new_cap < target_size)
                new_cap *= 2;
            reserve(new_cap);
        }
    }

    void ensure_index_to_key_capacity() {
        if (index_to_key_.size() < capacity_) {
            index_to_key_.resize(capacity_);
        }
    }

    // Make room according to LRU if max_items_ is set
    void ensure_space_by_lru(size_t count) {
        if (max_items_ == 0)
            return; // unlimited
        // Evict until we can accommodate `count`
        while (size_ + count > max_items_) {
            evict_one_lru();
        }
    }

    void enforce_limit_after_setting() {
        if (max_items_ == 0)
            return;
        while (size_ > max_items_) {
            evict_one_lru();
        }
        // Optionally reclaim memory that was over-reserved historically
        shrink_to_limit();
    }

    // Move key to front of LRU (most recently used)
    void touch_lru(unsigned int key) const {
        auto it = lru_pos_.find(key);
        if (it != lru_pos_.end()) {
            lru_list_.splice(lru_list_.begin(), lru_list_, it->second);
            it->second = lru_list_.begin();
        } else {
            lru_list_.push_front(key);
            lru_pos_[key] = lru_list_.begin();
        }
    }

    // Evict least recently used key
    void evict_one_lru() {
        if (size_ == 0 || lru_list_.empty())
            return;

        const unsigned int victim_key = lru_list_.back();
        lru_list_.pop_back();
        lru_pos_.erase(victim_key);

        auto it = id_data_map_.find(victim_key);
        if (it == id_data_map_.end())
            return; // should not happen

        const size_t victim_pos = it->second;
        id_data_map_.erase(it);

        const size_t last_pos = size_ - 1;
        if (victim_pos != last_pos) {
            // Move bytes
            std::memcpy(data_.data() + victim_pos * data_size_, data_.data() + last_pos * data_size_, data_size_);

            // Update mapping for the moved key
            const unsigned int moved_key = index_to_key_[last_pos];
            index_to_key_[victim_pos] = moved_key;
            id_data_map_[moved_key] = static_cast<unsigned int>(victim_pos);
        }

        // (Optional) clear the last slot key; not strictly necessary
        index_to_key_[last_pos] = 0;

        --size_;
    }

    // private:
    size_t data_size_;
    size_t size_;
    size_t capacity_;
    std::vector<char> data_;

    // external id -> physical position
    std::map<unsigned int, unsigned int> id_data_map_;

    // LRU bookkeeping (mutable to allow updates in const accessors)
    mutable std::list<unsigned int> lru_list_; // front = most recent, back = least recent
    mutable std::unordered_map<unsigned int, typename std::list<unsigned int>::iterator> lru_pos_;

    // physical position -> external id (for fast updates when swapping)
    std::vector<unsigned int> index_to_key_;

    // maximum number of items (0 = unlimited)
    size_t max_items_;
};

} // namespace hnswlib
