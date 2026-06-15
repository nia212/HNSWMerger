#pragma once
#include "hnswalg.h"
#include "VarArray.h"
#include <mutex>
#include <omp.h>
#include <random>
#include <iostream>
#include <vector>

#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <string>
#include <stdlib.h>
#include <algorithm>
#include <stdio.h>

#include "io_stats.hpp"

using namespace std;

extern float global_counter;

namespace hnswlib {

int read_from_disk(const std::string &location, size_t offset, char *dest, size_t data_size) {
    try {
        std::ifstream input(location, std::ios::binary);
        input.exceptions(std::ifstream::failbit | std::ifstream::badbit);

        input.seekg(offset, std::ios::beg);
        input.read(dest, data_size);

        std::streamsize got = input.gcount();
        if (got > 0)
            io_add_read(static_cast<std::size_t>(got));
        return 0;
    } catch (const std::ios_base::failure &e) {
        return 1;
    }
}

int write_to_disk(const std::string &location, size_t offset, const char *dest, size_t data_size) {
    try {
        std::fstream output(location, std::ios::in | std::ios::out | std::ios::binary);
        output.exceptions(std::fstream::failbit | std::fstream::badbit);

        output.seekp(offset, std::ios::beg);
        output.write(dest, data_size);

        io_add_write(data_size);
        return 0;
    } catch (const std::ios_base::failure &e) {
        return 1;
    }
}

template <typename dist_t>
class HierarchicalNSW_ME : public AlgorithmInterface<dist_t> {
public:
    size_t max_elements_{0};
    mutable std::atomic<size_t> cur_element_count{0};
    size_t size_data_per_element_{0};
    size_t size_dist_per_element_{0};
    size_t size_links_per_element_{0};
    size_t size_dist_links_per_element_{0};
    size_t M_{0};
    size_t maxM_{0};
    size_t maxM0_{0};
    size_t ef_construction_{0};

    double mult_{0.0}, revSize_{0.0};
    int maxlevel_{0};

    std::unique_ptr<VisitedListPool> visited_list_pool_{nullptr};

    tableint enterpoint_node_{0};

    size_t size_links_level0_{0};
    size_t offsetData_{0}, offsetLevel0_{0}, label_offset_{0};

    char *data_level0_memory_{nullptr};
    char *data_memory_{nullptr};
    char **linkLists_{nullptr};
    char *dist_level0_memory_{nullptr};
    char **dist_linkLists_{nullptr};
    std::vector<int> element_levels_;
    VarArray data_;
    size_t level0_in_memory_offset_{0};

    size_t data_size_{0};

    DISTFUNC<dist_t> fstdistfunc_;
    void *dist_func_param_{nullptr};

    // structure: list0, dist0, linklist+distlist (order by id)
    size_t index_offset_list0{0};
    size_t index_offset_dist0{0};
    size_t index_offset_linklist{0};
    size_t index_offset_distlist{0};
    size_t index_link_dist_length{0};
    std::vector<size_t> index_linklist_offset;
    std::vector<size_t> index_distlist_offset;
    std::string file_path;

    HierarchicalNSW_ME(SpaceInterface<dist_t> *s, size_t max_elements, size_t M = 16, size_t ef_construction = 200) :
        element_levels_(max_elements), index_linklist_offset(max_elements + 1), index_distlist_offset(max_elements + 1) {
        max_elements_ = max_elements;
        data_size_ = s->get_data_size();
        fstdistfunc_ = s->get_dist_func();
        dist_func_param_ = s->get_dist_func_param();
        if (M <= 10000) {
            M_ = M;
        } else {
            HNSWERR << "warning: M parameter exceeds 10000 which may lead to adverse effects." << std::endl;
            HNSWERR << "         Cap to 10000 will be applied for the rest of the processing." << std::endl;
            M_ = 10000;
        }
        maxM_ = M_;
        maxM0_ = M_ * 2;
        ef_construction_ = std::max(ef_construction, M_);

        size_links_level0_ = maxM0_ * sizeof(tableint) + sizeof(linklistsizeint);
        size_data_per_element_ = size_links_level0_ + data_size_ + sizeof(labeltype);
        size_dist_per_element_ = maxM0_ * sizeof(dist_t);
        offsetData_ = size_links_level0_;
        label_offset_ = size_links_level0_ + data_size_;
        offsetLevel0_ = 0;
        data_ = VarArray(data_size_);

        cur_element_count = 0;

        visited_list_pool_ = std::unique_ptr<VisitedListPool>(new VisitedListPool(1, max_elements));

        enterpoint_node_ = -1;
        maxlevel_ = -1;
        linkLists_ = (char **)malloc(max_elements_ * sizeof(void *));
        if (linkLists_ == nullptr)
            throw std::runtime_error("Not enough memory: HierarchicalNSW failed to allocate linklists");

        dist_linkLists_ = (char **)malloc(max_elements_ * sizeof(void *));
        if (dist_linkLists_ == nullptr)
            throw std::runtime_error("Not enough memory: HierarchicalNSW failed to allocate dist_linkLists_");

        for (int i = 0; i < max_elements_; i++) {
            linkLists_[i] = nullptr;
            dist_linkLists_[i] = nullptr;
        }

        size_links_per_element_ = maxM_ * sizeof(tableint) + sizeof(linklistsizeint);
        size_dist_links_per_element_ = maxM_ * sizeof(dist_t);
    }

    static inline uint8_t clamp_cast_u8(float x) noexcept {
        if (std::isnan(x))
            return 0;
        if (x < 0.f)
            return 0;
        if (x > 255.f)
            return 255;
        return static_cast<uint8_t>(x); // 向零截断
    }

    // void load_graph_data(const std::string &path,      // file path
    //                      std::uint64_t file_offset,    // offset to level 0 data in file
    //                      std::uint64_t segment_length, // size of one full record
    //                      std::uint64_t data_size,      // size of the slice to extract
    //                      std::uint64_t data_offset,    // offset of the slice within each record
    //                      std::uint64_t count,          // number of records to read
    //                      char *&data_memory            // output buffer
    // ) {
    //     // Preconditions
    //     if (data_offset + data_size > segment_length) {
    //         throw std::invalid_argument("slice range exceeds element size");
    //     }

    //     // Free previous buffer if any and handle empty request
    //     if (data_memory) {
    //         delete[] data_memory;
    //         data_memory = nullptr;
    //     }
    //     if (count == 0)
    //         return;

    //     // Allocate output buffer: contiguous [rec0_slice][rec1_slice]...
    //     const std::uint64_t total_out_bytes = count * data_size;
    //     if (data_size != 0 && total_out_bytes / data_size != count) {
    //         throw std::overflow_error("output size overflow");
    //     }
    //     data_memory = new char[static_cast<std::size_t>(total_out_bytes)];

    //     // Open source file and seek to region start
    //     std::ifstream src(path, std::ios::binary);
    //     if (!src)
    //         throw std::runtime_error("failed to open file: " + path);
    //     src.exceptions(std::ios::badbit); // throw on fatal I/O errors

    //     src.seekg(static_cast<std::streamoff>(file_offset), std::ios::beg);
    //     if (!src)
    //         throw std::runtime_error("seekg failed (offset past EOF?)");

    //     // Choose an input chunk size that is a multiple of record size (8 MiB baseline)
    //     constexpr std::uint64_t kBaseChunk = 8ULL * 1024 * 1024;
    //     std::uint64_t records_per_chunk = std::max<std::uint64_t>(1, kBaseChunk / segment_length);
    //     std::uint64_t in_chunk_bytes = records_per_chunk * segment_length;
    //     if (in_chunk_bytes == 0) { // extremely large record size fallback
    //         records_per_chunk = 1;
    //         in_chunk_bytes = segment_length;
    //     }

    //     std::vector<char> inbuf(static_cast<std::size_t>(in_chunk_bytes));

    //     std::uint64_t remaining_records = count;
    //     std::uint64_t produced = 0; // number of slices written

    //     while (remaining_records > 0) {
    //         const std::uint64_t to_read_records = std::min<std::uint64_t>(remaining_records, records_per_chunk);
    //         const std::uint64_t to_read_bytes = to_read_records * segment_length;

    //         src.read(inbuf.data(), static_cast<std::streamsize>(to_read_bytes));
    //         const std::streamsize got = src.gcount();
    //         if (got == 0) {
    //             throw std::runtime_error("unexpected EOF before completing region");
    //         }

    //         const std::uint64_t full = static_cast<std::uint64_t>(got) / segment_length;
    //         if (full == 0) {
    //             throw std::runtime_error("partial record encountered; region not aligned or truncated");
    //         }

    //         // Copy the requested slice from each full record into data_memory
    //         const char *pin = inbuf.data();
    //         char *pout_base = data_memory + produced * data_size;
    //         for (std::uint64_t i = 0; i < full; ++i) {
    //             const char *src_rec = pin + i * segment_length + data_offset;
    //             char *dst_rec = pout_base + i * data_size;
    //             memcpy(dst_rec, src_rec, static_cast<std::size_t>(data_size));
    //         }

    //         produced += full;
    //         remaining_records -= full;

    //         // If fewer full records than requested were read, ensure we've actually finished
    //         if (full < to_read_records && remaining_records != 0) {
    //             throw std::runtime_error("unexpected EOF within region");
    //         }

    //         // Push back any tail bytes (< one record) so next read starts at record boundary
    //         const std::uint64_t consumed_bytes = full * segment_length;
    //         const std::uint64_t rem_bytes = static_cast<std::uint64_t>(got) - consumed_bytes;
    //         if (rem_bytes > 0) {
    //             src.clear();
    //             src.seekg(-static_cast<std::streamoff>(rem_bytes), std::ios::cur);
    //             if (!src)
    //                 throw std::runtime_error("seekg back failed");
    //         }
    //     }

    //     if (produced != count) {
    //         throw std::runtime_error("produced slice count mismatch");
    //     }
    // }

    void load_graph_data(const std::string &path,      // file path
                         std::uint64_t file_offset,    // offset to level 0 data in file
                         std::uint64_t segment_length, // size of one full record
                         std::uint64_t data_size,      // size of the slice to extract (here = 4 for float)
                         std::uint64_t data_offset,    // offset of the slice within each record
                         std::uint64_t count,          // number of records to read
                         char *&data_memory,           // output buffer
                         bool cast_float_to_u8 = false // NEW: true -> float->uint8_t and output 1B per record
    ) {
        if (data_offset + data_size > segment_length) {
            throw std::invalid_argument("slice range exceeds element size");
        }

        // 释放旧缓冲
        delete[] data_memory;
        data_memory = nullptr;
        if (count == 0)
            return;

        // 计算输出单条大小
        const std::uint64_t out_item_size = cast_float_to_u8 ? data_size / 4 : data_size;
        const std::uint64_t total_out_bytes = count * out_item_size;
        if (out_item_size != 0 && total_out_bytes / out_item_size != count) {
            throw std::overflow_error("output size overflow");
        }
        data_memory = new char[static_cast<std::size_t>(total_out_bytes)];

        // 打开文件
        std::ifstream src(path, std::ios::binary);
        if (!src)
            throw std::runtime_error("failed to open file: " + path);
        src.exceptions(std::ios::badbit);

        src.seekg(static_cast<std::streamoff>(file_offset), std::ios::beg);
        if (!src)
            throw std::runtime_error("seekg failed (offset past EOF?)");

        // 选择chunk
        constexpr std::uint64_t kBaseChunk = 8ULL * 1024 * 1024;
        std::uint64_t records_per_chunk = std::max<std::uint64_t>(1, kBaseChunk / segment_length);
        std::uint64_t in_chunk_bytes = records_per_chunk * segment_length;
        if (in_chunk_bytes == 0) {
            records_per_chunk = 1;
            in_chunk_bytes = segment_length;
        }

        std::vector<char> inbuf(static_cast<std::size_t>(in_chunk_bytes));

        std::uint64_t remaining_records = count;
        std::uint64_t produced = 0;

        while (remaining_records > 0) {
            const std::uint64_t to_read_records = std::min<std::uint64_t>(remaining_records, records_per_chunk);
            const std::uint64_t to_read_bytes = to_read_records * segment_length;

            src.read(inbuf.data(), static_cast<std::streamsize>(to_read_bytes));
            const std::streamsize got = src.gcount();
            if (got == 0)
                throw std::runtime_error("unexpected EOF before completing region");
            if (got > 0)
                io_add_read(static_cast<std::size_t>(got));

            const std::uint64_t full = static_cast<std::uint64_t>(got) / segment_length;
            if (full == 0)
                throw std::runtime_error("partial record encountered; region not aligned or truncated");

            const char *pin = inbuf.data();

            if (cast_float_to_u8) {
                uint8_t *pout = reinterpret_cast<uint8_t *>(data_memory) + produced * (data_size / 4);
                for (std::uint64_t i = 0; i < full; ++i) {
                    for (std::uint64_t j = 0; j < data_size; j += 4) {
                        float f;
                        std::memcpy(&f, pin + i * segment_length + data_offset + j, sizeof(float));
                        pout[i * (data_size / 4) + (j / 4)] = clamp_cast_u8(f);
                    }
                }
            } else {
                // 原样拷贝 data_size 字节/条
                char *pout_base = data_memory + produced * data_size;
                for (std::uint64_t i = 0; i < full; ++i) {
                    const char *src_rec = pin + i * segment_length + data_offset;
                    char *dst_rec = pout_base + i * data_size;
                    std::memcpy(dst_rec, src_rec, static_cast<std::size_t>(data_size));
                }
            }

            produced += full;
            remaining_records -= full;

            if (full < to_read_records && remaining_records != 0) {
                throw std::runtime_error("unexpected EOF within region");
            }

            const std::uint64_t consumed_bytes = full * segment_length;
            const std::uint64_t rem_bytes = static_cast<std::uint64_t>(got) - consumed_bytes;
            if (rem_bytes > 0) {
                src.clear();
                src.seekg(-static_cast<std::streamoff>(rem_bytes), std::ios::cur);
                if (!src)
                    throw std::runtime_error("seekg back failed");
            }
        }

        if (produced != count)
            throw std::runtime_error("produced slice count mismatch");
    }

    void load_graph_data_filter(const std::string &path,      // file path
                                std::uint64_t file_offset,    // offset to level 0 data in file
                                std::uint64_t segment_length, // size of one full record
                                std::uint64_t data_size,      // size of the slice to extract
                                std::uint64_t data_offset,    // offset of the slice within each record
                                vector<int> &cnt              // number of records to read
    ) {
        set<int> id_set = set<int>(cnt.begin(), cnt.end());
        std::ifstream src(path, std::ios::binary);

        src.seekg(static_cast<std::streamoff>(file_offset), std::ios::beg);

        constexpr std::uint64_t kBaseChunk = 8ULL * 1024 * 1024;
        std::uint64_t records_per_chunk = std::max<std::uint64_t>(1, kBaseChunk / segment_length);
        std::uint64_t in_chunk_bytes = records_per_chunk * segment_length;
        if (in_chunk_bytes == 0) {
            records_per_chunk = 1;
            in_chunk_bytes = segment_length;
        }

        std::vector<char> inbuf(static_cast<std::size_t>(in_chunk_bytes));

        std::uint64_t remaining_records = getCurrentElementCount();
        std::uint64_t produced = 0; // number of slices written

        while (remaining_records > 0) {
            const std::uint64_t to_read_records = std::min<std::uint64_t>(remaining_records, records_per_chunk);
            const std::uint64_t to_read_bytes = to_read_records * segment_length;

            src.read(inbuf.data(), static_cast<std::streamsize>(to_read_bytes));
            const std::streamsize got = src.gcount();
            if (got > 0)
                io_add_read(static_cast<std::size_t>(got));
            const std::uint64_t full = static_cast<std::uint64_t>(got) / segment_length;

            const char *pin = inbuf.data();
            for (std::uint64_t i = 0; i < full; ++i) {
                if (id_set.find(i + produced) == id_set.end())
                    continue;

                const char *src_rec = pin + i * segment_length + data_offset;
                data_.push_back(src_rec, i + produced);
            }
            produced += full;
            remaining_records -= full;

            const std::uint64_t consumed_bytes = full * segment_length;
            const std::uint64_t rem_bytes = static_cast<std::uint64_t>(got) - consumed_bytes;
            if (rem_bytes > 0) {
                src.clear();
                src.seekg(-static_cast<std::streamoff>(rem_bytes), std::ios::cur);
            }
        }

        if (cnt.size() != data_.size()) {
            throw std::runtime_error("produced slice count mismatch");
        }
    }

    void clear_linklist(int level) {
        if (level > 0) {
            for (tableint i = 0; i < cur_element_count; i++) {
                if (linkLists_[i] != nullptr) {
                    free(linkLists_[i]);
                    linkLists_[i] = nullptr;
                }
                if (dist_linkLists_[i]) {
                    free(dist_linkLists_[i]);
                    dist_linkLists_[i] = nullptr;
                }
            }
        }
        if (level == 0) {
            if (data_level0_memory_)
                free(data_level0_memory_);
            data_level0_memory_ = nullptr;
            if (dist_level0_memory_)
                free(dist_level0_memory_);
            dist_level0_memory_ = nullptr;
            level0_in_memory_offset_ = 0;
        }
    }

    void clear_distlist(int level) {
        if (level != 0) {
            for (tableint i = 0; i < cur_element_count; i++) {
                if (dist_linkLists_[i] != nullptr) {
                    free(dist_linkLists_[i]);
                    dist_linkLists_[i] = nullptr;
                }
            }
        } else {
            free(dist_level0_memory_);
            dist_level0_memory_ = nullptr;
        }
    }

    void clear_level0() {
        if (linkLists_) {
            free(linkLists_);
            linkLists_ = nullptr;
        }
        if (dist_linkLists_) {
            free(dist_linkLists_);
            dist_linkLists_ = nullptr;
        }
        std::vector<size_t>().swap(index_linklist_offset);
        std::vector<size_t>().swap(index_distlist_offset);
        std::vector<int>().swap(element_levels_);
    }

    void clear() {
        clear_linklist(1);
        // clear_distlist();

        free(linkLists_);
        free(dist_linkLists_);
        linkLists_ = nullptr;
        dist_linkLists_ = nullptr;

        free(data_level0_memory_);
        data_level0_memory_ = nullptr;
        free(dist_level0_memory_);
        dist_level0_memory_ = nullptr;
        free(data_memory_);
        data_memory_ = nullptr;
        data_.release();

        level0_in_memory_offset_ = 0;
        cur_element_count = 0;
        visited_list_pool_.reset(nullptr);
    }

    ~HierarchicalNSW_ME() { clear(); }

    void addPoint(const void *datapoint, labeltype label, bool replace_deleted = false) override {}

    std::priority_queue<std::pair<dist_t, labeltype>> searchKnn(const void *, size_t, BaseFilterFunctor *isIdAllowed = nullptr) const override { return {}; }

    std::vector<std::pair<dist_t, labeltype>> searchKnnCloserFirst(const void *query_data, size_t k, BaseFilterFunctor *isIdAllowed = nullptr) const override { return {}; }

    void saveIndex(const std::string &location) override {}

    size_t getCurrentElementCount() { return cur_element_count; }

    size_t loadMetaData(const std::string &location, SpaceInterface<dist_t> *s, size_t max_elements_i = 0) {
        std::ifstream input(location, std::ios::binary);
        size_t data_size = 0;
        if (!input.is_open())
            throw std::runtime_error("Cannot open file");
        file_path = location;
        input.seekg(0, input.end);
        std::streampos total_filesize = input.tellg();
        input.seekg(0, input.beg);

        readBinaryPOD(input, offsetLevel0_);
        data_size += sizeof(offsetLevel0_);

        readBinaryPOD(input, max_elements_);
        data_size += sizeof(max_elements_);
        readBinaryPOD(input, cur_element_count);
        data_size += sizeof(cur_element_count);

        size_t max_elements = max_elements_i;
        if (max_elements < cur_element_count)
            max_elements = max_elements_;
        max_elements_ = max_elements;
        readBinaryPOD(input, size_data_per_element_);
        data_size += sizeof(size_data_per_element_);
        readBinaryPOD(input, size_dist_per_element_);
        data_size += sizeof(size_dist_per_element_);
        readBinaryPOD(input, label_offset_);
        data_size += sizeof(label_offset_);
        readBinaryPOD(input, offsetData_);
        data_size += sizeof(offsetData_);
        readBinaryPOD(input, maxlevel_);
        data_size += sizeof(maxlevel_);
        readBinaryPOD(input, enterpoint_node_);
        data_size += sizeof(enterpoint_node_);

        readBinaryPOD(input, maxM_);
        data_size += sizeof(maxM_);
        readBinaryPOD(input, maxM0_);
        data_size += sizeof(maxM0_);
        readBinaryPOD(input, M_);
        data_size += sizeof(M_);
        readBinaryPOD(input, mult_);
        data_size += sizeof(mult_);
        readBinaryPOD(input, ef_construction_);
        data_size += sizeof(ef_construction_);

        io_add_read(static_cast<std::size_t>(data_size));

        data_size_ = s->get_data_size();
        fstdistfunc_ = s->get_dist_func();
        dist_func_param_ = s->get_dist_func_param();

        auto pos = input.tellg();

        input.seekg(cur_element_count * size_data_per_element_, input.cur);
        input.seekg(cur_element_count * size_dist_per_element_, input.cur);
        for (size_t i = 0; i < cur_element_count; i++) {
            if (input.tellg() < 0 || input.tellg() >= total_filesize) {
                throw std::runtime_error("Index seems to be corrupted or unsupported");
            }

            unsigned int linkListSize;
            readBinaryPOD(input, linkListSize);
            io_add_read(static_cast<std::size_t>(linkListSize));
            if (linkListSize != 0) {
                input.seekg(linkListSize, input.cur);
            }

            unsigned int distLinkListSize;
            readBinaryPOD(input, distLinkListSize);
            io_add_read(static_cast<std::size_t>(distLinkListSize));
            if (distLinkListSize != 0) {
                input.seekg(distLinkListSize, input.cur);
            }
        }

        if (input.tellg() != total_filesize)
            throw std::runtime_error("Index seems to be corrupted or unsupported");

        input.clear();

        input.seekg(pos, input.beg);
        input.seekg(cur_element_count * size_data_per_element_, input.cur);
        input.seekg(cur_element_count * size_dist_per_element_, input.cur);

        size_links_per_element_ = maxM_ * sizeof(tableint) + sizeof(linklistsizeint);
        size_dist_links_per_element_ = maxM_ * sizeof(dist_t);

        size_links_level0_ = maxM0_ * sizeof(tableint) + sizeof(linklistsizeint);

        element_levels_ = std::vector<int>(max_elements);
        size_t offset = 0;
        for (size_t i = 0; i < cur_element_count; i++) {
            unsigned int linkListSize;

            readBinaryPOD(input, linkListSize);
            io_add_read(static_cast<std::size_t>(linkListSize));
            offset += linkListSize + sizeof(linkListSize);
            if (linkListSize == 0) {
                element_levels_[i] = 0;
            } else {
                element_levels_[i] = linkListSize / size_links_per_element_;
                input.seekg(linkListSize, input.cur);
            }
            linkLists_[i] = nullptr;
            index_distlist_offset[i] = offset;

            unsigned int distLinkListSize;
            readBinaryPOD(input, distLinkListSize);
            io_add_read(static_cast<std::size_t>(distLinkListSize));
            offset += distLinkListSize + sizeof(distLinkListSize);
            if (distLinkListSize > 0) {
                input.seekg(distLinkListSize, input.cur);
            }
            dist_linkLists_[i] = nullptr;
            index_linklist_offset[i + 1] = offset;
        }

        input.close();
        index_offset_list0 = data_size;
        index_offset_dist0 = data_size + cur_element_count * size_data_per_element_;
        index_offset_linklist = index_offset_dist0 + cur_element_count * size_dist_per_element_;
        index_offset_distlist = index_offset_linklist;
        index_link_dist_length = offset;
        return 0;
    }

    size_t get_data_position(tableint internalId) const { return index_offset_list0 + internalId * size_data_per_element_ + offsetData_; }

    size_t get_linklist_level_position(tableint internalId, int layer) const {
        return index_offset_linklist + index_linklist_offset[internalId] + sizeof(unsigned int) + (layer - 1) * size_links_per_element_;
    }

    size_t get_link0_position(tableint internalId) const { return index_offset_list0 + internalId * size_data_per_element_; }

    size_t get_link_position(tableint internal_id, int level) const { return level == 0 ? get_link0_position(internal_id) : get_linklist_level_position(internal_id, level); }

    size_t get_disklink_level_position(tableint internalId, int layer) const {
        return index_offset_distlist + index_distlist_offset[internalId] + sizeof(unsigned int) + (layer - 1) * size_dist_links_per_element_;
    }

    size_t get_dist0_position(tableint internalId) const { return index_offset_dist0 + internalId * size_dist_per_element_; }

    size_t get_dist_position(tableint internal_id, int level) const { return level == 0 ? get_dist0_position(internal_id) : get_disklink_level_position(internal_id, level); }

    void get_data_from_disk_batch(tableint start, int count, char *&dest) {
        load_graph_data(file_path, index_offset_list0 + start * size_data_per_element_, size_data_per_element_, data_size_, offsetData_, count, dest);
    }

    char *get_data(tableint internal_id, int level = 1) {
        if (internal_id > getCurrentElementCount())
            return nullptr;
        if (data_.find(internal_id))
            return data_[internal_id];
        else if (level == 0) {
            // printf("%d, ", internal_id);
            return getDataByInternalId_from_memory(internal_id);
        } else {
            char *temp_data = (char *)malloc(data_size_);
            load_graph_data(file_path, index_offset_list0 + internal_id * size_data_per_element_, size_data_per_element_, data_size_, offsetData_, 1, temp_data);
            data_.push_back(temp_data, internal_id);
            free(temp_data);
            return data_[internal_id];
        }
    }

    inline char *getDataByInternalId_from_memory(tableint internal_id) {
        if (!data_.find(internal_id)) {
            const uint8_t *src = reinterpret_cast<const uint8_t *>(data_memory_ + internal_id * (data_size_ / 4));
            std::vector<char> buf(data_size_);
            for (std::size_t i = 0; i < data_size_ / 4; ++i) {
                float f = static_cast<float>(src[i]);
                std::memcpy(buf.data() + i * sizeof(float), &f, sizeof(float));
            }
            data_.push_back(buf.data(), static_cast<unsigned int>(internal_id));
        }
        return data_[internal_id];
    }

    void get_list_at_level_from_disk_batch(int level, int start, int count, bool with_link, bool with_dist, bool stat_flag = false) {
        size_t stat = 0;
        if (level > 0) {
            for (int i = start; i < start + count; i++) {
                if (element_levels_[i] < level)
                    continue;
                if (with_link) {
                    linkLists_[i] = (char *)malloc(size_links_per_element_);
                    read_from_disk(file_path, get_link_position(i, level), linkLists_[i], size_links_per_element_);
                    stat += size_links_per_element_;
                }
                if (with_dist) {
                    dist_linkLists_[i] = (char *)malloc(size_dist_links_per_element_);
                    read_from_disk(file_path, get_dist_position(i, level), dist_linkLists_[i], size_dist_links_per_element_);
                    stat += size_dist_links_per_element_;
                }
            }
        } else {
            if (with_link) {
                data_level0_memory_ = (char *)malloc(count * size_links_level0_);
                if (data_level0_memory_ == nullptr)
                    throw std::runtime_error("Not enough memory");
                load_graph_data(file_path, index_offset_list0 + start * size_data_per_element_, size_data_per_element_, size_links_level0_, 0, count, data_level0_memory_);
                stat += count * size_links_level0_;
            }

            if (with_dist) {
                dist_level0_memory_ = (char *)malloc(count * size_dist_per_element_);
                if (dist_level0_memory_ == nullptr)
                    throw std::runtime_error("Not enough memory");
                load_graph_data(file_path, index_offset_dist0 + start * size_dist_per_element_, size_dist_per_element_, size_dist_per_element_, 0, count, dist_level0_memory_);
                stat += count * size_dist_per_element_;
            }
            level0_in_memory_offset_ = start;
        }
        if (stat_flag)
            printf("load level %d from disk, start=%d, count=%d, stat=%.2f MB\n", level, start, count, (double)stat / (1024.0 * 1024.0));
    }

    // void get_list0_from_disk_point(int start, bool with_link, bool with_dist) {
    //     size_t stat = 0;
    //     if (with_link) {
    //         linkLists_[start] = nullptr;
    //         load_graph_data(file_path, index_offset_list0 + start * size_data_per_element_, size_data_per_element_, size_links_level0_, 0, 1, linkLists_[start]);
    //     }
    //     if (with_dist) {
    //         dist_linkLists_[start] = nullptr;
    //         load_graph_data(file_path, index_offset_dist0 + start * size_dist_per_element_, size_dist_per_element_, size_dist_per_element_, 0, 1, linkLists_[start]);
    //         stat += size_dist_links_per_element_;
    //     }
    // }

    // linklistsizeint *get_linklist0_point(tableint internal_id) {
    //     if (linkLists_[internal_id] == nullptr)
    //         get_list0_from_disk_point(internal_id, true, false);
    //     return (linklistsizeint *)linkLists_[internal_id];
    // }

    // float *get_dist0_point(tableint internal_id) {
    //     if (dist_linkLists_[internal_id] == nullptr)
    //         get_list0_from_disk_point(internal_id, false, true);
    //     return (float *)dist_linkLists_[internal_id];
    // }

    linklistsizeint *get_linklist0(tableint internal_id) const {
        return (linklistsizeint *)(data_level0_memory_ + (internal_id - level0_in_memory_offset_) * size_links_level0_); // offsetLevel0_ always be 0 here
    }

    linklistsizeint *get_linklist(tableint internal_id) const { return (linklistsizeint *)(linkLists_[internal_id]); }

    linklistsizeint *get_linklist_at_level_from_memory(tableint internal_id, int level) const { return level == 0 ? get_linklist0(internal_id) : get_linklist(internal_id); }

    void get_linklist_at_level_from_disk(tableint internal_id, int level, char *dest, size_t size) const {
        size_t offset = get_link_position(internal_id, level);
        read_from_disk(file_path, offset, dest, size);
    }

    float *get_dist0(tableint internal_id) const { return (float *)(dist_level0_memory_ + (internal_id - level0_in_memory_offset_) * size_dist_per_element_); }

    float *get_dist(tableint internal_id) const { return (float *)(dist_linkLists_[internal_id]); }

    float *get_dist_at_level_from_memory(tableint internal_id, int level) const { return level == 0 ? get_dist0(internal_id) : get_dist(internal_id); }

    void get_dist_at_level_from_disk(tableint internal_id, int level, char *dest, size_t size) const {
        size_t offset = get_dist_position(internal_id, level);
        read_from_disk(file_path, offset, dest, size);
    }

    unsigned short int getListCount(linklistsizeint *ptr) const { return *((unsigned short int *)ptr); }

    void setListCount(linklistsizeint *ptr, unsigned short int size) const { *((unsigned short int *)(ptr)) = *((unsigned short int *)&size); }

    void setLinkList(tableint internal_id, int level, linklistsizeint *source, size_t size) {
        size_t offset = get_link_position(internal_id, level);
        write_to_disk(file_path, offset, reinterpret_cast<char *>(source), size);
    }

    void setDistList(tableint internal_id, int level, dist_t *source, size_t size) {
        size_t offset = get_dist_position(internal_id, level);
        write_to_disk(file_path, offset, reinterpret_cast<char *>(source), size);
    }

    struct CompareByFirst {
        constexpr bool operator()(std::pair<dist_t, tableint> const &a, std::pair<dist_t, tableint> const &b) const noexcept { return a.first < b.first; }
    };

    void getNeighborsByHeuristic2(std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> &top_candidates, const size_t M,
                                  bool collect_metrics = true, tableint cur_c = -1, float alpha = 1.0) {
        if (collect_metrics && top_candidates.size() < M) {
            return;
        }

        std::priority_queue<std::pair<dist_t, tableint>> queue_closest;
        std::vector<std::pair<dist_t, tableint>> return_list;
        while (top_candidates.size() > 0) {
            if (top_candidates.top().second != cur_c)
                queue_closest.emplace(-top_candidates.top().first, top_candidates.top().second);
            top_candidates.pop();
        }
        while (queue_closest.size()) {
            if (return_list.size() >= M)
                break;
            std::pair<dist_t, tableint> curent_pair = queue_closest.top();
            dist_t dist_to_query = -curent_pair.first;
            queue_closest.pop();
            bool good = true;

            for (std::pair<dist_t, tableint> second_pair : return_list) {
                dist_t curdist = fstdistfunc_(get_data(second_pair.second), get_data(curent_pair.second), dist_func_param_);
                if (curdist * alpha < dist_to_query) {
                    good = false;
                    break;
                }
            }
            if (good) {
                return_list.push_back(curent_pair);
            }
        }

        for (std::pair<dist_t, tableint> curent_pair : return_list) {
            top_candidates.emplace(-curent_pair.first, curent_pair.second);
        }
        // delete temp1;
        // delete temp2;
    }

    void searchNodeOnEachLayer(std::vector<std::vector<int>> &resultVector, bool combine = false);

    void search2Layer(const void *query_data, tableint &enterpoint_node, int level_higher, int level_lower, int lambda,
                      std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> *top_candidates = nullptr, int offset = 0);
};

template <typename dist_t>
void HierarchicalNSW_ME<dist_t>::searchNodeOnEachLayer(std::vector<std::vector<int>> &resultVector, bool combine) {
    size_t elementCount = getCurrentElementCount();
    for (size_t iter = 0; iter < elementCount; iter++) {
        if (element_levels_[iter])
            resultVector[element_levels_[iter]].push_back(iter);
    }
    if (combine) {
        for (int i = resultVector.size() - 2; i > 0; --i) {
            resultVector[i].insert(resultVector[i].end(), resultVector[i + 1].begin(), resultVector[i + 1].end());
            sort(resultVector[i].begin(), resultVector[i].end());
        }
    }
}

template <typename dist_t>
void deepCopyOneLayerOnIndex(HierarchicalNSW_ME<dist_t> *index, HierarchicalNSW_ME<dist_t> *alg_hnsw, int level, int offset, std::vector<std::vector<int>> &layer_node_for_index) {
    index->get_list_at_level_from_disk_batch(level, 0, index->getCurrentElementCount(), true, true);
    int link_size = (level == 0) ? alg_hnsw->size_links_level0_ : alg_hnsw->size_links_per_element_;
    int dist_size = (level == 0) ? alg_hnsw->size_dist_per_element_ : alg_hnsw->size_dist_links_per_element_;
    for (int iter = 0; iter < layer_node_for_index[level].size(); iter++) {
        tableint cur_c = layer_node_for_index[level][iter];
        tableint new_c = cur_c + offset;

        linklistsizeint *ll_cur = index->get_linklist_at_level_from_memory(cur_c, level);

        int linklistCount = index->getListCount(ll_cur);
        tableint *data_cur = (tableint *)(ll_cur + 1);
        dist_t *ll_cur_dist = index->get_dist_at_level_from_memory(cur_c, level);
        for (int i = 0; i < linklistCount; i++) {
            data_cur[i] = data_cur[i] + offset;
        }
        alg_hnsw->setLinkList(new_c, level, data_cur, link_size);
        alg_hnsw->setDistList(new_c, level, ll_cur_dist, dist_size);
    }
    index->clear_linklist(level);
}
/*
 * This is the specific function for ME version of algorithm.
 *
 * This function, which is used in search2Layer, searching the graph from the
 * top layer to the level_low layer. If the enterpoint_node is -1, it will be
 * set to enterpoint_node_. Otherwise, the enterpoint_node will be used as the
 * starting point, as it's the closest point to the query on layer level_low
 * + 1. The function will also return the `top_candidates` which contains the
 * `lambda`-closest points to the query on layer level_low.
 */
template <typename dist_t>
void HierarchicalNSW_ME<dist_t>::search2Layer(const void *query_data, tableint &enterpoint_node, int level_higher, int level_lower, int lambda,
                                              std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> *top_candidates,
                                              int offset) {
    tableint currObj = enterpoint_node == -1 ? enterpoint_node_ : enterpoint_node;
    dist_t curdist = fstdistfunc_(query_data, get_data(currObj, level_higher), dist_func_param_);
    for (int level = level_higher; level > level_lower; level--) {
        bool changed = true;
        while (changed) {
            changed = false;
            unsigned int *data = (unsigned int *)get_linklist_at_level_from_memory(currObj, level);
            int size = getListCount(data);

            tableint *datal = (tableint *)(data + 1);

            for (int i = 0; i < size; i++) {
                tableint cand = datal[i];
                if (cand < 0 || cand > max_elements_) {
                    throw std::runtime_error("cand error");
                }
                dist_t d = fstdistfunc_(query_data, get_data(cand, level_higher), dist_func_param_);

                if (d < curdist) {
                    curdist = d;
                    currObj = cand;
                    changed = true;
                }
            }
        }
    }
    if (top_candidates == nullptr) {
        enterpoint_node = currObj;
        return;
    }

    VisitedList *vl = visited_list_pool_->getFreeVisitedList();
    vl_type *visited_array = vl->mass;
    vl_type visited_array_tag = vl->curV;

    std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> candidateSet;

    dist_t lowerBound;
    dist_t entry_dist;
    dist_t dist = fstdistfunc_(query_data, get_data(currObj, level_lower), dist_func_param_);
    top_candidates->emplace(dist, currObj + offset);
    enterpoint_node = currObj;
    entry_dist = dist;
    lowerBound = dist;
    candidateSet.emplace(-dist, currObj);
    visited_array[currObj] = visited_array_tag;

    while (!candidateSet.empty()) {
        std::pair<dist_t, tableint> curr_el_pair = candidateSet.top();
        if ((-curr_el_pair.first) > lowerBound && top_candidates->size() == lambda) {
            break;
        }
        candidateSet.pop();

        tableint curNodeNum = curr_el_pair.second;

        int *data = (int *)get_linklist_at_level_from_memory(curNodeNum, level_lower);
        size_t size = getListCount((linklistsizeint *)data);
        tableint *datal = (tableint *)(data + 1);
#ifdef USE_SSE
        _mm_prefetch((char *)(visited_array + *(data + 1)), _MM_HINT_T0);
        _mm_prefetch((char *)(visited_array + *(data + 1) + 64), _MM_HINT_T0);
        _mm_prefetch(get_data(*datal, level_lower), _MM_HINT_T0);
        _mm_prefetch(get_data(*(datal + 1), level_lower), _MM_HINT_T0);
#endif

        for (size_t j = 0; j < size; j++) {
            tableint candidate_id = *(datal + j);
#ifdef USE_SSE
            _mm_prefetch((char *)(visited_array + *(datal + j + 1)), _MM_HINT_T0);
            _mm_prefetch(get_data(*(datal + j + 1), level_lower), _MM_HINT_T0);
#endif
            if (visited_array[candidate_id] == visited_array_tag)
                continue;
            visited_array[candidate_id] = visited_array_tag;

            char *currObj1 = get_data(candidate_id, level_lower);

            dist_t dist1 = fstdistfunc_(query_data, currObj1, dist_func_param_);
            if (top_candidates->size() < lambda || lowerBound > dist1) {
                candidateSet.emplace(-dist1, candidate_id);
#ifdef USE_SSE
                _mm_prefetch(get_data(candidateSet.top().second, level_lower), _MM_HINT_T0);
#endif

                top_candidates->emplace(dist1, candidate_id + offset);
                if (entry_dist > dist1) {
                    enterpoint_node = candidate_id;
                    entry_dist = dist1;
                }

                while (top_candidates->size() > lambda)
                    top_candidates->pop();

                if (!top_candidates->empty())
                    lowerBound = top_candidates->top().first;
            }
        }
    }
    visited_list_pool_->releaseVisitedList(vl);
}

void copy_file_region(const std::string &src_path, size_t offset, size_t length, std::ofstream &dst) {
    size_t CHUNK_SIZE = 1ULL * 256 * 1024 * 1024; // 64MB per chunk

    std::ifstream src(src_path, std::ios::binary);
    if (!src) {
        throw std::runtime_error("Failed to open input file: " + src_path);
    }

    // 跳到 offset
    src.seekg(offset, std::ios::beg);
    if (!src) {
        throw std::runtime_error("seekg failed");
    }

    std::vector<char> buffer(std::min(CHUNK_SIZE, length));
    size_t remaining = length;

    while (remaining > 0) {
        size_t chunk = std::min(buffer.size(), remaining);

        src.read(buffer.data(), chunk);
        if (!src) {
            throw std::runtime_error("read failed");
        }
        io_add_read(static_cast<std::size_t>(chunk));

        dst.write(buffer.data(), chunk);
        if (!dst) {
            throw std::runtime_error("write failed");
        }
        io_add_write(static_cast<std::size_t>(chunk));

        remaining -= chunk;
        // printf("remaining: %ld\n", remaining);
    }
}

template <typename dist_t>
void HNSWMerger_ME(const std::string &location_index1, const std::string &location_index2, const std::string &location_merged, L2Space *space, size_t max_elements, size_t M = -1,
                   size_t ef_construction = -1, int lambda = 4, float alpha = 1.05) {
    size_t GB_size = 1 << 30;
    size_t memory_size = 10 * GB_size;

    double s0 = elapsed();
    HierarchicalNSW_ME<dist_t> *index1 = new hnswlib::HierarchicalNSW_ME<dist_t>(space, max_elements, M, ef_construction);
    index1->loadMetaData(location_index1, space);
    HierarchicalNSW_ME<dist_t> *index2 = new hnswlib::HierarchicalNSW_ME<dist_t>(space, max_elements, M, ef_construction);
    index2->loadMetaData(location_index2, space);
    printf("Load index time: %.3f s\n", elapsed() - s0);

    // HierarchicalNSW<dist_t> *index5 = new hnswlib::HierarchicalNSW<dist_t>(space, max_elements, M, ef_construction);
    // index5->loadIndex(location_index1, space);
    // uint8_t *data = (uint8_t *)(index5->getDataByInternalId(0));
    // for (int i = 0; i < 10; i++) {
    //     printf("(");
    //     for (int j = 0; j < 4; j++)
    //         printf("%d, ", data[i*4+j]);
    //     printf(")");
    // }
    // printf("\n");
    // float *data2 = (float *)(index5->getDataByInternalId(0));
    // for (int i = 0; i < 10; i++)
    //     printf("%u, ", static_cast<uint8_t>(data2[i]));
    // printf("\n");
    // printf("%d\n",index5->data_size_);
    // exit(0);

    if (index1->cur_element_count > index2->cur_element_count || index1->maxlevel_ > index2->maxlevel_) {
        std::swap(index1, index2);
        printf("Swap index with more current element count to index1.\n");
    }

    size_t element_count_for_index1 = index1->getCurrentElementCount();
    size_t element_count_for_index2 = index2->getCurrentElementCount();
    size_t index1_offset = 0;
    size_t index2_offset = element_count_for_index1;

    M = M == -1 ? std::max(index1->M_, index2->M_) : M;
    ef_construction = ef_construction == -1 ? std::max(index1->ef_construction_, index2->ef_construction_) : ef_construction;
    size_t new_max_elements = index1->max_elements_ + index2->max_elements_;
    size_t maxLevel = std::max(index1->maxlevel_, index2->maxlevel_);

    // Create merged index and write original data into merged index
    {
        HierarchicalNSW<dist_t> *alg_hnsw = new HierarchicalNSW<dist_t>(space, new_max_elements, M, ef_construction);
        alg_hnsw->cur_element_count.store(index1->cur_element_count + index2->cur_element_count);
        if (index2->maxlevel_ < index1->maxlevel_) {
            throw std::logic_error("index1 is higher than index2, wrong case");
        }
        maxLevel = index2->maxlevel_;
        alg_hnsw->setMaxLevel(maxLevel);
        alg_hnsw->enterpoint_node_ = index2->enterpoint_node_ + index2_offset;

        std::ofstream output(location_merged, std::ios::binary);
        std::streampos position;

        writeBinaryPOD(output, alg_hnsw->offsetLevel0_);
        writeBinaryPOD(output, alg_hnsw->max_elements_);
        writeBinaryPOD(output, alg_hnsw->cur_element_count);
        writeBinaryPOD(output, alg_hnsw->size_data_per_element_);
        writeBinaryPOD(output, alg_hnsw->size_dist_per_element_);
        writeBinaryPOD(output, alg_hnsw->label_offset_);
        writeBinaryPOD(output, alg_hnsw->offsetData_);
        writeBinaryPOD(output, alg_hnsw->maxlevel_);
        writeBinaryPOD(output, alg_hnsw->enterpoint_node_);
        writeBinaryPOD(output, alg_hnsw->maxM_);
        writeBinaryPOD(output, alg_hnsw->maxM0_);
        writeBinaryPOD(output, alg_hnsw->M_);
        writeBinaryPOD(output, alg_hnsw->mult_);
        writeBinaryPOD(output, alg_hnsw->ef_construction_);

        copy_file_region(location_index1, index1->index_offset_list0, element_count_for_index1 * alg_hnsw->size_data_per_element_, output);
        copy_file_region(location_index2, index2->index_offset_list0, element_count_for_index2 * alg_hnsw->size_data_per_element_, output);
        copy_file_region(location_index1, index1->index_offset_dist0, element_count_for_index1 * alg_hnsw->size_dist_per_element_, output);
        copy_file_region(location_index2, index2->index_offset_dist0, element_count_for_index2 * alg_hnsw->size_dist_per_element_, output);
        printf("ckpt1, time: %.3f s\n", elapsed() - s0);
        copy_file_region(location_index1, index1->index_offset_linklist, index1->index_link_dist_length, output);
        copy_file_region(location_index2, index2->index_offset_linklist, index2->index_link_dist_length, output);
        delete alg_hnsw;
    }
    printf("ckpt2, time: %.3f s\n", elapsed() - s0);
    HierarchicalNSW_ME<dist_t> *alg_hnsw = new hnswlib::HierarchicalNSW_ME<float>(space, max_elements, M, ef_construction);
    alg_hnsw->loadMetaData(location_merged, space);

    tableint *entry_point_collect_index1_on_index2 = new tableint[element_count_for_index1];
    memset(entry_point_collect_index1_on_index2, -1, element_count_for_index1 * sizeof(int));

    std::vector<std::vector<int>> layer_node_for_index1(maxLevel + 1);
    index2->searchNodeOnEachLayer(layer_node_for_index1, true);
    std::vector<std::vector<int>> layer_node_for_index2(maxLevel + 1);
    index2->searchNodeOnEachLayer(layer_node_for_index2, true);
    printf("ckpt3, time: %.3f s\n", elapsed() - s0);

    size_t alg_hnsw_maxM_ = M;
    size_t alg_hnsw_maxM0_ = M * 2;
    printf("initial ends, time: %.3f s\n", elapsed() - s0);

    StopW sw2 = StopW();
    StopW sw3 = StopW();
    StopW sw1 = StopW();

    double sum0 = 0.0, sum1 = 0.0;
    double sum2 = 0.0;
    char *data_point;

    //! This is for level > 0
    for (int level = maxLevel; level > 0; level -= 1) {
        printf("Merging layer %d\n", level);

        if (layer_node_for_index1[level].size() > 0 && layer_node_for_index2[level].size() == 0) {
            throw std::logic_error("index1 is higher than index2, wrong case");
        }
        if (layer_node_for_index1[level].size() == 0 && layer_node_for_index2[level].size() > 0) {
            deepCopyOneLayerOnIndex(index2, alg_hnsw, level, index2_offset, layer_node_for_index2);
            continue;
        }

        index1->data_.set_max_items(1 << 19);
        index2->data_.set_max_items(524288);

        index2->load_graph_data_filter(index2->file_path, index2->index_offset_list0, index2->size_data_per_element_, index2->data_size_, index2->offsetData_,
                                       layer_node_for_index2[level]);

        index2->get_list_at_level_from_disk_batch(level, 0, index2->cur_element_count, true, false, true); //

        sw2.reset();
        sw3.reset();

        size_t Mcurmax = level ? alg_hnsw_maxM_ : alg_hnsw_maxM0_;
        std::vector<std::vector<std::pair<dist_t, tableint>>> candidateSetIndex2(element_count_for_index2);

        size_t ll_cur_size = (level == 0) ? alg_hnsw->size_links_level0_ : alg_hnsw->size_links_per_element_;
        linklistsizeint *ll_cur;
        linklistsizeint *ll_new = (linklistsizeint *)malloc(ll_cur_size);

        size_t dist_size = (level == 0) ? alg_hnsw->size_dist_per_element_ : alg_hnsw->size_dist_links_per_element_;
        dist_t *dist;
        dist_t *dist_new = (dist_t *)malloc(dist_size);

        int batch_size = 4096;
        char *batch_data = nullptr;
        int cnt = 0;
        for (int b = 0; b < index1->cur_element_count; b += batch_size) {
            index1->get_list_at_level_from_disk_batch(level, b, batch_size, true, true); //
            index1->get_data_from_disk_batch(b, batch_size, batch_data);
            index1->data_.push_back(batch_data, b, batch_size);

            for (int i = b; i < b + batch_size && i < index1->cur_element_count; ++i) {
                if (i % 200000 == 0 || (level == 0 && i == b)) {
                    printf("log[%d] time %f, cnt = %d\n", i, sw2.getElapsedTimeMicro() / 1e6, cnt);
                    printf("count1: %d, count2: %d, count_m: %d\n", index1->data_.size(), index2->data_.size(), alg_hnsw->data_.size());
                }
                tableint cur_c = i;
                tableint new_c = cur_c + index1_offset;

                data_point = index1->data_[cur_c];

                tableint &last_entry_point = entry_point_collect_index1_on_index2[cur_c];
                if (index1->element_levels_[i] < level) {
                    index2->search2Layer(data_point, last_entry_point, level, level - 1, lambda, nullptr, index2_offset);
                    continue;
                }
                cnt++;
                std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, typename HierarchicalNSW_ME<dist_t>::CompareByFirst> top_candidates;

                index2->search2Layer(data_point, last_entry_point, level, level, lambda, &top_candidates, index2_offset);

                auto temp_queue = top_candidates;
                while (!temp_queue.empty()) {
                    std::pair<dist_t, tableint> current = temp_queue.top();
                    temp_queue.pop();
                    candidateSetIndex2[current.second - index2_offset].push_back(std::make_pair(current.first, cur_c + index1_offset));

                    if (!alg_hnsw->data_.find(current.second)) {
                        alg_hnsw->data_.push_back(index2->get_data(current.second - index2_offset, level), current.second);
                    }
                }

                ll_cur = index1->get_linklist_at_level_from_memory(cur_c, level);
                int linklistCount = index1->getListCount(ll_cur);
                tableint *data = (tableint *)(ll_cur + 1);
                dist = index1->get_dist_at_level_from_memory(cur_c, level);

                tableint *data_new = (tableint *)(ll_new + 1);

                if (top_candidates.size() + linklistCount > Mcurmax) {
                    for (size_t iter = 0; iter < linklistCount; iter++) {
                        top_candidates.emplace(dist[iter], data[iter] + index1_offset);
                        data_point = index1->get_data(data[iter], level);

                        if (!alg_hnsw->data_.find(data[iter] + index1_offset))
                            alg_hnsw->data_.push_back(data_point, data[iter] + index1_offset);
                    }
                    alg_hnsw->getNeighborsByHeuristic2(top_candidates, Mcurmax, false, -1, alpha);

                    alg_hnsw->setListCount(ll_new, top_candidates.size());
                    for (size_t idx = 0; top_candidates.size() > 0; idx++) {
                        data_new[idx] = top_candidates.top().second;
                        dist_new[idx] = top_candidates.top().first;
                        top_candidates.pop();
                    }
                } else {
                    int length = top_candidates.size() + linklistCount;
                    alg_hnsw->setListCount(ll_new, length);

                    for (size_t idx = 0; idx < linklistCount; idx++) {
                        data_new[idx] = data[idx] + index1_offset;
                        dist_new[idx] = dist[idx];
                    }
                    size_t offset = linklistCount;
                    for (size_t idx = 0; top_candidates.size() > 0; idx++) {
                        data_new[idx + offset] = top_candidates.top().second;
                        dist_new[idx + offset] = top_candidates.top().first;
                        top_candidates.pop();
                    }
                }
                alg_hnsw->setLinkList(new_c, level, ll_new, ll_cur_size);
                alg_hnsw->setDistList(new_c, level, dist_new, dist_size);
                // TODO: batch write back
            }
            alg_hnsw->data_.release();
            index1->clear_linklist(level);
            index1->data_.release();
        }
        // index1->data_.release();
        index2->data_.release();
        delete batch_data;
        sum0 += sw2.getElapsedTimeMicro();
        printf("time for searching index2 on index1: %.3f s [total cost %.3f s]\n", sw2.getElapsedTimeMicro() / 1000000, elapsed() - s0);

        sw2.reset();
        sw3.reset();
        index2->get_list_at_level_from_disk_batch(level, 0, index2->cur_element_count, false, true, true);

        index1->load_graph_data_filter(index1->file_path, index1->index_offset_list0, index1->size_data_per_element_, index1->data_size_, index1->offsetData_,
                                       layer_node_for_index1[level]);

        int batch_size_2 = 4096;
        char *batch_data_2 = (char *)malloc(batch_size_2 * index2->data_size_);
        for (int b = 0; b < layer_node_for_index2[level].size(); b += batch_size_2) {
            // index2->get_data_from_disk_batch(layer_node_for_index2[level][b], batch_size_2, batch_data_2);
            for (int i = b; i < b + batch_size_2 && i < layer_node_for_index2[level].size(); ++i) {
                tableint cur_c = layer_node_for_index2[level][i];
                index2->get_data(cur_c, level);
            }
            printf("log[%d] time %f, cnt = %d\n", b, sw3.getElapsedTimeMicro() / 1e6, cnt);
            for (int i = b; i < b + batch_size_2 && i < layer_node_for_index2[level].size(); ++i) {
                tableint cur_c = layer_node_for_index2[level][i];
                tableint new_c = cur_c + index2_offset;
                std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, typename HierarchicalNSW_ME<dist_t>::CompareByFirst> top_candidates;

                ll_cur = index2->get_linklist_at_level_from_memory(cur_c, level);
                int linklistCount = index2->getListCount(ll_cur);
                tableint *data = (tableint *)(ll_cur + 1);
                dist = index2->get_dist_at_level_from_memory(cur_c, level);

                tableint *data_new = (tableint *)(ll_new + 1);

                if (candidateSetIndex2[cur_c].size() > 0) {
                    const std::vector<std::pair<dist_t, tableint>> &candidates = candidateSetIndex2[cur_c];
                    int candidate_size = candidates.size();

                    if (candidate_size + linklistCount > Mcurmax) {
                        for (const auto &p : candidates) {
                            tableint neighbor_id = p.second;
                            dist_t distance = p.first;
                            top_candidates.emplace(distance, neighbor_id);

                            if (!alg_hnsw->data_.find(neighbor_id)) {
                                alg_hnsw->data_.push_back(index1->get_data(neighbor_id - index1_offset, level), neighbor_id);
                            }
                        }

                        for (size_t iter = 0; iter < linklistCount; iter++) {
                            top_candidates.emplace(dist[iter], data[iter] + index2_offset);

                            if (!alg_hnsw->data_.find(data[iter] + index2_offset))
                                alg_hnsw->data_.push_back(index2->get_data(data[iter]), data[iter] + index2_offset, level);
                        }
                        alg_hnsw->getNeighborsByHeuristic2(top_candidates, Mcurmax, false, -1, alpha);

                        alg_hnsw->setListCount(ll_new, top_candidates.size());
                        for (size_t idx = 0; top_candidates.size() > 0; idx++) {
                            data_new[idx] = top_candidates.top().second;
                            dist_new[idx] = top_candidates.top().first;
                            top_candidates.pop();
                        }

                    } else {
                        int length = candidates.size() + linklistCount;
                        alg_hnsw->setListCount(ll_new, length);

                        for (size_t idx = 0; idx < linklistCount; idx++) {
                            data_new[idx] = data[idx] + index2_offset;
                            dist_new[idx] = dist[idx];
                        }
                        size_t offset = linklistCount;
                        for (size_t idx = 0; idx < candidate_size; idx++) {
                            data_new[idx + offset] = candidates[idx].second;
                            dist_new[idx + offset] = candidates[idx].first;
                        }
                    }
                } else {
                    alg_hnsw->setListCount(ll_new, linklistCount);
                    for (size_t idx = 0; idx < linklistCount; idx++) {
                        data_new[idx] = data[idx] + index2_offset;
                        dist_new[idx] = dist[idx];
                    }
                }
                alg_hnsw->setLinkList(new_c, level, ll_new, ll_cur_size);
                alg_hnsw->setDistList(new_c, level, dist_new, dist_size);
            }
            alg_hnsw->data_.release();
        }
        index1->data_.release();
        index2->data_.release();
        delete batch_data_2;
        delete ll_new;
        delete dist_new;
        index2->clear_linklist(level);
        printf("level=%d, time for setting links: %.3f s [total cost %.3f s]\n", level, sw3.getElapsedTimeMicro() / 1000000, elapsed() - s0);
    }

    //! LEVEL = 0, too much random access
    // index2->load_graph_data(index2->file_path, index2->index_offset_list0, index2->size_data_per_element_, index2->data_size_, index2->offsetData_, index2->cur_element_count,
    //                         index2->data_memory_);
    for (int level = 0; level >= 0; level -= 1) {
        printf("Merging layer %d\n", level);

        // if (layer_node_for_index1[level].size() > 0 && layer_node_for_index2[level].size() == 0) {
        //     throw std::logic_error("index1 is higher than index2, wrong case");
        // }
        // if (layer_node_for_index1[level].size() == 0 && layer_node_for_index2[level].size() > 0) {
        //     deepCopyOneLayerOnIndex(index2, alg_hnsw, level, index2_offset, layer_node_for_index2);
        //     continue;
        // }

        index1->data_.set_max_items(1 << 20); // larger than batch_size * 5
        index2->data_.set_max_items(1 << 20);
        index1->clear_level0();
        index2->clear_level0();
        alg_hnsw->clear_level0();
        layer_node_for_index1.assign(maxLevel + 1, {});
        layer_node_for_index2.assign(maxLevel + 1, {});

        index2->get_list_at_level_from_disk_batch(level, 0, index2->cur_element_count, true, false, true); //
        // printf("ckpt5\n");
        // sleep(5);
        index2->load_graph_data(index2->file_path, index2->index_offset_list0, index2->size_data_per_element_, index2->data_size_, index2->offsetData_, index2->cur_element_count,
                                index2->data_memory_, true);
        sw2.reset();
        sw3.reset();
        // printf("ckpt6\n");
        // sleep(5);
        // exit(0);

        size_t Mcurmax = level ? alg_hnsw_maxM_ : alg_hnsw_maxM0_;
        std::vector<std::vector<std::pair<dist_t, tableint>>> candidateSetIndex2(element_count_for_index2);

        size_t ll_cur_size = (level == 0) ? alg_hnsw->size_links_level0_ : alg_hnsw->size_links_per_element_;
        linklistsizeint *ll_cur;
        linklistsizeint *ll_new = (linklistsizeint *)malloc(ll_cur_size);

        size_t dist_size = (level == 0) ? alg_hnsw->size_dist_per_element_ : alg_hnsw->size_dist_links_per_element_;
        dist_t *dist;
        dist_t *dist_new = (dist_t *)malloc(dist_size);

        int batch_size = 4096;
        char *batch_data = nullptr;
        for (int b = 0; b < index1->cur_element_count; b += batch_size) {
            index1->get_list_at_level_from_disk_batch(level, b, batch_size, true, true); //
            index1->get_data_from_disk_batch(b, batch_size, batch_data);

            // float *xxxx = (float *)batch_data;
            // printf("%p\n", batch_data);
            // for (int xx = 0; xx < 10; xx++)
            //     printf("%lf, ", xxxx[xx]);
            // printf("\n");
            // sleep(5);

            index1->data_.push_back(batch_data, b, batch_size);
            // printf("log[%d] time %f\n", b, sw2.getElapsedTimeMicro() / 1e6);

            for (int i = b; i < b + batch_size && i < index1->cur_element_count; ++i) {
                tableint cur_c = i;
                tableint new_c = cur_c + index1_offset;

                data_point = index1->data_[cur_c];
                // float *xxxx = (float *)data_point;
                // for (int xx = 0; xx < 10; xx++)
                //     printf("%lf, ", xxxx[xx]);
                // printf("\n");
                // sleep(5);

                tableint &last_entry_point = entry_point_collect_index1_on_index2[cur_c];
                // if (index1->element_levels_[i] < level) {
                //     index2->search2Layer(data_point, last_entry_point, level, level - 1, lambda, nullptr, index2_offset);
                //     continue;
                // }
                std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, typename HierarchicalNSW_ME<dist_t>::CompareByFirst> top_candidates;
                index2->search2Layer(data_point, last_entry_point, level, level, lambda, &top_candidates, index2_offset);
                // printf("(%d, %d), ", last_entry_point, index2->data_.cnt);

                auto temp_queue = top_candidates;
                while (!temp_queue.empty()) {
                    std::pair<dist_t, tableint> current = temp_queue.top();
                    temp_queue.pop();
                    candidateSetIndex2[current.second - index2_offset].push_back(std::make_pair(current.first, cur_c + index1_offset));

                    if (!alg_hnsw->data_.find(current.second)) {
                        alg_hnsw->data_.push_back(index2->get_data(current.second - index2_offset, level), current.second);
                    }
                }

                ll_cur = index1->get_linklist_at_level_from_memory(cur_c, level);
                int linklistCount = index1->getListCount(ll_cur);
                tableint *data = (tableint *)(ll_cur + 1);
                dist = index1->get_dist_at_level_from_memory(cur_c, level);

                tableint *data_new = (tableint *)(ll_new + 1);

                if (top_candidates.size() + linklistCount > Mcurmax) {
                    // printf("level: %d, id: %d, size: %d\n", level, cur_c, linklistCount);
                    for (size_t iter = 0; iter < linklistCount; iter++) {
                        top_candidates.emplace(dist[iter], data[iter] + index1_offset);
                        data_point = index1->get_data(data[iter]);

                        if (!alg_hnsw->data_.find(data[iter] + index1_offset))
                            alg_hnsw->data_.push_back(data_point, data[iter] + index1_offset);
                    }
                    alg_hnsw->getNeighborsByHeuristic2(top_candidates, Mcurmax, false, -1, alpha);

                    alg_hnsw->setListCount(ll_new, top_candidates.size());
                    for (size_t idx = 0; top_candidates.size() > 0; idx++) {
                        data_new[idx] = top_candidates.top().second;
                        dist_new[idx] = top_candidates.top().first;
                        top_candidates.pop();
                    }
                } else {
                    int length = top_candidates.size() + linklistCount;
                    // printf("length[%d]: %d, %d\n", cur_c, top_candidates.size() , linklistCount);
                    alg_hnsw->setListCount(ll_new, length);

                    for (size_t idx = 0; idx < linklistCount; idx++) {
                        data_new[idx] = data[idx] + index1_offset;
                        dist_new[idx] = dist[idx];
                    }
                    size_t offset = linklistCount;
                    for (size_t idx = 0; top_candidates.size() > 0; idx++) {
                        data_new[idx + offset] = top_candidates.top().second;
                        dist_new[idx + offset] = top_candidates.top().first;
                        top_candidates.pop();
                    }
                }
                alg_hnsw->setLinkList(new_c, level, ll_new, ll_cur_size);
                alg_hnsw->setDistList(new_c, level, dist_new, dist_size);
                // TODO: batch write back
                if (i % 200000 == 0 || (level == 0 && (i + 1) % batch_size == 0)) {
                    printf("log[%d] time %f\n", i, sw2.getElapsedTimeMicro() / 1e6);
                    printf("count1: %d, count2: %d, count_m: %d\n", index1->data_.size(), index2->data_.size(), alg_hnsw->data_.size());
                    printf("index1 random access: %d, index2 random access: %d\n", index1->data_.cnt, index2->data_.cnt);
                }
                index2->data_.release();
                alg_hnsw->data_.release();
            }
            index1->clear_linklist(level);
            index1->data_.release();
        }
        // delete batch_data;
        sum0 += sw2.getElapsedTimeMicro();
        printf("time for searching index2 on index1: %.3f s [total cost %.3f s]\n", sw2.getElapsedTimeMicro() / 1000000, elapsed() - s0);

        sw2.reset();
        sw3.reset();
        index2->clear_linklist(level);

        int batch_size_2 = 4096;
        char *batch_data_2 = nullptr;

        index1->load_graph_data(index1->file_path, index1->index_offset_list0, index1->size_data_per_element_, index1->data_size_, index1->offsetData_, index1->cur_element_count,
                                index1->data_memory_, true);

        for (int b = 0; b < index2->cur_element_count; b += batch_size_2) {
            index2->get_list_at_level_from_disk_batch(level, b, batch_size_2, true, true, false); //
            // index2->get_data_from_disk_batch(b, batch_size_2, batch_data_2);
            // index2->data_.push_back(batch_data_2, b, batch_size_2);
            printf("log[%d] time %f\n", b, sw3.getElapsedTimeMicro() / 1e6);
            for (int i = b; i < b + batch_size_2 && i < index2->cur_element_count; ++i) {
                tableint cur_c = i;
                tableint new_c = cur_c + index2_offset;
                std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, typename HierarchicalNSW_ME<dist_t>::CompareByFirst> top_candidates;

                ll_cur = index2->get_linklist_at_level_from_memory(cur_c, level);
                int linklistCount = index2->getListCount(ll_cur);
                tableint *data = (tableint *)(ll_cur + 1);
                dist = index2->get_dist_at_level_from_memory(cur_c, level);

                tableint *data_new = (tableint *)(ll_new + 1);

                if (candidateSetIndex2[cur_c].size() > 0) {
                    const std::vector<std::pair<dist_t, tableint>> &candidates = candidateSetIndex2[cur_c];
                    int candidate_size = candidates.size();

                    if (candidate_size + linklistCount > Mcurmax) {
                        for (const auto &p : candidates) {
                            tableint neighbor_id = p.second;
                            dist_t distance = p.first;
                            top_candidates.emplace(distance, neighbor_id);

                            if (!alg_hnsw->data_.find(neighbor_id)) {
                                alg_hnsw->data_.push_back(index1->get_data(neighbor_id - index1_offset, level), neighbor_id);
                            }
                        }

                        for (size_t iter = 0; iter < linklistCount; iter++) {
                            top_candidates.emplace(dist[iter], data[iter] + index2_offset);

                            if (!alg_hnsw->data_.find(data[iter] + index2_offset))
                                alg_hnsw->data_.push_back(index2->get_data(data[iter], level), data[iter] + index2_offset);
                        }
                        alg_hnsw->getNeighborsByHeuristic2(top_candidates, Mcurmax, false, -1, alpha);

                        alg_hnsw->setListCount(ll_new, top_candidates.size());
                        for (size_t idx = 0; top_candidates.size() > 0; idx++) {
                            data_new[idx] = top_candidates.top().second;
                            dist_new[idx] = top_candidates.top().first;
                            top_candidates.pop();
                        }

                    } else {
                        int length = candidates.size() + linklistCount;
                        alg_hnsw->setListCount(ll_new, length);

                        for (size_t idx = 0; idx < linklistCount; idx++) {
                            data_new[idx] = data[idx] + index2_offset;
                            dist_new[idx] = dist[idx];
                        }
                        size_t offset = linklistCount;
                        for (size_t idx = 0; idx < candidate_size; idx++) {
                            data_new[idx + offset] = candidates[idx].second;
                            dist_new[idx + offset] = candidates[idx].first;
                        }
                    }
                } else {
                    alg_hnsw->setListCount(ll_new, linklistCount);
                    for (size_t idx = 0; idx < linklistCount; idx++) {
                        data_new[idx] = data[idx] + index2_offset;
                        dist_new[idx] = dist[idx];
                    }
                }
                alg_hnsw->setLinkList(new_c, level, ll_new, ll_cur_size);
                alg_hnsw->setDistList(new_c, level, dist_new, dist_size);
            }
            alg_hnsw->data_.release();
            index2->clear_linklist(level);
            index1->data_.release();
            index2->data_.release();
        }
        free(ll_new);
        free(dist_new);
        index2->clear_linklist(level);
        printf("level=%d, time for setting links: %.3f s [total cost %.3f s]\n", level, sw3.getElapsedTimeMicro() / 1000000, elapsed() - s0);
    }
    io_print_stats("HNSW merger iostream I/O");
    delete alg_hnsw;
    delete index1;
    delete index2;
    
}

} // namespace hnswlib