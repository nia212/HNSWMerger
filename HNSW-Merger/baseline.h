/*
 * This baseline algorithm is based on https://www.elastic.co/search-labs/blog/hnsw-graphs-speed-up-merging
 * Code is converted from Apache Lucene repository https://github.com/apache/lucene
 * For fair comparing with our algorithm, we only convert key function from java code,
 * and use the original HNSWlib code as much as possible.
 */
#include "extension.h"
#include <mutex>
#include <omp.h>

namespace hnswlib {
FILE *f;

inline int ceilDiv(int a, int b) {
    return (a + b - 1) / b;
}

inline long encode(int value1, int value2) {
    return (static_cast<long>(-value1) << 32) | (static_cast<long>(value2) & 0xFFFFFFFFL);
}

inline int decodeValue1(long encoded) {
    return static_cast<int>(-(encoded >> 32));
}

inline int decodeValue2(long encoded) {
    return static_cast<int>(encoded & 0xFFFFFFFFL);
}

// based on apache /lucene/lucene/core/src/java/org/apache/lucene/util/hnsw/UpdateGraphsUtils.java
template <typename dist_t>
std::unordered_set<int> computeJoinSet(HierarchicalNSW<dist_t> *index) {
    int N = index->cur_element_count;
    std::priority_queue<long, std::vector<long>, std::greater<long>> heap;

    std::unordered_set<int> joinSet;
    std::vector<bool> stale(N, false);
    std::vector<short> counts(N, 0);

    long gExit = 0L;
    for (int v = 0; v < N; ++v) {
        linklistsizeint *ll_cur = index->get_linklist_at_level(v, 0);
        size_t degree = index->getListCount(ll_cur);
        int k = (degree < 9 ? 2 : ceilDiv(degree, 4));
        gExit += k;
        int gain = k + degree;
        heap.push(encode(gain, v));
    }

    long gTot = 0L;
    while (gTot < gExit && !heap.empty()) {
        long topEncoded = heap.top();
        heap.pop();

        int gain = decodeValue1(topEncoded);
        int v = decodeValue2(topEncoded);

        linklistsizeint *ll_cur = index->get_linklist_at_level(v, 0);
        size_t degree = index->getListCount(ll_cur);

        std::vector<int> ns;
        ns.reserve(degree);
        tableint *data = (tableint *)(ll_cur + 1);
        for (size_t iter = 0; iter < degree; iter++) {
            ns.push_back(data[iter]);
        }

        int k = (degree < 9 ? 2 : ceilDiv(degree, 4));

        if (stale[v]) {
            int newGain = std::max(0, k - counts[v]);
            for (int i = 0; i < ns.size(); i++) {
                int u = ns[i];
                if (counts[u] < k && joinSet.find(u) == joinSet.end()) {
                    newGain += 1;
                }
            }
            if (newGain > 0) {
                heap.push(encode(newGain, v));
                stale[v] = false;
            }
        } else {
            joinSet.insert(v);
            gTot += gain;

            bool markNeighboursStale = (counts[v] < k);
            for (int i = 0; i < ns.size(); i++) {
                int u = ns[i];
                if (markNeighboursStale) {
                    stale[u] = true;
                }
                if (counts[u] < (k - 1)) {
                    linklistsizeint *ll_cur = index->get_linklist_at_level(u, 0);
                    size_t degree = index->getListCount(ll_cur);
                    tableint *data = (tableint *)(ll_cur + 1);
                    for (size_t iter = 0; iter < degree; iter++) {
                        stale[data[iter]] = true;
                    }
                }
                counts[u] += 1;
            }
        }
    }

    return joinSet;
}

// same logic with original HNSWlib addPoint function, except changing searching function in inserting phase to a beam search
template <typename dist_t>
tableint HierarchicalNSW<dist_t>::addPoint(labeltype label,
                                           const void *data_point,
                                           std::unordered_set<tableint> *eps0) {
    tableint cur_c = 0;
    {
        std::unique_lock<std::mutex> lock_table(label_lookup_lock);
        auto search = label_lookup_.find(label);
        if (search != label_lookup_.end()) {
            tableint existingInternalId = search->second;
            if (allow_replace_deleted_) {
                if (isMarkedDeleted(existingInternalId)) {
                    throw std::runtime_error("Can't use addPoint to update deleted elements if replacement of deleted elements is enabled.");
                }
            }
            lock_table.unlock();

            if (isMarkedDeleted(existingInternalId)) {
                unmarkDeletedInternal(existingInternalId);
            }
            updatePoint(data_point, existingInternalId, 1.0);

            return existingInternalId;
        }

        if (cur_element_count >= max_elements_) {
            throw std::runtime_error("The number of elements exceeds the specified limit");
        }

        cur_c = cur_element_count;
        cur_element_count++;
        label_lookup_[label] = cur_c;
    }

    std::unique_lock<std::mutex> lock_el(link_list_locks_[cur_c]);
    int curlevel = getRandomLevel(mult_);

    element_levels_[cur_c] = curlevel;

    std::unique_lock<std::mutex> templock(global);
    int maxlevelcopy = maxlevel_;
    if (curlevel <= maxlevelcopy)
        templock.unlock();
    tableint currObj = enterpoint_node_;
    tableint enterpoint_copy = enterpoint_node_;

    memset(data_level0_memory_ + cur_c * size_data_per_element_ + offsetLevel0_, 0, size_data_per_element_);
    memset(dist_level0_memory_ + cur_c * size_dist_per_element_ + offsetLevel0_, 0, size_dist_per_element_);

    memcpy(getExternalLabeLp(cur_c), &label, sizeof(labeltype));
    memcpy(getDataByInternalId(cur_c), data_point, data_size_);

    if (curlevel) {
        linkLists_[cur_c] = (char *)malloc(size_links_per_element_ * curlevel + 1);
        if (linkLists_[cur_c] == nullptr)
            throw std::runtime_error("Not enough memory: addPoint failed to allocate linklist");
        memset(linkLists_[cur_c], 0, size_links_per_element_ * curlevel + 1);
        dist_linkLists_[cur_c] = (char *)malloc(size_dist_links_per_element_ * curlevel + 1);
        if (dist_linkLists_[cur_c] == nullptr)
            throw std::runtime_error("Not enough memory: addPoint failed to allocate dist_linkLists_");
        memset(dist_linkLists_[cur_c], 0, size_dist_links_per_element_ * curlevel + 1);
    }

    if ((signed)currObj != -1) {
        if (curlevel < maxlevelcopy) {
            dist_t curdist = fstdistfunc_(data_point, getDataByInternalId(currObj), dist_func_param_);
            for (int level = maxlevelcopy; level > curlevel; level--) {
                bool changed = true;
                while (changed) {
                    changed = false;
                    unsigned int *data;
                    std::unique_lock<std::mutex> lock(link_list_locks_[currObj]);
                    data = get_linklist(currObj, level);
                    int size = getListCount(data);

                    tableint *datal = (tableint *)(data + 1);
                    for (int i = 0; i < size; i++) {
                        tableint cand = datal[i];
                        if (cand < 0 || cand > max_elements_)
                            throw std::runtime_error("cand error");
                        dist_t d = fstdistfunc_(data_point, getDataByInternalId(cand), dist_func_param_);
                        if (d < curdist) {
                            curdist = d;
                            currObj = cand;
                            changed = true;
                        }
                    }
                }
            }
        }
        //    beamCandidates0 = new GraphBuilderKnnCollector(Math.min(beamWidth / 2, M * 3));
        bool epDeleted = isMarkedDeleted(enterpoint_copy);
        std::unordered_set<tableint> eps;
        eps.insert(currObj);

        for (int level = std::min(curlevel, maxlevelcopy); level >= 0; level--) {
            if (level > maxlevelcopy || level < 0) // possible?
                throw std::runtime_error("Level error");
            size_t Mcurmax = level ? maxM_ : maxM0_;
            if (level == 0 && eps0 != nullptr && eps0->size() > 0) eps = *eps0;

            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates;
            top_candidates = ExtendSearchBaseLayer(data_point,
                                                   level,
                                                   &eps);
            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> cand(top_candidates);
            eps.clear();
            while (!cand.empty()) {
                if (cand.size() <= Mcurmax) {
                    eps.insert(cand.top().second);
                }
                cand.pop();
            }

            if (epDeleted) {
                top_candidates.emplace(fstdistfunc_(data_point, getDataByInternalId(enterpoint_copy), dist_func_param_), enterpoint_copy);
                if (top_candidates.size() > ef_construction_)
                    top_candidates.pop();
            }
            currObj = mutuallyConnectNewElement(data_point, cur_c, top_candidates, level, false);
        }
    } else {
        enterpoint_node_ = 0;
        maxlevel_ = curlevel;
    }

    if (curlevel > maxlevelcopy) {
        enterpoint_node_ = cur_c;
        maxlevel_ = curlevel;
    }
    return cur_c;
}

template <typename dist_t>
HierarchicalNSW<dist_t> *HNSWMerger_ES(HierarchicalNSW<dist_t> *index1, HierarchicalNSW<dist_t> *index2, L2Space *space, size_t M = -1, size_t ef_construction = -1) {

    double s0 = elapsed();
    if (index1->cur_element_count < index2->cur_element_count || index1->maxlevel_ < index2->maxlevel_) {
        std::swap(index1, index2);
        printf("Swap index with more current element count to index1.\n");
    }
    size_t element_count_for_index1 = index1->getCurrentElementCount();
    M = M == -1 ? std::max(index1->M_, index2->M_) : M;
    ef_construction = ef_construction == -1 ? std::max(index1->ef_construction_, index2->ef_construction_) : ef_construction;
    size_t new_max_elements = index1->max_elements_ + index2->max_elements_;
    HierarchicalNSW<dist_t> *alg_hnsw = new HierarchicalNSW<dist_t>(space, new_max_elements, M, ef_construction);
    size_t maxLevel = index1->maxlevel_;
    alg_hnsw->setMaxLevel(maxLevel);
    alg_hnsw->cur_element_count.store(index1->cur_element_count);

    alg_hnsw->data_level0_memory_ = (char *)malloc(alg_hnsw->max_elements_ * alg_hnsw->size_data_per_element_);
    if (alg_hnsw->data_level0_memory_ == nullptr)
        throw std::runtime_error("Not enough memory: loadIndex failed to allocate level0");

    alg_hnsw->dist_level0_memory_ = (char *)malloc(alg_hnsw->max_elements_ * alg_hnsw->size_dist_per_element_);
    if (alg_hnsw->dist_level0_memory_ == nullptr)
        throw std::runtime_error("Not enough memory: loadIndex failed to allocate level0");

    alg_hnsw->enterpoint_node_ = index1->enterpoint_node_;

    alg_hnsw->label_lookup_.clear();
    alg_hnsw->label_lookup_.reserve(alg_hnsw->max_elements_);

    for (int id = 0; id < element_count_for_index1; id++) {
        alg_hnsw->label_lookup_[index1->getExternalLabel(id)] = id;
    }
    for (int id = 0; id < index1->cur_element_count; id++) {
        if (index1->element_levels_[id] < 1)
            continue;
        int level = index1->element_levels_[id];
        alg_hnsw->element_levels_[id] = level;
        tableint new_c = id;
        alg_hnsw->linkLists_[new_c] = (char *)malloc(alg_hnsw->size_links_per_element_ * level + 1);
        if (alg_hnsw->linkLists_[new_c] == nullptr) {
            throw std::runtime_error("Not enough memory: addPoint failed to allocate linklist");
        }
        memcpy(alg_hnsw->linkLists_[new_c],
               index1->linkLists_[id],
               alg_hnsw->size_links_per_element_ * level + 1);
        alg_hnsw->dist_linkLists_[new_c] = (char *)malloc(alg_hnsw->size_dist_links_per_element_ * level + 1);
        if (alg_hnsw->dist_linkLists_[new_c] == nullptr) {
            throw std::runtime_error("Not enough memory: addPoint failed to allocate dist_linklist");
        }
        memcpy(alg_hnsw->dist_linkLists_[new_c],
               index1->dist_linkLists_[id],
               alg_hnsw->size_dist_links_per_element_ * level + 1);
    }
    memcpy(alg_hnsw->data_level0_memory_,
           index1->data_level0_memory_,
           element_count_for_index1 * alg_hnsw->size_data_per_element_);
    memcpy(alg_hnsw->dist_level0_memory_,
           index1->dist_level0_memory_,
           element_count_for_index1 * alg_hnsw->size_dist_per_element_);

    printf("time for new layer: %f\n", elapsed() - s0);

    s0 = elapsed();
    int size = index2->cur_element_count;
    std::unordered_set<int> j = computeJoinSet(index2);
    std::unordered_map<int, int> ordMapS;
    ordMapS.reserve(size);
    printf("Join set size: %d\n", j.size());
    printf("time for computeJoinSet: %f\n", elapsed() - s0);

    int count = 0;
    s0 = elapsed();
    for (int node : j) {
        tableint newId = alg_hnsw->addPoint(index2->getExternalLabel(node), index2->getDataByInternalId(node), nullptr);
        ordMapS[node] = newId;
        if ((++count) % 200000 == 0) {
            printf("Processing node %d, time: %f\n", count, elapsed() - s0);
        }
    }
    printf("time for addPoint in join set: %f\n", elapsed() - s0);

    s0 = elapsed();
    count = 0;
    for (int u = 0; u < size; ++u) {
        if ((u + 1) % 200000 == 0) {
            printf("Processing node %d(%d), time: %f\n", u + 1, count, elapsed() - s0);
        }
        if (j.find(u) != j.end()) {
            continue;
        }
        count++;
        std::unordered_set<tableint> eps;
        linklistsizeint *ll_cur = index2->get_linklist_at_level(u, 0);
        size_t degree = index2->getListCount(ll_cur);
        tableint *data = (tableint *)(ll_cur + 1);
        for (size_t iter = 0; iter < degree; iter++) {
            int v = data[iter];
            if (v < u || j.find(v) != j.end()) {
                tableint newv = ordMapS[v];
                eps.insert(newv);

                linklistsizeint *ll_cur2 = alg_hnsw->get_linklist_at_level(newv, 0);
                size_t degree2 = alg_hnsw->getListCount(ll_cur2);
                tableint *data2 = (tableint *)(ll_cur2 + 1);
                for (size_t iter2 = 0; iter2 < degree2; iter2++) {
                    tableint friendOrd = data2[iter2];
                    eps.insert(friendOrd);
                }
            }
        }
        tableint newId = alg_hnsw->addPoint(index2->getExternalLabel(u), index2->getDataByInternalId(u), &eps);
        ordMapS[u] = newId;
    }
    printf("time for addPoint in non-join set: %f\n", elapsed() - s0);
    return alg_hnsw;
}

} // namespace hnswlib