#pragma once
#include "hnswalg.h"
#include <mutex>
#include <omp.h>
#include <random>

using namespace std;

extern float global_counter;

namespace hnswlib {

template <typename dist_t>
void HierarchicalNSW<dist_t>::searchNodeOnEachLayer(
    std::vector<std::vector<int>> &resultVector,
    bool combine) {
    size_t elementCount = getCurrentElementCount();
    for (size_t iter = 0; iter < elementCount; iter++) {
        if (!isMarkedDeleted(iter)) {
            resultVector[element_levels_[iter]].push_back(iter);
        }
    }
    if (combine) {
        for (int i = resultVector.size() - 2; i >= 0; --i) {
            resultVector[i].insert(resultVector[i].end(), resultVector[i + 1].begin(), resultVector[i + 1].end());
        }
    }
}

/*
 * This function, which is used in search2Layer, searching the graph from the
 * top layer to the level_low layer. If the enterpoint_node is -1, it will be
 * set to enterpoint_node_. Otherwise, the enterpoint_node will be used as the
 * starting point, as it's the closest point to the query on layer level_low
 * + 1. The function will also return the `top_candidates` which contains the
 * `lambda`-closest points to the query on layer level_low.
 */
template <typename dist_t>
void HierarchicalNSW<dist_t>::search2Layer(const void *query_data,
                                           tableint &enterpoint_node,
                                           int level_higher,
                                           int level_lower,
                                           int lambda,
                                           std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> *top_candidates,
                                           int offset) {
    tableint currObj = enterpoint_node == -1 ? enterpoint_node_ : enterpoint_node;
    dist_t curdist = fstdistfunc_(query_data, getDataByInternalId(currObj), dist_func_param_);
    for (int level = level_higher; level > level_lower; level--) {
        bool changed = true;
        while (changed) {
            changed = false;
            unsigned int *data = (unsigned int *)get_linklist_at_level(currObj, level);
            int size = getListCount(data);

            tableint *datal = (tableint *)(data + 1);

            for (int i = 0; i < size; i++) {
                tableint cand = datal[i];
                if (cand < 0 || cand > max_elements_) {
                    throw std::runtime_error("cand error");
                }
                dist_t d = fstdistfunc_(query_data, getDataByInternalId(cand), dist_func_param_);

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
    if (!isMarkedDeleted(currObj)) {
        dist_t dist = fstdistfunc_(query_data, getDataByInternalId(currObj), dist_func_param_);
        top_candidates->emplace(dist, currObj + offset);
        enterpoint_node = currObj;
        entry_dist = dist;
        lowerBound = dist;
        candidateSet.emplace(-dist, currObj);
    } else {
        entry_dist = std::numeric_limits<dist_t>::max();
        lowerBound = std::numeric_limits<dist_t>::max();
        candidateSet.emplace(-lowerBound, currObj);
    }
    visited_array[currObj] = visited_array_tag;

    while (!candidateSet.empty()) {
        std::pair<dist_t, tableint> curr_el_pair = candidateSet.top();
        if ((-curr_el_pair.first) > lowerBound && top_candidates->size() == lambda) {
            break;
        }
        candidateSet.pop();

        tableint curNodeNum = curr_el_pair.second;

        int *data = (int *)get_linklist_at_level(curNodeNum, level_lower);
        size_t size = getListCount((linklistsizeint *)data);
        tableint *datal = (tableint *)(data + 1);
#ifdef USE_SSE
        _mm_prefetch((char *)(visited_array + *(data + 1)), _MM_HINT_T0);
        _mm_prefetch((char *)(visited_array + *(data + 1) + 64), _MM_HINT_T0);
        _mm_prefetch(getDataByInternalId(*datal), _MM_HINT_T0);
        _mm_prefetch(getDataByInternalId(*(datal + 1)), _MM_HINT_T0);
#endif

        for (size_t j = 0; j < size; j++) {
            tableint candidate_id = *(datal + j);
#ifdef USE_SSE
            _mm_prefetch((char *)(visited_array + *(datal + j + 1)), _MM_HINT_T0);
            _mm_prefetch(getDataByInternalId(*(datal + j + 1)), _MM_HINT_T0);
#endif
            if (visited_array[candidate_id] == visited_array_tag)
                continue;
            visited_array[candidate_id] = visited_array_tag;
            char *currObj1 = (getDataByInternalId(candidate_id));

            dist_t dist1 = fstdistfunc_(query_data, currObj1, dist_func_param_);
            if (top_candidates->size() < lambda || lowerBound > dist1) {
                candidateSet.emplace(-dist1, candidate_id);
#ifdef USE_SSE
                _mm_prefetch(getDataByInternalId(candidateSet.top().second), _MM_HINT_T0);
#endif

                if (!isMarkedDeleted(candidate_id)) {
                    top_candidates->emplace(dist1, candidate_id + offset);
                    if (entry_dist > dist1) {
                        enterpoint_node = candidate_id;
                        entry_dist = dist1;
                    }
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

// This function extends the search to BEAM SEARCH from the base layer to the next layer
// using the entry point and the eps set. It returns a priority queue of
// the top candidates found in the search, sorted by distance.
// Compared to the original SearchBaseLayer, this function uses a
// beam search approach rather than a simple entry point search.
template <typename dist_t>
std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, typename HierarchicalNSW<dist_t>::CompareByFirst>
HierarchicalNSW<dist_t>::ExtendSearchBaseLayer(const void *query_data,
                                               int level,
                                               std::unordered_set<tableint> *eps,
                                               size_t local_ef) {
    size_t beamWidth;
    if (local_ef == -1) { 
        local_ef = ef_construction_;
        beamWidth = level == 0 ? (size_t)(maxM_ * 2 / 3) : (size_t)(maxM_); 
    }
    else{
        beamWidth = level == 0 ? (size_t)(local_ef / 3) : (size_t)(local_ef/2); 
    }
    VisitedList *vl = visited_list_pool_->getFreeVisitedList();
    vl_type *visited_array = vl->mass;
    vl_type visited_array_tag = vl->curV;

    std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates;
    std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates2;
    std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> candidateSet;

    dist_t lowerBound = std::numeric_limits<dist_t>::max();

    for (const tableint &val : *eps) {
        dist_t dist = fstdistfunc_(query_data, getDataByInternalId(val), dist_func_param_);
        top_candidates.emplace(dist, val);
        candidateSet.emplace(-dist, val);
        visited_array[val] = visited_array_tag;
        if (top_candidates.size() > beamWidth) {
            top_candidates2.emplace(top_candidates.top());
            top_candidates.pop();
            if (top_candidates2.size() > local_ef) {
                top_candidates2.pop();
            }
        }
    }
    if (top_candidates.empty()) {
        lowerBound = std::numeric_limits<dist_t>::max();
    } else {
        if (top_candidates.size() + top_candidates2.size() >= local_ef || top_candidates2.size() == 0) {
            lowerBound = top_candidates.top().first;
        } else {
            lowerBound = top_candidates2.top().first;
        }
    }

    while (!candidateSet.empty()) {
        std::pair<dist_t, tableint> curr_el_pair = candidateSet.top();
        if ((-curr_el_pair.first) > lowerBound && (top_candidates.size() + top_candidates2.size() >= local_ef)) {
            break;
        }
        candidateSet.pop();

        tableint curNodeNum = curr_el_pair.second;

        std::unique_lock<std::mutex> lock(link_list_locks_[curNodeNum]);

        int *data = (int *)get_linklist_at_level(curNodeNum, level);
        size_t size = getListCount((linklistsizeint *)data);
        tableint *datal = (tableint *)(data + 1);
#ifdef USE_SSE
        _mm_prefetch((char *)(visited_array + *(data + 1)), _MM_HINT_T0);
        _mm_prefetch((char *)(visited_array + *(data + 1) + 64), _MM_HINT_T0);
        _mm_prefetch(getDataByInternalId(*datal), _MM_HINT_T0);
        _mm_prefetch(getDataByInternalId(*(datal + 1)), _MM_HINT_T0);
#endif

        for (size_t j = 0; j < size; j++) {
            tableint candidate_id = *(datal + j);
#ifdef USE_SSE
            _mm_prefetch((char *)(visited_array + *(datal + j + 1)), _MM_HINT_T0);
            _mm_prefetch(getDataByInternalId(*(datal + j + 1)), _MM_HINT_T0);
#endif
            if (visited_array[candidate_id] == visited_array_tag)
                continue;
            visited_array[candidate_id] = visited_array_tag;
            char *currObj1 = (getDataByInternalId(candidate_id));

            dist_t dist1 = fstdistfunc_(query_data, currObj1, dist_func_param_);
            if ((top_candidates.size() + top_candidates2.size() < local_ef) || lowerBound > dist1) {
                candidateSet.emplace(-dist1, candidate_id);
#ifdef USE_SSE
                _mm_prefetch(getDataByInternalId(candidateSet.top().second), _MM_HINT_T0);
#endif

                if (!isMarkedDeleted(candidate_id)) {
                    top_candidates.emplace(dist1, candidate_id);
                }

                while (top_candidates.size() > beamWidth) {
                    top_candidates2.emplace(top_candidates.top());
                    top_candidates.pop();
                    if (top_candidates2.size() > local_ef - beamWidth) {
                        top_candidates2.pop();
                    }
                }

                if (!top_candidates.empty()) {
                    if (top_candidates.size() + top_candidates2.size() >= local_ef || top_candidates2.size() == 0) {
                        lowerBound = top_candidates.top().first;
                    } else {
                        lowerBound = top_candidates2.top().first;
                    }
                }
            }
        }
    }
    visited_list_pool_->releaseVisitedList(vl);
    while (!top_candidates2.empty()) {
        top_candidates.emplace(top_candidates2.top());
        top_candidates2.pop();
    }
    return top_candidates;
}

template <typename dist_t>
std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, typename HierarchicalNSW<dist_t>::CompareByFirst>
HierarchicalNSW<dist_t>::ExtendSearchBaseLayer2(const void *query_data,
                                               int level,
                                               std::unordered_set<tableint> *eps,
                                               size_t local_ef) {
    size_t beamWidth;
    if (local_ef == -1) { 
        local_ef = ef_construction_;
        // beamWidth = level == 0 ? (size_t)(maxM_ * 2 / 3) : (size_t)(maxM_); 
    }
    else{
        // beamWidth = level == 0 ? (size_t)(local_ef / 3) : (size_t)(local_ef/2); 
    }
    VisitedList *vl = visited_list_pool_->getFreeVisitedList();
    vl_type *visited_array = vl->mass;
    vl_type visited_array_tag = vl->curV;

    std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates;
    // std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates2;
    std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> candidateSet;

    dist_t lowerBound = std::numeric_limits<dist_t>::max();

    for (const tableint &val : *eps) {
        if(visited_array[val] == visited_array_tag) continue;

        dist_t dist = fstdistfunc_(query_data, getDataByInternalId(val), dist_func_param_);
        top_candidates.emplace(dist, val);
        candidateSet.emplace(-dist, val);
        visited_array[val] = visited_array_tag;
        if(top_candidates.size() > local_ef) {
            top_candidates.pop();
        }
        // if (top_candidates.size() > beamWidth) {
        //     top_candidates2.emplace(top_candidates.top());
        //     top_candidates.pop();
        //     if (top_candidates2.size() > local_ef) {
        //         top_candidates2.pop();
        //     } 
        // }
    }
    lower_bound = top_candidates.empty() ? std::numeric_limits<dist_t>::max() : top_candidates.top().first;
    // if (top_candidates.empty()) {
    //     lowerBound = std::numeric_limits<dist_t>::max();
    // } else {
    //     if (top_candidates.size() + top_candidates2.size() >= local_ef || top_candidates2.size() == 0) {
    //         lowerBound = top_candidates.top().first;
    //     } else {
    //         lowerBound = top_candidates2.top().first;
    //     }
    // }

    while (!candidateSet.empty()) {
        std::pair<dist_t, tableint> curr_el_pair = candidateSet.top();
        // if ((-curr_el_pair.first) > lowerBound && (top_candidates.size() + top_candidates2.size() >= local_ef)) {
        if ((-curr_el_pair.first) > lowerBound && (top_candidates.size() >= local_ef)) {
            break;
        }
        candidateSet.pop();

        tableint curNodeNum = curr_el_pair.second;

        std::unique_lock<std::mutex> lock(link_list_locks_[curNodeNum]);

        int *data = (int *)get_linklist_at_level(curNodeNum, level);
        size_t size = getListCount((linklistsizeint *)data);
        tableint *datal = (tableint *)(data + 1);
#ifdef USE_SSE
        _mm_prefetch((char *)(visited_array + *(data + 1)), _MM_HINT_T0);
        _mm_prefetch((char *)(visited_array + *(data + 1) + 64), _MM_HINT_T0);
        _mm_prefetch(getDataByInternalId(*datal), _MM_HINT_T0);
        _mm_prefetch(getDataByInternalId(*(datal + 1)), _MM_HINT_T0);
#endif

        for (size_t j = 0; j < size; j++) {
            tableint candidate_id = *(datal + j);
#ifdef USE_SSE
            _mm_prefetch((char *)(visited_array + *(datal + j + 1)), _MM_HINT_T0);
            _mm_prefetch(getDataByInternalId(*(datal + j + 1)), _MM_HINT_T0);
#endif
            if (visited_array[candidate_id] == visited_array_tag)
                continue;
            visited_array[candidate_id] = visited_array_tag;
            char *currObj1 = (getDataByInternalId(candidate_id));

            dist_t dist1 = fstdistfunc_(query_data, currObj1, dist_func_param_);
            if ((top_candidates.size() < local_ef) || lowerBound > dist1) {
                candidateSet.emplace(-dist1, candidate_id);
#ifdef USE_SSE
                _mm_prefetch(getDataByInternalId(candidateSet.top().second), _MM_HINT_T0);
#endif

                
                top_candidates.emplace(dist1, candidate_id);

                while (top_candidates.size() > local_ef) {
                    top_candidates.pop();
                }

                if (!top_candidates.empty())
                    lowerBound = top_candidates.top().first;
            }
        }
    }
    visited_list_pool_->releaseVisitedList(vl);
    return top_candidates;
}

template <typename dist_t>
void HierarchicalNSW<dist_t>::deepCopyOneLayerOnIndex(HierarchicalNSW<dist_t> *index,
                                                      int level,
                                                      int offset,
                                                      std::vector<std::vector<int>> &layer_node_for_index) {
    for (int iter = 0; iter < layer_node_for_index[level].size(); iter++) {
        tableint cur_c = layer_node_for_index[level][iter];
        tableint new_c = cur_c + offset;

        linklistsizeint *ll_cur = index->get_linklist(cur_c, level);
        linklistsizeint *ll_new = get_linklist(new_c, level);
        size_t linklistCount = index->getListCount(ll_cur);
        setListCount(ll_new, linklistCount);
        tableint *data_cur = (tableint *)(ll_cur + 1);
        tableint *data_new = (tableint *)(ll_new + 1);
        dist_t *ll_cur_dist = index->get_dist_at_level(cur_c, level);
        dist_t *ll_new_dist = get_dist_at_level(new_c, level);

        if (level == 0) {
            for (int i = 0; i < linklistCount; i++) {
                data_new[i] = data_cur[i] + offset;
            }
            memcpy(ll_new_dist, ll_cur_dist, index->size_dist_per_element_);
        } else {
            for (int i = 0; i < linklistCount; i++) {
                data_new[i] = data_cur[i] + offset;
            }
            memcpy(ll_new_dist, ll_cur_dist, index->size_dist_links_per_element_);
        }
    }
}

class StopW {
    std::chrono::steady_clock::time_point time_begin;

public:
    StopW() {
        time_begin = std::chrono::steady_clock::now();
    }

    float getElapsedTimeMicro() {
        std::chrono::steady_clock::time_point time_end = std::chrono::steady_clock::now();
        return (std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_begin).count());
    }

    void reset() {
        time_begin = std::chrono::steady_clock::now();
    }
};

template <typename dist_t>
HierarchicalNSW<dist_t> *HNSWMerger(HierarchicalNSW<dist_t> *index1, HierarchicalNSW<dist_t> *index2, L2Space *space, int lambda = 4, float alpha = 1.05, size_t M = -1, size_t ef_construction = -1) {

    double s0 = elapsed();
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
    HierarchicalNSW<dist_t> *alg_hnsw = new HierarchicalNSW<dist_t>(space, new_max_elements, M, ef_construction);
    size_t maxLevel = std::max(index1->maxlevel_, index2->maxlevel_);
    alg_hnsw->setMaxLevel(maxLevel);
    alg_hnsw->cur_element_count.store(index1->cur_element_count + index2->cur_element_count);

    alg_hnsw->data_level0_memory_ = (char *)malloc(alg_hnsw->max_elements_ * alg_hnsw->size_data_per_element_);
    if (alg_hnsw->data_level0_memory_ == nullptr)
        throw std::runtime_error("Not enough memory: loadIndex failed to allocate level0");

    alg_hnsw->dist_level0_memory_ = (char *)malloc(alg_hnsw->max_elements_ * alg_hnsw->size_dist_per_element_);
    if (alg_hnsw->dist_level0_memory_ == nullptr)
        throw std::runtime_error("Not enough memory: loadIndex failed to allocate level0");

    std::vector<std::vector<int>> layer_node_for_index1(maxLevel + 2);
    index1->searchNodeOnEachLayer(layer_node_for_index1, true);
    std::vector<std::vector<int>> layer_node_for_index2(maxLevel + 2);
    index2->searchNodeOnEachLayer(layer_node_for_index2, true);

    tableint *entry_point_collect_index1_on_index2 = new tableint[element_count_for_index1];
    memset(entry_point_collect_index1_on_index2, -1, element_count_for_index1 * sizeof(int));

    int entry_point_index1 = index2->enterpoint_node_;

    memcpy(alg_hnsw->data_level0_memory_,
           index1->data_level0_memory_,
           element_count_for_index1 * alg_hnsw->size_data_per_element_);
    memcpy(alg_hnsw->data_level0_memory_ + element_count_for_index1 * alg_hnsw->size_data_per_element_,
           index2->data_level0_memory_,
           (element_count_for_index2)*alg_hnsw->size_data_per_element_);
    memcpy(alg_hnsw->dist_level0_memory_,
           index1->dist_level0_memory_,
           element_count_for_index1 * alg_hnsw->size_dist_per_element_);
    memcpy(alg_hnsw->dist_level0_memory_ + element_count_for_index1 * alg_hnsw->size_dist_per_element_,
           index2->dist_level0_memory_,
           (element_count_for_index2)*alg_hnsw->size_dist_per_element_);

    int layer_max_min = std::min(index1->maxlevel_, index2->maxlevel_);
    int bound = layer_node_for_index1[layer_max_min + 1].size() + layer_node_for_index2[layer_max_min + 1].size();
    if (bound == 0) {
        bound = 1;
    }
    std::unordered_set<int> newLayerSet;
    if (layer_node_for_index1[layer_max_min].size() + layer_node_for_index2[layer_max_min].size() > M * bound) {
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        int lambda = 0;
        auto filter_and_remove = [&](std::vector<int> &layer_nodes, std::vector<int> &upper_layer_nodes, int offset) {
            std::vector<int> newLayer;
            auto it = layer_nodes.begin();
            while (it != layer_nodes.end()) {
                if (distribution(alg_hnsw->level_generator_) < 1.0 / M || lambda == M) {
                    tableint old_c = *it;
                    tableint new_c = old_c + offset;
                    newLayer.push_back(new_c);
                    upper_layer_nodes.push_back(old_c);
                    newLayerSet.insert(new_c);
                    lambda = 0;

                    alg_hnsw->linkLists_[new_c] = (char *)malloc(alg_hnsw->size_links_per_element_ * (layer_max_min + 1) + 1);
                    if (alg_hnsw->linkLists_[new_c] == nullptr)
                        throw std::runtime_error("Not enough memory: addPoint failed to allocate linklist");
                    memset(alg_hnsw->linkLists_[new_c], 0, alg_hnsw->size_links_per_element_ * (layer_max_min + 1) + 1);
                    alg_hnsw->element_levels_[new_c] = (layer_max_min + 1);

                    alg_hnsw->dist_linkLists_[new_c] = (char *)malloc(alg_hnsw->size_dist_links_per_element_ * (layer_max_min + 1) + 1);
                    if (alg_hnsw->dist_linkLists_[new_c] == nullptr)
                        throw std::runtime_error("Not enough memory: addPoint failed to allocate dist_linkLists");
                    memset(alg_hnsw->dist_linkLists_[new_c], 0, alg_hnsw->size_dist_links_per_element_ * (layer_max_min + 1) + 1);
                } else {
                    lambda++;
                    ++it;
                }
            }
            if (newLayer.size() > 0) {
                if (newLayer.size() > M) {
                    throw std::logic_error("too many points on new-top here, debug for less neighbors");
                }
                for (int newLayerIter = 0; newLayerIter < newLayer.size(); newLayerIter++) {
                    linklistsizeint *ll_cur = alg_hnsw->get_linklist_at_level(newLayer[newLayerIter], layer_max_min + 1);
                    alg_hnsw->setListCount(ll_cur, newLayer.size() - 1);
                    tableint *data = (tableint *)(ll_cur + 1);
                    dist_t *dist = alg_hnsw->get_dist_at_level(newLayer[newLayerIter], layer_max_min + 1);
                    for (size_t it = 0, idx = 0; it < newLayer.size(); it++) {
                        if (it == newLayerIter)
                            continue;
                        dist_t dist0 =
                            alg_hnsw->fstdistfunc_(alg_hnsw->getDataByInternalId(newLayer[newLayerIter]), alg_hnsw->getDataByInternalId(newLayer[it]), alg_hnsw->dist_func_param_);
                        data[idx] = newLayer[it];
                        dist[idx] = dist0;
                        idx++;
                    }
                }
            }
        };
        if (layer_node_for_index1[layer_max_min + 1].size() == 0 && layer_node_for_index1[layer_max_min].size() > M) {
            filter_and_remove(layer_node_for_index1[layer_max_min], layer_node_for_index1[layer_max_min + 1], index1_offset);
            index1->setMaxLevel(index1->maxlevel_ + 1);
            index1->enterpoint_node_ = layer_node_for_index1[layer_max_min + 1][0];
        }
        if (layer_node_for_index2[layer_max_min + 1].size() == 0 && layer_node_for_index2[layer_max_min].size() > M) {
            filter_and_remove(layer_node_for_index2[layer_max_min], layer_node_for_index2[layer_max_min + 1], index2_offset);
            index2->setMaxLevel(index2->maxlevel_ + 1);
            index2->enterpoint_node_ = layer_node_for_index2[layer_max_min + 1][0];
        }
    }
    if (index2->maxlevel_ < index1->maxlevel_) {
        throw std::logic_error("index1 is higher than index2, wrong case");
    }
    maxLevel = index2->maxlevel_;
    alg_hnsw->setMaxLevel(maxLevel);
    alg_hnsw->enterpoint_node_ = index2->enterpoint_node_ + index2_offset;

    alg_hnsw->label_lookup_.clear();
    alg_hnsw->label_lookup_.reserve(alg_hnsw->max_elements_);

    for (int id = 0; id < element_count_for_index1; id++) {
        alg_hnsw->label_lookup_[index1->getExternalLabel(id)] = id + index1_offset;
    }
    for (int id = 0; id < element_count_for_index2; id++) {
        alg_hnsw->label_lookup_[index2->getExternalLabel(id)] = id + index2_offset;
    }

    auto allocateMemory = [&](hnswlib::HierarchicalNSW<dist_t> *index, int element_count_offset) {
        for (int id = 0; id < index->cur_element_count; id++) {
            if (index->element_levels_[id] < 1 || newLayerSet.find(id + element_count_offset) != newLayerSet.end())
                continue;
            int level = index->element_levels_[id];
            alg_hnsw->element_levels_[id + element_count_offset] = level;
            tableint new_c = id + element_count_offset;
            alg_hnsw->linkLists_[new_c] = (char *)malloc(alg_hnsw->size_links_per_element_ * level + 1);
            if (alg_hnsw->linkLists_[new_c] == nullptr) {
                throw std::runtime_error("Not enough memory: addPoint failed to allocate linklist");
            }
            memset(alg_hnsw->linkLists_[new_c], 0, alg_hnsw->size_links_per_element_ * level + 1);
            alg_hnsw->dist_linkLists_[new_c] = (char *)malloc(alg_hnsw->size_dist_links_per_element_ * level + 1);
            if (alg_hnsw->dist_linkLists_[new_c] == nullptr) {
                throw std::runtime_error("Not enough memory: addPoint failed to allocate dist_linklist");
            }
            memset(alg_hnsw->dist_linkLists_[new_c], 0, alg_hnsw->size_dist_links_per_element_ * level + 1);
        }
    };
    allocateMemory(index1, index1_offset);
    allocateMemory(index2, index2_offset);

    printf("time for new layer: %f\n", elapsed() - s0);

    StopW sw2 = StopW();
    StopW sw3 = StopW();
    StopW sw1 = StopW();

    double sum0 = 0.0, sum1 = 0.0;
    double sum2 = 0.0;
    for (int level = maxLevel; level >= 0; level -= 1) {
        if (layer_node_for_index1[level].size() > 0 && layer_node_for_index2[level].size() == 0) // copy all data for index 1 on this layer to new index
        {
            alg_hnsw->deepCopyOneLayerOnIndex(index1,
                                              level, index1_offset, layer_node_for_index1);
            continue;
        }
        if (layer_node_for_index1[level].size() == 0 && layer_node_for_index2[level].size() > 0) // copy all data for index 2 on this layer to new index
        {
            alg_hnsw->deepCopyOneLayerOnIndex(index2,
                                              level, index2_offset, layer_node_for_index2);
            continue;
        }
        sw2.reset();
        sw3.reset();

        std::vector<std::vector<std::pair<dist_t, tableint>>> candidateSetIndex2(element_count_for_index2);
        std::vector<std::mutex> mtx(element_count_for_index2);

        int batch_size = 64;
#pragma omp parallel for schedule(dynamic, batch_size)
        for (int b = 0; b < layer_node_for_index1[level].size(); b += batch_size) {
            for (int i = b; i < b + batch_size && i < layer_node_for_index1[level].size(); ++i) {
                tableint cur_c = layer_node_for_index1[level][i];
                char *data_point = index1->getDataByInternalId(cur_c);
                tableint &last_entry_point = entry_point_collect_index1_on_index2[cur_c];
                int higherLevel = last_entry_point == -1 ? index2->maxlevel_ : level;
                std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, typename HierarchicalNSW<dist_t>::CompareByFirst> top_candidates;

                linklistsizeint *ll_cur;
                size_t linklistCount;
                tableint *data;
                dist_t *dist;
                tableint candidate_id;
                dist_t dist1;
                dist_t lowerBound;
                size_t Mcurmax = level ? alg_hnsw->maxM_ : alg_hnsw->maxM0_;

                index2->search2Layer(data_point,
                                     last_entry_point,
                                     higherLevel,
                                     level,
                                     lambda,
                                     &top_candidates,
                                     index2_offset);

                auto temp_queue = top_candidates;
                while (!temp_queue.empty()) {
                    std::pair<dist_t, tableint> current = temp_queue.top();
                    temp_queue.pop();
                    std::lock_guard<std::mutex> lk(mtx[current.second - index2_offset]);
                    candidateSetIndex2[current.second - index2_offset].push_back(std::make_pair(current.first, cur_c + index1_offset));
                }

                if (top_candidates.size() + index1->getListCount(index1->get_linklist_at_level(cur_c, level)) > Mcurmax) {
                    ll_cur = index1->get_linklist_at_level(cur_c, level);
                    linklistCount = index1->getListCount(ll_cur);
                    data = (tableint *)(ll_cur + 1);
                    dist = index1->get_dist_at_level(cur_c, level);
                    for (size_t iter = 0; iter < linklistCount; iter++) {
                        candidate_id = data[iter] + index1_offset;
                        dist1 = dist[iter];
                        top_candidates.emplace(dist1, candidate_id);
                    }
                    alg_hnsw->getNeighborsByHeuristic2(top_candidates, Mcurmax, false, -1, alpha);
                    ll_cur = alg_hnsw->get_linklist_at_level(cur_c + index1_offset, level);
                    alg_hnsw->setListCount(ll_cur, top_candidates.size());
                    data = (tableint *)(ll_cur + 1);
                    dist_t *distData = (dist_t *)alg_hnsw->get_dist_at_level(cur_c + index1_offset, level);
                    for (size_t idx = 0; top_candidates.size() > 0; idx++) {
                        data[idx] = top_candidates.top().second;
                        distData[idx] = top_candidates.top().first;
                        top_candidates.pop();
                    }
                } else {
                    linklistsizeint *ll_cur_index1 = index1->get_linklist_at_level(cur_c, level);
                    size_t linklistCount_index1 = index1->getListCount(ll_cur_index1);
                    tableint *data_index1 = (tableint *)(ll_cur_index1 + 1);
                    dist_t *dist_index1 = index1->get_dist_at_level(cur_c, level);

                    ll_cur = alg_hnsw->get_linklist_at_level(cur_c + index1_offset, level);
                    alg_hnsw->setListCount(ll_cur, top_candidates.size() + linklistCount_index1);
                    data = (tableint *)(ll_cur + 1);
                    dist_t *distData = (dist_t *)alg_hnsw->get_dist_at_level(cur_c + index1_offset, level);
                    size_t offset = top_candidates.size();
                    for (size_t idx = 0; top_candidates.size() > 0; idx++) {
                        data[idx] = top_candidates.top().second;
                        distData[idx] = top_candidates.top().first;
                        top_candidates.pop();
                    }
                    for (size_t idx = 0; idx < linklistCount_index1; idx++) {
                        data[offset + idx] = data_index1[idx] + index1_offset;
                        distData[offset + idx] = dist_index1[idx];
                    }
                }
            }
        }

        sum0 += sw2.getElapsedTimeMicro();
        sw2.reset();

#pragma omp parallel for
        for (int b = 0; b < layer_node_for_index2[level].size(); b += batch_size) {
            for (int i = b; i < b + batch_size && i < layer_node_for_index2[level].size(); ++i) {
                tableint cur_c = layer_node_for_index2[level][i];
                char *data_point = index2->getDataByInternalId(cur_c);

                std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, typename HierarchicalNSW<dist_t>::CompareByFirst> top_candidates;
                if (candidateSetIndex2[cur_c].size() > 0) {
                    const std::vector<std::pair<dist_t, tableint>> &candidates = candidateSetIndex2[cur_c];
                    int candidate_size = candidates.size();
                    linklistsizeint *ll_cur;
                    size_t linklistCount;
                    tableint *data;
                    dist_t *dist;
                    tableint candidate_id;
                    dist_t dist1;
                    dist_t lowerBound;
                    size_t Mcurmax = level ? alg_hnsw->maxM_ : alg_hnsw->maxM0_;

                    if (candidate_size + index2->getListCount(index2->get_linklist_at_level(cur_c, level)) > Mcurmax) {
                        for (const auto &p : candidates) {
                            tableint neighbor_id = p.second;
                            dist_t distance = p.first;
                            top_candidates.emplace(distance, neighbor_id);
                        }

                        ll_cur = index2->get_linklist_at_level(cur_c, level);
                        linklistCount = index2->getListCount(ll_cur);
                        data = (tableint *)(ll_cur + 1);
                        dist = index2->get_dist_at_level(cur_c, level);
                        for (size_t iter = 0; iter < linklistCount; iter++) {
                            candidate_id = data[iter] + index2_offset;
                            dist1 = dist[iter];
                            top_candidates.emplace(dist1, candidate_id);
                        }
                        alg_hnsw->getNeighborsByHeuristic2(top_candidates, Mcurmax, false, -1, alpha);
                        ll_cur = alg_hnsw->get_linklist_at_level(cur_c + index2_offset, level);
                        alg_hnsw->setListCount(ll_cur, top_candidates.size());
                        data = (tableint *)(ll_cur + 1);
                        dist_t *distData = (dist_t *)alg_hnsw->get_dist_at_level(cur_c + index2_offset, level);
                        for (size_t idx = 0; top_candidates.size() > 0; idx++) {
                            data[idx] = top_candidates.top().second;
                            distData[idx] = top_candidates.top().first;
                            top_candidates.pop();
                        }
                    } else {
                        linklistsizeint *ll_cur_index2 = index2->get_linklist_at_level(cur_c, level);
                        size_t linklistCount_index2 = index2->getListCount(ll_cur_index2);
                        tableint *data_index2 = (tableint *)(ll_cur_index2 + 1);
                        dist_t *dist_index2 = index2->get_dist_at_level(cur_c, level);

                        ll_cur = alg_hnsw->get_linklist_at_level(cur_c + index2_offset, level);
                        alg_hnsw->setListCount(ll_cur, candidate_size + linklistCount_index2);
                        data = (tableint *)(ll_cur + 1);
                        dist_t *distData = (dist_t *)alg_hnsw->get_dist_at_level(cur_c + index2_offset, level);
                        size_t offset = candidate_size;
                        for (size_t idx = 0; idx < candidate_size; idx++) {
                            data[idx] = candidates[idx].second;
                            distData[idx] = candidates[idx].first;
                        }
                        for (size_t idx = 0; idx < linklistCount_index2; idx++) {
                            data[offset + idx] = data_index2[idx] + index2_offset;
                            distData[offset + idx] = dist_index2[idx];
                        }
                    }
                } else {
                    linklistsizeint *ll_cur_index2 = index2->get_linklist_at_level(cur_c, level);
                    size_t linklistCount_index2 = index2->getListCount(ll_cur_index2);
                    tableint *data_index2 = (tableint *)(ll_cur_index2 + 1);
                    dist_t *dist_index2 = index2->get_dist_at_level(cur_c, level);

                    linklistsizeint *ll_cur;
                    size_t linklistCount;
                    tableint *data;
                    dist_t *dist;
                    tableint candidate_id;
                    dist_t dist1;
                    dist_t lowerBound;
                    size_t Mcurmax = level ? alg_hnsw->maxM_ : alg_hnsw->maxM0_;

                    ll_cur = alg_hnsw->get_linklist_at_level(cur_c + index2_offset, level);
                    alg_hnsw->setListCount(ll_cur, linklistCount_index2);
                    data = (tableint *)(ll_cur + 1);
                    dist_t *distData = (dist_t *)alg_hnsw->get_dist_at_level(cur_c + index2_offset, level);

                    for (size_t idx = 0; idx < linklistCount_index2; idx++) {
                        data[idx] = data_index2[idx] + index2_offset;
                        distData[idx] = dist_index2[idx];
                    }
                }
            }
        }

        sum1 += sw2.getElapsedTimeMicro();
        sum2 += sw3.getElapsedTimeMicro();
    }
    printf("mergeIndex1BasedOnIndex2Connection sum0: %f\n", sum0 / 1000000);
    printf("mergeIndex1BasedOnIndex2Connection sum1: %f\n", sum1 / 1000000);
    printf("mergeIndex1BasedOnIndex2Connection sum2: %f\n", sum2 / 1000000);
    global_counter += sum2 / 1000000;
    return alg_hnsw;
}

// Backward Search comparison algorithm for ablation study
template <typename dist_t>
HierarchicalNSW<dist_t> *HNSWMerger_BS(HierarchicalNSW<dist_t> *index1, HierarchicalNSW<dist_t> *index2, L2Space *space, int lambda = 4, float alpha = 1.05, size_t M = -1, size_t ef_construction = -1) {

    double s0 = elapsed();
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
    HierarchicalNSW<dist_t> *alg_hnsw = new HierarchicalNSW<dist_t>(space, new_max_elements, M, ef_construction);
    size_t maxLevel = std::max(index1->maxlevel_, index2->maxlevel_);
    alg_hnsw->setMaxLevel(maxLevel);
    alg_hnsw->cur_element_count.store(index1->cur_element_count + index2->cur_element_count);

    alg_hnsw->data_level0_memory_ = (char *)malloc(alg_hnsw->max_elements_ * alg_hnsw->size_data_per_element_);
    if (alg_hnsw->data_level0_memory_ == nullptr)
        throw std::runtime_error("Not enough memory: loadIndex failed to allocate level0");

    alg_hnsw->dist_level0_memory_ = (char *)malloc(alg_hnsw->max_elements_ * alg_hnsw->size_dist_per_element_);
    if (alg_hnsw->dist_level0_memory_ == nullptr)
        throw std::runtime_error("Not enough memory: loadIndex failed to allocate level0");


    std::vector<std::vector<int>> layer_node_for_index1(maxLevel + 2);
    index1->searchNodeOnEachLayer(layer_node_for_index1, true);
    std::vector<std::vector<int>> layer_node_for_index2(maxLevel + 2);
    index2->searchNodeOnEachLayer(layer_node_for_index2, true);

    tableint *entry_point_collect_index1_on_index2 = new tableint[element_count_for_index1];
    tableint *entry_point_collect_index2_on_index1 = new tableint[element_count_for_index2];
    memset(entry_point_collect_index1_on_index2, -1, element_count_for_index1 * sizeof(int));
    memset(entry_point_collect_index2_on_index1, -1, element_count_for_index2 * sizeof(int));

    int entry_point_index1 = index2->enterpoint_node_;
    int entry_point_index2 = index1->enterpoint_node_;

    memcpy(alg_hnsw->data_level0_memory_,
           index1->data_level0_memory_,
           element_count_for_index1 * alg_hnsw->size_data_per_element_);
    memcpy(alg_hnsw->data_level0_memory_ + element_count_for_index1 * alg_hnsw->size_data_per_element_,
           index2->data_level0_memory_,
           (element_count_for_index2)*alg_hnsw->size_data_per_element_);
    memcpy(alg_hnsw->dist_level0_memory_,
           index1->dist_level0_memory_,
           element_count_for_index1 * alg_hnsw->size_dist_per_element_);
    memcpy(alg_hnsw->dist_level0_memory_ + element_count_for_index1 * alg_hnsw->size_dist_per_element_,
           index2->dist_level0_memory_,
           (element_count_for_index2)*alg_hnsw->size_dist_per_element_);

    std::vector<int> newLayer;
    if (layer_node_for_index1[maxLevel].size() + layer_node_for_index2[maxLevel].size() > M) {
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        int lambda = 0;
        auto filter_and_remove = [&](std::vector<int> &layer_nodes, HierarchicalNSW<dist_t> *index,
                                     int offset) {
            auto it = layer_nodes.begin();
            while (it != layer_nodes.end()) {
                if (distribution(alg_hnsw->level_generator_) < 1.0 / M || lambda == M) {
                    newLayer.push_back((*it) + offset);
                    lambda = 0;
                    tableint new_c = *it;
                    alg_hnsw->linkLists_[new_c + offset] = (char *)malloc(alg_hnsw->size_links_per_element_ * (maxLevel + 1) + 1);
                    if (alg_hnsw->linkLists_[new_c + offset] == nullptr)
                        throw std::runtime_error("Not enough memory: addPoint failed to allocate linklist");
                    memset(alg_hnsw->linkLists_[new_c + offset], 0, alg_hnsw->size_links_per_element_ * (maxLevel + 1) + 1);
                    alg_hnsw->element_levels_[new_c + offset] = (maxLevel + 1);

                    alg_hnsw->dist_linkLists_[new_c + offset] = (char *)malloc(alg_hnsw->size_dist_links_per_element_ * (maxLevel + 1) + 1);
                    if (alg_hnsw->dist_linkLists_[new_c + offset] == nullptr)
                        throw std::runtime_error("Not enough memory: addPoint failed to "
                                                 "allocate dist_linkLists");
                    memset(alg_hnsw->dist_linkLists_[new_c + offset], 0, alg_hnsw->size_dist_links_per_element_ * (maxLevel + 1) + 1);

                    memcpy(alg_hnsw->getDataByInternalId(new_c + offset), index->getDataByInternalId(new_c), alg_hnsw->data_size_);
                    alg_hnsw->setExternalLabel(new_c + offset, index->getExternalLabel(new_c));

                    it = layer_nodes.erase(it);
                } else {
                    lambda++;
                    ++it;
                }
            }
        };
        filter_and_remove(layer_node_for_index1[maxLevel], index1,
                          index1_offset);
        filter_and_remove(layer_node_for_index2[maxLevel], index2,
                          index2_offset);
        if (newLayer.size() > 0) { 
            for (int newLayerIter = 0; newLayerIter < newLayer.size(); newLayerIter++) {
                linklistsizeint *ll_cur = alg_hnsw->get_linklist_at_level(newLayer[newLayerIter], maxLevel + 1);
                alg_hnsw->setListCount(ll_cur, newLayer.size() - 1);
                tableint *data = (tableint *)(ll_cur + 1);
                dist_t *dist = alg_hnsw->get_dist_at_level(newLayer[newLayerIter], maxLevel + 1);
                for (size_t it = 0, idx = 0; it < newLayer.size(); it++) {
                    if (it == newLayerIter)
                        continue;
                    dist_t dist0 =
                        alg_hnsw->fstdistfunc_(alg_hnsw->getDataByInternalId(newLayer[newLayerIter]), alg_hnsw->getDataByInternalId(newLayer[it]), alg_hnsw->dist_func_param_);
                    data[idx] = newLayer[it];
                    dist[idx] = dist0;
                    idx++;
                }
            }
            alg_hnsw->setMaxLevel(maxLevel + 1);
            alg_hnsw->enterpoint_node_ = newLayer[0];
        } else {
            alg_hnsw->enterpoint_node_ = layer_node_for_index2[maxLevel][0] + index2_offset;
        }
    } 
    else {
        alg_hnsw->enterpoint_node_ = layer_node_for_index2[maxLevel][0] + index2_offset;
    }

    alg_hnsw->label_lookup_.clear();
    alg_hnsw->label_lookup_.reserve(alg_hnsw->max_elements_);

    for (int id = 0; id < element_count_for_index1; id++) {
        alg_hnsw->label_lookup_[index1->getExternalLabel(id)] = id + index1_offset;
    }
    for (int id = 0; id < element_count_for_index2; id++) {
        alg_hnsw->label_lookup_[index2->getExternalLabel(id)] = id + index2_offset;
    }

    std::unordered_set<int> newLayerSet(newLayer.begin(), newLayer.end());
    auto allocateMemory = [&](hnswlib::HierarchicalNSW<dist_t> *index, int element_count_offset) {
        for (int id = 0; id < index->cur_element_count; id++) {
            if (index->element_levels_[id] < 1 || newLayerSet.find(id + element_count_offset) != newLayerSet.end())
                continue;
            int level = index->element_levels_[id];
            alg_hnsw->element_levels_[id + element_count_offset] = level;
            tableint new_c = id + element_count_offset;
            alg_hnsw->linkLists_[new_c] = (char *)malloc(alg_hnsw->size_links_per_element_ * level + 1);
            if (alg_hnsw->linkLists_[new_c] == nullptr) {
                throw std::runtime_error("Not enough memory: addPoint failed to allocate linklist");
            }
            memset(alg_hnsw->linkLists_[new_c], 0, alg_hnsw->size_links_per_element_ * level + 1);
            alg_hnsw->dist_linkLists_[new_c] = (char *)malloc(alg_hnsw->size_dist_links_per_element_ * level + 1);
            if (alg_hnsw->dist_linkLists_[new_c] == nullptr) {
                throw std::runtime_error("Not enough memory: addPoint failed to allocate dist_linklist");
            }
            memset(alg_hnsw->dist_linkLists_[new_c], 0, alg_hnsw->size_dist_links_per_element_ * level + 1);
        }
    };
    allocateMemory(index1, index1_offset);
    allocateMemory(index2, index2_offset);

    printf("time for new layer: %f\n", elapsed() - s0);

    StopW sw2 = StopW();
    StopW sw3 = StopW();
    StopW sw1 = StopW();

    double sum0 = 0.0, sum1 = 0.0;
    double sum2 = 0.0;
    for (int level = maxLevel; level >= 0; level -= 1) {
        if (layer_node_for_index1[level].size() > 0 && layer_node_for_index2[level].size() == 0) // copy all data for index 1 on this layer to new index
        {
            alg_hnsw->deepCopyOneLayerOnIndex(index1,
                                              level, index1_offset, layer_node_for_index1);
            continue;
        }
        if (layer_node_for_index1[level].size() == 0 && layer_node_for_index2[level].size() > 0) // copy all data for index 2 on this layer to new index
        {
            alg_hnsw->deepCopyOneLayerOnIndex(index2,
                                              level, index2_offset, layer_node_for_index2);
            continue;
        }
        sw2.reset();
        sw3.reset();

        int batch_size = 64;
#pragma omp parallel for schedule(dynamic, batch_size)
        for (int b = 0; b < layer_node_for_index1[level].size(); b += batch_size) {
            for (int i = b; i < b + batch_size && i < layer_node_for_index1[level].size(); ++i) {
                int tid = omp_get_thread_num();
                tableint cur_c = layer_node_for_index1[level][i];
                char *data_point = index1->getDataByInternalId(cur_c);
                tableint &last_entry_point = entry_point_collect_index1_on_index2[cur_c];
                int higherLevel = last_entry_point == -1 ? index2->maxlevel_ : level;
                std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, typename HierarchicalNSW<dist_t>::CompareByFirst> top_candidates;

                linklistsizeint *ll_cur;
                size_t linklistCount;
                tableint *data;
                dist_t *dist;
                tableint candidate_id;
                dist_t dist1;
                dist_t lowerBound;
                size_t Mcurmax = level ? alg_hnsw->maxM_ : alg_hnsw->maxM0_;

                index2->search2Layer(data_point,
                                     last_entry_point,
                                     higherLevel,
                                     level,
                                     lambda,
                                     &top_candidates,
                                     index2_offset);

                if (top_candidates.size() + index1->getListCount(index1->get_linklist_at_level(cur_c, level)) > Mcurmax) {
                    ll_cur = index1->get_linklist_at_level(cur_c, level);
                    linklistCount = index1->getListCount(ll_cur);
                    data = (tableint *)(ll_cur + 1);
                    dist = index1->get_dist_at_level(cur_c, level);
                    for (size_t iter = 0; iter < linklistCount; iter++) {
                        candidate_id = data[iter] + index1_offset;
                        dist1 = dist[iter];
                        top_candidates.emplace(dist1, candidate_id);
                    }
                    alg_hnsw->getNeighborsByHeuristic2(top_candidates, Mcurmax, false, -1, alpha);
                    ll_cur = alg_hnsw->get_linklist_at_level(cur_c + index1_offset, level);
                    alg_hnsw->setListCount(ll_cur, top_candidates.size());
                    data = (tableint *)(ll_cur + 1);
                    dist_t *distData = (dist_t *)alg_hnsw->get_dist_at_level(cur_c + index1_offset, level);
                    for (size_t idx = 0; top_candidates.size() > 0; idx++) {
                        data[idx] = top_candidates.top().second;
                        distData[idx] = top_candidates.top().first;
                        top_candidates.pop();
                    }
                } else {
                    linklistsizeint *ll_cur_index1 = index1->get_linklist_at_level(cur_c, level);
                    size_t linklistCount_index1 = index1->getListCount(ll_cur_index1);
                    tableint *data_index1 = (tableint *)(ll_cur_index1 + 1);
                    dist_t *dist_index1 = index1->get_dist_at_level(cur_c, level);

                    ll_cur = alg_hnsw->get_linklist_at_level(cur_c + index1_offset, level);
                    alg_hnsw->setListCount(ll_cur, top_candidates.size() + linklistCount_index1);
                    data = (tableint *)(ll_cur + 1);
                    dist_t *distData = (dist_t *)alg_hnsw->get_dist_at_level(cur_c + index1_offset, level);
                    size_t offset = top_candidates.size();
                    for (size_t idx = 0; top_candidates.size() > 0; idx++) {
                        data[idx] = top_candidates.top().second;
                        distData[idx] = top_candidates.top().first;
                        top_candidates.pop();
                    }
                    for (size_t idx = 0; idx < linklistCount_index1; idx++) {
                        data[offset + idx] = data_index1[idx] + index1_offset;
                        distData[offset + idx] = dist_index1[idx];
                    }
                }
            }
        }

        sum0 += sw2.getElapsedTimeMicro();
        sw2.reset();

#pragma omp parallel for schedule(dynamic, batch_size)
        for (int b = 0; b < layer_node_for_index2[level].size(); b += batch_size) {
            for (int i = b; i < b + batch_size && i < layer_node_for_index2[level].size(); ++i) {
                int tid = omp_get_thread_num();
                tableint cur_c = layer_node_for_index2[level][i];
                char *data_point = index2->getDataByInternalId(cur_c);
                tableint &last_entry_point = entry_point_collect_index2_on_index1[cur_c];
                int higherLevel = last_entry_point == -1 ? index1->maxlevel_ : level;
                std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, typename HierarchicalNSW<dist_t>::CompareByFirst> top_candidates;

                linklistsizeint *ll_cur;
                size_t linklistCount;
                tableint *data;
                dist_t *dist;
                tableint candidate_id;
                dist_t dist1;
                dist_t lowerBound;
                size_t Mcurmax = level ? alg_hnsw->maxM_ : alg_hnsw->maxM0_;

                index1->search2Layer(data_point,
                                     last_entry_point,
                                     higherLevel,
                                     level,
                                     lambda,
                                     &top_candidates,
                                     index1_offset);

                if (top_candidates.size() + index2->getListCount(index2->get_linklist_at_level(cur_c, level)) > Mcurmax) {
                    ll_cur = index2->get_linklist_at_level(cur_c, level);
                    linklistCount = index2->getListCount(ll_cur);
                    data = (tableint *)(ll_cur + 1);
                    dist = index2->get_dist_at_level(cur_c, level);
                    for (size_t iter = 0; iter < linklistCount; iter++) {
                        candidate_id = data[iter] + index2_offset;
                        dist1 = dist[iter];
                        top_candidates.emplace(dist1, candidate_id);
                    }
                    alg_hnsw->getNeighborsByHeuristic2(top_candidates, Mcurmax, false, -1, alpha);
                    ll_cur = alg_hnsw->get_linklist_at_level(cur_c + index2_offset, level);
                    alg_hnsw->setListCount(ll_cur, top_candidates.size());
                    data = (tableint *)(ll_cur + 1);
                    dist_t *distData = (dist_t *)alg_hnsw->get_dist_at_level(cur_c + index2_offset, level);
                    for (size_t idx = 0; top_candidates.size() > 0; idx++) {
                        data[idx] = top_candidates.top().second;
                        distData[idx] = top_candidates.top().first;
                        top_candidates.pop();
                    }
                } else {
                    linklistsizeint *ll_cur_index1 = index2->get_linklist_at_level(cur_c, level);
                    size_t linklistCount_index1 = index2->getListCount(ll_cur_index1);
                    tableint *data_index1 = (tableint *)(ll_cur_index1 + 1);
                    dist_t *dist_index1 = index2->get_dist_at_level(cur_c, level);

                    ll_cur = alg_hnsw->get_linklist_at_level(cur_c + index2_offset, level);
                    alg_hnsw->setListCount(ll_cur, top_candidates.size() + linklistCount_index1);
                    data = (tableint *)(ll_cur + 1);
                    dist_t *distData = (dist_t *)alg_hnsw->get_dist_at_level(cur_c + index2_offset, level);
                    size_t offset = top_candidates.size();
                    for (size_t idx = 0; top_candidates.size() > 0; idx++) {
                        data[idx] = top_candidates.top().second;
                        distData[idx] = top_candidates.top().first;
                        top_candidates.pop();
                    }
                    for (size_t idx = 0; idx < linklistCount_index1; idx++) {
                        data[offset + idx] = data_index1[idx] + index2_offset;
                        distData[offset + idx] = dist_index1[idx];
                    }
                }
            }
        }

        sum1 += sw2.getElapsedTimeMicro();
        sum2 += sw3.getElapsedTimeMicro();
    }
    printf("mergeIndex1BasedOnIndex2Connection sum0: %f\n", sum0 / 1000000);
    printf("mergeIndex1BasedOnIndex2Connection sum1: %f\n", sum1 / 1000000);
    printf("mergeIndex1BasedOnIndex2Connection sum2: %f\n", sum2 / 1000000);
    global_counter += sum2 / 1000000;
    return alg_hnsw;
}

} // namespace hnswlib