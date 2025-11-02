/*
 * This baseline algorithm is based on https://arxiv.org/abs/2505.16064
 * Code is converted from authors' Github repository https://github.com/aponom84/merging-navigable-graphs
 * For fair comparing with our algorithm, we only convert key function from python code,
 * and use the original HNSWlib code as much as possible.
 */
#include "extension.h"
#include <mutex>
#include <omp.h>

using namespace std;

namespace hnswlib {

template <typename dist_t>
HierarchicalNSW<dist_t> *BasicMerge(HierarchicalNSW<dist_t> *index1, HierarchicalNSW<dist_t> *index2, L2Space *space, size_t M = -1, size_t ef_construction = -1) {
    double s0 = elapsed();
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
    index1->searchNodeOnEachLayer(layer_node_for_index1);
    std::vector<std::vector<int>> layer_node_for_index2(maxLevel + 2);
    index2->searchNodeOnEachLayer(layer_node_for_index2);

    alg_hnsw->enterpoint_node_ = layer_node_for_index2[maxLevel][0] + index2_offset;

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
            if (index->element_levels_[id] < 1)
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
    printf("time for new layer: %f\n", elapsed() - s0);

    return alg_hnsw;
}

template <typename dist_t>
HierarchicalNSW<dist_t> *HNSWMerger_Naive(HierarchicalNSW<dist_t> *index1, HierarchicalNSW<dist_t> *index2, L2Space *space, size_t search_ef = 40, size_t M = -1, size_t ef_construction = -1) {
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
    alg_hnsw = BasicMerge(index1, index2, space, M, ef_construction);

    size_t maxLevel = alg_hnsw->maxlevel_;

    std::vector<std::vector<int>> layer_node_for_index1(maxLevel + 2);
    index1->searchNodeOnEachLayer(layer_node_for_index1, true);
    std::vector<std::vector<int>> layer_node_for_index2(maxLevel + 2);
    index2->searchNodeOnEachLayer(layer_node_for_index2, true);


    for (int level = maxLevel; level >= 0; level -= 1) {
        if (layer_node_for_index1[level].size() > 0 && layer_node_for_index2[level].size() == 0) // copy all data for index 1 on this layer to new index
        {
            alg_hnsw->deepCopyOneLayerOnIndex(index1, level, index1_offset, layer_node_for_index1);
            continue;
        }
        if (layer_node_for_index1[level].size() == 0 && layer_node_for_index2[level].size() > 0) // copy all data for index 2 on this layer to new index
        {
            alg_hnsw->deepCopyOneLayerOnIndex(index2, level, index2_offset, layer_node_for_index2);
            continue;
        }

        linklistsizeint *ll_cur;
        size_t linklistCount;
        tableint *data;
        dist_t *dist;
        size_t Mcurmax = level ? alg_hnsw->maxM_ : alg_hnsw->maxM0_;

        for (int iter = 0; iter < layer_node_for_index1[level].size(); iter++) {
            tableint cur_c = layer_node_for_index1[level][iter];
            char *data_point = index1->getDataByInternalId(cur_c);

            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, typename HierarchicalNSW<dist_t>::CompareByFirst> top_candidates;

            tableint entry_point = -1;

            // We follow the idea in their paper, using HNSW search in index construction phase, not limiting the search to M neighbors
            index2->search2Layer(data_point, entry_point, maxLevel, level, search_ef, &top_candidates, index2_offset);
            ll_cur = index1->get_linklist_at_level(cur_c, level);
            linklistCount = index1->getListCount(ll_cur);
            data = (tableint *)(ll_cur + 1);
            dist = index1->get_dist_at_level(cur_c, level);
            for (size_t iter = 0; iter < linklistCount; iter++) {
                top_candidates.emplace(dist[iter], data[iter] + index1_offset);
            }

            alg_hnsw->getNeighborsByHeuristic2(top_candidates, Mcurmax, false);
            ll_cur = alg_hnsw->get_linklist_at_level(cur_c + index1_offset, level);
            alg_hnsw->setListCount(ll_cur, top_candidates.size());
            data = (tableint *)(ll_cur + 1);
            dist = (dist_t *)alg_hnsw->get_dist_at_level(cur_c + index1_offset, level);
            for (size_t idx = 0; top_candidates.size() > 0; idx++) {
                data[idx] = top_candidates.top().second;
                dist[idx] = top_candidates.top().first;
                top_candidates.pop();
            }
        }

        for (int iter = 0; iter < layer_node_for_index2[level].size(); iter++) {
            tableint cur_c = layer_node_for_index2[level][iter];
            char *data_point = index2->getDataByInternalId(cur_c);

            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, typename HierarchicalNSW<dist_t>::CompareByFirst> top_candidates;

            tableint entry_point = -1;

            // We follow the idea in their paper, using HNSW search in index construction phase, not limiting the search to M neighbors
            index1->search2Layer(data_point, entry_point, maxLevel, level, search_ef, &top_candidates, index1_offset);
            ll_cur = index2->get_linklist_at_level(cur_c, level);
            linklistCount = index2->getListCount(ll_cur);
            data = (tableint *)(ll_cur + 1);
            dist = index2->get_dist_at_level(cur_c, level);
            for (size_t iter = 0; iter < linklistCount; iter++) {
                top_candidates.emplace(dist[iter], data[iter] + index2_offset);
            }

            alg_hnsw->getNeighborsByHeuristic2(top_candidates, Mcurmax, false);
            ll_cur = alg_hnsw->get_linklist_at_level(cur_c + index2_offset, level);
            alg_hnsw->setListCount(ll_cur, top_candidates.size());
            data = (tableint *)(ll_cur + 1);
            dist = (dist_t *)alg_hnsw->get_dist_at_level(cur_c + index2_offset, level);
            for (size_t idx = 0; top_candidates.size() > 0; idx++) {
                data[idx] = top_candidates.top().second;
                dist[idx] = top_candidates.top().first;
                top_candidates.pop();
            }
        }
    }
    return alg_hnsw;
}

template <typename dist_t>
HierarchicalNSW<dist_t> *HNSWMerger_IGTM(HierarchicalNSW<dist_t> *index1, HierarchicalNSW<dist_t> *index2, L2Space *space,
                                         size_t jump_ef = 40, size_t local_ef = 10, size_t next_step_k = 6, size_t next_step_ef = 6, size_t search_M = 40,
                                         size_t M = -1, size_t ef_construction = -1) {
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
    alg_hnsw = BasicMerge(index1, index2, space, M, ef_construction);

    size_t maxLevel = alg_hnsw->maxlevel_;

    std::vector<std::vector<int>> layer_node_for_index1(maxLevel + 2);
    index1->searchNodeOnEachLayer(layer_node_for_index1, true);
    std::vector<std::vector<int>> layer_node_for_index2(maxLevel + 2);
    index2->searchNodeOnEachLayer(layer_node_for_index2, true);


    for (int level = maxLevel; level >= 0; level -= 1) {
        if (layer_node_for_index1[level].size() > 0 && layer_node_for_index2[level].size() == 0) // copy all data for index 1 on this layer to new index
        {
            alg_hnsw->deepCopyOneLayerOnIndex(index1, level, index1_offset, layer_node_for_index1);
            continue;
        }
        if (layer_node_for_index1[level].size() == 0 && layer_node_for_index2[level].size() > 0) // copy all data for index 2 on this layer to new index
        {
            alg_hnsw->deepCopyOneLayerOnIndex(index2, level, index2_offset, layer_node_for_index2);
            continue;
        }

        linklistsizeint *ll_cur;
        size_t linklistCount;
        tableint *data;
        dist_t *dist;
        size_t Mcurmax = level ? alg_hnsw->maxM_ : alg_hnsw->maxM0_;

        auto merge_all = [&](
                             HierarchicalNSW<dist_t> *index1,
                             HierarchicalNSW<dist_t> *index2,
                             const std::vector<std::vector<int>> &layer_node_for_index1,
                             int index1_offset,
                             int index2_offset) {
            std::unordered_set<tableint> not_done1(layer_node_for_index1[level].begin(), layer_node_for_index1[level].end());
            int global_count = 0;

            while (!not_done1.empty()) {
                tableint cur_c = *not_done1.begin();
                not_done1.erase(not_done1.begin());

                char *data_point = index1->getDataByInternalId(cur_c);

                std::priority_queue<
                    std::pair<dist_t, tableint>,
                    std::vector<std::pair<dist_t, tableint>>,
                    typename HierarchicalNSW<dist_t>::CompareByFirst>
                    starting_points;
                std::unordered_set<tableint> starting_eps;

                tableint entry_point = -1;
                index2->search2Layer(
                    data_point,
                    entry_point,
                    maxLevel, 
                    level, 
                    jump_ef,
                    &starting_points,
                    0);
                while (starting_points.size() > (size_t)search_M)
                    starting_points.pop();
                while (!starting_points.empty()) {
                    starting_eps.emplace(starting_points.top().second);
                    starting_points.pop();
                }

                while (true) {
                    if ((++global_count % 200000) == 0) {
                        printf(
                            "[log] Index merge process %d / %zu, time: %f \n",
                            global_count,
                            layer_node_for_index1[level].size(),
                            elapsed() - s0);
                    }

                    std::priority_queue<
                        std::pair<dist_t, tableint>,
                        std::vector<std::pair<dist_t, tableint>>,
                        typename HierarchicalNSW<dist_t>::CompareByFirst>
                        origin = index2->ExtendSearchBaseLayer(
                            data_point,
                            level,
                            &starting_eps,
                            local_ef);
                    starting_eps.clear();

                    while (origin.size() > (size_t)Mcurmax) {
                        origin.pop();
                    }

                    std::priority_queue<
                        std::pair<dist_t, tableint>,
                        std::vector<std::pair<dist_t, tableint>>,
                        typename HierarchicalNSW<dist_t>::CompareByFirst>
                        top_candidates;
                    while (!origin.empty()) {
                        top_candidates.emplace(
                            origin.top().first,
                            origin.top().second + index2_offset);
                        if (origin.size() <= (size_t)search_M)
                            starting_eps.emplace(origin.top().second);
                        origin.pop();
                    }

                    ll_cur = index1->get_linklist_at_level(cur_c, level);
                    linklistCount = index1->getListCount(ll_cur);
                    data = reinterpret_cast<tableint *>(ll_cur + 1);
                    dist = reinterpret_cast<dist_t *>(
                        index1->get_dist_at_level(cur_c, level));
                    for (size_t i = 0; i < linklistCount; ++i) {
                        top_candidates.emplace(
                            dist[i], data[i] + index1_offset);
                    }

                    alg_hnsw->getNeighborsByHeuristic2(
                        top_candidates, Mcurmax, false);
                    ll_cur = alg_hnsw->get_linklist_at_level(
                        cur_c + index1_offset, level);
                    alg_hnsw->setListCount(ll_cur, top_candidates.size());
                    data = reinterpret_cast<tableint *>(ll_cur + 1);
                    dist = reinterpret_cast<dist_t *>(
                        alg_hnsw->get_dist_at_level(
                            cur_c + index1_offset, level));
                    for (size_t i = 0; !top_candidates.empty(); ++i) {
                        data[i] = top_candidates.top().second;
                        dist[i] = top_candidates.top().first;
                        top_candidates.pop();
                    }

                    std::priority_queue<
                        std::pair<dist_t, tableint>,
                        std::vector<std::pair<dist_t, tableint>>,
                        typename HierarchicalNSW<dist_t>::CompareByFirst>
                        next_candidate;
                    index1->search2Layer(
                        data_point, cur_c, level, level,
                        next_step_ef, &next_candidate, 0);

                    while (next_candidate.size() > next_step_k)
                        next_candidate.pop();

                    tableint new_c = -1;
                    while (!next_candidate.empty()) {
                        auto top = next_candidate.top();
                        if (not_done1.find(top.second) != not_done1.end())
                            new_c = top.second;
                        next_candidate.pop();
                    }
                    if (new_c == -1) break;
                    cur_c = new_c;
                    not_done1.erase(cur_c);
                    data_point = index1->getDataByInternalId(cur_c);
                }
            }
        };
        merge_all(index1, index2, layer_node_for_index1, index1_offset, index2_offset);
        merge_all(index2, index1, layer_node_for_index2, index2_offset, index1_offset);
    }
    return alg_hnsw;
}

template <typename dist_t>
HierarchicalNSW<dist_t> *HNSWMerger_CGTM(HierarchicalNSW<dist_t> *index1, HierarchicalNSW<dist_t> *index2, L2Space *space,
                                         size_t jump_ef = 40, size_t local_ef = 10, size_t next_step_k = 6, size_t search_M = 40,
                                         size_t M = -1, size_t ef_construction = -1) {
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
    alg_hnsw = BasicMerge(index1, index2, space, M, ef_construction);

    size_t maxLevel = alg_hnsw->maxlevel_;

    std::vector<std::vector<int>> layer_node_for_index1(maxLevel + 2);
    index1->searchNodeOnEachLayer(layer_node_for_index1, true);
    std::vector<std::vector<int>> layer_node_for_index2(maxLevel + 2);
    index2->searchNodeOnEachLayer(layer_node_for_index2, true);

    for (int level = maxLevel; level >= 0; level -= 1) {
        // printf("[log] layer: %d\n", level);
        if (layer_node_for_index1[level].size() > 0 && layer_node_for_index2[level].size() == 0) // copy all data for index 1 on this layer to new index
        {
            alg_hnsw->deepCopyOneLayerOnIndex(index1, level, index1_offset, layer_node_for_index1);
            continue;
        }
        if (layer_node_for_index1[level].size() == 0 && layer_node_for_index2[level].size() > 0) // copy all data for index 2 on this layer to new index
        {
            alg_hnsw->deepCopyOneLayerOnIndex(index2, level, index2_offset, layer_node_for_index2);
            continue;
        }

        linklistsizeint *ll_cur;
        size_t linklistCount;
        tableint *data;
        dist_t *dist;
        size_t Mcurmax = level ? alg_hnsw->maxM_ : alg_hnsw->maxM0_;

        std::unordered_set<tableint> not_done(layer_node_for_index1[level].begin(), layer_node_for_index1[level].end());
        for (int iter = 0; iter < layer_node_for_index2[level].size(); iter++) {
            not_done.emplace(layer_node_for_index2[level][iter] + index2_offset);
        }

        int global_count = 0;
        while (!not_done.empty()) {
            tableint cur_c = *not_done.begin();

            char *data_point = alg_hnsw->getDataByInternalId(cur_c);
            tableint entry_point = -1;
            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, typename HierarchicalNSW<dist_t>::CompareByFirst> starting_points1;
            std::unordered_set<tableint> starting_eps1;
            index1->search2Layer(data_point, entry_point, maxLevel, level, jump_ef, &starting_points1, 0);
            while (starting_points1.size() > 0) {
                if (starting_points1.size() <= search_M)
                    starting_eps1.emplace(starting_points1.top().second);
                starting_points1.pop();
            }

            entry_point = -1;
            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, typename HierarchicalNSW<dist_t>::CompareByFirst> starting_points2;
            std::unordered_set<tableint> starting_eps2;
            index2->search2Layer(data_point, entry_point, maxLevel, level, jump_ef, &starting_points2, 0);
            while (starting_points2.size() > 0) {
                if (starting_points2.size() <= search_M)
                    starting_eps2.emplace(starting_points2.top().second);
                starting_points2.pop();
            }

            while (true) {
                not_done.erase(cur_c);
                data_point = alg_hnsw->getDataByInternalId(cur_c);

                if ((++global_count) % 200000 == 0) {
                    printf("[log] Index merge process %d / %d, time %f\n", global_count, layer_node_for_index1[level].size() + layer_node_for_index2[level].size(), elapsed() - s0);
                }
                tableint new_c = -1;
                dist_t distBound = std::numeric_limits<dist_t>::max();

                std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, typename HierarchicalNSW<dist_t>::CompareByFirst>
                    origin1;
                std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, typename HierarchicalNSW<dist_t>::CompareByFirst> top_candidates1;
                origin1 = index1->ExtendSearchBaseLayer(data_point, level, &starting_eps1, local_ef);
                starting_eps1.clear();
                while (origin1.size() > (size_t)Mcurmax) {
                    origin1.pop();
                }
                while (!origin1.empty()) {
                    top_candidates1.emplace(origin1.top().first, origin1.top().second + index1_offset);
                    starting_eps1.emplace(origin1.top().second);
                    if (origin1.size() <= next_step_k && not_done.find(origin1.top().second + index1_offset) != not_done.end() && origin1.top().first < distBound) {
                        distBound = origin1.top().first;
                        new_c = origin1.top().second + index1_offset;
                    }
                    origin1.pop();
                }

                std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, typename HierarchicalNSW<dist_t>::CompareByFirst> origin2;
                std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, typename HierarchicalNSW<dist_t>::CompareByFirst> top_candidates2;
                origin2 = index2->ExtendSearchBaseLayer(data_point, level, &starting_eps2, local_ef);
                starting_eps2.clear();
                while (origin2.size() > (size_t)Mcurmax) {
                    origin2.pop();
                }
                while (!origin2.empty()) {
                    top_candidates2.emplace(origin2.top().first, origin2.top().second + index2_offset);
                    starting_eps2.emplace(origin2.top().second);
                    if (origin2.size() <= next_step_k && not_done.find(origin2.top().second + index2_offset) != not_done.end() && origin2.top().first < distBound) {
                        distBound = origin2.top().first;
                        new_c = origin2.top().second + index2_offset;
                    }
                    origin2.pop();
                }

                std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, typename HierarchicalNSW<dist_t>::CompareByFirst> top_candidates;

                if (cur_c < index2_offset) {
                    top_candidates = top_candidates2;
                    ll_cur = index1->get_linklist_at_level(cur_c - index1_offset, level);
                    linklistCount = index1->getListCount(ll_cur);
                    data = (tableint *)(ll_cur + 1);
                    dist = index1->get_dist_at_level(cur_c - index1_offset, level);
                    for (size_t iter = 0; iter < linklistCount; iter++) {
                        top_candidates.emplace(dist[iter], data[iter] + index1_offset);
                    }
                } else {
                    top_candidates = top_candidates1;
                    ll_cur = index2->get_linklist_at_level(cur_c - index2_offset, level);
                    linklistCount = index2->getListCount(ll_cur);
                    data = (tableint *)(ll_cur + 1);
                    dist = index2->get_dist_at_level(cur_c - index2_offset, level);
                    for (size_t iter = 0; iter < linklistCount; iter++) {
                        top_candidates.emplace(dist[iter], data[iter] + index2_offset);
                    }
                }
                alg_hnsw->getNeighborsByHeuristic2(top_candidates, Mcurmax, false);

                ll_cur = alg_hnsw->get_linklist_at_level(cur_c, level);
                alg_hnsw->setListCount(ll_cur, top_candidates.size());
                data = (tableint *)(ll_cur + 1);
                dist = (dist_t *)alg_hnsw->get_dist_at_level(cur_c, level);
                for (size_t idx = 0; top_candidates.size() > 0; idx++) {
                    data[idx] = top_candidates.top().second;
                    dist[idx] = top_candidates.top().first;
                    top_candidates.pop();
                }

                if (new_c == -1) {
                    break;
                }
                cur_c = new_c;
            }
        }
    }
    return alg_hnsw;
}

} // namespace hnswlib