#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <algorithm>

// Enum definitions
enum WorkloadType {
    SIFT1M,
    SIFT10M,
    DEEP10M,
    TURING10M,
    SIFT100M,
    DEEP100M,
    SIFT1B,
    GLOVE25,
    COHERE10M
};

enum MergeMethod {
    REBUILD,
    INSERT,
    TWO_MERGE,
    MULTI_TWO_MERGE,
    ES,
    NGM,
    IGTM,
    CGTM,
    ABLATION_C,
    BACKWARD_SEARCH,
    MEMORY_EFFICIENCY,
    MEMORY_EFFICIENCY_LIMIT,
};

enum MultiTestMethod {
    LARGE_FIRST,
    SMALL_FIRST,
    RANDOM
};

// Config struct
struct Config {
    WorkloadType workload_type = SIFT10M;
    MergeMethod merge_method = TWO_MERGE;
    MultiTestMethod multi_test_method = LARGE_FIRST;
    int dim=128;
    long max_elements=1e7;
    int nb=1e7;
    int M=32;
    int ef_construction=64;
    int k=100;
    int kk=100;
    int nq=10000;
    int iterations=1;
    int lrange=0;
    int rrange=0;
    bool rerun = false;
    int thread = 1;
    int lambda = 4;
    float alpha = 1.05;
    bool save_index = true;
    std::string base_filepath;
    std::string query_filepath;
    std::string groundtruth_filepath;
    std::vector<std::string> index_path;
    std::string save_path = "./";
    std::vector<int> efs_array;
};

// Parse string to WorkloadType enum
WorkloadType parseWorkloadType(const std::string &s) {
    if (s == "SIFT1M") return SIFT1M;
    if (s == "SIFT10M") return SIFT10M;
    if (s == "DEEP10M") return DEEP10M;
    if (s == "TURING10M") return TURING10M;
    if (s == "SIFT100M") return SIFT100M;
    if (s == "DEEP100M") return DEEP100M;
    if (s == "SIFT1B") return SIFT1B;
    if (s == "GLOVE25") return GLOVE25;
    if (s == "COHERE10M") return COHERE10M;
    std::cerr << "Unknown WorkloadType: " << s << std::endl;
    std::exit(1);
}

std::string workloadTypeToString(WorkloadType type) {
    switch (type) {
    case SIFT1M: return "SIFT1M";
    case SIFT10M: return "SIFT10M";
    case DEEP10M: return "DEEP10M";
    case TURING10M: return "TURING10M";
    case SIFT100M: return "SIFT100M";
    case DEEP100M: return "DEEP100M";
    case SIFT1B: return "SIFT1B";
    case GLOVE25: return "GLOVE25";
    case COHERE10M: return "COHERE10M";
    default: return "UNKNOWN";
    }
}

// Parse string to MergeMethod enum
MergeMethod parseMergeMethod(const std::string &s) {
    if (s == "REBUILD") return REBUILD;
    if (s == "INSERT") return INSERT;
    if (s == "TWO_MERGE") return TWO_MERGE;
    if (s == "MULTI_TWO_MERGE") return MULTI_TWO_MERGE;
    if (s == "ES") return ES;
    if (s == "NGM") return NGM;
    if (s == "IGTM") return IGTM;
    if (s == "CGTM") return CGTM;
    if (s == "ABLATION_C") return ABLATION_C;
    if (s == "BACKWARD_SEARCH") return BACKWARD_SEARCH;
    if (s == "MEMORY_EFFICIENCY") return MEMORY_EFFICIENCY;
    if (s == "MEMORY_EFFICIENCY_LIMIT") return MEMORY_EFFICIENCY_LIMIT;
    std::cerr << "Unknown MergeMethod: " << s << std::endl;
    std::exit(1);
}

std::string mergeMethodToString(MergeMethod method) {
    switch (method) {
    case REBUILD: return "REBUILD";
    case INSERT: return "INSERT";
    case TWO_MERGE: return "TWO_MERGE";
    case MULTI_TWO_MERGE: return "MULTI_TWO_MERGE";
    case ES: return "ES";
    case NGM: return "NGM";
    case IGTM: return "IGTM";
    case CGTM: return "CGTM";
    case ABLATION_C: return "ABLATION_C";
    case BACKWARD_SEARCH: return "BACKWARD_SEARCH";
    case MEMORY_EFFICIENCY: return "MEMORY_EFFICIENCY";
    case MEMORY_EFFICIENCY_LIMIT: return "MEMORY_EFFICIENCY_LIMIT";
    default: return "UNKNOWN";
    }
}

// Parse string to MultiTestMethod enum
MultiTestMethod parseMultiTestMethod(const std::string &s) {
    if (s == "LARGE_FIRST") return LARGE_FIRST;
    if (s == "SMALL_FIRST") return SMALL_FIRST;
    if (s == "RANDOM") return RANDOM;
    std::cerr << "Unknown MultiTestMethod: " << s << std::endl;
    std::exit(1);
}

std::string multiTestMethodToString(MultiTestMethod method) {
    switch (method) {
    case LARGE_FIRST: return "LARGE_FIRST";
    case SMALL_FIRST: return "SMALL_FIRST";
    case RANDOM: return "RANDOM";
    default: return "UNKNOWN";
    }
}

// Parse comma-separated integers into vector<int>
std::vector<int> parseIntList(const std::string &s) {
    std::vector<int> result;
    std::istringstream ss(s);
    std::string token;
    while (std::getline(ss, token, ',')) {
        size_t start = token.find_first_not_of(" \t\r\n");
        size_t end = token.find_last_not_of(" \t\r\n");
        if (start == std::string::npos) continue;
        std::string num_str = token.substr(start, end - start + 1);
        result.push_back(std::stoi(num_str));
    }
    return result;
}

// Parse comma-separated strings into vector<string>
std::vector<std::string> parseStringList(const std::string &s) {
    std::vector<std::string> result;
    std::istringstream ss(s);
    std::string token;
    while (std::getline(ss, token, ',')) {
        size_t start = token.find_first_not_of(" \t\r\n");
        size_t end = token.find_last_not_of(" \t\r\n");
        if (start == std::string::npos) continue;
        std::string path_str = token.substr(start, end - start + 1);
        result.push_back(path_str);
    }
    return result;
}

// Set default values for dim, max_elements, nb, k, kk, nq based on workload_type
void setDefaultsByWorkload(Config &cfg) {
    switch (cfg.workload_type) {
    case SIFT1M:
        cfg.dim = 128;
        cfg.max_elements = 1e6;
        cfg.nb = 1e6;
        cfg.k = 100;
        cfg.kk = 1000;
        cfg.nq = 10000;
        break;
    case SIFT10M:
        cfg.dim = 128;
        cfg.max_elements = 10e6;
        cfg.nb = 10e6;
        cfg.k = 100;
        cfg.kk = 1000;
        cfg.nq = 10000;
        break;
    case DEEP10M:
        cfg.dim = 96;
        cfg.max_elements = 10e6;
        cfg.nb = 10e6;
        cfg.k = 100;
        cfg.kk = 100;
        cfg.nq = 10000;
        break;
    case TURING10M:
        cfg.dim = 100;
        cfg.max_elements = 10e6;
        cfg.nb = 10e6;
        cfg.k = 100;
        cfg.kk = 100;
        cfg.nq = 10000;
        break;
    case SIFT100M:
        cfg.dim = 128;
        cfg.max_elements = 100e6;
        cfg.nb = 100e6;
        cfg.k = 1000;
        cfg.kk = 100;
        cfg.nq = 10000;
        break;
    case DEEP100M:
        cfg.dim = 96;
        cfg.max_elements = 100e6;
        cfg.nb = 100e6;
        cfg.k = 100;
        cfg.kk = 100;
        cfg.nq = 10000;
        break;
    default:
        std::cerr << "Unhandled WorkloadType for defaults" << std::endl;
        std::exit(1);
    }
}

// Main function: load configuration from file into Config object
Config loadConfig(const std::string &filename) {
    Config cfg;
    setDefaultsByWorkload(cfg);
    std::map<std::string, std::string> kv;
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Failed to open config file: " << filename << std::endl;
        std::exit(1);
    }

    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#') continue;
        auto pos_eq = line.find('=');
        if (pos_eq == std::string::npos) continue;

        std::string key = line.substr(0, pos_eq);
        std::string value = line.substr(pos_eq + 1);

        auto trim = [](std::string &s) {
            size_t start = s.find_first_not_of(" \t\r\n");
            size_t end = s.find_last_not_of(" \t\r\n");
            if (start == std::string::npos) {
                s.clear();
            } else {
                s = s.substr(start, end - start + 1);
            }
        };
        trim(key);
        trim(value);

        if (!key.empty()) {
            kv[key] = value;
        }
    }
    infile.close();

    // Parse workload_type first and set default values
    cfg.workload_type = parseWorkloadType(kv.at("workload_type"));

    // Parse merge_method
    cfg.merge_method = parseMergeMethod(kv.at("merge_method"));

    // If config file provides explicit overrides for dim, max_elements, nb, k, kk, nq, apply them
    if (kv.count("dim")) cfg.dim = std::stoi(kv.at("dim"));
    if (kv.count("max_elements")) cfg.max_elements = std::stol(kv.at("max_elements"));
    if (kv.count("nb")) cfg.nb = std::stoi(kv.at("nb"));
    if (kv.count("k")) cfg.k = std::stoi(kv.at("k"));
    if (kv.count("kk")) cfg.kk = std::stoi(kv.at("kk"));
    if (kv.count("nq")) cfg.nq = std::stoi(kv.at("nq"));
    if (kv.count("rerun")) cfg.rerun = (kv.at("rerun") == "true");
    if (kv.count("thread")) cfg.thread = std::stoi(kv.at("thread"));
    if (kv.count("lambda")) cfg.lambda = std::stoi(kv.at("lambda"));
    if (kv.count("alpha")) cfg.alpha = std::stof(kv.at("alpha"));
    if (kv.count("save_index")) cfg.save_index = (kv.at("save_index") == "true");

    if (kv.count("multi_test_method")) cfg.multi_test_method = parseMultiTestMethod(kv.at("multi_test_method"));

    // Parse remaining numeric fields
    if (kv.count("M")) cfg.M = std::stoi(kv.at("M"));
    if (kv.count("ef_construction")) cfg.ef_construction = std::stoi(kv.at("ef_construction"));
    if (kv.count("iterations")) cfg.iterations = std::stoi(kv.at("iterations"));
    if (kv.count("lrange")) cfg.lrange = std::stoi(kv.at("lrange"));
    if (kv.count("rrange")) cfg.rrange = std::stoi(kv.at("rrange"));

    // Parse file paths
    if (kv.count("base_filepath")) cfg.base_filepath = kv.at("base_filepath");
    if (kv.count("query_filepath")) cfg.query_filepath = kv.at("query_filepath");
    if (kv.count("groundtruth_filepath")) cfg.groundtruth_filepath = kv.at("groundtruth_filepath");
    if (kv.count("save_path")) cfg.save_path = kv.at("save_path");

    // Parse list fields
    cfg.index_path = parseStringList(kv.at("index_path"));
    cfg.efs_array = parseIntList(kv.at("efs_array"));

    return cfg;
}
