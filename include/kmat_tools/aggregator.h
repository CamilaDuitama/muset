#ifndef AGGREGATOR_H
#define AGGREGATOR_H

#include <vector>
#include <stdint.h>
#include <iostream>

#include <fmt/format.h>

#include <kmat_tools/utils.h>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

// This interface is used to abstract the possible aggregating statistics for per-sample k-mer count in unitigs
class Aggregator {
public:
    virtual ~Aggregator() = default;
    // Process the counts for a single k-mer belonging to a unitig
    virtual void process_kmer(size_t unitig_id, const std::vector<uint32_t>& kmer_counts) = 0;
    // Calculate the final abundance for a given unitig and sample
    virtual std::pair<double, double> get_abundance_fraction(size_t unitig_id, size_t sample_id, size_t unitig_num_kmers) const = 0;
};

// Reimplementation of Riccardo's mean computation logic
class MeanAggregator: public Aggregator {
  private:
  std::vector<std::vector<std::pair<uint32_t,uint32_t>>> m_samples_count;
  size_t m_num_samples;
  double m_min_fraction;
  // the pair contains <nb_present, abundance_sum>

  public:

  MeanAggregator(size_t num_samples, size_t num_utgs, double min_fraction):
  m_num_samples(num_samples), m_min_fraction(min_fraction) {
    m_samples_count.resize(num_utgs);
    for(size_t i {0}; i < num_utgs; i++){
      m_samples_count[i].resize(m_num_samples, {0,0});
    }
  }

  ~MeanAggregator() override = default;

  void process_kmer(size_t unitig_id, const std::vector<uint32_t>& kmer_counts) override {
    if (unitig_id >= m_samples_count.size()) {
      spdlog::debug(fmt::format("ERROR: GOT UTG ID {} BUT MAX POSSIBLE IS {}", unitig_id, m_samples_count.size()));
      return;
    }
    if (kmer_counts.size() != m_num_samples) {
      spdlog::debug(fmt::format("ERROR: GOT THESE NUMBER OF SAMPLES {}. EXPECTED {}", kmer_counts.size(), m_num_samples));
        return;
    }

    auto& samples = m_samples_count[unitig_id];
    size_t idx{0};
    for(auto& [nb_present,abundance_sum]: samples) {
        uint32_t num = kmer_counts[idx];
        // std::cout << "[INSIDE-" << unitig_id << "]: I see " << num << " and report " << uint32_t{num > 0} << ". old sum: " << nb_present << "\t";
        nb_present = kmat::add_sat(nb_present, uint32_t{num > 0});
        // std::cout << "New sum:" << nb_present << std::endl;
        abundance_sum = kmat::add_sat(abundance_sum, num);
        idx++;
    }
  }

  std::pair<double,double> get_abundance_fraction(size_t unitig_id, size_t sample_id, size_t unitig_num_kmers) const override{
    if (unitig_id >= m_samples_count.size() || sample_id >= m_num_samples || unitig_num_kmers == 0) {
            return std::pair(0.0,0.0);
        }

    const auto& [nb_present, abundance_sum] = m_samples_count[unitig_id][sample_id];
    // std::cout << "[INSIDE] GOT " << nb_present << " present out of " << unitig_num_kmers << std::endl;
    // return std::pair(static_cast<double>(abundance_sum) / unitig_num_kmers,static_cast<double>(nb_present)/ unitig_num_kmers);
    double fraction = static_cast<double>(nb_present) / unitig_num_kmers;
    double abundance = fraction >= m_min_fraction ? static_cast<double>(abundance_sum) / unitig_num_kmers : 0.0;
    return std::pair(abundance, fraction);
  }

};


// Median computation logic
class MedianAggregator: public Aggregator {
  private:
  std::vector<std::vector<std::map<uint32_t, uint32_t>>> m_samples_count;
  size_t m_num_samples;
  double m_min_fraction;
  // trying to reduce memory footprint, instead of having a vector of counts, I store the possible counts in a dictionary
  // I hope that most of the counts per sample will be the same so that I can save space.

  public:

  MedianAggregator(size_t num_samples, size_t num_utgs, double min_fraction):
  m_num_samples(num_samples), m_min_fraction(min_fraction) {
    m_samples_count.resize(num_utgs);
    for(int i {0}; i < num_utgs; i++){
      m_samples_count[i].resize(m_num_samples);
    }
  }

  ~MedianAggregator() override = default;

  void process_kmer(size_t unitig_id, const std::vector<uint32_t>& kmer_counts) override {
    if (unitig_id >= m_samples_count.size()) {
      spdlog::debug(fmt::format("ERROR: GOT UTG ID {} BUT MAX POSSIBLE IS {}", unitig_id, m_samples_count.size()));
      return;
    }
    if (kmer_counts.size() != m_num_samples) {
      spdlog::debug(fmt::format("ERROR: GOT THESE NUMBER OF SAMPLES {}. EXPECTED {}", kmer_counts.size(), m_num_samples));
        return;
    }

    auto& samples = m_samples_count[unitig_id];
    for(size_t idx {0}; idx < m_samples_count[unitig_id].size(); idx++) {
        samples[idx][kmer_counts[idx]]++;
    }
  }

  std::pair<double, double> get_abundance_fraction(size_t unitig_id, size_t sample_id, size_t unitig_num_kmers) const override{

    const auto& sample_map = m_samples_count[unitig_id][sample_id];

    // 1 - computing fraction by counting number of occurrences.
    size_t nb_present {0};
    for (const auto& [key, value] : sample_map) {
        if (key > 0) {
            nb_present += value;
        }
    }

    // 2 - computing abundance by finding the median value from the dictionoary of occurencies.
    // Since I am using a std::map, keys are already sorted in ascending order. I just need to traverse the map
    // and stop at the middle value(s). In case odd, take the value, in case even, take
    double abundance {0.0};

    if (nb_present != 0)
    {
        // 0-based indices
        const std::size_t idx1 = (nb_present - 1) / 2;
        const std::size_t idx2 =  nb_present      / 2;

        std::size_t running = 0;
        int         m1 {};              // value at idx1
        int         m2 {};              // value at idx2
        bool        have_m1 = false;

        // std::cout << "IDX1: " << idx1 << "; IDX2: " << idx2 << std::endl;
        // for (auto const& [key, cnt] : sample_map)
        // {
        //   if (key <= 0)
        //     continue;
        //   std::cout << "key: " << key << ", running: " << running << ", count: " << cnt << std::endl;
        //   running += cnt;
        //   if (!have_m1 && running > idx1) { m1 = key; have_m1 = true; std::cout << "m1 is: " << key << ", ";}
        //   if (running > idx2) { m2 = key; std::cout << "m2 is: " << key << std::endl; break; }
        // }
        abundance = (idx1 == idx2) ? static_cast<double>(m2) : (static_cast<double>(m1) + m2) / 2.0;
    }

    // returning
    const double fraction = static_cast<double>(nb_present) / static_cast<double>(unitig_num_kmers);
    if (fraction < m_min_fraction) abundance = 0.0;
    return {abundance, fraction};
  }

};

#endif