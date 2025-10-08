#include <kmat_tools/aggregator.h>
#include <gtest/gtest.h>
#include <random>
#include <vector>
#include <algorithm>

// Random number generator
std::mt19937& get_rng() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return gen;
}

// Function to generate random count values between 0 and 100 per k-mer
std::vector<uint32_t> random_counts(size_t num_samples, uint32_t max_val = 100) {
    std::vector<uint32_t> counts(num_samples);
    std::uniform_int_distribution<uint32_t> dist(1, max_val);
    for (auto& c : counts) {
        c = dist(get_rng());
    }
    return counts;
}

// Helper function to compute the expected fraction
double expected_fract(const std::vector<std::vector<uint32_t>>& kmers, size_t sample_idx) {
    if (kmers.empty()) return 0.0;

    double sum = 0.0;
    for (const auto& kmer : kmers) {
        if (kmer[sample_idx] != 0) sum += 1;
    }
    return sum / kmers.size();
}

// Helper function to compute the expected mean
double expected_mean(const std::vector<std::vector<uint32_t>>& kmers, size_t sample_idx) {
    if (kmers.empty()) return 0.0;

    double sum = 0.0;
    for (const auto& kmer : kmers) {
        sum += kmer[sample_idx];
    }
    return sum / kmers.size();
}

// Helper function to compute the expected median
double expected_median(const std::vector<uint32_t>& values) {
    if (values.empty()) return 0.0;

    std::vector<uint32_t> sorted = values;
    std::sort(sorted.begin(), sorted.end());

    size_t n = sorted.size();
    if (n % 2 == 0) {
        return (sorted[n/2 - 1] + sorted[n/2]) / 2.0;
    } else {
        return sorted[n/2];
    }
}


// Test fixture for randomized tests
class AggregatorRandomTest : public ::testing::Test {
protected:
    void SetUp() override {
        get_rng().seed(std::random_device()());
    }
};

// MeanAggregator without filtering
TEST_F(AggregatorRandomTest, MeanAggregator_1000RandomTests) {
    const size_t num_tests = 1000;
    const size_t max_samples = 5;
    const size_t max_kmers = 20;
    const double min_fraction = 0.5;

    for (size_t test_num = 0; test_num < num_tests; ++test_num) {
        // Random configuration
        size_t num_samples = 1 + get_rng()() % max_samples;
        size_t num_kmers = 1 + get_rng()() % max_kmers;

        MeanAggregator agg(num_samples, 1, min_fraction);

        // Filling k-mers with randomized count to the unitig
        std::vector<std::vector<uint32_t>> kmers;
        for (size_t i = 0; i < num_kmers; ++i) {
            auto counts = random_counts(num_samples);
            kmers.push_back(counts);
            agg.process_kmer(0, counts);
        }

        // Verification for each sample
        for (size_t sample = 0; sample < num_samples; ++sample) {
            auto [computed_abundance, computed_fraction] = agg.get_abundance_fraction(0, sample, num_kmers);

            // Verify fraction
            double expected_fraction = expected_fract(kmers, sample);
            EXPECT_DOUBLE_EQ(computed_fraction, expected_fraction);

            // Verify mean
            double expected_abundance = expected_mean(kmers, sample);
            EXPECT_NEAR(computed_abundance, expected_abundance, 0.001)
                << "Test " << test_num << ", Sample " << sample
                << ": Expected [Fraction] " << expected_fraction << ", got " << computed_fraction << "\t"
                << ": Expected [Abundance] " << expected_abundance << ", got " << computed_abundance;
        }
    }
}

// MeanAggregator with min fraction filtering
TEST_F(AggregatorRandomTest, MeanAggregator_MinFraction_1000RandomTests) {
    const size_t num_tests = 1000;
    const size_t max_samples = 3;
    const size_t max_kmers = 15;

    for (size_t test_num = 0; test_num < num_tests; ++test_num) {
      // Random configuration with random min_fraction
      size_t num_samples = 1 + get_rng()() % max_samples;
      size_t num_kmers = 5 + get_rng()() % max_kmers;
      double min_fraction = 0.1 + (get_rng()() % 90) / 100.0; // 0.1 to 0.99

      MeanAggregator agg(num_samples, 1, min_fraction);

      // Generate random kmers
      std::vector<std::vector<uint32_t>> all_kmers;
      for (size_t i = 0; i < num_kmers; ++i) {
          all_kmers.push_back(random_counts(num_samples));
      }

        // Randomly process some kmers (simulating missing data)
        std::bernoulli_distribution process_dist(0.7); // 70% chance to process each kmer
        std::vector<bool> processed(num_kmers, false);

        for (size_t i = 0; i < num_kmers; ++i) {
            if (process_dist(get_rng())) {
                agg.process_kmer(0, all_kmers[i]);
                processed[i] = true;
            }
            else {
              agg.process_kmer(0,std::vector<uint32_t>(num_samples,0));
            }
        }

      // Verify results for each sample
      for (size_t sample = 0; sample < num_samples; ++sample) {
          auto [computed_abundance, computed_fraction] = agg.get_abundance_fraction(0, sample, num_kmers);
          double count {0.0};
          for (int i {0}; i < num_kmers; i++){
            if (processed[i]){
              if (all_kmers[i][sample] > 0) count+=1;
            }
          }
          double expected_fraction {count / num_kmers};
          // Verify fraction
          EXPECT_NEAR(computed_fraction, expected_fraction, 0.001);
          // Verify filtering
          if (expected_fraction < min_fraction) {
              EXPECT_DOUBLE_EQ(computed_abundance, 0.0);
          } else {
              // Calculate expected mean from processed kmers
              double sum = 0.0;
              for (size_t i = 0; i < num_kmers; ++i) {
                  if (processed[i]) {
                      sum += all_kmers[i][sample];
                  }
              }
              double expected_abundance = expected_fraction >= min_fraction ? sum / num_kmers : 0.0;
              EXPECT_NEAR(computed_abundance, expected_abundance, 0.001)
                  << "Test " << test_num << ", Sample " << sample
              << ": Expected [Fraction] " << expected_fraction << ", got " << computed_fraction << "\t"
              << ": Expected [Abundance] " << expected_abundance << ", got " << computed_abundance;
            }
        }
    }
}


// Run 1000 randomized tests for MedianAggregator
TEST_F(AggregatorRandomTest, MedianAggregator_1000RandomTests) {
    const size_t num_tests = 1000;
    const size_t max_samples = 5;
    const size_t max_kmers = 25;
    const double min_fraction = 0.5;

    for (size_t test_num = 0; test_num < num_tests; ++test_num) {
        // Random configuration
        size_t num_samples = 1 + get_rng()() % max_samples;
        size_t num_kmers = 1 + get_rng()() % max_kmers;

        MedianAggregator agg(num_samples, 1, min_fraction);

        // Generate and process random kmers
        std::vector<std::vector<uint32_t>> kmers;
        for (size_t i = 0; i < num_kmers; ++i) {
            auto counts = random_counts(num_samples, 20);
            kmers.push_back(counts);
            agg.process_kmer(0, counts);
        }

        // Verify results for each sample
        for (size_t sample = 0; sample < num_samples; ++sample) {
            auto [computed_abundance, computed_fraction] = agg.get_abundance_fraction(0, sample, num_kmers);
            // Verify fraction
            EXPECT_DOUBLE_EQ(computed_fraction, 1.0);

            // Collect all counts for this sample
            std::vector<uint32_t> sample_counts;
            for (const auto& kmer : kmers) {
                sample_counts.push_back(kmer[sample]);
            }

            // Verify median
            double expected_abundance = expected_median(sample_counts);
            EXPECT_NEAR(computed_abundance, expected_abundance, 0.001)
                << "Test " << test_num << ", Sample " << sample
                << ": Expected " << expected_abundance << ", got " << computed_abundance;
        }
    }
}

// MeanAggregator with min fraction filtering
TEST_F(AggregatorRandomTest, MedianAggregator_MinFraction_1000RandomTests) {
    const size_t num_tests = 1000;
    const size_t max_samples = 3;
    const size_t max_kmers = 15;

    for (size_t test_num = 0; test_num < num_tests; ++test_num) {
      // Random configuration with random min_fraction
      size_t num_samples = 1 + get_rng()() % max_samples;
      size_t num_kmers = 5 + get_rng()() % max_kmers;
      double min_fraction = 0.1 + (get_rng()() % 90) / 100.0; // 0.1 to 0.99

      MedianAggregator agg(num_samples, 1, min_fraction);

      // Generate random kmers
      std::vector<std::vector<uint32_t>> all_kmers;
      for (size_t i = 0; i < num_kmers; ++i) {
          all_kmers.push_back(random_counts(num_samples));
      }

        // Randomly process some kmers (simulating missing data)
        std::bernoulli_distribution process_dist(0.7); // 70% chance to process each kmer
        std::vector<bool> processed(num_kmers, false);

        for (size_t i = 0; i < num_kmers; ++i) {
            if (process_dist(get_rng())) {
                agg.process_kmer(0, all_kmers[i]);
                processed[i] = true;
            }
            else {
              agg.process_kmer(0,std::vector<uint32_t>(num_samples,0));
            }
        }

      // Verify results for each sample
      for (size_t sample = 0; sample < num_samples; ++sample) {
          auto [computed_abundance, computed_fraction] = agg.get_abundance_fraction(0, sample, num_kmers);
          double count {0.0};
          for (int i {0}; i < num_kmers; i++){
            if (processed[i]){
              if (all_kmers[i][sample] > 0) count+=1;
            }
          }
          double expected_fraction {count / num_kmers};
          // Verify fraction
          EXPECT_NEAR(computed_fraction, expected_fraction, 0.001);
          // Verify filtering
          if (expected_fraction < min_fraction) {
              EXPECT_DOUBLE_EQ(computed_abundance, 0.0);
          } else {
              // Calculate expected median from processed kmers
              std::vector<uint32_t> sample_counts;
              for (size_t i = 0; i < num_kmers; ++i) {
                  if (processed[i]) {
                      sample_counts.push_back(all_kmers[i][sample]);
                  }
              }
              double expected_abundance = expected_median(sample_counts);
              EXPECT_NEAR(computed_abundance, expected_abundance, 0.001)
                  << "Test " << test_num << ", Sample " << sample
              << ": Expected [Fraction] " << expected_fraction << ", got " << computed_fraction << "\t"
              << ": Expected [Abundance] " << expected_abundance << ", got " << computed_abundance;
            }
        }
    }
}
