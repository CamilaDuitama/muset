#pragma once

#include <algorithm>

#define WITH_KM_IO
#include <kmtricks/public.hpp>

#include <kmat_tools/cli/filter.h>

namespace kmat {

template<size_t MAX_K>
class FilterTask : public km::ITask
{
    using count_type = typename km::selectC<DMAX_C>::type;

public:
    
    FilterTask(std::string &input, std::string &output, std::size_t &nb_kmers, std::size_t &nb_retained, filter_opt_t opts, bool compress = true)
        : km::ITask(4, false), m_input(input), m_output(output), m_nb_kmers(nb_kmers), m_nb_retained(nb_retained), m_opts(opts), m_compress(compress)
    {}

    void preprocess() {}
    void postprocess() {}

    void exec()
    {
        km::MatrixReader reader(m_input);
        km::Kmer<MAX_K> kmer; kmer.set_k(m_opts->kmer_size);
        
        std::size_t nb_samples{reader.infos().nb_counts};
        std::vector<count_type> counts(nb_samples);
        
        km::MatrixWriter<8192> writer(m_output,
            m_opts->kmer_size,
            reader.infos().count_slots,
            nb_samples,
            reader.infos().id,
            reader.infos().partition,
            m_compress);

        while (reader.template read<MAX_K, DMAX_C>(kmer, counts)) {
            m_nb_kmers++;
            std::size_t nb_absent{0};
            std::size_t nb_present{0};
            for (auto c : counts) {
                if(c >= m_opts->min_abundance){ 
                    nb_present++;
                } else { 
                    nb_absent++; 
                }
            }

            bool enough_absent = (!m_opts->min_nb_absent_set && nb_absent >= m_opts->min_frac_absent * nb_samples)
                || (m_opts->min_nb_absent_set && nb_absent >= m_opts->min_nb_absent);    
            bool enough_present = (!m_opts->min_nb_present_set && nb_present >= m_opts->min_frac_present * nb_samples)
                || (m_opts->min_nb_present_set && nb_present >= m_opts->min_nb_present);
            
            if(enough_absent && enough_present) {
                m_nb_retained++;
                writer.template write<MAX_K, DMAX_C>(kmer,counts);
            }
        }
    }

private:

    std::string& m_input;
    std::string& m_output;
    std::size_t& m_nb_kmers;
    std::size_t& m_nb_retained;
    filter_opt_t m_opts;
    bool m_compress;
};

};