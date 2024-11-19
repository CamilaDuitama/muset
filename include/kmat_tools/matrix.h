#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <exception>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include <fmt/format.h>

#include <kmat_tools/utils.h>

namespace kmat {

template<size_t buf_size = 16384>
class TextMatrixReader {

    using stream_t = std::unique_ptr<std::istream>;
    using buffer_t = std::array<char, buf_size>;

  public:

    TextMatrixReader () : m_stream(new std::fstream{}) {}

    TextMatrixReader (const std::string& path)
      : m_path(path), m_stream(new std::fstream{path, std::ios::in})
    {
      if (!m_stream->good()) {
        throw std::runtime_error(fmt::format("cannot open {}", path));
      }

      m_stream->rdbuf()->pubsetbuf(m_buf.data(), m_buf.size());
    }

    TextMatrixReader (TextMatrixReader const &) = delete;
    TextMatrixReader (TextMatrixReader&&) = delete;
    TextMatrixReader& operator= (TextMatrixReader const &) = delete;
    TextMatrixReader& operator= (TextMatrixReader&&) = delete;

    // return reference to last read line
    inline const std::string& line() const {
      return m_line;
    }

    inline size_t line_count() const {
      return m_line_count;
    }

    // read kmer and discard the rest of the line
    inline bool read_kmer(std::string &kmer) {

      if (!this->get_nonempty_line()) { return false; }

      auto idx = m_line.find_first_of(" \t");
      kmer = idx == std::string::npos ? 
        std::string_view{ m_line } : 
        std::string_view{ m_line.data(), idx };

      if(!is_valid_kmer(kmer)) {
        throw std::runtime_error(
          fmt::format("bad character found in k-mer \"{}\" at line {}", kmer, this->line_count())
        );
      }
      
      return kmer.length() > 0;
    }

    // read kmer and the remaing part of the line
    inline bool read_kmer_and_line(std::string &kmer, std::string &line) {
  
      if (!this->get_nonempty_line()) { return false; }

      line.clear();
      
      auto idx = m_line.find_first_of(" \t");
      if (idx == std::string::npos) {
        kmer = m_line;
        if(!is_valid_kmer(kmer)) {
          throw std::runtime_error(
            fmt::format("bad character found in k-mer \"{}\" at line {}", kmer, this->line_count())
          );
        }
        return kmer.length() > 0;
      }

      kmer = std::string_view{ m_line.data(), idx };
      if(!is_valid_kmer(kmer)) {
        throw std::runtime_error(
          fmt::format("bad character found in k-mer \"{}\" at line {}", kmer, this->line_count())
        );
      }

      idx = m_line.find_first_not_of(" \t", idx);
      if (idx != std::string::npos) {
        line = std::string_view{m_line.data()+idx, m_line.size()-idx};
      }

      return true;
    }

    inline bool read_kmer_counts(std::string &kmer, std::vector<size_t> &counts) {

      if (!this->get_nonempty_line()) { return false; }

      counts.clear();

      auto idx = m_line.find_first_of(" \t");
      kmer = idx == std::string::npos ? 
        std::string_view{ m_line } : 
        std::string_view{ m_line.data(), idx };

      if(!is_valid_kmer(kmer)) {
        throw std::runtime_error(
          fmt::format("bad character found in k-mer \"{}\" at line {}", kmer, this->line_count())
        );
      }

      while (idx != std::string::npos) {

        idx = m_line.find_first_not_of(" \t", idx); // skip whitespaces
        if (idx == std::string::npos) { break; }

        size_t value{0};
        auto [ptr, ec] = std::from_chars(m_line.data()+idx, m_line.data()+m_line.size(), value);
        if (ec != std::errc()) {
          throw std::runtime_error(fmt::format("{}: error loading counts at line {}", this->m_path, this->line_count()));
        }

        counts.push_back(value);
        idx = ptr - m_line.data();
      }

      return kmer.length() > 0;
    }

    

  private:

    inline bool get_nonempty_line() {

      do {
        if (!m_stream->good()) {
          return false;
        }

        std::getline(*m_stream, m_line);

        if (m_stream->eof()) {
          return false;
        } else if (m_stream->bad() || m_stream->fail()) {
          throw std::runtime_error(fmt::format("{}: unexpected read error", m_path));
        }
        
        m_line_count++;

      } while (m_line.find_first_not_of(" \t") == std::string::npos);

      return true;
    }

    std::string m_path;
    stream_t    m_stream;
    buffer_t    m_buf;

    std::string m_line;
    size_t      m_line_count{0};
};

};