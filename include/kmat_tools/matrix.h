#pragma once

#include <array>
#include <cstddef>
#include <exception>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <string_view>

#include <fmt/format.h>

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

    // read kmer and discard the rest of the line
    inline bool read_kmer(std::string &kmer) {

      if (!this->getline()) { return false; }
      kmer = std::string_view{ m_line, m_line.find_first_of(" \t") };
      return true;
    }

    // read kmer and the remaing part of the line
    inline bool read_kmer_and_line(std::string &kmer, std::string &line) {
  
      if (!this->getline()) { return false; }
      
      auto data = m_line.data();
      
      auto idx = m_line.find_first_of(" \t");
      kmer = std::string_view{data,idx};

      idx = m_line.find_first_not_of(" \t", idx+1);
      line = std::string_view{data+idx, m_line.size()-idx};

      return true;
    }

  private:

    inline bool getline() {
      
      if (!m_stream->good()) {
        return false;
      }

      std::getline(*m_stream, m_line);

      if (m_stream->eof() && m_line.empty()) {
        return false;
      } else if (m_stream->bad() || m_stream->fail()) {
        throw std::runtime_error(fmt::format("bad unexpected error reading from {}", m_path));
      }

      m_line_count++;
      return true;
    }

    std::string m_path;
    stream_t    m_stream;
    buffer_t    m_buf;

    std::string m_line;
    size_t      m_line_count{0};
};

};