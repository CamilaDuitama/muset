#ifndef MATRIX_WRITER_H
#define MATRIX_WRITER_H

#include <vector>
#include <string>
#include <fstream>

#include <zlib.h>

#include <fmt/format.h>
#include <kmat_tools/utils.h>

namespace fs = std::filesystem;

// Interface to write the matrix into a file
class MatrixWriter {
  public:
    virtual ~MatrixWriter() = default;

    virtual void write_row(
        const std::string& unitig_identifier,
        const std::vector<double>& abundances,
        const std::vector<double>& fractions) = 0;
};


// Text matrix output
class TextMatrixWriter: public MatrixWriter {
  public:
    TextMatrixWriter(const std::string& prefix, bool write_fraction_matrix)
    : m_write_fraction(write_fraction_matrix) {

      // ABUNDANCE MATRIX STREAM
      std::string abundance_file = fmt::format("{}.abundance.mat", prefix);
      m_abundance_stream.open(abundance_file);
      if (!m_abundance_stream.good()) { throw std::runtime_error(fmt::format("cannot open output file \"{}\" for writing", abundance_file)); }
      m_abundance_stream << std::fixed << std::setprecision(2);

      // FRACTION MATRIX STREAM
      if (m_write_fraction) {
        std::string frac_file = fmt::format("{}.frac.mat", prefix);
        m_fraction_stream.open(frac_file);
        if (!m_fraction_stream.good()) { throw std::runtime_error(fmt::format("cannot open output file \"{}\" for writing", frac_file)); }
        m_fraction_stream << std::fixed << std::setprecision(2);
      }
    }

    ~TextMatrixWriter() override{
      if (m_abundance_stream.is_open()) {
        m_abundance_stream.close();
      }
      if (m_fraction_stream.is_open()) {
        m_fraction_stream.close();
      }
    }



    void write_row(const std::string& unitig_identifier, const std::vector<double>& abundances, const std::vector<double>& fractions) override{
      // ABUNDANCE MATRIX STREAM
      m_abundance_stream << unitig_identifier;
      for (const auto& val : abundances) {
        m_abundance_stream << " " << val;
      }
      m_abundance_stream << "\n";

      // FRACTION MATRIX STREAM
      if (m_write_fraction) {
        m_fraction_stream << unitig_identifier;
        for (const auto& val : fractions) {
          m_fraction_stream << " " << val;
        }
        m_fraction_stream << "\n";
      }
    }


  private:
    std::ofstream m_abundance_stream;
    std::ofstream m_fraction_stream;
    bool m_write_fraction;

};

class CompressedTSVMatrixWriter : public MatrixWriter {
  public:
    CompressedTSVMatrixWriter(const std::string& prefix, bool write_fraction_matrix, size_t number_samples):
    m_write_fraction(write_fraction_matrix), m_number_samples(number_samples) {
      // ABUNDANCE FILE
      std::string abundance_file = fmt::format("{}.abundance.tsv.gz", prefix);

      m_abundance_stream = gzopen(abundance_file.c_str(), "wb");
      if (!m_abundance_stream) {
        throw std::runtime_error(fmt::format("Cannot open abundance file: {}", abundance_file));
      }

      gzputs(m_abundance_stream, "unitig_id\t");
      for (size_t i {0}; i < m_number_samples; ++i) {
          gzprintf(m_abundance_stream, "sample_%zu\t", i);
      }
      gzputs(m_abundance_stream, "\n");

      // FRACTION FILE
      if (m_write_fraction) {
        std::string fraction_file = fmt::format("{}.frac.tsv.gz", prefix);
        m_fraction_stream = gzopen(fraction_file.c_str(), "wb");
        if (!m_fraction_stream) {
            throw std::runtime_error(fmt::format("Cannot open fraction file: {}", fraction_file));
        }

        // Write header for fraction file
        gzputs(m_fraction_stream, "unitig_id\t");
        for (size_t i = 0; i < m_number_samples; ++i) {  // Same number as above
            gzprintf(m_fraction_stream, "sample_%zu\t", i);
        }
        gzputs(m_fraction_stream, "\n");
      }

    }

    ~CompressedTSVMatrixWriter() override {
      if (m_abundance_stream) {
          gzclose(m_abundance_stream);
      }
      if (m_fraction_stream) {
          gzclose(m_fraction_stream);
      }
    }

    void write_row(const std::string& unitig_identifier, const std::vector<double>& abundances, const std::vector<double>& fractions) override{
      // ABUNDANCE MATRIX STREAM
      gzprintf(m_abundance_stream, "%s", unitig_identifier.c_str());
      for (const auto& val : abundances) {
        gzprintf(m_abundance_stream, "\t%.2f", val);
      }
      gzputs(m_abundance_stream, "\n");

      // FRACTION MATRIX STREAM
      if (m_write_fraction) {
        gzprintf(m_fraction_stream, "%s", unitig_identifier.c_str());
        for (const auto& val : fractions) {
          gzprintf(m_fraction_stream, "\t%.2f", val);
        }
        gzputs(m_fraction_stream, "\n");
      }
    }

  private:
    gzFile m_abundance_stream = nullptr;
    gzFile m_fraction_stream = nullptr;
    bool m_write_fraction;
    size_t m_number_samples;
};

// class HDF5MatrixWriter : public MatrixWriter {
// public:
//     HDF5MatrixWriter(const std::string& prefix, bool write_fraction_matrix, size_t number_samples, size_t number_unitigs)
//         : m_write_fraction(write_fraction_matrix), m_num_cols(number_samples), m_current_row(0), m_num_rows(number_unitigs) {

//         std::string h5_file = fmt::format("{}.h5", prefix);
//         m_file = H5::H5File(h5_file, H5F_ACC_TRUNC);

//         // --- Create Dataspace and Properties for Extendable Datasets ---
//         hsize_t final_dims[2] = {m_num_rows, m_num_cols};
//         H5::DataSpace matrix_dataspace(2, final_dims);
//         H5::DataSpace id_dataspace(1,m_num_rows);

//         // --- Create the Datasets ---
//         // Variable-length string type for identifiers
//         H5::StrType id_type(H5::PredType::C_S1, H5T_VARIABLE);

//         // Create dataset for unitig IDs
//         m_id_dataset = m_file.createDataSet("unitig_ids", id_type, id_dataspace);

//         // Create dataset for abundances
//         m_abundance_dataset = m_file.createDataSet("abundances", H5::PredType::NATIVE_DOUBLE, matrix_dataspace);

//         if (m_write_fraction) {
//             m_fraction_dataset = m_file.createDataSet("fractions", H5::PredType::NATIVE_DOUBLE, matrix_dataspace);
//         }
//     }

//     ~HDF5MatrixWriter() override {
//         // H5:: objects manage their own resources (RAII)
//     }

//     void write_row(const std::string& unitig_identifier, const std::vector<double>& abundances, const std::vector<double>& fractions) override {
//         if (abundances.size() != m_num_cols || (m_write_fraction && fractions.size() != m_num_cols)) {
//             throw std::invalid_argument("Row size does not match number of samples.");
//         }

//         // --- Select Hyperslab (the new row) in the file ---
//         hsize_t count[2] = {1, m_num_cols};
//         hsize_t offset[2] = {m_current_row, 0};
//         H5::DataSpace file_dataspace = m_abundance_dataset.getSpace();
//         file_dataspace.selectHyperslab(H5S_SELECT_SET, count, offset);

//         hsize_t count_1d[1] = {1};
//         hsize_t offset_1d[1] = {m_current_row};
//         H5::DataSpace file_dataspace_1d = m_id_dataset.getSpace();
//         file_dataspace_1d.selectHyperslab(H5S_SELECT_SET, count_1d, offset_1d);

//         // --- Define Memory Dataspace for the data being written ---
//         H5::DataSpace mem_dataspace(2, count);
//         H5::DataSpace mem_dataspace_1d(1, count_1d);

//         // --- Write Data ---
//         m_abundance_dataset.write(abundances.data(), H5::PredType::NATIVE_DOUBLE, mem_dataspace, file_dataspace);
//         if (m_write_fraction) {
//             m_fraction_dataset.write(fractions.data(), H5::PredType::NATIVE_DOUBLE, mem_dataspace, file_dataspace);
//         }

//         const char* id_c_str = unitig_identifier.c_str();
//         H5::StrType id_type(H5::PredType::C_S1, H5T_VARIABLE);
//         m_id_dataset.write(&id_c_str, id_type, mem_dataspace_1d, file_dataspace_1d);

//         m_current_row++;
//     }

// private:
//     H5::H5File m_file;
//     H5::DataSet m_id_dataset;
//     H5::DataSet m_abundance_dataset;
//     H5::DataSet m_fraction_dataset;
//     bool m_write_fraction;
//     hsize_t m_num_cols;
//     hsize_t m_current_row;
//     hsize_t m_num_rows;
// };

#endif // MATRIX_WRITER_H