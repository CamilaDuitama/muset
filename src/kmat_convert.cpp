#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <fmt/format.h>
#include <fmt/ranges.h>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include <kseq++/seqio.hpp>
#include <nlohmann/json.hpp>

#include <kmat_tools/cli/cli_common.h>
#include <kmat_tools/cli/convert.h>


namespace kmat {

int main_convert(convert_opt_t opt)
{
    namespace fs = std::filesystem;
    using json = nlohmann::json;

    if (opt->min_frac_set && !opt->ap_flag){
        spdlog::warn("Minimum fraction (-f) was set without option -p. Option -f will not be used.");
    } 

    // Get the input and output filenames from the arguments and check that the input file exist.
    std::string unitigs_filename = opt->inputs[0];
    if(!fs::is_regular_file(unitigs_filename)){
        throw std::runtime_error(fmt::format("input unitig file \"{}\" does not exist", unitigs_filename));
    }

    std::string color_dump_Filename = opt->inputs[1];
    if(!fs::is_regular_file(color_dump_Filename)) {
        throw std::runtime_error(fmt::format("color dump file \"{}\" does not exist", color_dump_Filename));
    }

    std::string color_query_Filename = opt->inputs[2];
    if(!fs::is_regular_file(color_query_Filename)) {
        throw std::runtime_error(fmt::format("query output file \"{}\" does not exist", color_query_Filename));
    }

    spdlog::info("reading color names");
    
    std::vector<std::string> color_names;  // Vector to store the values of "x"
    std::string line;
    color_names.push_back("UnitigID");
    
    // Parse the color dump file and get the names of the colors (to output in the heades of the csv)
    std::ifstream colorDumpFile(color_dump_Filename);
    while (std::getline(colorDumpFile, line)) {
        try {
            // Parse the line as a JSON object
            json jsonObj = json::parse(line);

            // Check if the key "x" exists in the JSON object
            if (jsonObj.contains("color_name")) {
                // Extract the value of "x" and store it in the vector
                color_names.push_back(jsonObj["color_name"].get<std::string>());
            }
        } catch (json::parse_error& e) {
            spdlog::error("Cannot find field color_index in a line of the color dump file.");
            spdlog::error(e.what());
            std::exit(EXIT_FAILURE);
        }
    }
    colorDumpFile.close();

    spdlog::info("building unitig matrix");

    uint64_t num_colors {color_names.size()-1};
    std::vector<float> presence_values(num_colors, 0);  // Vector to store the values of "x"

    std::ostream* fpout = &std::cout;
    std::ofstream ofs;
    if(!(opt->out_fname).empty()) {
        ofs.open((opt->out_fname).c_str());
        if(!ofs.good()) {
            spdlog::error(fmt::format("cannot open output file \"{}\"", opt->out_fname));
            std::exit(EXIT_FAILURE);
        }
        fpout = &ofs;
    }

    //*fpout << std::fixed << std::setprecision(2);
    std::string separator{opt->out_csv ? "," : " "};

    if (!opt->no_header){
        *fpout << color_names[0];
        for (size_t i = 1; i < color_names.size(); i++) {
            *fpout << separator << color_names[i];
        }
        *fpout << "\n";
    }

    int line_counter{0};
    float curr_value;
    std::ifstream colorQueryFile(color_query_Filename);

    klibpp::KSeq unitig;
    klibpp::SeqStreamIn utg_ssi(unitigs_filename.c_str());
    while (std::getline(colorQueryFile, line)) {
        try {
            // Parse the line as a JSON object
            json jsonObj = json::parse(line);

            std::fill(presence_values.begin(), presence_values.end(), 0);  // Vector to store the values of "x"
            // Check if the key "x" exists and is an object
            if (jsonObj.contains("matches") && jsonObj["matches"].is_object()) {
                json nestedObj = jsonObj["matches"];
                for (auto it = nestedObj.begin(); it != nestedObj.end(); ++it) {
                    auto key = it.key();
                    auto value = it.value();
                    if (!key.empty() && value.is_number()) {  // Check if key is not empty and value is a number
                        uint64_t index = std::stoi(key);  // Convert key to integer index
                        curr_value = value.get<float>();
                        if (index < presence_values.size()) {
                            if (opt->ap_flag) {
                                presence_values[index] = (curr_value > 0 && curr_value >= opt->min_frac);
                            } else {
                                presence_values[index] = curr_value;
                            }
                        }
                    }
                }
                utg_ssi >> unitig;
                // std::cerr << unitig.seq << " " << unitig.name << std::endl;
                *fpout << (opt->out_write_seq ? unitig.seq : unitig.name);
                for (size_t i=0; i < presence_values.size(); i++) {
                    *fpout << separator << presence_values[i];
                }
                *fpout << "\n";
            }
        } catch (json::parse_error& e) {
            spdlog::error("Error reading ggcat files. Check that the format is correct.");
            spdlog::error(e.what());
            std::exit(EXIT_FAILURE);
        }
    }

    // closing files (output depends)
    colorQueryFile.close();
    if(!(opt->out_fname).empty()) {
        ofs.close();
        spdlog::info(fmt::format("presence-absence unitig matrix written to \"{}\"", (opt->out_fname).c_str()));
    }

    return 0;
}

};