#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <filesystem>

#include "../external/json/json.hpp"
#include "common.h"

using json = nlohmann::json;

int main_convert(int argc, char* argv[]) {

  bool ap_flag = false;
  std::string out_fname;
  bool help_opt = false;
  double min = 0.8;

  int c;
  while ((c = getopt(argc, argv, "o:m:hp")) != -1) {
    switch (c) {
      case 'p':
        ap_flag = true;
        break;
      case 'm':
        min = std::atof(optarg);
        break;
      case 'o':
        out_fname = optarg;
        break;
      case 'h':
        help_opt = true;
        break;
      case '?':
        return 1;
      default:
        abort();
    }
  }

  if(argc-optind != 2 || help_opt) {
    std::cout << "Usage: kmat_tools convert [options] <color_names_dump.jsonl> <query_output.jsonl>\n\n";
    std::cout << "Converts the jsonl output of ggcat into an unitig matrix in csv format.\n\n";
    std::cout << "Options:\n";
    std::cout << "  -p      if you want the matrix to be absence/presence (i.e. 0/1) and not\n";
    std::cout << "           with k-mer presence ratios.\n";
    std::cout << "  -m float minimum value to set the presence to 1 [0.8]\n";
    std::cout << "  -o FILE  write unitig matrix to FILE [stdout]\n";
    std::cout << "  -h       print this help message\n";
    return 0;
  }

    // Get the input and output filenames from the arguments and check that the input file exist.
    
    std::string color_dump_Filename = argv[optind];
    if(!std::filesystem::exists(color_dump_Filename.c_str())) {
        std::cerr << "[error] color dump file \"" << color_dump_Filename << "\" does not exist" << std::endl;
        return 1;
    }

    std::string color_query_Filename = argv[optind+1];
    if(!std::filesystem::exists(color_query_Filename.c_str())) {
        std::cerr << "[error] query output file \"" << color_query_Filename << "\" does not exist" << std::endl;
        return 1;
    }

    std::vector<std::string> color_names;  // Vector to store the values of "x"
    std::string line;
    color_names.push_back("Unitigs_id");
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
            std::cerr << "Parse error. Cannot find field color_index in a line of the color dump file. " << e.what() << std::endl;
            return 1;
        }
    }

    colorDumpFile.close();

    // Print the values stored in the vector
    //std::cout << "Values of key 'color_name': " << std::endl;
    //for (const std::string& name : color_names) {
    //    std::cout << name << std::endl;
    //}

    u_int64_t num_colors {color_names.size() - 1};
    std::vector<int> keys(num_colors) ;
    std::iota (std::begin(keys), std::end(keys), 0); // Fill with 0, 1, ..., color_names.size().

    std::vector<float> presence_values(num_colors, 0);  // Vector to store the values of "x"
    std::string csv_line;

    //std::ofstream outputFile(out_fname);
    std::ostream* fpout = &std::cout;
    std::ofstream ofs;
    if(!out_fname.empty()) {
        ofs.open(out_fname.c_str());
        if(!ofs.good()) {
            std::cerr << "[error] cannot open output file \"" << out_fname << "\"\n";
            return 1;
        }
        fpout = &ofs;
    }

   //*fpout << std::fixed << std::setprecision(2);

    for (uint64_t i = 0; i < color_names.size()-1; i++) {
        *fpout << color_names[i] << ",";
    }

    *fpout << color_names[static_cast<int>(color_names.size()) - 1] << std::endl;
    int line_counter {0};
    float curr_value;
    std::ifstream colorQueryFile(color_query_Filename);
    while (std::getline(colorQueryFile, line)) {
        try {
            // Parse the line as a JSON object
            json jsonObj = json::parse(line);

            std::fill(presence_values.begin(), presence_values.end(), 0);;  // Vector to store the values of "x"
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
                            if (ap_flag) {
                                if (curr_value > min){
                                    presence_values[index] = 1;
                                }
                                else {
                                    presence_values[index] = 0;
                                }
                            } else {
                                presence_values[index] = curr_value;
                            }
                        }
                    }
                }
            *fpout << line_counter++ << ",";
            for (uint64_t i = 0; i < presence_values.size()-1; i++) {
                *fpout << presence_values[i] << ",";
                }
            *fpout << presence_values[presence_values.size() - 1] << std::endl;
            
            }
        } catch (json::parse_error& e) {
            std::cerr << "Parse error. Error while readins ggcat output. Check that the format is correct." << e.what() << std::endl;
            return 1;
        }
    }

    // closing files (output depends)
    colorQueryFile.close();
    if(!out_fname.empty()) {
        ofs.close();
    }

    return 0;
}