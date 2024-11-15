#include <iostream>
#include <string>
#include <chrono>
#include "kff_io.hpp"


/* This main provide a way to test the reading speed of a file.
 * It uses the high level API.
*/


using namespace std;


// --- Test functions ---
void high_level_api_kmer_read(const string & filename);


/** Parse args and perform the right tests
 **/
int main(int argc, char const *argv[])
{
	if (argc != 2) {
		cerr << "to perform this benchmark: reader <file.kff>" << endl;
		exit(1);
	}

	cout << "KFF version 1.0" << endl;
	cout << "API version " << KFF_VERSION_MAJOR << "." << KFF_VERSION_MINOR << endl;

	// Open the file and read metadata
	string filename = argv[1];
	high_level_api_kmer_read(filename);

	return 0;
}


using namespace std::chrono;
int64_t current_time_micros(){
    return (std::chrono::duration_cast< microseconds >(high_resolution_clock::now().time_since_epoch())).count();
}


/** Uses the high level API to read the file. The kmers are enumerated one by one.
 **/
void high_level_api_kmer_read(const string & filename) {
	cout << "Test high_level_api_kmer_read" << endl;

	Kff_reader reader(filename);

	uint8_t * kmer;
	uint8_t * data;

	uint8_t xor_test = 0;

	int64_t read_time = current_time_micros();
	while (reader.has_next()) {
		// read the kmer and its data
		reader.next_kmer(kmer, data);
		// Get the kmer/data caracteristics
		uint64_t k = reader.k;
		uint64_t data_size = reader.data_size;

		// read the values to be sure that the compiler do not over optimize the previous lines
		// kmer
		for (size_t byte_idx=0 ; byte_idx<(k+3)/4 ; byte_idx++) {
			xor_test ^= kmer[byte_idx];
		}
		// data
		for (size_t byte_idx=0 ; byte_idx<data_size ; byte_idx++) {
			xor_test ^= data[byte_idx];
		}
	}

	read_time = current_time_micros() - read_time;
	cout << "Time to read: " << static_cast<double>(read_time)/1000000.0 << endl;
	cout << "Useless xor result: " << (uint64_t)xor_test << endl;
}