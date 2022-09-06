#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
namespace Pulsejet::Shims {
	inline float CosF(float x) {
		return cosf(x);
	}

	inline float Exp2f(float x) {
		return exp2f(x);
	}

	inline float SinF(float x) {
		return sinf(x);
	}

	inline float SqrtF(float x) {
		return sqrtf(x);
	}
}
#include <Pulsejet/Pulsejet.hpp>
#define FMT_HEADER_ONLY
#include <fmt/format.h>

using namespace std;

static vector<uint8_t> ReadFile(const char *fileName) {
	ifstream inputFile(fileName, ios::binary | ios::ate);
	const size_t inputFileSize = inputFile.tellg();
	inputFile.seekg(0, ios::beg);

	vector<uint8_t> ret(inputFileSize);
	inputFile.read(reinterpret_cast<char *>(ret.data()), inputFileSize);

	return ret;
}

static vector<float> ReadFileFloats(const char *fileName) {
	ifstream inputFile(fileName, ios::binary | ios::ate);
	const size_t inputFileSize = inputFile.tellg();
	inputFile.seekg(0, ios::beg);

	vector<float> ret(inputFileSize);
	inputFile.read(reinterpret_cast<char *>(ret.data()), inputFileSize);

	return ret;
}

static void SaveFile(const char *fileName, const uint8_t *data, size_t size) {
	ofstream outputFile(fileName, ios::binary | ios::ate);
	outputFile.write(reinterpret_cast<const char *>(data), size);
}

static void SaveFile(const char *fileName, const float *data, size_t sample_count) {
	auto ptr = reinterpret_cast<const uint8_t *>(data);
	SaveFile(fileName, ptr, sample_count * sizeof(float));
}

static void SaveFile(const char *fileName, vector<float> &data) {
	SaveFile(fileName, data.data(), data.size());
}

static void SaveFile(const char *fileName, vector<uint8_t> &data) {
	SaveFile(fileName, data.data(), data.size());
}

int main() {
	fmt::print("hello\n");
	auto input = ReadFileFloats("orig.raw");
	// auto numSamples = input.size();
	auto numSamples = 441000;
	fmt::print("numSamples: {}\n", numSamples);
	double totalBitsEstimate;
	auto encoded_stream = Pulsejet::Encode(
		reinterpret_cast<const float *>(input.data()),
		numSamples,
		44100,
		16,
		totalBitsEstimate);
	fmt::print("totalBitsEstimate: {:.2f}\n", totalBitsEstimate);
	SaveFile("encoded.pj", encoded_stream);
	fmt::print("Saved encoded\n");

	uint32_t out_num_samples = 0;
	auto samples = Pulsejet::Decode(encoded_stream.data(), &out_num_samples);
	SaveFile("native.raw", samples, out_num_samples);
	fmt::print("Saved decoded\n");
    return 0;
}
