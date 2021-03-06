#include <iostream>
#include <fstream>
#include <thread>
#include <chrono>
#include <boost/filesystem.hpp>

uint32_t get_temp(std::string device_path) {
	std::ifstream ifstr(device_path + "/w1_slave");
	size_t pos;
	std::string res;
	do {
		getline(ifstr, res);
		pos = res.find("t=");
	} while (pos == std::string::npos);
	return std::stoi(res.substr(pos + 2));
}

bool is_suitable(std::string dir) {
	return dir.substr(dir.rfind('/') + 1, 3) == "28-";
}

std::vector<std::string> get_devices() {
	std::vector<std::string> res;
	boost::filesystem::path p("test");
	//boost::filesystem::path p("/sys/bus/w1/devices");

	boost::filesystem::directory_iterator end_itr;

	// cycle through the directory
	for (boost::filesystem::directory_iterator itr(p); itr != end_itr; itr++) {
		std::string path = itr->path().string();
		if (!is_regular_file(itr->path()) && is_suitable(path))
			res.emplace_back(path);
	}
	return res;
}

int main() {
	std::vector<std::string> devices = get_devices();
	int tmp;
	std::ofstream of;
	while (1) {
		of.open("temps", std::ios::trunc);
		for (std::string s : devices) {
			tmp = get_temp(s);
			std::cout << "device: '" << s << "'; temperature: '" << tmp << "'" << std::endl;
			of << tmp << std::endl;
		}
		of.close();
		std::this_thread::sleep_for(std::chrono::seconds(10));
	}
}
