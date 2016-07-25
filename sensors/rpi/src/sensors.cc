#include <iostream>
#include <fstream>
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
	boost::filesystem::path p("../sensors/test");
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
	for (std::string s : devices)
		std::cout << "device: '" << s << "'; temperature: '" << get_temp(s) << "'" << std::endl;
}
