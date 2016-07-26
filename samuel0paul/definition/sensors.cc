#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>

std::vector<int> get_temperatures() {
	system("./get_temperatures.sh");
	std::vector<int> res;
	std::string str;
	std::ifstream ifstr("temps");
	while (getline(ifstr, str))
		res.emplace_back(std::stoi(str));
	return res;
}

int main() {
	std::vector<int> temps = get_temperatures();
	for (int i : temps)
		std::cout << "temp: '" << i << "';"<< std::endl;
}
