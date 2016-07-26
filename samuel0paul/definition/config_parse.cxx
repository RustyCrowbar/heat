#include <boost/lexical_cast.hpp>
#include <fstream>
#include <iostream>

#include "../header/config_parse.h"
#include "../header/sensors.hh"

std::vector<int> temps;

prec_t read_temp(std::string str)
{
	if (str.compare(0, 6, "sensor") == 0) {
		uint32_t nb = std::stoi(str.substr(6, 42));
		std::cout << "using sensor '" << nb << "', reading '" << temps[nb] / 1000 << "'C" << std::endl;
		return boost::lexical_cast<prec_t>(temps[nb]) / 1000;
	} else
		return boost::lexical_cast<prec_t>(str) / 1000;
}

bool parse_config(struct Config& config, const char* filename)
{
	std::string str;
	std::ifstream confstr(filename);
	temps = get_temperatures();
	std::cout << "Got " << temps.size() << " values from sensors." << std::endl;

	if (!confstr.is_open())
	{
		std::cerr << "Error opening file " << filename << std::endl;
		return false;
	}

	while (confstr >> str)
	{
		if (str.compare("north") == 0)
		{
			// North wall source
			config.north_source = true;
			confstr >> str;
			config.north_temp = read_temp(str);
		}
		else if (str.compare("east") == 0)
		{
			// East wall source
			config.east_source = true;
			confstr >> str;
			config.east_temp = read_temp(str);
		}
		else if (str.compare("south") == 0)
		{
			// South wall source
			config.south_source = true;
			confstr >> str;
			config.south_temp = read_temp(str);
		}
		else if (str.compare("west") == 0)
		{
			// West wall source
			config.west_source = true;
			confstr >> str;
			config.west_temp = read_temp(str);
		}
		else if (str.compare("point_src") == 0)
		{
			// Punctual source
			confstr >> str;
			dim_t x = boost::lexical_cast<dim_t>(str);
			confstr >> str;
			dim_t y = boost::lexical_cast<dim_t>(str);
			confstr >> str;
			prec_t temp = read_temp(str);
			config.point_sources.emplace_back(x, y, temp);
		}
		else if (str.compare("point") == 0)
		{
			// Punctual initial temperature
			confstr >> str;
			dim_t x = boost::lexical_cast<dim_t>(str);
			confstr >> str;
			dim_t y = boost::lexical_cast<dim_t>(str);
			confstr >> str;
			prec_t temp = read_temp(str);
			config.point_initial_temps.emplace_back(x, y, temp);

		}
		else if (str.compare("initial") == 0)
		{
			// Initialization temperature for all points
			confstr >> str;
			config.initial_temp = read_temp(str);
		}
		else
		{
			std::cerr << "Error parsing configuration file" << std::endl
				<< "unexpected token: " << str
				<< std::endl;
			return false;
		}
	}

	return true;
}

void print_config(const struct Config& config)
{
	std::cout << std::boolalpha
		<< "Configuration:"
		<< std::endl;

	std::cout << "\tInitial temperature: " << config.initial_temp
		<< std::endl;
	std::cout << "\tIndividual initial temperatures:" << std::endl;
	for (auto& point : config.point_initial_temps)
	{
		std::cout << "\t" << std::get<0>(point) << ", "
			<< std::get<1>(point) << "\t"
			<< std::get<2>(point)
			<< std::endl;
	}

	std::cout << "\tNorth source: " << config.north_source;
	if (config.north_source)
		std::cout << ", " << config.north_temp;
	std::cout << std::endl;

	std::cout << "\tEast source: " << config.east_source;
	if (config.east_source)
		std::cout << ", " << config.east_temp;
	std::cout << std::endl;

	std::cout << "\tSouth source: " << config.south_source;
	if (config.south_source)
		std::cout << ", " << config.south_temp;
	std::cout << std::endl;

	std::cout << "\tWest source: " << config.west_source;
	if (config.west_source)
		std::cout << ", " << config.west_temp;
	std::cout << std::endl;

	std::cout << "\tIndividual sources:" << std::endl;
	for (auto& point : config.point_sources)
	{
		std::cout << "\t" << std::get<0>(point) << ", "
			<< std::get<1>(point) << "\t"
			<< std::get<2>(point)
			<< std::endl;
	}
}
