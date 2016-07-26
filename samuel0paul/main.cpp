/**
*	2-D Steady State Conduction without Heat Generation | main.cpp
*
*	@libs	[Nodes.a, NodesHelper.a]
*	@header	[Nodes.h, NodesHelper.h]
*
*	@author Samuel0Paul <paulsamuelvishesh@live.com>
**/

#include <atomic>
#include <boost/lexical_cast.hpp>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <thread>
#include <unistd.h>

#include "header/config_parse.h"
#include "header/defs.h"
#include "header/Nodes.h"
#include "header/NodesHelper.h"

using std::cout;
using std::endl;
using std::clog;
using std::cin;
using std::make_pair;


static void print_usage(const char* program_name)
{
	std::cout << "Usage: " << program_name << " [ " <<
		"-x width | " <<
		"-y height | " <<
		"-e epsilon | " <<
		"-p[grain_size] ]" <<
		std::endl;
}

static void print_info(dim_t x_len, dim_t y_len, prec_t epsilon,
		       bool using_threads, size_t grain_size)
{
	cout << std::nounitbuf;
	cout << std::setprecision(4) << std::fixed << std::boolalpha;
	cout << "############################ HMT Assignment ##########################" << endl
		 << " 2-D Steady State Conduction with and without Heat Generation" << endl << endl
		 << " Author: Samuel Paul Vishesh [UR11ME145] <paulsamuelvishesh@live.com>" << endl
		 << "######################################################################" << endl
		 << " Parallelized and enhanced by:" << endl
		 << " Valentin Grimaldi, Mathieu Corre, Geoffrey Le Gourriérec" << endl
		 << "######################################################################" << endl
		 << endl; 

	std::cout << "Computation model:" << endl
		<< "\t" << y_len << " lines" << endl
		<< "\t" << x_len << " columns" << endl
		<< "\t" << "Precision: " << epsilon << endl
		<< "\t" << "Using parallel version: " << using_threads;
	if (using_threads)
		std::cout << " [Granularity:" << grain_size << "]";
	std::cout << std::endl;
}

static bool check_config(struct Config& config, dim_t x_len, dim_t y_len)
{
	for (auto& point : config.point_initial_temps)
		if (std::get<POINT_X>(point) >= x_len || std::get<POINT_Y>(point) >= y_len)
			return false;
	for (auto& point : config.point_sources)
		if (std::get<POINT_X>(point) >= x_len || std::get<POINT_Y>(point) >= y_len)
			return false;
	return true;
}

static void apply_config(HMT::Nodes<prec_t>& nodes, struct Config& config)
{
	for (auto& point : config.point_initial_temps)
		nodes.setTemperature(std::get<POINT_X>(point), std::get<POINT_Y>(point),
				     std::get<POINT_TEMP>(point));

	if (config.north_source)
		nodes.setWallSource(HMT::Nodes<prec_t>::NORTH, config.north_temp);
	if (config.east_source)
		nodes.setWallSource(HMT::Nodes<prec_t>::EAST, config.east_temp);
	if (config.south_source)
		nodes.setWallSource(HMT::Nodes<prec_t>::SOUTH, config.south_temp);
	if (config.west_source)
		nodes.setWallSource(HMT::Nodes<prec_t>::WEST, config.west_temp);

	for (auto& point : config.point_sources)
		nodes.setHeatSource(std::get<POINT_X>(point), std::get<POINT_Y>(point),
				    std::get<POINT_TEMP>(point));
}

int main(int argc, char *argv[])
{
	char filename[30] = { '\0' };
	int iter = 0;
	// Computation parameters and their default values
	dim_t x_len = 30;
	dim_t y_len = 30;
	prec_t epsilon = 1.0;
	bool using_threads = false;
	size_t grain_size = 2;
	unsigned long long ticks = 0;
	int opt;
	struct Config conf;

	while ((opt = getopt(argc, argv, "x:y:e:p::h")) != -1)
	{
		switch (opt)
		{
		case 'x':
			x_len = boost::lexical_cast<dim_t>(optarg);
			break;
		case 'y':
			y_len = boost::lexical_cast<dim_t>(optarg);
			break;
		case 'e':
			epsilon = boost::lexical_cast<prec_t>(optarg);
			break;
		case 'p':
			using_threads = true;
			if (optarg)
				grain_size = boost::lexical_cast<size_t>(optarg);
			break;
		case 'h':
			print_usage(argv[0]);
			printf("\t-x width  : set number of columns\n");
			printf("\t-y height : set number of lines\n");
			printf("\t-e epsilon: set precision\n");
			printf("\t-p[arg]   : use parallel version with granularity=arg (default:2)\n");
			printf("\t-h        : show this help\n");
			return 0;
		default:
			print_usage(argv[0]);
			return 1;
		}
	}
	if (!parse_config(conf, "heat.config"))
		return 1;
	if (!check_config(conf, x_len, y_len))
	{
		std::cerr << "Error: some heat sources are out of bounds."
			<< std::endl
			<< "(" << x_len << " columns, "
			<< y_len << " lines)"
			<< std::endl;
		return 2;
	}

	print_info(x_len, y_len, epsilon, using_threads, grain_size);
	print_config(conf);
	std::cout << "######" << std::endl << std::endl;

	HMT::Nodes<prec_t> nodes(x_len, y_len, conf.initial_temp, grain_size);
	nodes.canUseThreads(using_threads);
	apply_config(nodes, conf);
	while (!nodes.hasCalculated())
	{
		nodes.calculate(epsilon);
		ticks += nodes.getDuration().count();

		sprintf(filename, "temperature_%04d.txt", iter++);
		std::ofstream outfile(filename);

		outfile << nodes;

		outfile.close();
	}
	cout << "Time taken: " << ticks << "ns" << endl;
	cout << "Iterations: " << nodes.getItterCount() << endl;

	return 0;
}
