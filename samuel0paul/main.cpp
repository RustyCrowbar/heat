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

/*
void runTest(void)
{
	test::NodesWithoutHeatSrc<prec_t> testNodesWOHSrcTE{12, 30,
		500.0f, 100.0f, 100.0f, 100.0f, 0.0000001f, true};
	testNodesWOHSrcTE.test();
	test::NodesWithoutHeatSrc<prec_t> testNodesWOHSrcWTE{12, 30,
		500.0f, 100.0f, 100.0f, 100.0f, 0.0000001f, false};
	testNodesWOHSrcWTE.test();

	std::vector<std::pair<std::pair<uint64_t, uint64_t>, prec_t>> heatSrcs = {
		make_pair(make_pair(2, 2), 300.0f),
		make_pair(make_pair(5, 5), -1000.0f)
	};
	test::NodesWithHeatSrc<prec_t> testNodesWHSrcTE{12, 30,
		500.0f, 100.0f, 100.0f, 100.0f, 0.0000001f, false,
		heatSrcs};
	testNodesWHSrcTE.test();
	test::NodesWithHeatSrc<prec_t> testNodesWHSrcWTE{12, 30,
		500.0f, 100.0f, 100.0f, 100.0f, 0.0000001f, true,
		heatSrcs};
	testNodesWHSrcWTE.test();
}
*/

static void print_info(dim_t x_len, dim_t y_len, prec_t epsilon,
		       bool using_threads)
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
		<< "\t" << "Using parallel version: " << using_threads << endl;
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
	unsigned long long ticks = 0;
	int opt;
	struct Config conf;

	while ((opt = getopt(argc, argv, "w:h:e:p")) != -1)
	{
		switch (opt)
		{
		case 'w':
			x_len = boost::lexical_cast<dim_t>(optarg);
			break;
		case 'h':
			y_len = boost::lexical_cast<dim_t>(optarg);
			break;
		case 'e':
			epsilon = boost::lexical_cast<prec_t>(optarg);
			break;
		case 'p':
			using_threads = true;
			break;
		default:
			printf("Usage: %s [-w width | -h height | -e epsilon | -p]\n", argv[0]);
			return 1;
		}
	}
	if (!parse_config(conf, "heat.config"))
		return 1;

	print_info(x_len, y_len, epsilon, using_threads);
	print_config(conf);
	std::cout << "######" << std::endl << std::endl;

	HMT::Nodes<prec_t> nodes(x_len, y_len, 100.0);
	nodes.canUseThreads(using_threads);
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
	cout << "Iterations: " << iter << endl;

	//runTest();
	
	return 0;
}
