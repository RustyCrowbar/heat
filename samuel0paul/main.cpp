/**
*	2-D Steady State Conduction without Heat Generation | main.cpp
*
*	@libs	[Nodes.a, NodesHelper.a]
*	@header	[Nodes.h, NodesHelper.h]
*
*	@author Samuel0Paul <paulsamuelvishesh@live.com>
**/

#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <thread>

#include "header/Nodes.h"
#include "header/NodesHelper.h"

#include "test/test.cpp"

using std::cout;	using std::endl;
using std::clog;
using std::cin;
using std::make_pair;

using prec_t = long double;

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

int main(int argc, char const *argv[])
{
	cout << std::nounitbuf;
	cout << std::setprecision(4) << std::fixed << std::boolalpha;
	cout << "############################ HMT Assignment ##########################" << endl
		 << " 2-D Steady State Conduction with and without Heat Generation" << endl << endl
		 << " Author: Samuel Paul Vishesh [UR11ME145] <paulsamuelvishesh@live.com>" << endl
		 << "######################################################################" << endl << endl;

	char filename[30] = { '\0' };
	int iter = 0;
	const uint64_t x_len = 12;
	const uint64_t y_len = 12;
	const prec_t epsilon = 1.0;
	unsigned long long ticks = 0;

	HMT::Nodes<prec_t> nodes(x_len, y_len, 100.0);
	nodes.setHeatSource(7, 9, 400.0);
	nodes.canUseThreads(false);
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
