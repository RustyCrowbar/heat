/**
The MIT License (MIT)

Copyright (c) 2014 Samuel Vishesh Paul

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
**/

#ifndef NODES_CXX
#define NODES_CXX

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <chrono>
#include <thread>
#include <mutex>
#include <tbb/tbb.h>

#include "../header/Nodes.h"
#include "../header/defs.h"

using std::cout;	using std::endl;
using std::clog;
using std::cin;
using std::make_pair;

using namespace tbb::flow;
typedef continue_node<continue_msg> node_t;
typedef const continue_msg msg_t;

namespace HMT
{

/* Functor given to TBB to calculate a new iteration for all inner nodes.
 * This CANNOT check against an epsilon if it should stop
 * afterwards, because the parallel_for functor must be const.
 */
template <typename T>
class TBB_Calculator
{
public:
	TBB_Calculator(T nextNodes, T oldNodes)
		: nodes_{nextNodes}
		, nodesOld_{oldNodes}
	{}

	void operator()(const tbb::blocked_range2d<size_t>& range) const
	{
		dim_t row_begin = range.rows().begin();
		dim_t row_end = range.rows().end();
		dim_t col_begin = range.cols().begin();
		dim_t col_end = range.cols().end();

		for (dim_t i = row_begin; i < row_end; ++i) {
			for (dim_t j = col_begin; j < col_end; ++j) {
				if (!nodes_[i][j].second) {
					nodes_[i][j].first =
						(nodesOld_[i - 1][j].first +
						 nodesOld_[i + 1][j].first +
						 nodesOld_[i][j - 1].first +
						 nodesOld_[i][j + 1].first)
						/ 4;
				}
			}
		}
	}

public:
	T nodes_;
	T nodesOld_; 
};

/* Functor given to TBB to find the largest value difference between
 * the previous iteration and the current one. */
template <typename T>
class TBB_Reductor
{
public:
	TBB_Reductor(T curNodes, T oldNodes)
		: diff{0}
		, nodes_{curNodes}
		, nodesOld_{oldNodes}
	{}

	TBB_Reductor(TBB_Reductor& other, tbb::split)
		: diff(other.diff)
		, nodes_{other.nodes_}
		, nodesOld_{other.nodesOld_}
	{}

	void operator()(const tbb::blocked_range2d<size_t>& range)
	{
		dim_t row_begin = range.rows().begin();
		dim_t row_end = range.rows().end();
		dim_t col_begin = range.cols().begin();
		dim_t col_end = range.cols().end();
		prec_t local_diff = diff;

		for (dim_t i = row_begin; i < row_end; ++i)
			for (dim_t j = col_begin; j < col_end; ++j) {
				if (local_diff < std::fabs(nodesOld_[i][j].first - nodes_[i][j].first))
					local_diff = std::fabs(nodesOld_[i][j].first - nodes_[i][j].first);
			}
		diff = local_diff;
	}
	
	void join(const TBB_Reductor& other)
	{
		diff = std::max(diff, other.diff);
	}

public:
	prec_t diff;

private:
	T nodes_;
	T nodesOld_;
};

template<typename T>
std::chrono::nanoseconds Nodes<T>::getDuration(void) const
{
	if (_hasCalculated) {
		return std::chrono::duration_cast<std::chrono::nanoseconds>(_endTime - _startTime);
	} else {
		return std::chrono::nanoseconds(0);
	}
}

template<typename T>
uint64_t Nodes<T>::getItterCount(void) const
{
	if (_hasCalculated)
		return _itterCnt;
	else {
		return 0;
	}
}

template<typename T>
void Nodes<T>::testBuffers(void) const
{
	for (dim_t i = 0; i < _nodeY; ++i) {
		for (dim_t j = 0; j < _nodeX; ++j) {
			cout << _nodes[i][j].first << ", ";
		}
		cout << endl;
	}
	for (dim_t i = 0; i < _nodeY; ++i) {
		for (dim_t j = 0; j < _nodeX; ++j) {
			cout << _nodesOld[i][j].first << ", ";
		}
		cout << endl;
	}
}

template<typename T>
Nodes<T>::Nodes(const dim_t nodeX, const dim_t nodeY,
		const T initial_temp)
	: _hasHeatSource{false}
	, _hasCalculated{false}
	, _canUseThreads{false}
	, _nodeX{nodeX}
	, _nodeY{nodeY}
	, _itterCnt{0}
{
	initBuffer(initial_temp);
}

template<typename T>
void Nodes<T>::initBuffer(const T initial_temp)
{
	_nodes.reserve(_nodeY);
	_nodesOld.reserve(_nodeY);
	for (dim_t i = 0; i < _nodeY; ++i) {
		std::vector<std::pair<T, bool>> tmp;
		tmp.reserve(_nodeX);
		for (dim_t j = 0; j < _nodeX; ++j) {
			tmp.push_back(make_pair(initial_temp, false));
		}
		_nodes.push_back(tmp);
		_nodesOld.push_back(std::move(tmp));
	}
	for (dim_t i = 0; i < _nodeY; i++) {
		_nodes[i].shrink_to_fit();
		_nodesOld[i].shrink_to_fit();
	}
	_nodes.shrink_to_fit();
	_nodesOld.shrink_to_fit();
}

template<typename T>
void Nodes<T>::setWallSource(const enum Wall wall, const T& temp)
{
	switch (wall)
	{
	case NORTH:
		for (dim_t i = 0; i < _nodeX; ++i)
			setHeatSource(i, 0, temp);
		break;

	case EAST:
		for (dim_t i = 0; i < _nodeY; ++i)
			setHeatSource(_nodeX - 1, i, temp);
		break;
	case SOUTH:
		for (dim_t i = 0; i < _nodeX; ++i)
			setHeatSource(i, _nodeY - 1, temp);
		break;
	case WEST:
		for (dim_t i = 0; i < _nodeY; ++i)
			setHeatSource(0, i, temp);
		break;
	default:
		std::cerr << "Error: unknown wall" << std::endl;
		break;
	}
}

template<typename T>
void Nodes<T>::setHeatSource(const dim_t& posX, const dim_t& posY, const T& temp)
{
	_hasHeatSource = true;
	_nodes[posY][posX].first = temp;
	_nodes[posY][posX].second = true;
}

template <typename T>
void Nodes<T>::setTemperature(const dim_t posX, const dim_t posY, const T temp)
{
	_nodes[posY][posX].first = temp;
}

template<typename T>
void Nodes<T>::canUseThreads(const bool choice) noexcept(true)
{
	_canUseThreads = choice;
}

template<typename T>
bool Nodes<T>::canUseThreads(void) const noexcept(true)
{
	return _canUseThreads;
}

template<typename T>
bool Nodes<T>::hasCalculated(void) const
{
	return _hasCalculated;
}

template <typename T>
void Nodes<T>::calculate(const prec_t epsilon)
{
	if (_canUseThreads)
		calculateWThread(epsilon);
	else
		calculateWoutThread(epsilon);
}

template<typename T>
void Nodes<T>::calculateWoutThread(const prec_t epsilon)
{
	if (_hasCalculated)
		return;

	_startTime = std::chrono::high_resolution_clock::now();
	++(_itterCnt);
	for (dim_t i = 0; i < _nodeY; ++i)
		for (dim_t j = 0; j < _nodeX; ++j)
			_nodesOld[i][j] = _nodes[i][j];

	calculateOuterNodes();

	prec_t diff = 0.0f;

	// Every non-border node
	for (dim_t i = 1; i < _nodeY - 1; ++i) {
		for (dim_t j = 1; j < _nodeX - 1; ++j) {
			if (_nodes[i][j].second != true) {
				_nodes[i][j].first =
					(_nodesOld[i - 1][j].first +
					 _nodesOld[i + 1][j].first +
					 _nodesOld[i][j - 1].first +
					 _nodesOld[i][j + 1].first)
					/ 4;
				if (diff < std::fabs(_nodesOld[i][j].first - _nodes[i][j].first)) {
					diff = std::fabs(_nodesOld[i][j].first - _nodes[i][j].first);
				}
			}
		}
	}
	_endTime = std::chrono::high_resolution_clock::now();
	if (diff <= epsilon)
		_hasCalculated = true;
}

template <typename T>
void Nodes<T>::calculateOuterNodes(void)
{
	// Corner nodes
	if (!_nodes[0][0].second)
	{
		_nodes[0][0].first =
			(_nodesOld[0][1].first +
			 _nodesOld[1][0].first)
			 / 2;
	}
	if (!_nodes[0][_nodeX - 1].second)
	{
		_nodes[0][_nodeX - 1].first =
			(_nodesOld[0][_nodeX - 2].first +
			 _nodesOld[1][_nodeX - 1].first)
			 / 2;
	}
	if (!_nodes[_nodeY - 1][0].second)
	{
		_nodes[_nodeY - 1][0].first =
			(_nodesOld[_nodeY - 1][1].first +
			 _nodesOld[_nodeY - 2][0].first)
			 / 2;
	}
	if (!_nodes[_nodeY - 1][_nodeX - 1].second)
	{
		_nodes[_nodeY - 1][_nodeX - 1].first =
			(_nodesOld[_nodeY - 1][_nodeX - 2].first +
			 _nodesOld[_nodeY - 2][_nodeX - 1].first)
			 / 2;
	}

	// Border nodes
	for (dim_t i = 1; i < _nodeY - 1; ++i)
	{
		if (!_nodes[i][0].second)
		{
			_nodes[i][0].first =
				(_nodesOld[i][1].first +
				 _nodesOld[i - 1][0].first +
				 _nodesOld[i + 1][0].first)
				 / 3;
		}
		if (!_nodes[i][_nodeX - 1].second)
		{
			_nodes[i][_nodeX - 1].first =
				(_nodesOld[i][_nodeX - 2].first +
				 _nodesOld[i - 1][_nodeX - 1].first +
				 _nodesOld[i + 1][_nodeX - 1].first)
				 / 3;
		}
	}
	for (dim_t j = 1; j < _nodeX - 1; ++j)
	{
		if (!_nodes[0][j].second)
		{
			_nodes[0][j].first =
				(_nodesOld[1][j].first +
				 _nodesOld[0][j - 1].first +
				 _nodesOld[0][j + 1].first)
				 / 3;
		}
		if (!_nodes[_nodeY - 1][j].second)
		{
			_nodes[_nodeY - 1][j].first =
				(_nodesOld[_nodeY - 2][j].first +
				 _nodesOld[_nodeY - 1][j - 1].first +
				 _nodesOld[_nodeY - 1][j + 1].first)
				 / 3;
		}

	}
}

template <typename T>
void Nodes<T>::clear(const T temp)
{
	for (size_t i = 0; i < _nodeX; ++i)
	{
		for (size_t j = 0; j < _nodeY; ++j)
		{
			_nodes[i][j].first = temp;
			_nodes[i][j].second = false;
		}
	}

	_hasCalculated = false;
	_itterCnt = 0;
}

template<typename T>
void Nodes<T>::calculateWThread(const prec_t epsilon)
{
	if (_hasCalculated)
		return;

	_startTime = std::chrono::high_resolution_clock::now();

	size_t iter_cnt = 2;
	tbb::task_group_context tbb_context;
	graph graph(tbb_context);

	node_t check_stop(graph, [&](msg_t&) {

		if (! iter_cnt--)
			tbb_context.cancel_group_execution();
	});

	// Node to copy last iteration nodes to "old" matrix
	node_t copy_to_old(graph, [this](msg_t&) {

		// Increment iteration counter
		++(_itterCnt);

		for (dim_t i = 0; i < _nodeY; ++i)
			for (dim_t j = 0; j < _nodeX; ++j)
				_nodesOld[i][j] = _nodes[i][j];
	});

	// Node to compute temperature on all inner nodes
	node_t inner_nodes(graph, [this](msg_t&) {

		TBB_Calculator<decltype(_nodes.begin())>
			calculator(_nodes.begin(), _nodesOld.begin());

		tbb::parallel_for(tbb::blocked_range2d<size_t>(1, _nodeY - 1,
							       1, _nodeX - 1),
				  calculator);
	});

	// Node to compute temperature on all outer nodes
	node_t outer_nodes(graph, [this](msg_t&) {

		// Corner nodes
		if (!_nodes[0][0].second)
		{
			_nodes[0][0].first =
				(_nodesOld[0][1].first +
				 _nodesOld[1][0].first)
				/ 2;
		}
		if (!_nodes[0][_nodeX - 1].second)
		{
			_nodes[0][_nodeX - 1].first =
				(_nodesOld[0][_nodeX - 2].first +
				 _nodesOld[1][_nodeX - 1].first)
				/ 2;
		}
		if (!_nodes[_nodeY - 1][0].second)
		{
			_nodes[_nodeY - 1][0].first =
				(_nodesOld[_nodeY - 1][1].first +
				 _nodesOld[_nodeY - 2][0].first)
				/ 2;
		}
		if (!_nodes[_nodeY - 1][_nodeX - 1].second)
		{
			_nodes[_nodeY - 1][_nodeX - 1].first =
				(_nodesOld[_nodeY - 1][_nodeX - 2].first +
				 _nodesOld[_nodeY - 2][_nodeX - 1].first)
				/ 2;
		}

		// Border nodes
		for (dim_t i = 1; i < _nodeY - 1; ++i)
		{
			if (!_nodes[i][0].second)
			{
				_nodes[i][0].first =
					(_nodesOld[i][1].first +
					 _nodesOld[i - 1][0].first +
					 _nodesOld[i + 1][0].first)
					/ 3;
			}
			if (!_nodes[i][_nodeX - 1].second)
			{
				_nodes[i][_nodeX - 1].first =
					(_nodesOld[i][_nodeX - 2].first +
					 _nodesOld[i - 1][_nodeX - 1].first +
					 _nodesOld[i + 1][_nodeX - 1].first)
					/ 3;
			}
		}
		for (dim_t j = 1; j < _nodeX - 1; ++j)
		{
			if (!_nodes[0][j].second)
			{
				_nodes[0][j].first =
					(_nodesOld[1][j].first +
					 _nodesOld[0][j - 1].first +
					 _nodesOld[0][j + 1].first)
					/ 3;
			}
			if (!_nodes[_nodeY - 1][j].second)
			{
				_nodes[_nodeY - 1][j].first =
					(_nodesOld[_nodeY - 2][j].first +
					 _nodesOld[_nodeY - 1][j - 1].first +
					 _nodesOld[_nodeY - 1][j + 1].first)
					/ 3;
			}

		}
	});

	// Node to find if we reached required precision level
	node_t check_precision(graph, [this](msg_t&) {

		// This is ugly, but I can't manage another solution, sorry.
		TBB_Reductor<decltype(_nodes.begin())>
			reductor(_nodes.begin(), _nodesOld.begin());

		tbb::parallel_reduce(tbb::blocked_range2d<size_t>(1, _nodeY - 1,
								  1, _nodeX - 1),
				     reductor);
		diff_ = reductor.diff;
	});

	/* Link all nodes together. check_stop and check_precision are
	 * differentiated because if we hold a tbb::task_group_context in
	 * the class, we cannot restart the graph. And we can only capture
	 * the "this" pointer in check_precision. */
	make_edge(check_stop, copy_to_old);
	make_edge(copy_to_old, inner_nodes);
	make_edge(copy_to_old, outer_nodes);
	make_edge(inner_nodes, check_precision);
	make_edge(outer_nodes, check_precision);
	make_edge(check_precision, check_stop);

	check_stop.try_put(continue_msg());
	graph.wait_for_all();

	_endTime = std::chrono::high_resolution_clock::now();
	if (diff_ <= epsilon)
		_hasCalculated = true;
}

template<typename T>
T Nodes<T>::getTemp(const dim_t& posX, const dim_t& posY) const
{
	if (_hasCalculated)
		return _nodes[posY][posX].first;
	else {
		return 0;
	}
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const Nodes<T>& obj)
{
	for (const auto& i : obj._nodes) {
		for (const auto& node : i)
			os << node.first << " ";
		os << endl;
	}
	return os;
}

}

#endif
