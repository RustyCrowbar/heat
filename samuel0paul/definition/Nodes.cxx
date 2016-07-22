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

using std::cout;	using std::endl;
using std::clog;
using std::cin;
using std::make_pair;

using prec_t = long double;

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
		uint64_t row_begin = range.rows().begin();
		uint64_t row_end = range.rows().end();
		uint64_t col_begin = range.cols().begin();
		uint64_t col_end = range.cols().end();

		for (uint64_t i = row_begin; i < row_end; ++i) {
			for (uint64_t j = col_begin; j < col_end; ++j) {
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
		uint64_t row_begin = range.rows().begin();
		uint64_t row_end = range.rows().end();
		uint64_t col_begin = range.cols().begin();
		uint64_t col_end = range.cols().end();
		prec_t local_diff = diff;

		for (uint64_t i = row_begin; i < row_end; ++i)
			for (uint64_t j = col_begin; j < col_end; ++j) {
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
	if (this->_hasCalculated) {
		return std::chrono::duration_cast<std::chrono::nanoseconds>(this->_endTime - this->_startTime);
	} else {
		//!TODO	implement error handling and notification
		return std::chrono::nanoseconds(0);
	}
}

template<typename T>
uint64_t Nodes<T>::getItterCount(void) const
{
	if (this->_hasCalculated)
		return this->_itterCnt;
	else {
		//!TODO	implement error handling and notification
		return 0;
	}
}

template<typename T>
void Nodes<T>::testBuffers(void) const
{
	for (uint64_t i = 0; i < this->_nodeY; ++i) {
		for (uint64_t j = 0; j < this->_nodeX; ++j) {
			cout << this->_nodes[i][j].first << ", ";
		}
		cout << endl;
	}
	for (uint64_t i = 0; i < this->_nodeY; ++i) {
		for (uint64_t j = 0; j < this->_nodeX; ++j) {
			cout << this->_nodesOld[i][j].first << ", ";
		}
		cout << endl;
	}
}

template<typename T>
Nodes<T>::Nodes(const uint64_t nodeX, const uint64_t nodeY,
		const T initial_temp)
{
	this->_nodeX = nodeX;
	this->_nodeY = nodeY;
	this->initBuffer(initial_temp);
	this->_hasHeatSource = false;
	this->_hasCalculated = false;
	this->_canUseThreads = false;
}

template<typename T>
void Nodes<T>::initBuffer(const T initial_temp)
{
	this->_nodes.reserve(this->_nodeY);
	this->_nodesOld.reserve(this->_nodeY);
	for (uint64_t i = 0; i < this->_nodeY; ++i) {
		std::vector<std::pair<T, bool>> tmp;
		tmp.reserve(this->_nodeX);
		for (uint64_t j = 0; j < this->_nodeX; ++j) {
			tmp.push_back(make_pair(initial_temp, false));
		}
		this->_nodes.push_back(tmp);
		this->_nodesOld.push_back(std::move(tmp));
	}
	for (uint64_t i = 0; i < this->_nodeY; i++) {
		this->_nodes[i].shrink_to_fit();
		this->_nodesOld[i].shrink_to_fit();
	}
	this->_nodes.shrink_to_fit();
	this->_nodesOld.shrink_to_fit();
}

template<typename T>
void Nodes<T>::setWallSources(const T& northTemp, const T& eastTemp, const T& southTemp, const T& westTemp)
{
	for (uint64_t i = 0; i < this->_nodeY; ++i) {
		setHeatSource(0, i, westTemp);
		setHeatSource(this->_nodeX - 1, i, eastTemp);
	}
	for (uint64_t i = 0; i < this->_nodeX; ++i) {
		setHeatSource(i, 0, northTemp);
		setHeatSource(i, this->_nodeY - 1, southTemp);
	}
}

template<typename T>
void Nodes<T>::setHeatSource(const uint64_t& posX, const uint64_t& posY, const T& temp)
{
	this->_hasHeatSource = true;
	this->_nodes[posY][posX].first = temp;
	this->_nodes[posY][posX].second = true;
}

template <typename T>
void Nodes<T>::setTemperature(const uint64_t posX, const uint64_t posY, const T temp)
{
	this->_nodes[posY][posX].first = temp;
}

template<typename T>
void Nodes<T>::canUseThreads(const bool choice) noexcept(true)
{
	this->_canUseThreads = choice;
}

template<typename T>
bool Nodes<T>::canUseThreads(void) const noexcept(true)
{
	return this->_canUseThreads;
}

template<typename T>
bool Nodes<T>::hasCalculated(void) const
{
	return this->_hasCalculated;
}

template <typename T>
void Nodes<T>::calculate(const prec_t epsilon)
{
	// This is an ugly solution. Change it sometime.
	if (this->_canUseThreads)
		calculateWThread(epsilon);
	else
		calculateWoutThread(epsilon);
}

template<typename T>
void Nodes<T>::calculateWoutThread(const prec_t epsilon)
{
	if (!this->_hasCalculated) {
		this->_startTime = std::chrono::high_resolution_clock::now();
		++(this->_itterCnt);
		for (uint64_t i = 0; i < this->_nodeY; ++i)
			for (uint64_t j = 0; j < this->_nodeX; ++j)
				this->_nodesOld[i][j] = this->_nodes[i][j];

		prec_t diff = 0.0f;

		// Corner nodes
		if (!this->_nodes[0][0].second)
		{
			this->_nodes[0][0].first =
				(this->_nodesOld[0][1].first +
				 this->_nodesOld[1][0].first)
				/ 2;
		}
		if (!this->_nodes[0][this->_nodeX - 1].second)
		{
			this->_nodes[0][this->_nodeX - 1].first =
				(this->_nodesOld[0][this->_nodeX - 2].first +
				 this->_nodesOld[1][this->_nodeX - 1].first)
				/ 2;
		}
		if (!this->_nodes[this->_nodeY - 1][0].second)
		{
			this->_nodes[this->_nodeY - 1][0].first =
				(this->_nodesOld[this->_nodeY - 1][1].first +
				 this->_nodesOld[this->_nodeY - 2][0].first)
				/ 2;
		}
		if (!this->_nodes[this->_nodeY - 1][this->_nodeX - 1].second)
		{
			this->_nodes[this->_nodeY - 1][this->_nodeX - 1].first =
				(this->_nodesOld[this->_nodeY - 1][this->_nodeX - 2].first +
				 this->_nodesOld[this->_nodeY - 2][this->_nodeX - 1].first)
				/ 2;
		}

		// Border nodes
		for (uint64_t i = 1; i < this->_nodeY - 1; ++i)
		{
			if (!this->_nodes[i][0].second)
			{
				this->_nodes[i][0].first =
					(this->_nodesOld[i][1].first +
					 this->_nodesOld[i - 1][0].first +
					 this->_nodesOld[i + 1][0].first)
					/ 3;
			}
			if (!this->_nodes[i][this->_nodeX - 1].second)
			{
				this->_nodes[i][this->_nodeX - 1].first =
					(this->_nodesOld[i][this->_nodeX - 2].first +
					 this->_nodesOld[i - 1][this->_nodeX - 1].first +
					 this->_nodesOld[i + 1][this->_nodeX - 1].first)
					/ 3;
			}
		}
		for (uint64_t j = 1; j < this->_nodeX - 1; ++j)
		{
			if (!this->_nodes[0][j].second)
			{
				this->_nodes[0][j].first =
					(this->_nodesOld[1][j].first +
					 this->_nodesOld[0][j - 1].first +
					 this->_nodesOld[0][j + 1].first)
					/ 3;
			}
			if (!this->_nodes[this->_nodeY - 1][j].second)
			{
				this->_nodes[this->_nodeY - 1][j].first =
					(this->_nodesOld[this->_nodeY - 2][j].first +
					 this->_nodesOld[this->_nodeY - 1][j - 1].first +
					 this->_nodesOld[this->_nodeY - 1][j + 1].first)
					/ 3;
			}

		}

		// Every non-border node
		for (uint64_t i = 1; i < this->_nodeY - 1; ++i) {
			for (uint64_t j = 1; j < this->_nodeX - 1; ++j) {
				if (this->_nodes[i][j].second != true) {
					this->_nodes[i][j].first =
						(this->_nodesOld[i - 1][j].first +
						 this->_nodesOld[i + 1][j].first +
						 this->_nodesOld[i][j - 1].first +
						 this->_nodesOld[i][j + 1].first)
						/ 4;
					if (diff < std::fabs(this->_nodesOld[i][j].first - this->_nodes[i][j].first)) {
						diff = std::fabs(this->_nodesOld[i][j].first - this->_nodes[i][j].first);
					}
				}
			}
		}
		this->_endTime = std::chrono::high_resolution_clock::now();
		if (diff <= epsilon)
			this->_hasCalculated = true;
	}
	else
		this->_itterCnt = 0; // Ugly as hell
}

template <typename T>
void Nodes<T>::clear(const T temp)
{
	for (size_t i = 0; i < this->_nodeX; ++i)
	{
		for (size_t j = 0; j < this->_nodeY; ++j)
		{
			this->nodes_[i][j].first = temp;
			this->nodes_[i][j].second = false;
		}
	}

	this->_hasCalculated = false;
}

template<typename T>
void Nodes<T>::calculateWThread(const prec_t epsilon)
{
	this->_startTime = std::chrono::high_resolution_clock::now();
	++(this->_itterCnt);
	for (uint64_t i = 0; i < this->_nodeY; ++i)
		for (uint64_t j = 0; j < this->_nodeX; ++j)
			this->_nodesOld[i][j] = this->_nodes[i][j];

	TBB_Calculator<decltype(this->_nodes.begin())>
		calculator(this->_nodes.begin(), this->_nodesOld.begin());
	TBB_Reductor<decltype(this->_nodes.begin())>
		reductor(this->_nodes.begin(), this->_nodesOld.begin());

	/* dbg
	for (auto i = 1; i < _nodeY - 1; ++i)
		for (auto j = 1; j < _nodeX - 1; ++j)
			if (_nodes[i][j].first != _nodesOld[i][j].first)
				std::cout << "WTF" << std::endl;
	*/
	// Compute a new iteration
	tbb::parallel_for(tbb::blocked_range2d<size_t>(1, this->_nodeY - 1,
						       1, this->_nodeX - 1),
			  calculator);

	// Find out if we should stop now
	tbb::parallel_reduce(tbb::blocked_range2d<size_t>(1, this->_nodeY - 1,
						          1, this->_nodeX - 1),
			     reductor);
	this->_endTime = std::chrono::high_resolution_clock::now();
	// dbg
	std::cout << "diff = " << reductor.diff << std::endl;
	if (reductor.diff <= epsilon)
		this->_hasCalculated = true;
}

template<typename T>
T Nodes<T>::getTemp(const uint64_t& posX, const uint64_t& posY) const
{
	if (this->_hasCalculated)
		return this->_nodes[posY][posX].first;
	else {
		//!TODO implement error handling or error throw mechanism
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
