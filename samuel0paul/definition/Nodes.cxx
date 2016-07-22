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
using std::make_tuple;

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
				if (!std::get<1>(nodes_[i][j])) {
					std::get<0>(nodes_[i][j]) =
						(std::get<0>(nodesOld_[i - 1][j]) +
						 std::get<0>(nodesOld_[i + 1][j]) +
						 std::get<0>(nodesOld_[i][j - 1]) +
						 std::get<0>(nodesOld_[i][j + 1]))
						/ std::get<2>(nodes_[i][j]);
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
				if (local_diff < std::fabs(std::get<0>(nodesOld_[i][j]) - std::get<0>(nodes_[i][j])))
					local_diff = std::fabs(std::get<0>(nodesOld_[i][j]) - std::get<0>(nodes_[i][j]));
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
	for (uint64_t i = 1; i < this->_nodeY - 1; ++i) {
		for (uint64_t j = 1; j < this->_nodeX - 1; ++j) {
			cout << std::get<0>(this->_nodes[i][j]) << ", ";
		}
		cout << endl;
	}
	for (uint64_t i = 1; i < this->_nodeY - 1; ++i) {
		for (uint64_t j = 1; j < this->_nodeX - 1; ++j) {
			cout << std::get<0>(this->_nodesOld[i][j]) << ", ";
		}
		cout << endl;
	}
}

template<typename T>
Nodes<T>::Nodes(const uint64_t nodeX, const uint64_t nodeY,
		const T initial_temp)
{
	this->_nodeX = nodeX + 2;
	this->_nodeY = nodeY + 2;
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

	std::vector<std::tuple<T, bool, int>> tmp;
	tmp.reserve(this->_nodeX);

	for (uint64_t j = 0; j < this->_nodeX; ++j) {
		tmp.push_back(make_tuple(0, false, 0));
	}

	this->_nodes.push_back(tmp);
	this->_nodesOld.push_back(std::move(tmp));

	for (uint64_t i = 0; i < this->_nodeY - 2; ++i) {
		std::vector<std::tuple<T, bool, int>> tmp;
		tmp.reserve(this->_nodeX);

		tmp.push_back(make_tuple(0, false, 0));
		for (uint64_t j = 0; j < this->_nodeX - 2; ++j) {
			tmp.push_back(make_tuple(initial_temp, false,
						 4 -
						 (i == 0 || i == this->_nodeX - 3) -
						 (j == 0 || j == this->_nodeY - 3)));
		}
		tmp.push_back(make_tuple(0, false, 0));

		this->_nodes.push_back(tmp);
		this->_nodesOld.push_back(std::move(tmp));
	}

	std::vector<std::tuple<T, bool, int>> tmp2;
	tmp2.reserve(this->_nodeX);

	for (uint64_t j = 0; j < this->_nodeX; ++j) {
		tmp2.push_back(make_tuple(0, false, 0));
	}

	this->_nodes.push_back(tmp2);
	this->_nodesOld.push_back(std::move(tmp2));

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
	for (uint64_t i = 0; i < this->_nodeY - 2; ++i) {
		setHeatSource(0, i, westTemp);
		setHeatSource(this->_nodeX - 3, i, eastTemp);
	}
	for (uint64_t i = 1; i < this->_nodeX - 2; ++i) {
		setHeatSource(i, 0, northTemp);
		setHeatSource(i, this->_nodeY - 3, southTemp);
	}
}

template<typename T>
void Nodes<T>::setHeatSource(const uint64_t& posX, const uint64_t& posY, const T& temp)
{
	this->_hasHeatSource = true;
	std::get<0>(this->_nodes[posY + 1][posX + 1]) = temp;
	std::get<1>(this->_nodes[posY + 1][posX + 1]) = true;
}

template <typename T>
void Nodes<T>::setTemperature(const uint64_t posX, const uint64_t posY, const T temp)
{
	std::get<0>(this->_nodes[posY + 1][posX + 1]) = temp;
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

		// Every non-border node
		for (uint64_t i = 1; i < this->_nodeY - 1; ++i) {
			for (uint64_t j = 1; j < this->_nodeX - 1; ++j) {
				if (std::get<1>(this->_nodes[i][j]) != true) {
					std::get<0>(this->_nodes[i][j]) =
						(std::get<0>(this->_nodesOld[i - 1][j]) +
						 std::get<0>(this->_nodesOld[i + 1][j]) +
						 std::get<0>(this->_nodesOld[i][j - 1]) +
						 std::get<0>(this->_nodesOld[i][j + 1])) /
						std::get<2>(this->_nodes[i][j]);
					if (diff < std::fabs(std::get<0>(this->_nodesOld[i][j]) - std::get<0>(this->_nodes[i][j]))) {
						diff = std::fabs(std::get<0>(this->_nodesOld[i][j]) - std::get<0>(this->_nodes[i][j]));
					}
				}
			}
		}
		this->_endTime = std::chrono::high_resolution_clock::now();
		if (diff <= epsilon)
			this->_hasCalculated = true;
	}
}

template <typename T>
void Nodes<T>::clear(const T temp)
{
	for (size_t i = 1; i < this->_nodeX - 1; ++i)
	{
		for (size_t j = 1; j < this->_nodeY - 1; ++j)
		{
			std::get<0>(this->nodes_[i][j]) = temp;
			std::get<1>(this->nodes_[i][j]) = false;
		}
	}

	this->_hasCalculated = false;
	this->_itterCnt = 0;
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

	// Compute a new iteration
	tbb::parallel_for(tbb::blocked_range2d<size_t>(1, this->_nodeY - 1,
						       1, this->_nodeX - 1),
			  calculator);

	// Find out if we should stop now
	tbb::parallel_reduce(tbb::blocked_range2d<size_t>(1, this->_nodeY - 1,
						          1, this->_nodeX - 1),
			     reductor);
	this->_endTime = std::chrono::high_resolution_clock::now();
	if (reductor.diff <= epsilon)
		this->_hasCalculated = true;
}

template<typename T>
T Nodes<T>::getTemp(const uint64_t& posX, const uint64_t& posY) const
{
	if (this->_hasCalculated)
		return std::get<0>(this->_nodes[posY + 1][posX + 1]);
	else {
		//!TODO implement error handling or error throw mechanism
		return 0;
	}
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const Nodes<T>& obj)
{
	for (auto i = obj._nodes.begin() + 1; i != obj._nodes.end() - 1; ++i) {
		for (auto j = i->begin() + 1; j != i->end() - 1; ++j)
			os << std::get<0>(*j) << " ";
		os << endl;
	}
	return os;
}

}

#endif
