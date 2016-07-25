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
	TBB_Calculator(T nextNodes, T oldNodes, dim_t x_len, dim_t y_len)
		: nodes_{nextNodes}
		, nodesOld_{oldNodes}
		, x_len_{x_len}
		, y_len_{y_len}
	{}

	void operator()(const tbb::blocked_range2d<size_t>& range) const
	{
		dim_t row_begin = range.rows().begin();
		dim_t row_end = range.rows().end();
		dim_t col_begin = range.cols().begin();
		dim_t col_end = range.cols().end();

		for (dim_t i = row_begin; i < row_end; ++i) {
			for (dim_t j = col_begin; j < col_end; ++j) {

				// Ignore heat sources, they don't change
				if (nodes_[i][j].second)
					continue;

				if (i == 0)
				{
					if (j == 0) // Upper-left corner
						nodes_[i][j].first =
							(nodesOld_[i][j + 1].first +
							 nodesOld_[i + 1][j].first)
							 / 2;
					else if (j == (x_len_ - 1)) // Upper-right corner
						nodes_[i][j].first =
							(nodesOld_[i][j - 1].first +
							 nodesOld_[i + 1][j].first)
							 / 2;
					else // North border
						nodes_[i][j].first =
							(nodesOld_[i + 1][j].first +
							 nodesOld_[i][j - 1].first +
							 nodesOld_[i][j + 1].first)
							 / 3;

				}
				else if (i == (y_len_ - 1))
				{
					if (j == 0) // Lower-left corner
						nodes_[i][j].first =
							(nodesOld_[i][j + 1].first +
							 nodesOld_[i - 1][j].first)
							 / 2;
					else if (j == (x_len_ - 1)) // Lower-right corner
						nodes_[i][j].first =
							(nodesOld_[i][j - 1].first +
							 nodesOld_[i - 1][j].first)
							 / 2;
					else // South border
						nodes_[i][j].first =
							(nodesOld_[i - 1][j].first +
							 nodesOld_[i][j - 1].first +
							 nodesOld_[i][j + 1].first)
							 / 3;
				}
				else if (j == 0) // West border
				{
					nodes_[i][j].first =
						(nodesOld_[i][j + 1].first +
						 nodesOld_[i - 1][j].first +
						 nodesOld_[i + 1][j].first)
						 / 3;
				}
				else if (j == (x_len_ - 1)) // East border
				{
					nodes_[i][j].first =
						(nodesOld_[i][j - 1].first +
						 nodesOld_[i - 1][j].first +
						 nodesOld_[i + 1][j].first)
						 / 3;
				}
				else //	Inner nodes
				{
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
	dim_t x_len_;
	dim_t y_len_;
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
	for (dim_t i = 0; i < this->_nodeY; ++i) {
		for (dim_t j = 0; j < this->_nodeX; ++j) {
			cout << this->_nodes[i][j].first << ", ";
		}
		cout << endl;
	}
	for (dim_t i = 0; i < this->_nodeY; ++i) {
		for (dim_t j = 0; j < this->_nodeX; ++j) {
			cout << this->_nodesOld[i][j].first << ", ";
		}
		cout << endl;
	}
}

template<typename T>
Nodes<T>::Nodes(const dim_t nodeX, const dim_t nodeY,
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
	for (dim_t i = 0; i < this->_nodeY; ++i) {
		std::vector<std::pair<T, bool>> tmp;
		tmp.reserve(this->_nodeX);
		for (dim_t j = 0; j < this->_nodeX; ++j) {
			tmp.push_back(make_pair(initial_temp, false));
		}
		this->_nodes.push_back(tmp);
		this->_nodesOld.push_back(std::move(tmp));
	}
	for (dim_t i = 0; i < this->_nodeY; i++) {
		this->_nodes[i].shrink_to_fit();
		this->_nodesOld[i].shrink_to_fit();
	}
	this->_nodes.shrink_to_fit();
	this->_nodesOld.shrink_to_fit();
}

template<typename T>
void Nodes<T>::setWallSource(const enum Wall wall, const T& temp)
{
	switch (wall)
	{
	case NORTH:
		for (dim_t i = 0; i < this->_nodeX; ++i)
			setHeatSource(i, 0, temp);
		break;

	case EAST:
		for (dim_t i = 0; i < this->_nodeY; ++i)
			setHeatSource(this->_nodeX - 1, i, temp);
		break;
	case SOUTH:
		for (dim_t i = 0; i < this->_nodeX; ++i)
			setHeatSource(i, this->_nodeY - 1, temp);
		break;
	case WEST:
		for (dim_t i = 0; i < this->_nodeY; ++i)
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
	this->_hasHeatSource = true;
	this->_nodes[posY][posX].first = temp;
	this->_nodes[posY][posX].second = true;
}

template <typename T>
void Nodes<T>::setTemperature(const dim_t posX, const dim_t posY, const T temp)
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
		for (dim_t i = 0; i < this->_nodeY; ++i)
			for (dim_t j = 0; j < this->_nodeX; ++j)
				this->_nodesOld[i][j] = this->_nodes[i][j];

		calculateOuterNodes();

		prec_t diff = 0.0f;

		// Every non-border node
		for (dim_t i = 1; i < this->_nodeY - 1; ++i) {
			for (dim_t j = 1; j < this->_nodeX - 1; ++j) {
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
}

template <typename T>
void Nodes<T>::calculateOuterNodes(void)
{
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
	for (dim_t i = 1; i < this->_nodeY - 1; ++i)
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
	for (dim_t j = 1; j < this->_nodeX - 1; ++j)
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
	this->_itterCnt = 0;
}

template<typename T>
void Nodes<T>::calculateWThread(const prec_t epsilon)
{
	this->_startTime = std::chrono::high_resolution_clock::now();
	++(this->_itterCnt);
	for (dim_t i = 0; i < this->_nodeY; ++i)
		for (dim_t j = 0; j < this->_nodeX; ++j)
			this->_nodesOld[i][j] = this->_nodes[i][j];

	TBB_Calculator<decltype(this->_nodes.begin())>
		calculator(this->_nodes.begin(), this->_nodesOld.begin(),
			   this->_nodeX, this->_nodeY);
	TBB_Reductor<decltype(this->_nodes.begin())>
		reductor(this->_nodes.begin(), this->_nodesOld.begin());

	// Compute a new iteration
	tbb::parallel_for(tbb::blocked_range2d<size_t>(0, this->_nodeY - 1,
						       0, this->_nodeX - 1),
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
T Nodes<T>::getTemp(const dim_t& posX, const dim_t& posY) const
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
