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

#ifndef NODES_H
#define NODES_H

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <chrono>
#include <tbb/tbb.h>
#include <thread>

#include "defs.h"


namespace HMT
{

template<typename T>
class Nodes
{
public:
	enum Wall
	{
		NORTH,
		EAST,
		SOUTH,
		WEST
	};

	Nodes() = default;
	Nodes(const uint64_t nodeX, const uint64_t nodeY, const T initial_temp,
	      const size_t grain_size = 2);
	virtual ~Nodes() = default;

	void setWallSource(const enum Wall wall, const T& temp);
	void setHeatSource(const uint64_t& posX, const uint64_t& posY, const T& temp);
	void setTemperature(const uint64_t posX, const uint64_t posY, const T temp);
	void canUseThreads(const bool choice) noexcept(true);
	bool canUseThreads(void) const noexcept(true);
	bool hasCalculated(void) const;
	void calculate(const prec_t epsilon);
	void clear(const T temp);
	bool hasHeatSource(void) const noexcept(true);
	T getTemp(const uint64_t& posX, const uint64_t& posY) const;
	std::chrono::nanoseconds getDuration(void) const;
	uint64_t getItterCount(void) const;

	template<typename T1> friend std::ostream& operator<<(std::ostream&, const Nodes<T1>&);

	void testBuffers(void) const;

protected:
	void initBuffer(const T initial_temp);
	void calculateWThread(const prec_t epsilon);
	void calculateWoutThread(const prec_t epsilon);
	void calculateOuterNodes(void);
		
private:
	bool _hasHeatSource, _hasCalculated, _canUseThreads;
	uint64_t _nodeX, _nodeY, _itterCnt;
	std::vector<std::vector<std::pair<T, bool>>> _nodes, _nodesOld;
	std::chrono::time_point<std::chrono::high_resolution_clock> _startTime, _endTime;
	const size_t grain_size_;
	prec_t diff_;
};

template<typename T1> std::ostream& operator<<(std::ostream&, const Nodes<T1>&);

}

#include "../definition/Nodes.cxx"

#endif
