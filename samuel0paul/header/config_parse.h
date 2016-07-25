#pragma once

#include <tuple>
#include <vector>

#include "defs.h"


/* Indexes of differents pieces of information in the
 * Config tuples representing points */
#define POINT_X		0
#define POINT_Y		1
#define POINT_TEMP	2

/* Configuration structure used in conjunction with
 * the configuration file. Please change this in accordance
 * with the chosen format.
 *
 * Current format:
 * [ [north | east | south | west] value | initial value | [point | point_src] x y value ]
 * */
struct Config
{
	Config()
		: initial_temp{static_cast<prec_t>(0)}
		, north_source{false}
		, north_temp{static_cast<prec_t>(0)}
		, east_source{false}
		, east_temp{static_cast<prec_t>(0)}
		, south_source{false}
		, south_temp{static_cast<prec_t>(0)}
		, west_source{false}
		, west_temp{static_cast<prec_t>(0)}
	{}

	prec_t initial_temp;

	bool north_source;
	prec_t north_temp;

	bool east_source;
	prec_t east_temp;

	bool south_source;
	prec_t south_temp;

	bool west_source;
	prec_t west_temp;

	std::vector<std::tuple<dim_t, dim_t, prec_t>> point_sources;
	std::vector<std::tuple<dim_t, dim_t, prec_t>> point_initial_temps;
};

bool parse_config(struct Config& config, const char* filename);

void print_config(const struct Config& config);
