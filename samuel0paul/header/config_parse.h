#pragma once

#include <tuple>
#include <vector>

#include "defs.h"

/* Configuration structure used in conjunction with
 * the configuration file. Please change this in accordance
 * to the chose format.
 *
 * Current format:
 * [ [north | east | south | west] value | initial value | point x y value ]
 * */
struct Config
{
	Config()
		: initial_temp{static_cast<prec_t>(100)}
		, north_source{false}
		, north_temp{static_cast<prec_t>(100)}
		, east_source{true}
		, east_temp{static_cast<prec_t>(400)}
		, south_source{true}
		, south_temp{static_cast<prec_t>(250)}
		, west_source{false}
		, west_temp{static_cast<prec_t>(100)}
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
};

bool parse_config(struct Config& config, const char* filename);

void print_config(const struct Config& config);
