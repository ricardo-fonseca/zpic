/*
Copyright (C) 2017 Instituto Superior Tecnico

This file is part of the ZPIC Educational code suite

The ZPIC Educational code suite is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

The ZPIC Educational code suite is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with the ZPIC Educational code suite. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __ZDF__
#define __ZDF__


#include <stdint.h>
#include <stdio.h>

#define zdf_max_dims 3

enum zdf_file_access_mode {ZDF_WRITE, ZDF_READ};

typedef struct {
	FILE *fp;
	enum zdf_file_access_mode mode;
} t_zdf_file;

enum zdf_axis_type{ zdf_linear, zdf_log10, zdf_log2 };

typedef struct {
	enum zdf_axis_type type;

	double min;
	double max;

	char* label;
	char* units;

} t_zdf_grid_axis;

typedef struct {
	uint32_t ndims;
	uint64_t nx[zdf_max_dims];
	char* label;
	char* units;

    t_zdf_grid_axis *axis;
} t_zdf_grid_info;

typedef struct {
	int32_t n;
	double t;
	char* time_units;
} t_zdf_iteration;


typedef struct {
	char* name;
	uint32_t nquants;
	char** quants;
	char** units;

	uint64_t np;
} t_zdf_part_info;


// Low level interface

//int zdf_open_file( t_zdf_file* zdf, const char * restrict filename, enum zdf_file_access_mode mode );

int zdf_open_file( t_zdf_file* zdf, char* filename, enum zdf_file_access_mode mode );

int zdf_close_file( t_zdf_file* zdf );


// High level interface

int zdf_save_grid( const float* data, const t_zdf_grid_info *info, 
	const t_zdf_iteration *iteration, char const path[] );

int zdf_part_file_open( t_zdf_file *file, t_zdf_part_info *info, 
	const t_zdf_iteration *iteration, char const path[] );

int zdf_part_file_add_quant( t_zdf_file *zdf, const char *name, const float* data, 
	const unsigned np );


#endif


