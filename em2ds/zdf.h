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

#if !(defined _WIN32 || defined _WIN64)

  #if __SIZEOF_FLOAT__ != 4
  #error This version of ZDF only support systems were sizeof(float) = 4
  #endif

  #if __SIZEOF_DOUBLE__ != 8
  #error This version of ZDF only support systems were sizeof(double) = 8
  #endif

#endif

#include <stdint.h>
#include <stdio.h>

/**
 * Number of bytes in each ZDF file unit
 */
#define BYTES_PER_ZDF_UNIT	(4)

/**
 * @brief Magic byte sequence length
 * 
 */
#define ZDF_MAGIC_LENGTH BYTES_PER_ZDF_UNIT

/**
 * @brief Magic byte sequence identifying ZDF files
 * 
 */
static const char zdf_magic[ZDF_MAGIC_LENGTH] = {'Z','D','F','1'};

/**
 * @brief Maximum array dimensions supported
 * 
 */
#define zdf_max_dims 3

/**
 * @brief ZDF data types
 * 
 */
enum zdf_data_type{ 
	zdf_null,		///< No type
	zdf_int8,		///< 8 bit signed intger
	zdf_uint8,		///< 8 bit unsigned integer
	zdf_int16,		///< 16 bit signed integer
	zdf_uint16,		///< 16 bit unsigned integer
	zdf_int32,		///< 32 bit signed integer
	zdf_uint32,		///< 32 bit unsigned integer
	zdf_int64,		///< 64 bit integer
	zdf_uint64,		///< 64 bit unsigned integer
	zdf_float32,	///< 32 bit floating point
	zdf_float64		///< 64 bit floating point
};

/**
 * @brief File access modes
 * 
 */
enum zdf_file_access_mode {
	ZDF_CREATE,		///< Create new file
	ZDF_READ,		///< Read file
	ZDF_UPDATE		///< Update (read/write) file
};

/**
 * @brief ZDF File
 * 
 */
typedef struct ZDF_File{
	FILE *fp;							///< File pointer
	enum zdf_file_access_mode mode;		///< Access mode
	uint32_t ndatasets;					///< Number of datasets in file
} t_zdf_file;

/**
 * @brief ZDF Dataset
 * 
 */
typedef struct ZDF_Dataset {
	char* name;							///< Dataset name
 	enum zdf_data_type data_type;		///< Data type
	uint32_t ndims;						///< Number of dimensions
	uint64_t count[zdf_max_dims];		///< Dimension for each direction
	void * data;						///< Pointer to data buffer
	uint64_t id;						///< Optional integer ID
	uint64_t offset;					///< File position for dataset header
} t_zdf_dataset;

/**
 * @brief ZDF data chunk
 * 
 */
typedef struct ZDF_Chunk {
	uint64_t count[zdf_max_dims];		///< Dimension of chunk
	uint64_t start[zdf_max_dims];		///< Start position of data chunk
	uint64_t stride[zdf_max_dims];		///< Chunk stride
	void * data;						///< Pointer to chunk data
} t_zdf_chunk;

/**
 * @brief Types of axis
 * 
 */
enum zdf_axis_type { 
	zdf_linear,	///< Linear axis
	zdf_log10,	///< Log10 axis
	zdf_log2	///< Log2 axis
};

/**
 * @brief Grid data axis information
 * 
 */
typedef struct ZDF_GridAxis {
	char* name;					///< Axis name
	enum zdf_axis_type type;	///< Axis type

	double min;					///< Minimum value
	double max;					///< Maximum value

	char* label;				///< Axis label
	char* units;				///< Axis units

} t_zdf_grid_axis;

/**
 * @brief Grid data information
 * 
 */
typedef struct ZDF_GridInfo {
	char* name;						///< Grid name
	uint32_t ndims;					///< Number of grid dimensions {1..zdf_max_dims}
	uint64_t count[zdf_max_dims];	///< Grid dimensions
	char* label;					///< Grid label
	char* units;					///< Grid units
    t_zdf_grid_axis *axis;			///< Grid axis information
} t_zdf_grid_info;

/**
 * @brief Iteration details
 * 
 */
typedef struct ZDF_Iteration{
	char* name;			///< Iteration name
	int32_t n;			///< Iteration number
	double t;			///< Time (phyical units)
	char* time_units;	///< Units for time
} t_zdf_iteration;

/**
 * @brief Particle data information
 * 
 */
typedef struct ZDF_PartInfo {
	char* name;			///< Particle species name
	char* label;		///< Particle species label
	uint64_t np;		///< Number of particles
	uint32_t nquants;	///< Number of quantities stored
	char** quants;		///< Names of quantities
	char** qlabels;		///< Labels for quantities
	char** qunits;		///< Units for quantities
} t_zdf_part_info;

/**
 * @brief Tracks data information
 * 
 */
typedef struct ZDF_TrackInfo {
	char* name;			///< Particle species name
	char* label;		///< Particle species label
	uint32_t ntracks;	///< Number of tracks in dataset
	uint32_t ndump;		///< Frequency at which the file was updated
	uint32_t niter;		///< Number of iterations between track points

	uint32_t nquants;	///< Number of quantities in tracks
	char** quants;		///< Names of quantities
	char** qlabels;		///< Labels for quantities
	char** qunits;		///< Units for quantities
} t_zdf_track_info;

// Low level interface

/**
 * Returns size of ZDF datatype
 * @param  data_type Data type id
 * @return           Returns type size in bytes or 0 for an invalid data type
 */
size_t zdf_sizeof( enum zdf_data_type data_type );

/**
 * Opens ZDF file
 * @param  zdf      ZDF file to open
 * @param  filename Filename of the ZDF file to open, including path
 * @param  mode     Can be one of ZDF_CREATE, ZDF_READ, or ZDF_UPDATE
 * @return          Returns 1 on success, 0 on error
 */
int zdf_open_file( t_zdf_file* zdf, const char* filename, enum zdf_file_access_mode mode );

/**
 * Closes ZDF file
 * @param  zdf ZDF file to close
 * @return     Returns 1 on success, 0 on error
 */
int zdf_close_file( t_zdf_file* zdf );

/**
 * Writes a vector of the specified datatype to file
 * @param  zdf       ZDF file descriptor
 * @param  data      Pointer to data
 * @param  data_type Data type descriptor
 * @param  len       Number of elements in vector
 * @return           Returns number of bytes written on success, 0 on error
 */
size_t zdf_vector_write( t_zdf_file* zdf, const void * data, enum zdf_data_type data_type, size_t len );

/**
 * Adds string element to ZDF file
 * @param  zdf  ZDF File handle
 * @param  name Element name
 * @param  str  String value
 * @return      Returns number of bytes written on success, 0 on error
 */
size_t zdf_add_string( t_zdf_file* zdf, const char* name, const char* str );

/**
 * Adds int32 element to ZDF file
 * @param  zdf    ZDF File handle
 * @param  name   Element name
 * @param  value  int32 value
 * @return        Returns number of bytes written on success, 0 on error
 */
size_t zdf_add_int32( t_zdf_file* zdf, const char* name, const int32_t value );

/**
 * Adds float64 element to ZDF file
 * @param  zdf    ZDF file handle
 * @param  name   Element name
 * @param  value  float64 value
 * @return        Returns number of bytes written on success, 0 on error
 */
size_t zdf_add_double( t_zdf_file* zdf, const char* name, const double value );

/**
 * Adds iteration metadata group to ZDF file
 * @param  zdf  ZDF file handle
 * @param  iter Iteration info
 * @return      Number of bytes written on success, 0 on error
 */
size_t zdf_add_iteration( t_zdf_file* zdf, const t_zdf_iteration* iter );

/**
 * Adds grid information metadata group to file
 * @param  zdf  File handle
 * @param  grid Grid information
 * @return      Number of bytes written on success, 0 on error
 */
size_t zdf_add_grid_info( t_zdf_file* zdf, const t_zdf_grid_info* grid );

/**
 * Adds particle information metadata group to file
 * @param  zdf  File handle
 * @param  part Particle information
 * @return      Number of bytes written on success, 0 on error
 */
size_t zdf_add_part_info( t_zdf_file* zdf, const t_zdf_part_info* part );

/**
 * Adds tracks information metadata group to file
 * @param  zdf    File handle
 * @param  tracks Tracks information
 * @return        Number of bytes written on success, 0 on error
 */
size_t zdf_add_track_info( t_zdf_file* zdf, const t_zdf_track_info* tracks );

/**
 * Adds dataset to file
 * @param  zdf     ZDF file handle
 * @param  dataset Dataset object (includes pointer to data)
 * @return         Number of bytes written on success, 0 on error. The dataset object is also
 *                 modified to include a unique id.
 */
size_t zdf_add_dataset( t_zdf_file* zdf, t_zdf_dataset* dataset );

// Chunked dataset interface

/**
 * Adds a chunked dataset header to file
 * @param  zdf     ZDF file handle
 * @param  dataset Dataset object
 * @return         Number of bytes written on success, 0 on error. The dataset object is also
 *                 modified to include a unique id.
 */
size_t zdf_start_cdset( t_zdf_file* zdf, t_zdf_dataset* dataset );

/**
 * Returns size of chunk header
 * @param  dataset Dataset object
 * @return         Chunk header size in bytes
 */
size_t size_zdf_chunk_header(const t_zdf_dataset* dataset);

/**
 * Adds a chunk header to file
 * @param  zdf     ZDF file handle
 * @param  dataset Dataset object
 * @param  chunk   Chunk object
 * @return         Number of bytes written on success, 0 on error.
 */
size_t zdf_write_chunk_header( t_zdf_file* zdf, t_zdf_dataset* dataset, t_zdf_chunk* chunk );

/**
 * Adds a chunk of dataset data to file
 * @param  zdf     ZDF file handle
 * @param  dataset Dataset object
 * @param  chunk   Chunk object (includes pointer to data). Data is expected to be contiguos in memory
 * @return         Number of bytes written on success, 0 on error.
 */
size_t zdf_write_cdset( t_zdf_file* zdf, t_zdf_dataset* dataset, t_zdf_chunk* chunk );

/**
 * Add an end marker for a chunked dataset to file
 * @param  zdf     ZDF file handle
 * @param  dataset Dataset object
 * @return         Number of bytes written on success, 0 on error.
 */
size_t zdf_end_cdset( t_zdf_file* zdf, t_zdf_dataset* dataset );

/**
 * Open a dataset and read the header. The dataset is selected by name.
 * @param  zdf     ZDF file handle
 * @param  dataset Dataset object
 * @return         Returns 1 if successful, 0 on error.
 */
size_t zdf_open_dataset( t_zdf_file* zdf, t_zdf_dataset* dataset );

/**
 * Extends dataset dimensions. After this call the file pointer is positioned at the
 * end of the file.
 *
 * @param  zdf       ZDF file handle
 * @param  dataset   Datasete object
 * @param  new_count New dimensions for the dataset, must be >= than the previous
 *                   dimensions.
 * @return           Returns 1 if successful, 0 otherwise.
 */
int zdf_extend_dataset( t_zdf_file* zdf, t_zdf_dataset* dataset, uint64_t* new_count );


// High level interface

/**
 * Opens ZDF file and adds TYPE, GRID, and ITERATION metadata
 * @param  zdf       File handle
 * @param  info      Grid information
 * @param  iteration Iteration information
 * @param  path      File path
 * @return           Returns 1 on success, 0 on error
 */
int zdf_open_grid_file( t_zdf_file *file, const t_zdf_grid_info *info,
	const t_zdf_iteration *iteration, char const path[] );

/**
 * Saves a ZDF grid file
 * @param  data       Pointer to grid data
 * @param  data_type  ZDF Data type of grid data
 * @param  info       Grid information
 * @param  iteration  Iteration information
 * @param  path       File path
 * @return            Returns 1 on success, 0 on error
 */
int zdf_save_grid( const void* data, enum zdf_data_type data_type, const t_zdf_grid_info *info,
	const t_zdf_iteration *iteration, char const path[] );

/**
 * Opens ZDF file and adds TYPE, PARTICLES, and ITERATION metadata
 * @param  zdf        File handle
 * @param  info       Particles information
 * @param  iteration  Iteration information
 * @param  path       File path
 * @return            1 on sucess, 0 on error
 */
int zdf_open_part_file( t_zdf_file *file, t_zdf_part_info *info,
	const t_zdf_iteration *iteration, char const path[] );

/**
 * Adds individual particle quantity to file
 * @param  zdf  File handle
 * @param  name Quantity name
 * @param  data Quantity data (float32)
 * @param  np   Number of particles
 * @return      Number of bytes written on success, 0 on error
 */
int zdf_add_quant_part_file( t_zdf_file *zdf, const char *name, const float* data,
	const uint64_t np );


#endif


