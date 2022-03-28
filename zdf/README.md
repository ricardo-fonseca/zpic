# ZDF data format

ZDF is a simple, self-describing, lightweight file format for writing scientific data. It was developed in the framework of the ZPIC educational code suite as an alternative to existing standards such as HDF5.

It has no external dependencies, and it consists of a minimum number (currently 1) of C files that can be easily included in any project. It is also portable, ensuring that it can be written/read on any system regardless of endianness, number of bits, etc.

All the relevant metadata can be included alongside with the data in the same file.

## Linking the ZDF library

The library is provided in source code format. To use the ZDF library in your program you only need 2 files: zdf.c and zdf.h. Just compile zdf.c alongside the rest of your code.

## Internal Format

* Files being with a magic number 0x5A, 0x44, 0x46, 0x31 (ASCII code of ‘ZDF1’)
* Data is written using __little-endian__ format, and always in units of 32 bits. When compiled in big-endian systems, the library will convert the data to little-endian before writing.
* Data is organized in the file as a sequence of records. Each record includes a type id + version tag and total size, that allows records to be easily skipped.

### Low-level data formats

The library provides routines to portably write the following scalar types. As mentioned before these are written to file in little-endian format.

* int32 (signed 32 bit integer)
* uint32 (unsigned 32 bit integer)
* int64 (signed 64 bit integer)
* uint64 (unsigned 64 bit integer)
* float64 (double)

The library also provides routines to write vectors of the following types:

* float32 (float)
* float64 (double)

Finally, the library provides routines to write strings (ASCII). The library first writes the string length (as an unsigned 32bit integer), then writes the string data padding with 0 values so that the final size is a multiple of 32 bits (4 bytes).

## ZDF records

The structure for every ZDF record is as follows:

* id_version (uint32)
  * This value will both identify the type of record, and the version number of the saving routine. The version number will be used for compatibility with older files: if a new feature has been introduced in version N, but the file was written before, the reading routine will not attempt to read in the new features.
* name (string)
  * This is a short name describing the record, encoded as described below. It can be set to an empty string.
* length (uint64)
  * The length in bytes of the record. If the program reading the data does not recognize the record type (or version) it can use this data to skip the record.

This can be followed by any combination of data.

The library currently supports the following record types:

* __ZDF\_INT32__
* __ZDF\_DOUBLE__
* __ZDF\_STRING__
* __ZDF\_DATASET__ 
* __ZDF\_ITERATION__ 
* __ZDF\_GRID_INFO__ 
* __ZDF\_PART_INFO__ 


## High-Level interfaces for writing ZDF files

For simplicity, the ZDF data format defines two high-level file formats for storing grid and particle data.

### Grid data

ZDF grid data files consist of a sequence of 4 records as folows:

* TYPE: string "grid"
* GRID: grid information record
	* Number of dimensions / number of cell in each dimension
	* Data label / units
	* Axis information (optional), including type, range, label and units
* ITERATION: iteration information record
	* Iteration (integer), simulation time (double) and units.
* DATASET: grid data

To write ZDF grid data files the library provides the `zdf_save_grid` routine. Currently this routine only supports single precision data:

```C
int zdf_save_grid( const float* data, const t_zdf_grid_info *info, 
	const t_zdf_iteration *iteration, char const path[] )
```

It takes the following parameters:

* data - pointer to n-D array to be saved
* info - pointer to `t_zdf_grid_info` structure describing the grid
* iter - pointer to `t_zdf_iteration` structure describing the current iteration
* path - directory where to write the ZDF file. The file name is constructed from the information supplied.


### Particle data

ZDF particle data files consist of a sequence of 3 + _Nq_ (with _Nq_ being the number of quantities saved) records as folows:

* TYPE: string "particles"
* PARTICLES: particle information record
	* Name of particle group
	* Number of particles in the group
	* Number of quantities (_Nq_) recorded
	* Name / units of each quantity
* ITERATION: iteration information record
	* Iteration (integer), simulation time (double) and units.
* N × DATASET: 1 dataset for each particle quantity

To write ZDF particle data files the library provides the `zdf_part_file_open`, `zdf_part_file_add_quant` and `zdf_close_file` routines. Writing a ZDF particle data file requires three steps:

1. Create the file with `zdf_part_file_open`
2. For each quantity that we want to add, write it to the file using `zdf_part_file_add_quant`
3. Close the file with `zdf_close_file`

```C
// Open ZDF particle file

int zdf_part_file_open( t_zdf_file *file, t_zdf_part_info *info, 
	const t_zdf_iteration *iteration, char const path[] );
```
* __file__ - ZDF file handle
* __info__ - pointer to `t_zdf_part_info` structure describing the particle file
	* This includes a list of the quantities to be saved
* __iter__ - pointer to `t_zdf_iteration` structure describing the current iteration
* __path__ - directory where to write the ZDF file. The file name is constructed from the information supplied.


```C
// Add quantity to ZDF particle file

int zdf_part_file_add_quant( t_zdf_file *zdf, const char *name, 
	const float*data, const unsigned np );
```

* __file__ - ZDF file handle
* __name__ - name of the quantity being saved
* __data__ - pointer to particle data being saved
* __np__ - number of points in data (must match the number of particles specified in info)


```C
// Close ZDF file

int zdf_close_file( t_zdf_file* zdf );

```
* __file__ - ZDF file handle


## Reading ZDF Files

