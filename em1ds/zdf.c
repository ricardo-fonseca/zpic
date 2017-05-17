/*
 *  zdf.c
 *  zpic
 *
 *  Created by Ricardo Fonseca on 4/1/17.
 *  Copyright 2010 Centro de Física dos Plasmas. All rights reserved.
 *
 */

#include "zdf.h"


#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <rpc/rpc.h>
#include <string.h>
#include <errno.h>

/*
void save_grid_visxd( const float* data, const int ndim, const int nx[], 
					  const float xmin[], const float xmax[],
					  const char *label, char* xlabel[],
					  const int n, const float t, char const path[])
*/

#define max_string_length 128

const unsigned int zdf_magic = 0x5A444630; // This corresponds to the ascii code of 'Z','D','F','0'

enum zdf_data_type{ zdf_null,
	                zdf_int8,  zdf_uint8,  zdf_int16, zdf_uint16, 
	                zdf_int32, zdf_uint32, zdf_int64, zdf_uint64,
	                zdf_float32, zdf_float64 };


const unsigned size_xdr_bool    = 4;
const unsigned size_xdr_int     = 4;
const unsigned size_xdr_u_int   = 4;
const unsigned size_xdr_u_int64 = 8;
const unsigned size_xdr_double  = 8;
const unsigned size_xdr_float   = 4;

#define ZDF_INT_ID       0x00010000
#define ZDF_DOUBLE_ID    0x00020000
#define ZDF_STRING_ID    0x00030000

#define ZDF_DATASET_ID   0x00100000

#define ZDF_ITERATION_ID 0x00200000
#define ZDF_GRID_INFO_ID 0x00210000
#define ZDF_PART_INFO_ID 0x00220000

/* -----------------------------------------------------------------------------------------------
  recursively create path if required
-------------------------------------------------------------------------------------------------- */

int create_path( const char path[] )
{
	char uppath[256], *p;
	int ierr = 0;
	
	if (mkdir(path,S_IRWXU | (S_IRGRP | S_IXGRP ) | (S_IROTH | S_IXOTH) )) {
		switch (errno) {
			case ENOENT : // A component of the path does not exist 
				
				// get upper path
				strncpy(uppath, path, 256);
				p = uppath + strlen(uppath);
				while(*p!='/') p--;
				*p=0;
				
				//recursively build the path
				if ( !create_path( uppath ) ) ierr = create_path( path );
				break;
				
			case EEXIST : /* if directory already exists ignore the error */
				ierr = 0;
				break;
			default: ierr = errno;
		}
	}
	return ierr;
}

/* -----------------------------------------------------------------------------------------------
  Open / Close ZDF file
-------------------------------------------------------------------------------------------------- */


int zdf_close_file( t_zdf_file* zdf ) {

	xdr_destroy( &zdf->xdrs );
	if ( fclose( zdf->fp ) ) {
		perror("(*error*) Unable to close ZDF file");
		return(-1);
	}
	zdf -> fp = NULL;
	return(0);
}

int zdf_open_file( t_zdf_file* zdf, char* filename, enum zdf_file_access_mode mode ){

//int zdf_open_file( t_zdf_file* zdf, const char *restrict filename, enum zdf_file_access_mode mode ) {

	zdf -> mode = mode;

	switch( mode ) {
		case ZDF_WRITE :
			// Open file for writing
			if (!(zdf->fp = fopen( filename, "w"))) {
				perror("(*error*) Unable to open ZDF file for writing");
				return(-1);
			}
			
			// Create XDR stream
			xdrstdio_create( &zdf->xdrs, zdf->fp, XDR_ENCODE );

			// Write magic number
			if (!xdr_u_long( &zdf->xdrs, (unsigned int *) &zdf_magic )) {
				fprintf(stderr, "(*error*) Unable to write magic number to ZDF file.");
				zdf_close_file( zdf );
				return(-1);
			}
			break;

		case ZDF_READ :
			// Open file for reading
			if (!(zdf->fp = fopen( filename, "r"))) {
				perror("(*error*) Unable to open ZDF file for reading");
				return(-1);
			}
			
			// Create XDR stream
			xdrstdio_create( &zdf->xdrs, zdf->fp, XDR_DECODE );

			unsigned int test_magic;

			// Read magic number
			if (!xdr_u_long( &zdf->xdrs, (unsigned int *) &test_magic )) {
				fprintf(stderr, "(*error*) Unable to read magic number from ZDF file.");
				zdf_close_file( zdf );
				return(-1);
			}

            // Check magic number
            if ( test_magic != zdf_magic ) {
				fprintf(stderr, "(*error*) Invalid magic number, file is not a proper ZDF file.");
				zdf_close_file( zdf );
				return(-1);
			}
			break;

		default:
			fprintf(stderr, "(*error*) zdf_open_file: unsupported mode\n");
			return(-1);
	}

	return(0);
}


/* -----------------------------------------------------------------------------------------------
  xdr_simple_string
-------------------------------------------------------------------------------------------------- */

bool_t xdr_simple_string( XDR* xdrs, char ** str ){

	unsigned len;
	char *p = *str;

	len = ( str ) ? strlen( *str ) : 0;
	if ( !xdr_u_int( xdrs, &len ) ) return(FALSE);

	if ( len > 0 ) {
		if ( !xdr_bytes( xdrs, &p, &len, max_string_length ) ) return(FALSE);
	}

	*str = p;

	return(TRUE);
}

uint64_t size_xdr_simple_string( const char *s )
{
	unsigned len = ( s ) ? strlen( s ) : 0;
	return size_xdr_u_int + ((len > 0) ? ( size_xdr_u_int + RNDUP(len) ) : 0);
}

/* -----------------------------------------------------------------------------------------------
  zdf records
-------------------------------------------------------------------------------------------------- */

typedef struct {
	uint32_t id_version;
	char*    name;
	uint64_t length;
} t_zdf_record;


bool_t xdr_zdf_record( XDR* xdrs, t_zdf_record* rec ){

    if ( !xdr_u_int( xdrs, &rec -> id_version ) ) return(FALSE);
 	if ( !xdr_simple_string( xdrs, &rec -> name ) ) return(FALSE);
	if ( !xdr_u_int64_t( xdrs, &rec -> length ) ) return(FALSE);

 	return(TRUE);
}

/* -----------------------------------------------------------------------------------------------
  zdf basic data tags
-------------------------------------------------------------------------------------------------- */
int zdf_add_string( t_zdf_file* zdf, char* name, char* str ){
	
	t_zdf_record rec;
	rec.id_version = ZDF_STRING_ID;
	rec.name = name;
	rec.length = size_xdr_simple_string( str );

    if ( !xdr_zdf_record( &zdf -> xdrs, &rec) ) return(-1);
    if ( !xdr_simple_string( &zdf -> xdrs, &str) ) return(-1);

    return(0);
}

int zdf_add_int( t_zdf_file* zdf, char* name, int* value ){
	
	t_zdf_record rec;
	rec.id_version = ZDF_INT_ID;
	rec.name = name;
	rec.length = size_xdr_int;

    if ( !xdr_zdf_record( &zdf -> xdrs, &rec) ) return(-1);
    if ( !xdr_int( &zdf -> xdrs, value) ) return(-1);

    return(0);
}

int zdf_add_double( t_zdf_file* zdf, char* name, double* value )
{
	t_zdf_record rec;
	rec.id_version = ZDF_DOUBLE_ID;
	rec.name = name;
	rec.length = size_xdr_double;

    if ( !xdr_zdf_record( &zdf -> xdrs, &rec) ) return(-1);
    if ( !xdr_double( &zdf -> xdrs, value) ) return(-1);

    return(0);
}

/* -----------------------------------------------------------------------------------------------
  zdf compound metadata tags
-------------------------------------------------------------------------------------------------- */

int zdf_add_iteration( t_zdf_file* zdf, char* name, t_zdf_iteration* iter ){

	t_zdf_record rec;
	rec.id_version = ZDF_ITERATION_ID;
	rec.name = name;
	rec.length = size_xdr_u_int + 
	             size_xdr_double +
	             size_xdr_simple_string( iter -> time_units );

    if ( !xdr_zdf_record( &zdf -> xdrs, &rec) ) return(-1);

 	if ( !xdr_u_int( &zdf -> xdrs, (unsigned int *) &iter ->n ) ) return(-1);
	if ( !xdr_double( &zdf -> xdrs, &iter->t ) ) return(-1);
 	if ( !xdr_simple_string( &zdf -> xdrs, &iter->time_units ) ) return(-1);

 	return(0);
}


uint64_t size_xdr_zdf_grid_info(const t_zdf_grid_info* grid) {
	uint64_t size;

	size = size_xdr_u_int + grid -> ndims * size_xdr_u_int64 +
	        size_xdr_simple_string(grid->label) +  size_xdr_simple_string(grid->units);

   	// Includes axis information
   	size += size_xdr_bool;
    if ( grid -> axis ) {

    	int i;
    	for(i=0; i<grid -> ndims; i++) 
    		size += size_xdr_int + 
    	            2 * size_xdr_double + 
    	            size_xdr_simple_string( grid -> axis[i].label ) +
    	            size_xdr_simple_string( grid -> axis[i].units );
    }

	return size;
}

int zdf_add_grid_info( t_zdf_file* zdf, char* name, t_zdf_grid_info* grid ){

	t_zdf_record rec;
	rec.id_version = ZDF_GRID_INFO_ID;
	rec.name = name;
	rec.length = size_xdr_zdf_grid_info( grid );

    if ( !xdr_zdf_record( &zdf -> xdrs, &rec) ) return(-1);

 	if ( !xdr_u_int( &zdf -> xdrs, (unsigned int *) &grid -> ndims ) ) return(-1);

    unsigned int i;
    for( i=0; i < grid -> ndims; i++) {
	 	if ( !xdr_u_int64_t( &zdf -> xdrs, &grid -> nx[i] ) ) return(-1);
    }
 	
 	if ( !xdr_simple_string( &zdf ->xdrs, &grid->label ) ) return(-1);
 	if ( !xdr_simple_string( &zdf ->xdrs, &grid->units ) ) return(-1);

    bool_t has_axis = ( grid -> axis != NULL);
 	if ( !xdr_bool( &zdf -> xdrs, &has_axis ) ) return(-1);

 	if ( has_axis ) {

	    for( i=0; i < grid -> ndims; i++) {

		    if ( !xdr_int( &zdf -> xdrs, (int *) &grid -> axis[i].type ) ) return(-1);
		    if ( !xdr_double( &zdf -> xdrs, &grid -> axis[i].min ) ) return(-1);
		    if ( !xdr_double( &zdf -> xdrs, &grid -> axis[i].max ) ) return(-1);

		    if ( !xdr_simple_string( &zdf -> xdrs, &grid -> axis[i].label ) ) return(-1);
		    if ( !xdr_simple_string( &zdf -> xdrs, &grid -> axis[i].units ) ) return(-1);
	    }
 	}

 	return(0);
}

uint64_t size_xdr_zdf_part_info(const t_zdf_part_info* part) {
	uint64_t size = size_xdr_simple_string(part->name) +	// name
			size_xdr_u_int; 								// nquants

	unsigned int i;
	for( i = 0; i < part -> nquants; i++) {
		size += size_xdr_simple_string( part -> quants[i] ) +
				size_xdr_simple_string( part -> units[i] );
	}

	size += size_xdr_u_int64; 					// np

	return size;
}

int zdf_add_part_info( t_zdf_file* zdf, char* name, t_zdf_part_info* part ){

	t_zdf_record rec = {
		.id_version = ZDF_PART_INFO_ID,
		.name = name,
		.length = size_xdr_zdf_part_info( part )
	};

    if ( !xdr_zdf_record( &zdf -> xdrs, &rec) ) return(-1);
	if ( !xdr_simple_string( &zdf ->xdrs, &part->name ) ) return(-1);
 	if ( !xdr_u_int( &zdf -> xdrs, (unsigned int *) &part -> nquants ) ) return(-1);

    unsigned int i;
    for( i=0; i < part -> nquants; i++) {
		if ( !xdr_simple_string( &zdf ->xdrs, &part->quants[i] ) ) return(-1);
    }
    for( i=0; i < part -> nquants; i++) {
		if ( !xdr_simple_string( &zdf ->xdrs, &part->units[i] ) ) return(-1);
    }
 	
 	if ( !xdr_u_int64_t( &zdf -> xdrs, &part -> np ) ) return(-1);

 	return(0);
}

/* -----------------------------------------------------------------------------------------------
  zdf dataset
-------------------------------------------------------------------------------------------------- */

typedef struct {
	enum zdf_data_type data_type;
	uint32_t ndims;
	uint64_t nx[zdf_max_dims];

	char* data;
} t_zdf_dataset;


uint64_t zdf_datatype_size( enum zdf_data_type data_type ) {
    uint64_t size;

    switch ( data_type ) {
    	case zdf_float32: size = 4; break;
    	case zdf_float64: size = 8; break;
    	default : size = 0;
    }

    return(size);
}

uint64_t size_xdr_zdf_dataset(const t_zdf_dataset* dataset) {
	
	unsigned int i;
    uint64_t data_size;

    data_size = zdf_datatype_size( dataset -> data_type );
    for( i = 0; i < dataset -> ndims; i++ ) {
    	data_size *= dataset -> nx[i];
    }

	return size_xdr_int + size_xdr_u_int + 
	       dataset -> ndims * size_xdr_u_int64 +
	       data_size ;
}

int zdf_add_dataset( t_zdf_file* zdf, char* name, t_zdf_dataset* dataset ){

	t_zdf_record rec;
	rec.id_version = ZDF_DATASET_ID;
	rec.name       = name;
	rec.length     = size_xdr_zdf_dataset( dataset );

    if ( !xdr_zdf_record( &zdf -> xdrs, &rec) ) return(-1);

 	if ( !xdr_int( &zdf -> xdrs, (int *) &dataset -> data_type ) ) return(-1);
	if ( !xdr_u_int( &zdf -> xdrs, &dataset -> ndims ) ) return(-1);

    unsigned int i;
 	unsigned int count;
    for( i=0, count = 1; i < dataset -> ndims; i++) {
    	count *= dataset -> nx[i];
	 	if ( !xdr_u_int64_t( &zdf -> xdrs, &dataset -> nx[i] ) ) return(-1);
    }

 	switch ( dataset -> data_type ) {
    	case zdf_float32: 
    	    if ( !xdr_vector( &zdf -> xdrs, dataset -> data, count, sizeof(float), 
    	    	            (xdrproc_t) xdr_float ) ) return(-1);
    	    break;
    	case zdf_float64:
    	    if ( !xdr_vector( &zdf -> xdrs, dataset -> data, count, sizeof(double), 
    	    				(xdrproc_t) xdr_double )) return(-1);
    	    break;
    	default: 
    		fprintf(stderr,"(*error*) zdf_add_dataset: Unsupported datatype.");
    		return(-1);
    		break;
 	}
 	
 	return(0);
}


/* -----------------------------------------------------------------------------------------------
  zdf high level interface
-------------------------------------------------------------------------------------------------- */
/*
int zdf_save_grid( const float* data, const int ndim, const int nx[], 
				    const float xmin[], const float xmax[],
				    const char* label, char* xlabel[],
				    const int n, const float t, char const path[])
*/

int zdf_save_grid( const float* data, const t_zdf_grid_info *_info, 
	const t_zdf_iteration *_iteration, char const path[] )
{

	int i;
	char filename[1024];

    // Set iteration info
    t_zdf_iteration iteration = { .n = _iteration->n, 
    	                          .t = _iteration->t, 
    	                          .time_units = _iteration->time_units};

    // Set axis info
    t_zdf_grid_axis axis[zdf_max_dims];

    if ( _info -> axis ) {
	    // Copy axis info
	    for( i = 0; i < _info->ndims; i ++) {
	    	axis[i] = (t_zdf_grid_axis) 
	    	          { .type = zdf_linear,
	    	            .min = _info -> axis[i].min,
	    	            .max = _info -> axis[i].max,
	    	            .label = _info -> axis[i].label ,
	    	            .units = _info -> axis[i].units };
	    	}

    } else {
    	// set default axis info
	    for( i = 0; i < _info->ndims; i ++) {
	    	axis[i] = (t_zdf_grid_axis) 
	    	          { .type = zdf_linear,
	    	            .min = 0,
	    	            .max = _info -> nx[i],
	    	            .label = "",
	    	            .units = "" };
	    	}
    }

    // Set grid info
    t_zdf_grid_info   grid_info = {
    	.ndims = _info->  ndims,
    	.label = _info -> label,
    	.units = _info -> units,
    	.axis  = axis
    };
    for( i = 0; i < _info->ndims; i ++) grid_info.nx[i] = _info->nx[i];

    // Set data
    t_zdf_dataset dataset = {
    	.data_type = zdf_float32,
    	.ndims = _info->ndims,
    	.data = (char *) data 
    };
    for( i = 0; i < _info->ndims; i ++) dataset.nx[i] = _info->nx[i];

	// Ensure that the path is available
	create_path( path );
	
	// Build filename
	sprintf( filename, "%s/%s-%06u.zdf", path, _info->label, _iteration -> n );
    //printf("Saving filename %s\n", filename );
	
    // Create ZDF file
    t_zdf_file zdf;
    if ( zdf_open_file( &zdf, filename, ZDF_WRITE ) ) {
    	fprintf(stderr,"(*error*) Unable to open ZDF file, aborting.");
    	return(-1);
    }

    // Add file type
    zdf_add_string( &zdf, "TYPE", "grid");

    // Add grid info
    zdf_add_grid_info( &zdf, "GRID", &grid_info );

    // Add iteration info
    zdf_add_iteration( &zdf, "ITERATION", &iteration );

    // Add dataset
    zdf_add_dataset( &zdf, "DATA", &dataset );

    // Close ZDF file and return
    return( zdf_close_file( &zdf ) );
}

int zdf_part_file_open( t_zdf_file *zdf, t_zdf_part_info *_info, 
	const t_zdf_iteration *_iteration, char const path[] ){

	char filename[1024];

    t_zdf_iteration iteration = *_iteration;
    t_zdf_part_info info = *_info;

	// Ensure that the path is available
	create_path( path );
	
	// Build filename
	sprintf( filename, "%s/%s-%s-%06u.zdf", path, "particles", _info->name, _iteration -> n );
    //printf("Saving filename %s\n", filename );
	
    // Create ZDF file
    if ( zdf_open_file( zdf, filename, ZDF_WRITE ) ) {
    	fprintf(stderr,"(*error*) Unable to open ZDF file, aborting.");
    	return(-1);
    }

    // Add file type
    zdf_add_string( zdf, "TYPE", "particles");

    // Add particle info
    zdf_add_part_info( zdf, "PARTICLES", &info );

    // Add iteration info
    zdf_add_iteration( zdf, "ITERATION", &iteration );

    return(0);
}

int zdf_part_file_add_quant( t_zdf_file *zdf, const char *name, const float* data, 
	const unsigned np ) {

    t_zdf_dataset dataset = {
    	.data_type = zdf_float32,
    	.ndims = 1,
    	.data = (char *) data
    };

    dataset.nx[0] = np;

    zdf_add_dataset( zdf, (char *) name, &dataset );

    return(0);

}
