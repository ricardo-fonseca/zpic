#ifndef __FILTER__
#define __FILTER__

enum filter_type { FILTER_NONE, FILTER_GAUSS, FILTER_SHARP };


typedef struct {

	enum filter_type type;

	int nk;
	float *Sk;

} t_filter;

void filter_new( t_filter *filter, int nk );
void filter_delete( t_filter *filter );

int filter_set( t_filter* const filter, enum filter_type const type, float const ck );

#endif
