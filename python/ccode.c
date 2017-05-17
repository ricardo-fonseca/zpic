#include "ccode.h"
#include <stdio.h>

int prints( t_struct * s) {
	printf("[%d,%g]\n", s->a, s->b);
	return(1);
}