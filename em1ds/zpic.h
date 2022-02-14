/*
 *  zpic.h
 *  zpic
 *
 *  Created by Ricardo Fonseca on 12/8/10.
 *  Copyright 2010 Centro de Física dos Plasmas. All rights reserved.
 *
 */

#ifndef __ZPIC__
#define __ZPIC__

/**
 * @brief Three component vector
 * 
 */
typedef struct Float3 {
	float x;	///< x vector component
	float y;	///< y vector component
	float z;	///< z vector component
} float3;

/* ANSI C does not define math constants */

#ifndef M_PI
#define M_PI        3.14159265358979323846264338327950288   ///< pi
#endif

#ifndef M_PI_2
#define M_PI_2      1.57079632679489661923132169163975144   ///< pi/2
#endif

#ifndef M_PI_4
#define M_PI_4      0.785398163397448309615660845819875721  ///< pi/4
#endif


#endif
