/*
 * This file is part of the Score-P software (http://www.score-p.org)
 *
 * Copyright (c) 2015-2016,
 * Technische Universitaet Muenchen, Germany
 *
 * This software may be modified and distributed under the terms of
 * a BSD-style license.  See the COPYING file in the package base
 * directory for details.
 *
 */

/* *INDENT-OFF* */
//#include <config.h>

#include <stdio.h>

#include <scorep/SCOREP_User.h>
#include <rrl/user_parameters.h>

int variable = 0;

static void
set_value( int value, int* old_value )
{
    old_value = &variable;

    variable = value;
}


static int
get_value( void )
{
    return variable;
}

int
main( int    argc,
      char** argv )
{
    int default_value = 1; /* default value */
    int retVal = 0; /* return value */
    int k, i;

    SCOREP_USER_REGION_DEFINE( mainRegion );

//    SCOREP_USER_REGION_DEFINE( region1 );

//    SCOREP_USER_METRIC_LOCAL( METRIC1 );
//    SCOREP_USER_METRIC_INIT( METRIC1, "METRIC1", "s", SCOREP_USER_METRIC_TYPE_INT64,
//                             SCOREP_USER_METRIC_CONTEXT_CALLPATH );

    for( k = 0; k < 30; k++ )
    {
        SCOREP_USER_OA_PHASE_BEGIN( mainRegion, "mainRegion", SCOREP_USER_REGION_TYPE_COMMON);

        ATP_PARAM_DECLARE("PARAMETER1", ATP_PARAM_TYPE_RANGE, default_value, "SCOREP_REGION", "Domain1");
        ATP_PARAM_GET("PARAMETER1", &retVal, "Domain1");
        printf("PARAMETER1: value = %d\n", retVal );

//        SCOREP_USER_REGION_BEGIN( region1, "region1", SCOREP_USER_REGION_TYPE_COMMON);

//        printf("tuning_test: value = %d\n", get_value() );

//        SCOREP_USER_METRIC_INT64( METRIC1, 1001 );

//        SCOREP_USER_REGION_END(region1);

        SCOREP_USER_OA_PHASE_END( mainRegion );
    }

    return retVal;
}
