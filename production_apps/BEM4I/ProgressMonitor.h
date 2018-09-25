/*!
 * @file    ProgresMonitor.h
 * @author  Michal Merta
 * @author  Jan Zapletal
 * @date    July 24, 2014
 * @brief   Monitoring and printing progress of current operation
 *
 */

#ifndef PROGRESSMONITOR_H
#define PROGRESSMONITOR_H

#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <chrono>
#include "omp.h"

namespace bem4i {

class ProgressMonitor {
  typedef std::chrono::high_resolution_clock clock;
  typedef std::chrono::milliseconds milliseconds;

public:

  static void init(
    const std::string & msg = "",
    unsigned long nSteps = 1,
    bool time = true
    ) {
    start = clock::now( );
    stepsTotal = nSteps;
    stepsDone = 0;
    percentDone = 0.0;
    currProgressBar = -1;
    keepTime = time;
    std::cout << msg << ", " << stepsTotal << std::endl;
    if ( stepsTotal <= 0 ) stepsTotal = 1;
  }

  static double step(
    unsigned long nSteps = 1
    ) {

    double ret = 0.0;
#pragma omp critical
    {
      //#pragma omp atomic update
      stepsDone += nSteps;
      //#pragma omp master
      //    {
      double newPercentDone = 100.0 * stepsDone / stepsTotal;

      if ( stepsDone >= stepsTotal ) {
        newPercentDone = 100.0;
      }

      if ( floor( newPercentDone ) > floor( percentDone ) ) {
        percentDone = newPercentDone;
        currProgressBar = std::floor( (double) percentDone / 5.0 );

        if ( stepsTotal > 1 ) {
          std::cout << "\r  " << progressBar[ currProgressBar ] << " " <<
            getPercentDone( );
          std::cout.flush( );
        }

        if ( percentDone == 100.0 ) {
          end = clock::now( );
          milliseconds ms =
            std::chrono::duration_cast<milliseconds>( end - start );
          ret = ms.count( ) / 1000.0;
          if ( keepTime ) {
            std::cout << "  Done in " << timeToString( ) << " seconds.";
          }
          std::cout << std::endl;
        }
      }
      //    } // end omp master
    } // end omp critical

    return ret;
  }

  static std::string getPercentDone( ) {
    return percentDoneToString( );
  }

private:

  static std::string percentDoneToString( ) {
    std::stringstream str;
    //str << std::fixed;
    //str << std::setprecision( 2 );
    str << std::setfill( ' ' );
    //str << std::setw( 6 );
    str << std::setw( 3 );
    str << floor( percentDone ) << " %";
    return str.str( );
  }

  static std::string timeToString( ) {
    milliseconds ms = std::chrono::duration_cast<milliseconds>( end - start );
    std::stringstream str;
    str << std::fixed;
    str << std::setprecision( 2 );
    str << ms.count( ) / 1000.0;
    return str.str( );
  }

  static unsigned long stepsTotal;
  static unsigned long stepsDone;
  static double percentDone;
  static int currProgressBar;
  static bool keepTime;
  static clock::time_point start;
  static clock::time_point end;

  static const std::vector< std::string > progressBar;

};

typedef std::chrono::high_resolution_clock clock;

// initialization of static members
unsigned long ProgressMonitor::stepsTotal = 0;
unsigned long ProgressMonitor::stepsDone = 0;
double ProgressMonitor::percentDone = 0.0;
int ProgressMonitor::currProgressBar = -1;
bool ProgressMonitor::keepTime = true;
clock::time_point ProgressMonitor::start = clock::now( );
clock::time_point ProgressMonitor::end = clock::now( );


const std::vector< std::string > ProgressMonitor::progressBar{
  "[                    ]",
  "[*                   ]",
  "[**                  ]",
  "[***                 ]",
  "[****                ]",
  "[*****               ]",
  "[******              ]",
  "[*******             ]",
  "[********            ]",
  "[*********           ]",
  "[**********          ]",
  "[***********         ]",
  "[************        ]",
  "[*************       ]",
  "[**************      ]",
  "[***************     ]",
  "[****************    ]",
  "[*****************   ]",
  "[******************  ]",
  "[******************* ]",
  "[********************]", };

}

#endif /* PROGRESSMONITOR_H */
