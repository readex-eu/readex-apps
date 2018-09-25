/* ----------------------------------------------------------------------
   miniMD is a simple, parallel molecular dynamics (MD) code.   miniMD is
   an MD microapplication in the Mantevo project at Sandia National
   Laboratories ( http://www.mantevo.org ). The primary
   authors of miniMD are Steve Plimpton (sjplimp@sandia.gov) , Paul Crozier
   (pscrozi@sandia.gov) and Christian Trott (crtrott@sandia.gov).

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This library is free software; you
   can redistribute it and/or modify it under the terms of the GNU Lesser
   General Public License as published by the Free Software Foundation;
   either version 3 of the License, or (at your option) any later
   version.

   This library is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this software; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA.  See also: http://www.gnu.org/licenses/lgpl.txt .

   For questions, contact Paul S. Crozier (pscrozi@sandia.gov) or
   Christian Trott (crtrott@sandia.gov).

   Please read the accompanying README and LICENSE files.
---------------------------------------------------------------------- */
//#define PRINTDEBUG(a) a
#define PRINTDEBUG(a)
#include "stdio.h"
#include "integrate.h"
#include "openmp.h"
#include "math.h"

#if defined (USE_SCOREP) || defined (USE_SCOREP_MANUAL)
#include "SCOREP_User.h"
#endif

Integrate::Integrate() {}
Integrate::~Integrate() {}

void Integrate::setup()
{
  dtforce = 0.5 * dt;
}

void Integrate::initialIntegrate()
{
  #pragma omp for schedule(static,CHUNKSIZE)

  for(MMD_int i = 0; i < nlocal; i++) {
    v[i * 3 + 0] += dtforce * f[i * 3 + 0];
    v[i * 3 + 1] += dtforce * f[i * 3 + 1];
    v[i * 3 + 2] += dtforce * f[i * 3 + 2];
    x[i * 3 + 0] += dt * v[i * 3 + 0];
    x[i * 3 + 1] += dt * v[i * 3 + 1];
    x[i * 3 + 2] += dt * v[i * 3 + 2];
  }
}

void Integrate::finalIntegrate()
{
  #pragma omp for schedule(static,CHUNKSIZE)

  for(MMD_int i = 0; i < nlocal; i++) {
    v[i * 3 + 0] += dtforce * f[i * 3 + 0];
    v[i * 3 + 1] += dtforce * f[i * 3 + 1];
    v[i * 3 + 2] += dtforce * f[i * 3 + 2];
  }

}

void Integrate::run(Atom &atom, Force* force, Neighbor &neighbor,
                    Comm &comm, Thermo &thermo, Timer &timer)
{
  int i, n;

  comm.timer = &timer;
  timer.array[TIME_TEST] = 0.0;

  int check_safeexchange = comm.check_safeexchange;

  mass = atom.mass;
  dtforce = dtforce / mass;
  //Use OpenMP threads only within the following loop containing the main loop.
  //Do not use OpenMP for setup and postprocessing.
  #pragma omp parallel private(i,n)
  {

#if defined (USE_SCOREP) || defined (USE_SCOREP_MANUAL)
    SCOREP_USER_REGION_DEFINE(R1)
#endif
    for(n = 0; n < ntimes; n++) {
#if defined (USE_SCOREP) || defined (USE_SCOREP_MANUAL)
      SCOREP_USER_OA_PHASE_BEGIN(R1, "INTEGRATE_RUN_LOOP", SCOREP_USER_REGION_TYPE_COMMON)
#endif

      #pragma omp barrier

      x = &atom.x[0][0];
      v = &atom.v[0][0];
      f = &atom.f[0][0];
      xold = &atom.xold[0][0];
      nlocal = atom.nlocal;

      initialIntegrate();
      #pragma omp barrier

      #pragma omp master
      timer.stamp();

      if((n + 1) % neighbor.every) {

        #pragma omp barrier
        comm.communicate(atom);
        #pragma omp master
        timer.stamp(TIME_COMM);
        #pragma omp barrier

      } else {
        //these routines are not yet ported to OpenMP
        {
          if(check_safeexchange) {
            #pragma omp master
            {
              double d_max = 0;

              for(i = 0; i < atom.nlocal; i++) {
                double dx = (x[3 * i + 0] - xold[3 * i + 0]);

                if(dx > atom.box.xprd) dx -= atom.box.xprd;

                if(dx < -atom.box.xprd) dx += atom.box.xprd;

                double dy = (x[3 * i + 1] - xold[3 * i + 1]);

                if(dy > atom.box.yprd) dy -= atom.box.yprd;

                if(dy < -atom.box.yprd) dy += atom.box.yprd;

                double dz = (x[3 * i + 2] - xold[3 * i + 2]);

                if(dz > atom.box.zprd) dz -= atom.box.zprd;

                if(dz < -atom.box.zprd) dz += atom.box.zprd;

                double d = dx * dx + dy * dy + dz * dz;

                if(d > d_max) d_max = d;
              }

              d_max = sqrt(d_max);

              if((d_max > atom.box.xhi - atom.box.xlo) || (d_max > atom.box.yhi - atom.box.ylo) || (d_max > atom.box.zhi - atom.box.zlo))
                printf("Warning: Atoms move further than your subdomain size, which will eventually cause lost atoms.\n"
                "Increase reneighboring frequency or choose a different processor grid\n"
                "Maximum move distance: %lf; Subdomain dimensions: %lf %lf %lf\n",
                d_max, atom.box.xhi - atom.box.xlo, atom.box.yhi - atom.box.ylo, atom.box.zhi - atom.box.zlo);

            }

          }


          //int tid = omp_get_thread_num();
          //printf("Check B: %i %i %i\n",comm.me,tid,n);
          #pragma omp master
          timer.stamp_extra_start();
          comm.exchange(atom);
          comm.borders(atom);
          #pragma omp master
          {
            timer.stamp_extra_stop(TIME_TEST);
            timer.stamp(TIME_COMM);
          }

          if(check_safeexchange)
            for(int i = 0; i < 3 * atom.nlocal; i++) atom.xold[i] = atom.x[i];
        }

        #pragma omp barrier

        neighbor.build(atom);
        #pragma omp barrier

        #pragma omp master
        timer.stamp(TIME_NEIGH);
      }

      force->evflag = (n + 1) % thermo.nstat == 0;
      force->compute(atom, neighbor, comm, comm.me);


      #pragma omp master
      timer.stamp(TIME_FORCE);

      if(neighbor.halfneigh && neighbor.ghost_newton) {
        comm.reverse_communicate(atom);

        #pragma omp master
        timer.stamp(TIME_COMM);
      }

      v = &atom.v[0][0];
      f = &atom.f[0][0];
      nlocal = atom.nlocal;

      #pragma omp barrier

      finalIntegrate();

      #pragma omp barrier

      if(thermo.nstat) thermo.compute(n + 1, atom, neighbor, force, timer, comm);

#if defined (USE_SCOREP) || defined (USE_SCOREP_MANUAL)
    SCOREP_USER_OA_PHASE_END(R1)
#endif
    } //end for-loop
  } //end OpenMP parallel
}
