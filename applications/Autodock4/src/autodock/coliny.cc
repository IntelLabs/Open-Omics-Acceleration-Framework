/*

 $Id: coliny.cc,v 1.17 2011/05/10 23:21:51 rhuey Exp $

 AutoDock

Copyright (C) 2009 The Scripps Research Institute. All rights reserved.

 AutoDock is a Trade Mark of The Scripps Research Institute.

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

 */

//
// coliny.cc
//
// This file provides the low-level hooks to the Coliny optimizers.
// Coliny optimizers can be accessed within AutoDock using the command
// line:
//
//	coliny <algname> <numtrials>
//
// If the "help" algname is used, this code prints the list of
// solvers that are supported within the current Coliny configuration.
//

#if USING_COLINY

#include "coliny.h"
#include <coliny/coliny.h>
#include <colin/Factory.h>
#include <colin/OptSolverHandle.h>
#include <colin/OptApplications.h>

using namespace std;
using namespace utilib;

//
// The AutoDock 'objective function' used within Coliny
//
//WARNING: obsolete- uses axis-angle but thinks it's quaternion
//         see eval.cc 2011 MP+RH
double ADEvalFn(/* not const */ double *const x, const int n);

#ifdef COLIN_3_0

//
// Global COLIN problem
//
colin::OptProblem<colin::NLP0_problem> coliny_problem;
//
// Global Coliny solver
//
colin::OptSolverHandle* handle = 0;
///coliny::ColinySolver<colin::OptProblem<BasicArray<double>, colin::AppResponse_Utilib>,BasicArray<double> > coliny_solver;

////
//// Initialize the "algname" optimizer over the given domain.  An initial
//// point is generate as the midpoint over the domain.
////
void coliny_init(const char *const algname, const char *const domain, const int num_vars)
{
//
// If 'algname' equals "help", then print the list of supported
// Coliny solvers and return.
//
   if (strcmp(algname, "help") == 0)
   {
      std::vector<std::string> names;
      std::vector<std::string> descriptions;
      std::vector<int> ndx;
      colin::get_solver_info(names, descriptions);
      utilib::order(ndx, names);
      if (names.size() == 0)
         ucout << "    None" << std::endl;
      else
      {
         ucout << std::endl;
         for (size_t i = 0; i < names.size(); i++)
         {
            ucout << "    " << names[ndx[i]] << std::endl;
            ucout << "         " << descriptions[ndx[i]] << std::endl;
         }
      }

      return;
   }
//
// If 'domain' equals "help", then return after print the
// solver options
//
   if (strcmp(domain, "help") == 0)
   {
      handle = colin::OptSolverFactory(algname);
      if (handle) {
	 handle->solver()->help_parameters(cout);
         cout << flush;
         }
      return;
   }
//
// Setup the OptProblem object
//
   coliny_problem.set_application( colin::new_application("",&ADEvalFn) );
   coliny_problem->config_real_domain(num_vars);
   coliny_problem->set_bounds(domain);
   for (unsigned int i = 6; i < coliny_problem->num_real_vars(); i++)
      coliny_problem->set_periodic_real_bound(i);
//
// Initialize the OptSolver object
//
   handle->solver()->set_problem(coliny_problem);
//
// Read in the OptSolver parameters
//
   ifstream ifstr;
   char fname[256];
   sprintf(fname, "%s.in", algname);
   ifstr.open(fname);
   if (ifstr)
      handle->solver()->read_parameter_values(ifstr);
//
// Create a default 'initial point'
//
   BasicArray<double> lower, upper;
   coliny_problem->get_real_bounds(lower, upper);
}


////
//// Perform minimize using a generic Coliny minimization utility
////
//// To turn on "full debugging", set the false flag to true.
////
void coliny_minimize(const int seed, const std::vector<double>& initpt,
                     /* not const */ std::vector<double>& finalpt,
                     /* unused */ const int& neval, /* unused */ const int& niters)
{
   coliny_problem->reset();

   BasicArray<double> initpt_;

#if 1
   initpt_ << initpt;
#else
   initpt_.resize(initpt.size());
   initpt_[0] = -9.151;
   initpt_[1] = 16.175;
   initpt_[2] = 28.005;
   initpt_[3] = 0.866;
   initpt_[4] = 0.496;
   initpt_[5] = 0.071;
   initpt_[6] = -0.847;
   initpt_[7] = 2 * 3.1416 * (-74.83 / 360);
   initpt_[8] = 2 * 3.1416 * (18.36 / 360);
   initpt_[9] = 2 * 3.1416 * (12.26 / 360);
   initpt_[10] = 2 * 3.1416 * (169.72 / 360);
   initpt_[11] = 2 * 3.1416 * (12.54 / 360);
   initpt_[12] = 2 * 3.1416 * (4.41 / 360);
   initpt_[13] = 2 * 3.1416 * (5.80 / 360);
   initpt_[14] = 2 * 3.1416 * (13.62 / 360);
   initpt_[15] = 2 * 3.1416 * (-0.83 / 360);
   initpt_[16] = 2 * 3.1416 * (-0.44 / 360);
#endif

#if 0
// TESTING STUFF
   cerr << "Bound ";
   for (unsigned int i = 0; i < coliny_problem.num_real_params(); i++)
      cerr << coliny_problem.periodic_real_bound(i);
   cerr << endl;
#endif

   handle->minimize(initpt_, seed, false, false);
   colin::ContextMngr().lexical_cast(handle->solver()->best().point,finalpt);
}


#else

//
// Global COLIN problem
//
colin::OptProblem<BasicArray<double>, colin::AppResponse_Utilib > coliny_problem;
//
// Global Coliny solver
//
coliny::ColinySolver<colin::OptProblem<BasicArray<double>, colin::AppResponse_Utilib > , BasicArray<double> > coliny_solver;




////
//// Initialize the "algname" optimizer over the given domain.  An initial
//// point is generate as the midpoint over the domain.
////
void coliny_init(char* algname, char* domain, int)
{
   //
   // If 'algname' equals "help", then return after calling
   // ColinySolver::initialize
   //
   if (strcmp(algname, "help") == 0)
   {
      coliny_solver.initialize(algname);
      return;
   }
   //
   // If 'domain' equals "help", then return after print the
   // solver options
   //
   if (strcmp(domain, "help") == 0)
   {
      coliny_solver.initialize(algname);
      coliny_solver.help_parameters(cout);
      cout << flush;
      return;
   }
   //
   // Setup the OptProblem object
   //
   colin::OptSetup(coliny_problem, &ADEvalFn, domain);
   for (unsigned int i = 6; i < coliny_problem.num_real_params(); i++)
      coliny_problem.set_periodic_real_bound(i);
   //
   // Initialize the OptSolver object
   //
   coliny_solver.initialize(algname);
   //
   // Read in the OptSolver parameters
   //
   ifstream ifstr;
   char fname[256];
   sprintf(fname, "%s.in", algname);
   ifstr.open(fname);
   if (ifstr)
      coliny_solver.read_parameter_values(ifstr);
   //
   // Create a default 'initial point'
   //
   BasicArray<double> lower, upper;
   coliny_problem.get_real_bounds(lower, upper);
}


////
//// Perform minimize using a generic Coliny minimization utility
////
//// To turn on "full debugging", set the false flag to true.
////
void coliny_minimize(const int seed, const std::vector<double>& initpt,
                     /* not const */ std::vector<double>& finalpt,
                     /* unused */ const int& neval, /* unused */ const int& niters)
{
   coliny_problem.reset_neval();

   BasicArray<double> initpt_;
#if 1
   initpt_ << initpt;
#else
initpt_.resize(initpt.size());
initpt_[0] = -9.151;
initpt_[1] = 16.175;
initpt_[2] = 28.005;
initpt_[3] = 0.866;
initpt_[4] = 0.496;
initpt_[5] = 0.071;
initpt_[6] = -0.847;
initpt_[7] = 2 * 3.1416 * (-74.83 / 360);
initpt_[8] = 2 * 3.1416 * (18.36 / 360);
initpt_[9] = 2 * 3.1416 * (12.26 / 360);
initpt_[10] = 2 * 3.1416 * (169.72 / 360);
initpt_[11] = 2 * 3.1416 * (12.54 / 360);
initpt_[12] = 2 * 3.1416 * (4.41 / 360);
initpt_[13] = 2 * 3.1416 * (5.80 / 360);
initpt_[14] = 2 * 3.1416 * (13.62 / 360);
initpt_[15] = 2 * 3.1416 * (-0.83 / 360);
initpt_[16] = 2 * 3.1416 * (-0.44 / 360);
#endif
   BasicArray<double> finalpt_;

#if 0
   // TESTING STUFF
   cerr << "Bound ";
   for (unsigned int i = 0; i < coliny_problem.num_real_params(); i++)
      cerr << coliny_problem.periodic_real_bound(i);
   cerr << endl;
#endif

   colin::real best_value;
   coliny_solver.minimize(coliny_problem, initpt_, seed, false, false, finalpt_, best_value);
   finalpt << finalpt_;
}

#endif
#endif
