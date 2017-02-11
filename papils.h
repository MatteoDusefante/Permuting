//
//  papils.h
//  project
//
//  Created by Matteo Dusefante on 15/02/16.
//  Copyright © 2016 Matteo Dusefante. All rights reserved.
//

#include <papi.h>
//#include "papi.h"
#include <cstring>
#include <omp.h>
#include <string>

#define MAXNE 7

namespace papils {

double start_time, time;

struct {
   int ne;
   const char **names;
   long long values[MAXNE];
   unsigned long real_time;
   unsigned long virtual_time;
   unsigned long cycles;
   int EventSet;
} papicounter;

/*
PAPI_Lx_TCM : Level x cache misses
PAPI_Lx_ICM : Level x instruction cache misses
PAPI_TLB_DM : Data translation lookaside buffer misses
PAPI_TLB_IM : Instruction translation lookaside buffer misses
*/

const char *names[] = {"PAPI_L1_TCM", "PAPI_L2_TCM", "PAPI_L3_TCM",
                       "PAPI_TLB_DM"};

const char *tlb[] = {"PAPI_TLB_DM", "PAPI_TLB_IM"};

#ifdef OLD
template <typename S> void critical(const S *s) { std::cout << s << std::endl; }
#endif

void papi_init() {

   papicounter.ne = 3;
   // papicounter.ne = 2;

   papicounter.names = names;
   // papicounter.names = icm;
   // papicounter.names = tlb;

   std::string str;

   int event;
   if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT)
      std::cout << "PAPI_library_init failed" << std::endl;
   papicounter.EventSet = PAPI_NULL;
   if (PAPI_create_eventset(&(papicounter.EventSet)) != PAPI_OK)
      std::cout << "PAPI_create_eventset failed" << std::endl;

   for (int i = 0; i < papicounter.ne; ++i) {
      if (PAPI_event_name_to_code((char *)papicounter.names[i], &event) !=
          PAPI_OK)
         std::cout << "PAPI_event_name_to_code failed" << std::endl;
      if (PAPI_add_event(papicounter.EventSet, event) != PAPI_OK)
         std::cout << "Failed to add event " << papicounter.names[i]
                   << std::endl;
   }
}

void papi_multiplex_init() {

   papicounter.ne = 4;
   papicounter.names = names;

   std::string str;
   int event;
   PAPI_option_t opt;

   if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT)
      std::cout << "PAPI_library_init failed" << std::endl;

   if (PAPI_thread_init(pthread_self) != PAPI_OK)
      std::cout << "PAPI_thread_init failed" << std::endl;

   if (PAPI_get_multiplex(papicounter.EventSet) != PAPI_OK)
      std::cout << "This event set is not enabled for multiplexing"
                << std::endl;

   if (PAPI_multiplex_init() != PAPI_OK)
      std::cout << "PAPI_init_multiplex failed" << std::endl;

   papicounter.EventSet = PAPI_NULL;
   if (PAPI_create_eventset(&(papicounter.EventSet)) != PAPI_OK)
      std::cout << "PAPI_create_eventset failed" << std::endl;

   if (PAPI_assign_eventset_component(papicounter.EventSet, 0) != PAPI_OK)
      std::cout << "PAPI_assign_event_component failed" << std::endl;

   if (PAPI_set_multiplex(papicounter.EventSet) != PAPI_OK)
      std::cout << "PAPI_set_multiplex failed" << std::endl;

   memset(&opt, 0x0, sizeof(PAPI_option_t));
   opt.inherit.inherit = PAPI_INHERIT_ALL;
   opt.inherit.eventset = papicounter.EventSet;
   if (PAPI_set_opt(PAPI_INHERIT, &opt) != PAPI_OK)
      std::cout << "PAPI_set_opt INHERIT failed" << std::endl;

   for (int i = 0; i < papicounter.ne; ++i) {
      if (PAPI_event_name_to_code((char *)papicounter.names[i], &event) !=
          PAPI_OK)
         std::cout << "PAPI_event_name_to_code failed" << std::endl;
      if (PAPI_add_event(papicounter.EventSet, event) != PAPI_OK)
         std::cout << "Failed to add event " << papicounter.names[i]
                   << std::endl;
      // std::cout << "Granularity( " << papicounter.names[i] << " ) : " <<
      // PAPI_get_opt("PAPI_GET_GRANUL", event ) <<
      // std::endl;
   }
}

void papi_start() {

   if (PAPI_reset(papicounter.EventSet) != PAPI_OK)
      std::cout << "PAPI reset failed" << std::endl;

   if (PAPI_start(papicounter.EventSet) != PAPI_OK)
      std::cout << "PAPI start failed" << std::endl;

   start_time = omp_get_wtime();

   papicounter.cycles = PAPI_get_virt_cyc();
   papicounter.real_time = PAPI_get_real_usec();
   // papicounter.virtual_time = PAPI_get_virt_usec();
}

void papi_stop() {
   if (PAPI_stop(papicounter.EventSet, papicounter.values) != PAPI_OK)
      std::cout << "PAPI stop failed" << std::endl;

   time = omp_get_wtime() - start_time;

   papicounter.cycles = PAPI_get_virt_cyc() - papicounter.cycles;
   papicounter.real_time = PAPI_get_real_usec() - papicounter.real_time;
   // papicounter.virtual_time = PAPI_get_virt_usec() -
   // papicounter.virtual_time;
}

void papi_print() {

   for (int i = 0; i < papicounter.ne; ++i)
      std::cout << papicounter.names[i] << " : " << papicounter.values[i] * 12
                << std::endl;

   std::cout << "Sec         : " << time << std::endl;

   // std::cout << "(Sec (real)    : " << papicounter.real_time << std::endl;
   // std::cout << "USec (virtual) : " << papicounter.virtual_time << std::endl;
   // std::cout << "Cycles : " << papicounter.cycles << std::endl;
}

template <typename T> void papi_collect(T *results) {

   // results[0] = papicounter.time;
   results->time = time;
   // results->time = papicounter.real_time;
   // results->time = papicounter.virtual_time;

   // std::cout << "Sec : " << time << std::endl;

   for (int i = 0; i < papicounter.ne; ++i)
      results->counters[i] = papicounter.values[i];
}
}
