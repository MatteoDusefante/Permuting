//
//  An Empirical Evaluation of Permuting in Parallel External Memory
//  utils.h
//
//  Created by Matteo Dusefante on 26/01/16.
//  Copyright © 2016 Matteo Dusefante. All rights reserved.
//

#include <chrono>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <string>
#include <sys/mman.h>
#include <thread>

namespace utils {

#define EVENTS 5
//#define DEBUG 1

/******************************************************************************************/

typedef struct data {
   double time;
   long long *counters;
   data() {
      time = 0.0;
      counters = new long long[EVENTS - 1]();
   };
   ~data() {
      if (counters)
         delete[] counters;
   };
} data;

/******************************************************************************************/

template <typename T>
void copy(T **source, T **destination, T *pointers_length, size_t layers) {

   for (size_t ly = 0; ly < layers; ++ly)
      for (size_t it = 0; it < pointers_length[ly] + 1; ++it)
         destination[ly][it] = source[ly][it];
}

template <typename T> void clean(T *array, size_t length) {

   for (size_t i = 0; i < length; ++i)
      array[i] = 0;
}

/******************************************************************************************/

template <typename T> void random_permutation(T *perm, size_t length) {

   std::srand(std::time(nullptr));

   // std::srand(1);

   for (size_t it = 0; it < length; ++it)
      perm[it] = it;

   for (size_t it = 0; it < length; ++it) {
      size_t index = (std::rand() % (length - it)) + it;
      std::swap(perm[index], perm[it]);
   }
}

/******************************************************************************************/

template <typename T> void identity_permutation(T *perm, size_t length) {

   for (size_t it = 0; it < length; ++it)
      perm[it] = it;
}

/******************************************************************************************/

template <typename T> void populate(T *in, size_t length) {

   for (size_t it = 0; it < length; ++it)
      in[it] = it + 1;
}

/******************************************************************************************/

#ifdef DEBUG

template <typename T, typename S>
void verify(T const *out_p, S const *out, size_t length, bool verbose) {

   for (size_t it = 0; it < length; ++it) {
      // std::cout << out_p[it] << " " << out[it] << std::endl;
      if (out_p[it] != out[it]) {
         std::cout << "ERROR " << out_p[it] << " " << out[it] << " at " << it
                   << std::endl;
         // return;
         exit(0);
      }
   }

   if (verbose)
      std::cout << "No Errors detected" << std::endl;
}

template <typename T, typename S> void verify(T const *out_p, S const *out, size_t length) {

   utils::verify(out_p, out, length, false);
}
#endif
}
