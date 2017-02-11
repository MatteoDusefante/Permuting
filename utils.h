//
//  utils.h
//  project
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

#define MAP_HUGETLB 0x40000
#define EVENTS 5
#define DEBUG 1

/******************************************************************************************/

typedef struct data {
   double time;
   // long long time;
   long long *counters;
   data() {
      time = 0.0;
      // time = 0;
      counters = new long long[EVENTS - 1]();
   };
   ~data() {
      if (counters)
         delete[] counters;
   };
} data;

/*
class Type {
public:
    Type();
    Type operator=(const double &rhs);
    Type operator=(const Type &rhs);
    Type operator+=(const Type &rhs);
    Type operator/=(const int &rhs);

    size_t x;
};

Type::Type() {
    x = 0;
}

Type Type::operator=(const Type &rhs) {

    x = rhs.x;
    return *this;
}

Type Type::operator=(const double &rhs) {

    x = rhs;
    return *this;
}

Type Type::operator+=(const Type &rhs) {

    x += rhs.x;
    return *this;
}

Type Type::operator/=(const int &rhs) {

    x /= rhs;
    return *this;
}
*/
/******************************************************************************************/

template <typename T>
void copy(T **source, T **destination, T *pointers_length, size_t layers) {

   for (size_t ly = 0; ly < layers; ++ly)
      for (auto it = 0; it < pointers_length[ly] + 1; ++it)
         destination[ly][it] = source[ly][it];
}

/******************************************************************************************/

void print_sep() { printf("***********************************\n"); }

/******************************************************************************************/

template <typename RandomAccessIterator>
void KFYShuffle(RandomAccessIterator begin, RandomAccessIterator end) {
   for (unsigned int n = end - begin - 1; n >= 1; --n) {
      unsigned int k = rand() % (n + 1);
      if (k != n) {
         std::iter_swap(begin + k, begin + n);
      }
   }
}

template <typename T> void random_permutation(T *perm, size_t length) {

   std::srand(std::time(nullptr));

   // std::srand(1);

   for (size_t it = 0; it < length; ++it)
      perm[it] = it;

   for (size_t it = 0; it < length; ++it) {
      T idx = (std::rand() % (length - it)) + it;
      std::swap(perm[idx], perm[it]);
   }
   return;
}

/******************************************************************************************/

template <typename T>
void random_permutation_table(T **table, size_t length, size_t layers) {

   // std::srand(std::time(nullptr));
   std::srand(1);

   for (size_t i = 0; i < length; ++i)
      table[0][i] = i;
   /*
      for(size_t i = 0; i < length; ++i) {
         size_t index = (std::rand() % (length - i)) + i;
         std::swap(table[0][index], table[0][i]);
      }*/

   for (size_t ly = 1; ly < layers; ++ly) {
      for (size_t i = 0; i < length; ++i)
         table[ly][i] = 0;
   }

   return;
}

/******************************************************************************************/

template <typename T> void identity_permutation(T *perm, size_t length) {

   for (size_t it = 0; it < length; ++it)
      perm[it] = it;

   return;
}

/******************************************************************************************/

template <typename T> void populate(T *in, size_t length) {

   for (size_t it = 0; it < length; ++it)
      in[it] = ++it;

   return;
}

/******************************************************************************************/

template <typename T>
void populate_table(T **table, size_t length, size_t layers) {

   for (size_t i = 0; i < length; ++i)
      table[0][i] = i + 1;

   for (size_t ly = 1; ly < layers; ++ly) {
      for (size_t i = 0; i < length; ++i)
         table[ly][i] = 0;
   }

   return;
}

/******************************************************************************************/

template <typename T, typename S>
void verify(T *out_p, S *out, size_t length, bool verbose) {

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
   return;
}

/******************************************************************************************/

template <typename T, typename S> void verify(T *out_p, S *out, size_t length) {

   utils::verify(out_p, out, length, false);
}

/******************************************************************************************/

#ifdef ALLOCATE
template <typename T> int *allocate(T name, size_t size) {

   const size_t prefix = 64;

   int *start = (int *)mmap(NULL, size + prefix, PROT_WRITE | PROT_READ,
                            MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB, -1, 0);

   if (start == (void *)-1) {
      std::cout << "failed to use mmap for " << name;
      exit(EXIT_FAILURE);
   }
   *((long long int *)start) = size;
   return start + prefix;
}
#endif

/******************************************************************************************/

template <typename T> void fill(T *start, T *end, size_t bucket_size) {

   *start++ = 0;

   for (; start != end; ++start)
      *start = *(start - 1) + bucket_size;
   return;
}
}
