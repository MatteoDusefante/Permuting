//
//  algorithms.h
//  project
//
//  Created by Matteo Dusefante on 26/01/16.
//  Copyright © 2016 Matteo Dusefante. All rights reserved.
//

#include <algorithm>

namespace algorithms {

template <typename T, typename P>
inline void permute(T *in, T *out, P *perm, size_t length) {

   for (size_t it = 0; it < length; ++it)
      out[perm[it]] = in[it];
   return;
}

/******************************************************************************************/

template <typename T, typename P>
inline void omp_permute(T *in, T *out, P *perm, size_t length, size_t cores) {

   size_t it;
   omp_set_dynamic(0);
#pragma omp parallel for private(it) num_threads(cores)
   for (it = 0; it < length; ++it)
      out[perm[it]] = in[it];
   return;
}

/******************************************************************************************/

template <typename T, typename S>
inline void fast_omp_bucketize(T **table_in, T **table_perm, T **buckets,
                               T *out, T *bucket_size, S **locks, size_t length,
                               size_t layers, size_t cores) {

   T bucket, index;
   size_t it;

   omp_set_dynamic(0);

   for (size_t ly = 0; ly < layers; ++ly) {
#pragma omp parallel for private(bucket, index) num_threads(cores)
      for (it = 0; it < length; ++it) {
         bucket = (T)floor(table_perm[ly][it] / bucket_size[ly]);
         omp_set_lock(&locks[ly][bucket]);
         index = buckets[ly][bucket]++;
         omp_unset_lock(&locks[ly][bucket]);
         table_in[ly + 1][index] = table_in[ly][it];
         table_perm[ly + 1][index] = table_perm[ly][it];
      }
   }

#pragma omp parallel for private(it) num_threads(cores)
   for (it = 0; it < length; ++it)
      out[table_perm[layers][it]] = table_in[layers][it];
}

/******************************************************************************************/

template <typename T, typename S>
inline void multi_layer_preprocessing(T **table, T **table_perm, T **buckets,
                                      T *bucket_size, S **locks, size_t length,
                                      size_t layers, size_t cores) {

   T index, bucket;
   size_t it;
   omp_set_dynamic(0);
   for (size_t ly = 0; ly < layers; ++ly) {
#pragma omp parallel for private(bucket, index) num_threads(cores)
      for (size_t it = 0; it < length; ++it) {
         bucket = (T)floor(table_perm[ly][it] / bucket_size[ly]);
         omp_set_lock(&locks[ly][bucket]);
         index = buckets[ly][bucket]++;
         omp_unset_lock(&locks[ly][bucket]);
         table[ly][it] = index;
         table_perm[ly + 1][index] = table_perm[ly][it];
      }
   }
#pragma omp parallel for private(it) num_threads(cores)
   for (it = 0; it < length; ++it)
      table[layers][it] = table_perm[layers][it];
}

/******************************************************************************************/

template <typename T>
inline void fast_omp_multi_layer_table(T **table_out, T **table, size_t length,
                                       size_t layers, size_t cores) {

   size_t it;
   omp_set_dynamic(0);
   for (size_t ly = 0; ly < layers + 1; ++ly)
#pragma omp parallel for private(it) num_threads(cores)
      for (it = 0; it < length; ++it)
         table_out[ly + 1][table[ly][it]] = table_out[ly][it];

   return;
}

/******************************************************************************************/

template <typename T>
inline void parallel_bucket_preprocessing(T **table_perm, T **table_pt,
                                          T *bucket_size, T *buckets,
                                          size_t cores, size_t length,
                                          size_t layers, size_t delta) {

   T bucket;
   size_t c, it, index;

   omp_set_dynamic(0);
   omp_set_num_threads(cores);

   for (size_t ly = 0; ly < layers; ++ly) {
#pragma omp parallel for private(c, it, bucket) shared(table_pt)               \
                                     schedule(static, delta)
      for (it = 0; it < length; ++it) {
         c = omp_get_thread_num();
         bucket = floor(table_perm[ly][it] / bucket_size[ly + 1]);
         table_pt[ly][c * buckets[ly] + bucket]++;
      }

      size_t prev = 0, curr = 0;
      for (size_t b = 0; b < (size_t)buckets[ly]; ++b) {
         for (size_t c = 0; c < cores; ++c) {
            curr = table_pt[ly][buckets[ly] * c + b];
            table_pt[ly][buckets[ly] * c + b] = prev;
            prev = table_pt[ly][buckets[ly] * c + b] + curr;

            // printf("%ld  %ld  %ld \n",curr,prev, table_pt[ly][buckets[ly] * c
            // + b]);

            // table_pt[ly][buckets[ly] * c + b] += table_pt[ly][buckets[ly] *
            // (c - 1) + b];
            // temp += table_pt[ly][buckets[ly] * (c - 1) + b];
         }
         // temp = table_pt[ly][buckets[ly] * (cores - 1) + b];
         // printf("-- %ld  %ld\n", prev, curr);
         // printf("\n");
      }

      T *tp = new T[buckets[ly] * cores];
      std::copy(&table_pt[ly][0], &table_pt[ly][0] + buckets[ly] * cores,
                &tp[0]);

#pragma omp parallel for private(c, it, index, bucket) shared(table_pt)        \
                                     schedule(static, delta)
      for (it = 0; it < length; ++it) {
         c = omp_get_thread_num();
         bucket = floor(table_perm[ly][it] / bucket_size[ly + 1]);
         index = tp[c * buckets[ly] + bucket]++;
         table_perm[ly + 1][index] = table_perm[ly][it];
      }

      delete[] tp;
   }
}

/******************************************************************************************/

template <typename T>
inline void parallel_bucket(T **table_in, T **table_perm, T **table_pt, T *out,
                            T *bucket_size, T *buckets, size_t cores,
                            size_t length, size_t layers, size_t delta) {

   size_t c, it, index;

   T bucket;

   omp_set_dynamic(0);
   omp_set_num_threads(cores);
   for (size_t ly = 0; ly < layers; ++ly) {
#pragma omp parallel for private(c, it, index, bucket) shared(table_pt)        \
                                     schedule(static, delta)
      for (it = 0; it < length; ++it) {
         c = omp_get_thread_num();
         bucket = floor(table_perm[ly][it] / bucket_size[ly + 1]);
         index = table_pt[ly][c * buckets[ly] + bucket]++;
         table_in[ly + 1][index] = table_in[ly][it];
         // table_perm[ly + 1][index] = table_perm[ly][it];
      }
   }

#pragma omp parallel for private(it)
   for (size_t it = 0; it < length; ++it)
      out[table_perm[layers][it]] = table_in[layers][it];
   return;
}

}
