//
//  An Empirical Evaluation of Permuting in Parallel External Memory
//  benchmark.h
//
//  Created by Matteo Dusefante on 31/05/16.
//  Copyright © 2016 Matteo Dusefante. All rights reserved.
//

#include "utils.h"
// dependencies
#include "algorithms.h"
#include "papils.h"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include "Eigen/Core"
#include "Eigen/SparseCore"
#pragma GCC diagnostic pop

namespace Eigen {
typedef Matrix<uint32_t, Dynamic, 1> VectorXuint32;
}

namespace benchmark {

template <typename T>
inline void omp_direct_algorithm(T **table_in, T **table_perm,
                                 __attribute__((unused)) T *control_sample,
                                 T *out, size_t length, size_t cores,
                                 struct utils::data *collection) {

   utils::data *temp = new utils::data[REP]();

   for (size_t rep = 0; rep < REP; ++rep) {

      papils::papi_start();

      algorithms::omp_permute(&table_in[0][0], &out[0], &table_perm[0][0],
                              length, cores);

      papils::papi_stop();
      papils::papi_collect(&temp[rep]);
   }

   utils::postprocess(&temp[0], &collection[0]);

   delete[] temp;

#ifdef DEBUG
   utils::verify(&control_sample[0], &out[0], length);
#endif
}

/******************************************************************************************/

template <typename T>
inline void eigen(T **table_in, T **table_perm,
                  __attribute__((unused)) T *control_sample, size_t length,
                  size_t cores, struct utils::data *collection) {

   utils::data *temp = new utils::data[REP]();

   Eigen::initParallel();
   Eigen::setNbThreads(cores);

   Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, uint32_t> perm(
       length);
   Eigen::VectorXuint32 in(length);
   Eigen::VectorXuint32 out(length);

   for (size_t i = 0; i < length; ++i) {
      in[i] = table_in[0][i];
      perm.indices()[i] = table_perm[0][i];
   }

   for (size_t rep = 0; rep < REP; ++rep) {

      papils::papi_start();

      out = perm * in;

      papils::papi_stop();
      papils::papi_collect(&temp[rep]);
   }

   utils::postprocess(&temp[0], &collection[0]);

#ifdef DEBUG
   utils::verify(&control_sample[0], &out[0], length);
#endif

   delete[] temp;
}

/******************************************************************************************/

template <typename T>
inline void omp_bucket_algorithm(T **table_in, T **table_perm, T *buckets,
                                 __attribute__((unused)) T *control_sample,
                                 T *out, size_t length, size_t layers,
                                 size_t cores, struct utils::data *collection) {

   T **b_table = new T *[layers];
   T *bucket_size = new T[layers]();
   omp_lock_t *locks = new omp_lock_t[buckets[layers]];
   for (size_t it = 0; it < buckets[layers]; ++it)
      omp_init_lock(&locks[it]);

   for (size_t ly = 0; ly < layers; ++ly) {
      bucket_size[ly] = length / buckets[ly] +
                        (((length / buckets[ly]) * buckets[ly]) != length);
      b_table[ly] = new T[buckets[ly]]();
      b_table[ly][0] = 0;
      for (size_t it = 1; it < buckets[ly]; ++it)
         b_table[ly][it] = b_table[ly][it - 1] + bucket_size[ly];
   }

   utils::data *temp = new utils::data[REP]();

   for (size_t rep = 0; rep < REP; ++rep) {

      papils::papi_start();

      algorithms::bucketize(&table_in[0], &table_perm[0], &b_table[0], &out[0],
                            &bucket_size[0], &locks[0], length, layers, cores);

      papils::papi_stop();
      papils::papi_collect(&temp[rep]);

      for (size_t ly = 0; ly < layers; ++ly) {
         b_table[ly][0] = 0;
         for (size_t it = 1; it < buckets[ly]; ++it) {
            b_table[ly][it] = b_table[ly][it - 1] + bucket_size[ly];
         }
      }
   }

   utils::postprocess(&temp[0], &collection[0]);

#ifdef DEBUG
   utils::verify(&control_sample[0], &out[0], length);
#endif

   for (size_t ly = 0; ly < layers; ++ly)
      delete[] b_table[ly];

   delete[] temp;
   delete[] locks;
   delete[] b_table;
}

/******************************************************************************************/

template <typename T>
inline void omp_table_algorithm(T **table_in, T **table_perm, T *buckets,
                                __attribute__((unused)) T *control_sample,
                                T *out, size_t length, size_t layers,
                                size_t cores, struct utils::data *collection) {

   T **table = new T *[layers + 1];
   T **b_table = new T *[layers];
   T *bucket_size = new T[layers]();
   omp_lock_t *locks = new omp_lock_t[buckets[layers - 1]];
   for (size_t it = 0; it < buckets[layers - 1]; ++it)
      omp_init_lock(&locks[it]);

   for (size_t ly = 0; ly < layers; ++ly) {
      table[ly] = new T[length]();
      bucket_size[ly] = length / buckets[ly] +
                        (((length / buckets[ly]) * buckets[ly]) != length);
      b_table[ly] = new T[buckets[ly]]();
      b_table[ly][0] = 0;
      for (size_t it = 1; it < buckets[ly]; ++it)
         b_table[ly][it] = b_table[ly][it - 1] + bucket_size[ly];
   }

   table[layers] = new T[length]();

   papils::papi_start();

   algorithms::multi_layer_table_preprocessing(
       &table[0], &table_perm[0], &b_table[0], &bucket_size[0], &locks[0],
       length, layers, cores);
   papils::papi_stop();
   // papils::papi_print();

   for (size_t ly = 0; ly < layers; ++ly)
      delete[] b_table[ly];

   delete[] locks;
   delete[] b_table;
   delete[] bucket_size;

   utils::data *temp = new utils::data[REP]();

   for (size_t rep = 0; rep < REP; ++rep) {

      papils::papi_start();

      algorithms::multi_layer_table(&table_in[0], &table[0], &out[0], length,
                                    layers, cores);

      papils::papi_stop();
      papils::papi_collect(&temp[rep]);
   }

#ifdef DEBUG
   for (size_t it = 0; it < length; ++it)
      out[it] = table_in[layers + 1][it];
   utils::verify(&control_sample[0], &out[0], length);
#endif

   utils::postprocess(&temp[0], &collection[0]);

   for (size_t ly = 0; ly < layers + 1; ++ly)
      delete[] table[ly];

   delete[] temp;
   delete[] table;
}

/******************************************************************************************/

template <typename T>
inline void parallel_algorithm(T **table_in, T **table_perm, T *buckets,
                               __attribute__((unused)) T *control_sample,
                               T *out, size_t length, size_t layers,
                               size_t cores, struct utils::data *collection) {

   T *pointers_length = new T[layers];
   T *bucket_size = new T[layers + 1];
   // T *nbuckets = new T[layers + 1];
   T **table_pt = new T *[layers];
   for (size_t ly = 0; ly < layers; ++ly)
      table_pt[ly] = new T[cores * buckets[ly]]();
   T **table_pt_copy = new T *[layers];
   for (size_t ly = 0; ly < layers; ++ly)
      table_pt_copy[ly] = new T[cores * buckets[ly]]();

   size_t delta =
       floor(length / cores) + (((length / cores) * cores) != length);

   size_t size = 0;

   bucket_size[0] = length;
   for (size_t ly = 0; ly < layers; ++ly) {
      // nbuckets[ly + 1] = buckets[ly];
      bucket_size[ly + 1] = floor(length / buckets[ly]) +
                            (((length / buckets[ly]) * buckets[ly]) != length);
      // bucket_size[ly + 1] = ceil(length / buckets[ly]);
      pointers_length[ly] = cores * buckets[ly];
      size += pointers_length[ly];
   }

   // algorithms::parallel_bucket_preprocessing(
   //     &table_perm[0], &table_pt[0], &nbuckets[0], &bucket_size[0],
   //     &pointers_length[0], cores, length, layers);

   papils::papi_start();

   algorithms::parallel_bucket_preprocessing(&table_perm[0], &table_pt[0],
                                             &bucket_size[0], &buckets[0],
                                             cores, length, layers, delta);

   papils::papi_stop();
   // papils::papi_print();

   utils::data *temp = new utils::data[REP]();

   for (size_t rep = 0; rep < REP; ++rep) {

      utils::copy(&table_pt[0], &table_pt_copy[0], &pointers_length[0], layers);

      papils::papi_start();

      // algorithms::parallel_bucket(&table_in[0], &table_perm[0], &table_pt[0],
      //                             &out[0], &nbuckets[0], &bucket_size[0],
      //                             cores,
      //                             length, layers);
      algorithms::parallel_bucket(&table_in[0], &table_perm[0],
                                  &table_pt_copy[0], &out[0], &bucket_size[0],
                                  &buckets[0], cores, length, layers, delta);

      papils::papi_stop();
      papils::papi_collect(&temp[rep]);
   }

#ifdef DEBUG
   utils::verify(&control_sample[0], &out[0], length);
#endif

   utils::postprocess(&temp[0], &collection[0]);

   delete[] temp;

   for (size_t ly = 0; ly < layers; ++ly) {
      delete[] table_pt[ly];
      delete[] table_pt_copy[ly];
   }

   delete[] table_pt;
   delete[] bucket_size;
   delete[] table_pt_copy;
   delete[] pointers_length;
}
}
