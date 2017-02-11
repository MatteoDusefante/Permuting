//
//  benchmark.h
//  project
//
//  Created by Matteo Dusefante on 31/05/16.
//  Copyright © 2016 Matteo Dusefante. All rights reserved.
//

#include "utils.h"
// dependecies
#include "algorithms.h"
#include "papils.h"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include "Eigen/Core"
#include "Eigen/SparseCore"
#pragma GCC diagnostic pop

namespace benchmark {

#define EIGEN 0
#define FUNC 0
#define REP 1

template <typename T>
inline void omp_direct_algorithm(T **table_in, T **table_perm,
                                 __attribute__((unused)) T *control_sample,
                                 T *out, size_t length, size_t events,
                                 size_t cores, struct utils::data *collection) {

   utils::data *temp = new utils::data[REP]();

   for (size_t rep = 0; rep < REP; ++rep) {

      papils::papi_start();

      algorithms::omp_permute(&table_in[0][0], &out[0], &table_perm[0][0],
                              length, cores);

      papils::papi_stop();
      papils::papi_collect(&temp[rep]);
   }

   for (size_t evnt = 0; evnt < events - 1; ++evnt) {
      for (size_t rep = 0; rep < REP; ++rep) {
         collection->counters[evnt] += temp[rep].counters[evnt];
      }
      collection->counters[evnt] /= REP;
   }
   for (size_t rep = 0; rep < REP; ++rep)
      collection->time += temp[rep].time;
   collection->time /= REP;

   delete[] temp;

#ifdef DEBUG
   utils::verify(&control_sample[0], &out[0], length);
#endif
}

/******************************************************************************************/

#ifdef EIGEN
template <typename T>
inline void eigen(T **table_in, T **table_perm,
                  __attribute__((unused)) T *control_sample, size_t length,
                  size_t events, size_t cores, struct utils::data *collection) {

   utils::data *temp = new utils::data[REP]();

   Eigen::initParallel();
   Eigen::setNbThreads(cores);

   Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(length);
   Eigen::VectorXi in(length);
   Eigen::VectorXi out(length);

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

#ifdef DEBUG
   utils::verify(&control_sample[0], &out[0], length);
#endif

   for (size_t evnt = 0; evnt < events - 1; ++evnt) {
      for (size_t rep = 0; rep < REP; ++rep) {
         collection->counters[evnt] += temp[rep].counters[evnt];
      }
      collection->counters[evnt] /= (REP * 2);
   }
   for (size_t rep = 0; rep < REP; ++rep)
      collection->time += temp[rep].time;
   collection->time /= (REP * 2);

   delete[] temp;
}
#endif

/******************************************************************************************/

template <typename T>
inline void omp_bucket_algorithm(T **table_in, T **table_perm, T *bucket_size,
                                 __attribute__((unused)) T *control_sample,
                                 T *out, size_t length, size_t layers,
                                 size_t events, size_t cores,
                                 struct utils::data *collection) {

   T **buckets = new T *[layers];
   omp_lock_t **locks = new omp_lock_t *[layers];

   for (size_t ly = 0; ly < layers; ++ly) {
      size_t nbuckets =
          length / bucket_size[ly] +
          (((length / bucket_size[ly]) * bucket_size[ly]) != length);
      buckets[ly] = new T[nbuckets]();
      locks[ly] = new omp_lock_t[nbuckets];
      buckets[ly][0] = 0;
      omp_init_lock(&locks[ly][0]);
      for (size_t it = 1; it < nbuckets; ++it) {
         buckets[ly][it] = buckets[ly][it - 1] + bucket_size[ly];
         omp_init_lock(&locks[ly][it]);
      }
   }

   utils::data *temp = new utils::data[REP]();

   for (size_t rep = 0; rep < REP; ++rep) {

      papils::papi_start();

      algorithms::fast_omp_bucketize(&table_in[0], &table_perm[0], &buckets[0],
                                     &out[0], &bucket_size[0], &locks[0],
                                     length, layers, cores);

      papils::papi_stop();
      papils::papi_collect(&temp[rep]);

      for (size_t ly = 0; ly < layers; ++ly) {
         size_t nbuckets =
             length / bucket_size[ly] +
             (((length / bucket_size[ly]) * bucket_size[ly]) != length);
         buckets[ly][0] = 0;
         for (size_t it = 1; it < nbuckets; ++it) {
            buckets[ly][it] = buckets[ly][it - 1] + bucket_size[ly];
         }
      }
   }

#ifdef DEBUG
   utils::verify(&control_sample[0], &out[0], length);
#endif

   for (size_t evnt = 0; evnt < events - 1; ++evnt) {
      for (size_t rep = 0; rep < REP; ++rep) {
         collection->counters[evnt] += temp[rep].counters[evnt];
      }
      collection->counters[evnt] /= REP;
   }
   for (size_t rep = 0; rep < REP; ++rep)
      collection->time += temp[rep].time;
   collection->time /= REP;

   for (size_t ly = 0; ly < layers; ++ly) {
      delete[] buckets[ly];
      delete[] locks[ly];
   }

   delete[] temp;
   delete[] locks;
   delete[] buckets;
}

/******************************************************************************************/

template <typename T>
inline void omp_table_algorithm(T **table_in, T **table_perm, T *bucket_size,
                                __attribute__((unused)) T *control_sample,
                                size_t length, size_t layers, size_t events,
                                size_t cores, struct utils::data *collection) {

   long **table = new T *[layers + 1];
   long **buckets = new T *[layers];
   omp_lock_t **locks = new omp_lock_t *[layers];

   for (size_t ly = 0; ly < layers; ++ly) {
      table[ly] = new T[length]();
      size_t nbuckets =
          length / bucket_size[ly] +
          (((length / bucket_size[ly]) * bucket_size[ly]) != length);
      buckets[ly] = new T[nbuckets]();
      locks[ly] = new omp_lock_t[nbuckets];
      buckets[ly][0] = 0;
      omp_init_lock(&locks[ly][0]);
      for (size_t it = 1; it < nbuckets; ++it) {
         buckets[ly][it] = buckets[ly][it - 1] + bucket_size[ly];
         omp_init_lock(&locks[ly][it]);
      }
   }

   table[layers] = new T[length]();

   papils::papi_start();

   algorithms::multi_layer_preprocessing(&table[0], &table_perm[0], &buckets[0],
                                         &bucket_size[0], &locks[0], length,
                                         layers, cores);

   papils::papi_stop();
   papils::papi_print();

   return;

   for (size_t ly = 0; ly < layers; ++ly) {
      delete[] buckets[ly];
      delete[] locks[ly];
   }

   delete[] locks;
   delete[] buckets;

   utils::data *temp = new utils::data[REP]();

   for (size_t rep = 0; rep < REP; ++rep) {

      papils::papi_start();

      algorithms::fast_omp_multi_layer_table(&table_in[0], &table[0], length,
                                             layers, cores);

      papils::papi_stop();
      papils::papi_collect(&temp[rep]);
   }

#ifdef DEBUG
   utils::verify(&control_sample[0], &table_in[layers + 1][0], length);
#endif

   for (size_t evnt = 0; evnt < events - 1; ++evnt) {
      for (size_t rep = 0; rep < REP; ++rep) {
         collection->counters[evnt] += temp[rep].counters[evnt];
      }
      collection->counters[evnt] /= REP;
   }
   for (size_t rep = 0; rep < REP; ++rep)
      collection->time += temp[rep].time;
   collection->time /= REP;

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
                               size_t events, size_t cores,
                               struct utils::data *collection) {

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

   // nbuckets[0] = 1;
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
   //    &table_perm[0], &table_pt[0], &nbuckets[0], &bucket_size[0],
   //    &pointers_length[0], cores, length, layers);

   papils::papi_start();

   algorithms::parallel_bucket_preprocessing(&table_perm[0], &table_pt[0],
                                             &bucket_size[0], &buckets[0],
                                             cores, length, layers, delta);

   papils::papi_stop();
   papils::papi_print();

   return;

   utils::data *temp = new utils::data[REP]();

   for (size_t rep = 0; rep < REP; ++rep) {

      utils::copy(&table_pt[0], &table_pt_copy[0], &pointers_length[0], layers);

      /*
              for(size_t ly = 0; ly < layers; ++ly) {
                  for(size_t i = 0; i < pointers_length[ly] + 1; ++i)
                      printf("%ld ", table_pt_copy[ly][i]);
                  printf(" -\n\n");
              }*/

      papils::papi_start();

      // algorithms::parallel_bucket(
      //    &table_in[0], &table_perm[0], &table_pt[0], &out[0], &nbuckets[0],
      //    &bucket_size[0], cores, length,
      //    layers);
      algorithms::parallel_bucket(&table_in[0], &table_perm[0],
                                  &table_pt_copy[0], &out[0], &bucket_size[0],
                                  &buckets[0], cores, length, layers, delta);

      papils::papi_stop();
      papils::papi_collect(&temp[rep]);
   }

#ifdef DEBUG
   utils::verify(&control_sample[0], &out[0], length);
#endif

   for (size_t evnt = 0; evnt < events - 1; ++evnt) {
      for (size_t rep = 0; rep < REP; ++rep) {
         collection->counters[evnt] += temp[rep].counters[evnt];
      }
      collection->counters[evnt] /= REP;
   }
   for (size_t rep = 0; rep < REP; ++rep)
      collection->time += temp[rep].time;
   collection->time /= REP;

   delete[] temp;

   for (size_t ly = 0; ly < layers; ++ly) {
      delete[] table_pt[ly];
      // delete[] table_pt_copy[ly];
   }

   delete[] pointers_length;
   delete[] table_pt;
   delete[] table_pt_copy;
}

/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/

template <typename T>
inline long long *direct_algorithm(T *in, T *perm, size_t length,
                                   size_t events) {

   long long *collection = new long long[events];

   T *out = new T[length];

   papils::papi_start();

   algorithms::permute(&in[0], &out[0], &perm[0], length);

   papils::papi_stop();
   papils::papi_collect(&collection[0]);

   delete[] out;

   return collection;
}

/******************************************************************************************/

template <typename T>
inline long long *bucket_algorithm(T *in, T *perm, T *bucket_size,
                                   T *control_sample, size_t length,
                                   size_t layers, size_t iter, size_t events) {

   int *out = new int[length];

   int **buckets = new int *[layers];
   int **table_in = new int *[layers + 1];
   int **table_perm = new int *[layers + 1];

   long long **collection = new long long *[iter];

   for (size_t i = 0; i < iter; ++i)
      collection[i] = new long long[events]();

   for (size_t i = 0; i < iter; ++i) {

      for (size_t ly = 0; ly < layers + 1; ++ly) {
         table_in[ly] = new int[length]();
         table_perm[ly] = new int[length]();
      }
      std::copy(in, in + length, &table_in[0][0]);
      std::copy(perm, perm + length, &table_perm[0][0]);

      for (int ly = 0; ly < layers; ++ly) {
         int nbuckets =
             length / bucket_size[ly] +
             (((length / bucket_size[ly]) * bucket_size[ly]) != length);
         buckets[ly] = new int[nbuckets]();
         for (int it = 1; it < nbuckets; ++it)
            buckets[ly][it] = buckets[ly][it - 1] + bucket_size[ly];
      }

      papils::papi_start();

      algorithms::fast_bucketize(&table_in[0], &table_perm[0], &buckets[0],
                                 &out[0], &bucket_size[0], length, layers);

      papils::papi_stop();
      papils::papi_collect(&collection[i][0]);

      utils::verify(&control_sample[0], &out[0], length);
   }

   delete[] out;
   delete[] buckets;
   delete[] table_in;
   delete[] table_perm;

   return &collection[0];
}

/******************************************************************************************/

template <typename T>
inline long long *table_algorithm(T *in, T *perm, T *bucket_size,
                                  T *control_sample, size_t length,
                                  size_t layers, size_t iter, size_t events) {

   T **table = new T *[layers + 1];
   T **table_out = new T *[layers + 2];

   for (size_t it = 0; it < layers + 1; ++it)
      table[it] = new T[length]();

   algorithms::multi_layer_preprocessing(&perm[0], &table[0], &bucket_size[0],
                                         length, 0, layers);

   long long **collection = new long long *[iter];

   for (size_t i = 0; i < iter; ++i)
      collection[i] = new long long[events]();

   for (size_t i = 0; i < iter; ++i) {

      for (int ly = 0; ly < layers + 2; ++ly)
         table_out[ly] = new T[length]();

      std::copy(in, in + length, &table_out[0][0]);

      papils::papi_start();

      algorithms::fast_multi_layer_table(&table_out[0], &table[0], length,
                                         layers);

      papils::papi_stop();
      papils::papi_collect(&collection[i][0]);

      utils::verify(&control_sample[0], &table_out[layers + 1][0], length);
   }

   delete[] table;
   delete[] table_out;

   return &collection[0];
}

/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/

template <typename T> inline T *min_value(T **collection, size_t iter) {

   T minv = collection[0][0];
   size_t index = 0;
   for (size_t i = 1; i < iter; ++i) {
      if (collection[i][0] < minv) {
         minv = collection[i][0];
         index = i;
      }
   }
   return collection[index];
}

/******************************************************************************************/

template <typename T>
inline T *filtered_average(T **collection, size_t events, size_t iter) {

   T *results = new T[events]();

   for (size_t j = 0; j < events; ++j) {
      T worst1, worst2;
      if (collection[0][j] > collection[1][j]) {
         worst1 = collection[0][j];
         worst2 = collection[1][j];
      } else {
         worst1 = collection[1][j];
         worst2 = collection[0][j];
      }
      for (size_t i = 0; i < iter; ++i) {
         if (collection[i][j] > worst2) {
            if (collection[i][j] > worst1) {
               worst2 = worst1;
               worst1 = collection[i][j];
            } else {
               worst2 = collection[i][j];
            }
         }
         results[j] += collection[i][j];
      }
      results[j] = (results[j] - worst1 - worst2) / (iter - 2);
      if (results[j] < 0)
         results[j] = 0;
   }

   for (size_t i = 0; i < iter; ++i)
      delete[] collection[i];

   delete[] collection;

   return results;
}

/******************************************************************************************/

template <typename T>
inline T *rms(T **collection, size_t events, size_t iter) {

   T *results = new T[events]();

   for (size_t j = 0; j < events; ++j) {
      for (size_t i = 0; i < iter; ++i)
         results[j] += (T)pow(collection[i][j], 2);
      results[j] = (T)sqrt(results[j] / iter);
   }

   for (size_t i = 0; i < iter; ++i)
      delete[] collection[i];

   delete[] collection;

   return results;
}

/******************************************************************************************/

template <typename T>
inline T *average(T **collection, size_t events, size_t iter) {

   T *results = new T[events]();

   for (size_t j = 0; j < events; ++j) {
      for (size_t i = 0; i < iter; ++i)
         results[j] += collection[i][j];
      results[j] = results[j] / iter;
   }

   for (size_t i = 0; i < iter; ++i)
      delete[] collection[i];

   delete[] collection;

   return results;
}

/******************************************************************************************/

template <typename T>
inline T *postprocess(T **collection, size_t events, size_t iter, size_t type) {

   T *results = new T[events]();

   if (iter == 1)
      return collection[0];

   switch (type) {
   case 0:
      return min_value(&collection[0], events);
      break;
   case 1:
      return average(&collection[0], events, iter);
      break;
   case 2:
      return filtered_average(&collection[0], events, iter);
      break;
   case 3:
      return rms(&collection[0], events, iter);
      break;
   }
   /*
   for(size_t i = 0; i < iter; ++i)
       delete[] collection[i];

   delete[] collection;*/

   return results;
}

/******************************************************************************************/

template <typename T> inline T compute_cost(size_t type, T *results) {

   size_t res;

   switch (type) {
   case 0:
      res = results[0]; // time
      break;
   case 1:
      res = results[0] * 100 + results[1] + results[2] + results[3]; // plain
      break;
   case 2:
      res = results[0] * 1000 + results[1] + 70 * results[2] +
            350 * results[3]; // cache latency
      break;
   }
   return res;
}
}
