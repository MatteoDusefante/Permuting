//
//  An Empirical Evaluation of Permuting in Parallel External Memory
//  benchmark.cpp
//
//  Created by Matteo Dusefante on 31/05/16.
//  Copyright © 2016 Matteo Dusefante. All rights reserved.
//

#include "benchmark.h"

#define start 500000000
#define end 500000000
#define delta 500000000
#define layers 1

int main(__attribute__((unused)) int argc,
         __attribute__((unused)) char *argv[]) {

   unsigned int events = EVENTS;
   unsigned int buckets = 1;

   uint32_t nbuckets[] = {10000};

   papils::papi_multiplex_init();
   std::ofstream outfile;
   outfile.open("output.txt", std::ios_base::app);

   std::string input = "";

   uint32_t *out_p = new uint32_t[end];
   uint32_t *out = new uint32_t[end];

   uint32_t **table_perm = new uint32_t *[layers + 1];
   uint32_t **table_in = new uint32_t *[layers + 1];

   for (size_t ly = 0; ly < layers + 1; ++ly)
      table_perm[ly] = new uint32_t[end]();

   for (size_t ly = 0; ly < layers + 2; ++ly)
      table_in[ly] = new uint32_t[end]();

   std::cout << "Start..." << std::endl;

   input += "#start=" + std::to_string(delta) + " delta=" +
            std::to_string(delta) + " end=" + std::to_string(end) + " layers=" +
            std::to_string(layers) + " nb=" + std::to_string(nbuckets[0]) +
            "\n";

   size_t cores = std::thread::hardware_concurrency();

#ifdef BUCKET
   results_bucket = new utils::data();
   utils::data *results_direct, *results_bucket, *results_eigen, *results_table,
       *results_parallel;
#else
   utils::data *results_direct, *results_eigen, *results_table,
       *results_parallel;
#endif

   for (size_t length = start; length < end + delta; length += delta) {

      results_direct = new utils::data();
      results_eigen = new utils::data();
      results_table = new utils::data();
      results_parallel = new utils::data();

      utils::random_permutation(&table_perm[0][0], length);
      utils::populate(&table_in[0][0], length);

#ifdef DEBUG
      algorithms::permute(&table_in[0][0], &out_p[0], &table_perm[0][0],
                          length);
#endif

      benchmark::omp_direct_algorithm(&table_in[0], &table_perm[0], &out_p[0],
                                      &out[0], length, events, cores,
                                      results_direct);
      utils::clean(&out[0], length);

      benchmark::eigen(&table_in[0], &table_perm[0], &out_p[0], length, events,
                       cores, results_eigen);

      utils::clean(&out[0], length);

      for (size_t it = 0; it < buckets; ++it) {

#ifdef BUCKET
         benchmark::omp_bucket_algorithm(
             &table_in[0], &table_perm[0], &nbuckets[0], &out_p[0], &out[0],
             length, layers, events, cores, results_bucket);

         utils::clean(&out[0], length);
#endif

         benchmark::omp_table_algorithm(
             &table_in[0], &table_perm[0], &nbuckets[0], &out_p[0], &out[0],
             length, layers, events, cores, results_table);

         utils::clean(&out[0], length);

         benchmark::parallel_algorithm(&table_in[0], &table_perm[0],
                                       &nbuckets[0], &out_p[0], &out[0], length,
                                       layers, events, cores, results_parallel);
         utils::clean(&out[0], length);
      }

      input += std::to_string(length) + "   ";
      input += std::to_string(results_direct->time) + "   ";
      for (size_t i = 0; i < events - 1; ++i)
         input += std::to_string(results_direct->counters[i]) + "   ";
#ifdef BUCKET
      input += std::to_string(results_bucket->time) + "   ";
      for (size_t i = 0; i < events; ++i)
         input += std::to_string(results_bucket->counters[i]) + "   ";
#endif
      input += std::to_string(results_table->time) + "   ";
      for (size_t i = 0; i < events - 1; ++i)
         input += std::to_string(results_table->counters[i]) + "   ";
      input += std::to_string(results_parallel->time) + "   ";
      for (size_t i = 0; i < events - 1; ++i)
         input += std::to_string(results_parallel->counters[i]) + "   ";
      input += std::to_string(results_eigen->time) + "   ";
      for (size_t i = 0; i < events - 1; ++i)
         input += std::to_string(results_eigen->counters[i]) + "   ";

      input += "\n";

      outfile << input;

      input = "";

      std::cout << length << " Completed" << std::endl;

      delete results_direct;
#ifdef BUCKET
      delete results_bucket;
#endif
      delete results_eigen;
      delete results_table;
      delete results_parallel;
   }

   for (int ly = 0; ly < layers + 1; ++ly)
      delete[] table_perm[ly];

   for (int ly = 0; ly < layers + 2; ++ly)
      delete[] table_in[ly];

   delete[] out;
   delete[] out_p;
   delete[] table_in;
   delete[] table_perm;

   outfile.close();

   return 0;
}
