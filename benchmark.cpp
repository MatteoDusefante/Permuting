//
//  benchmark.cpp
//  project
//
//  Created by Matteo Dusefante on 31/05/16.
//  Copyright © 2016 Matteo Dusefante. All rights reserved.
//

#include "benchmark.h"

int main(int argc, char *argv[]) {

   int start = 100000000;
   int end = 100000000;
   int delta = start; // end - start / 10;
   int events = EVENTS;
   int buckets = 1;
   int layers = 1;

   int nbuckets1[] = {10000};
   int nbuckets2[] = {10000};

   if (argc > 2) {
      start = atoi(argv[1]);
      end = atoi(argv[2]);
      delta = end - start / 10;
      nbuckets1[0] = nbuckets2[0] = atoi(argv[3]);
   }

   // papils::papi_init();
   papils::papi_multiplex_init();
   std::ofstream out;
   out.open("output.txt", std::ios_base::app);

   long *out_a = new long[end];
#ifdef BUCKET
   long *out_b = new long[end];
#endif
   long *out_d = new long[end];
   long *out_p = new long[end];
   long *tnb = new long[layers];
   long *pnb = new long[layers];

   std::string input = "";

   long **table_perm = new long *[layers + 1]();
   long **table_in = new long *[layers + 2]();

   for (int ly = 0; ly < layers + 1; ++ly)
      table_perm[ly] = new long[end]();

   for (int ly = 0; ly < layers + 2; ++ly)
      table_in[ly] = new long[end]();

   std::cout << "Start..." << std::endl;

   input += "#start=" + std::to_string(delta) + " delta=" +
            std::to_string(delta) + " end=" + std::to_string(end) + " layers=" +
            std::to_string(layers) + " bt=" + std::to_string(nbuckets1[0]) +
            " bp=" + std::to_string(nbuckets2[0]) + "\n";

   int cores = std::thread::hardware_concurrency();

#ifdef BUCKET
   utils::data *results_direct, *results_bucket, *results_eigen, *results_table,
       *results_parallel;
#else
   utils::data *results_direct, *results_eigen, *results_table,
       *results_parallel;
#endif

   for (int length = start; length < end + delta; length += delta) {

      results_direct = new utils::data();
#ifdef BUCKET
      results_bucket = new utils::data();
#endif
      results_eigen = new utils::data();
      results_table = new utils::data();
      results_parallel = new utils::data();

      utils::random_permutation_table(&table_perm[0], length, layers);
      utils::populate_table(&table_in[0], length, layers);

#ifdef DEBUG
      algorithms::permute(&table_in[0][0], &out_p[0], &table_perm[0][0],
                          length);
#endif

      benchmark::omp_direct_algorithm(&table_in[0], &table_perm[0], &out_p[0],
                                      &out_d[0], length, events, cores,
                                      results_direct);
      benchmark::eigen(&table_in[0], &table_perm[0], &out_p[0], length, events,
                       cores, results_eigen);

      for (int it = 0; it < buckets; ++it) {

         tnb[0] = nbuckets1[it];

         for (int ly = 1; ly < layers; ++ly)
            tnb[ly] = cores; // tnb[ly - 1] * 100;

#ifdef BUCKET
         benchmark::omp_bucket_algorithm(&table_in[0], &table_perm[0], &tnb[0],
                                         &out_p[0], &out_b[0], length, layers,
                                         events, cores, results_bucket);
#endif

         benchmark::omp_table_algorithm(&table_in[0], &table_perm[0], &tnb[0],
                                        &out_p[0], length, layers, events,
                                        cores, results_table);

         pnb[0] = nbuckets2[it];
         for (int ly = 1; ly < layers; ++ly)
            pnb[ly] = cores; // pnb[ly - 1] * 100;

         benchmark::parallel_algorithm(&table_in[0], &table_perm[0], &pnb[0],
                                       &out_p[0], &out_a[0], length, layers,
                                       events, cores, results_parallel);
      }

      input += std::to_string(length) + "   ";
      input += std::to_string(results_direct->time) + "   ";
      for (int i = 0; i < events - 1; ++i)
         input += std::to_string(results_direct->counters[i]) + "   ";
#ifdef BUCKET
      input += std::to_string(results_bucket->time) + "   ";
      for (int i = 0; i < events; ++i)
         input += std::to_string(results_bucket->counters[i]) + "   ";
#endif
      input += std::to_string(results_table->time) + "   ";
      for (int i = 0; i < events - 1; ++i)
         input += std::to_string(results_table->counters[i]) + "   ";
      input += std::to_string(results_parallel->time) + "   ";
      for (int i = 0; i < events - 1; ++i)
         input += std::to_string(results_parallel->counters[i]) + "   ";
      input += std::to_string(results_eigen->time) + "   ";
      for (int i = 0; i < events - 1; ++i)
         input += std::to_string(results_eigen->counters[i]) + "   ";

      input += "\n";

      out << input;

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

   delete[] out_a;
#ifdef BUCKET
   delete[] out_b;
#endif
   delete[] out_d;
   delete[] out_p;
   delete[] table_in;
   delete[] table_perm;

   out.close();

   return 0;
}
