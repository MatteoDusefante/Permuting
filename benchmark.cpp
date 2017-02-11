//
//  benchmark.cpp
//  project
//
//  Created by Matteo Dusefante on 31/05/16.
//  Copyright © 2016 Matteo Dusefante. All rights reserved.
//

// TO DO start=750000000 end=1000000000 delta=10000000 1000x1000000
// TO DO start=700000000 end=1000000000 delta=10000000 100x100000

#include "benchmark.h"

int main(int argc, char *argv[]) {

   // std::srand(std::time(nullptr));

   int delta = 100000000;
   int start = 100000000;
   // start = 900000000;
   // start = 750000000;
   int end = 100000000;
   // delta = 1000000;
   // end = 700000000;
   int events = EVENTS;
   int buckets = 1;
   int layers = 1;

   int nbuckets1[] = {100};
   int nbuckets2[] = {100};

   if (argc > 2) {
      start = delta = end = atoi(argv[1]);
      nbuckets1[0] = nbuckets2[0] = atoi(argv[2]);
      // delta = std::min(atoi(argv[1]), atoi(argv[2]));
      // end = std::max(atoi(argv[1]), atoi(argv[2]));
   }
   if (argc > 3)
      layers = atoi(argv[3]);

   // papils::papi_init();
   papils::papi_multiplex_init();
   std::ofstream out;
   out.open("output.txt", std::ios_base::app);

   long *out_a = new long[end];
   long *out_b = new long[end];
   long *out_d = new long[end];
   long *out_p = new long[end];
   long *tnb = new long[layers];
   long *pnb = new long[layers];

   std::string input = "";

   // for(int i = 0; i < 4 * events + 1; ++i)
   //    input += "0   ";
   // input += "\n";

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

   utils::data *results_direct, *results_bucket, *results_eigen, *results_table,
       *results_parallel;

   for (int length = start; length < end + delta; length += delta) {

      results_direct = new utils::data();
      results_bucket = new utils::data();
      results_eigen = new utils::data();
      results_table = new utils::data();
      results_parallel = new utils::data();

      utils::random_permutation_table(&table_perm[0], length, layers);
      utils::populate_table(&table_in[0], length, layers);

#ifdef DEBUG
      algorithms::permute(&table_in[0][0], &out_p[0], &table_perm[0][0],
                          length);
#endif

      // benchmark::omp_direct_algorithm(
      //    &table_in[0], &table_perm[0], &out_p[0], &out_d[0], length, events,
      //    cores, results_direct);
      // benchmark::eigen(&table_in[0], &table_perm[0], &out_p[0], length,
      // events, cores, results_direct);
      benchmark::eigen(&table_in[0], &table_perm[0], &out_p[0], length, events,
                       cores, results_eigen);
      /*
            for(int it = 0; it < buckets; ++it) {

               tnb[0] = nbuckets1[it];

               for(int ly = 1; ly < layers; ++ly)
                  tnb[ly] = cores; //tnb[ly - 1] * 100;

               //benchmark::omp_bucket_algorithm(&table_in[0], &table_perm[0],
         &tnb[0], &out_p[0], &out_b[0], length,
         layers,
               //    events, cores, results_bucket);

               benchmark::omp_table_algorithm(
                  &table_in[0], &table_perm[0], &tnb[0], &out_p[0], length,
         layers, events, cores, results_table);

               pnb[0] = nbuckets2[it];
               for(int ly = 1; ly < layers; ++ly)
                  pnb[ly] = cores; //pnb[ly - 1] * 100;

               benchmark::parallel_algorithm(&table_in[0], &table_perm[0],
         &pnb[0], &out_p[0], &out_a[0], length,
         layers,
                                             events, cores, results_parallel);
            }*/

      input += std::to_string(length) + "   ";
      input += std::to_string(results_direct->time) + "   ";
      for (int i = 0; i < events - 1; ++i)
         input += std::to_string(results_direct->counters[i]) + "   ";
      input += std::to_string(results_bucket->time) + "   ";
      for (int i = 0; i < events; ++i)
         input += std::to_string(results_bucket->counters[i]) + "   ";
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
      delete results_eigen;
      delete results_table;
      delete results_parallel;
   }

   // delete[] results_bucket;

   for (int ly = 0; ly < layers + 1; ++ly)
      delete[] table_perm[ly];

   for (int ly = 0; ly < layers + 2; ++ly)
      delete[] table_in[ly];

   delete[] out_a;
   delete[] out_b;
   delete[] out_d;
   delete[] out_p;
   delete[] table_in;
   delete[] table_perm;

   out.close();

   return 0;
}

/*
int main() {

    papils::papi_init();
    std::ofstream out;
    out.open("output.txt", std::ios_base::app);

    long long *results_naive = new long long[4];
    long long *results_bucket = new long long[4];
    long long *results_table = new long long[4];
    long long *results_preprocess = new long long[4];
    long long *results_b = new long long[4];
    long long *results_p = new long long[4];
    long long *results_t = new long long[4];

    // int bucketsize[] = { 3000, 4000, 4500, 5000, 6000, 9000, 11000, 40000,
45000, 50000 }; // plain
    // int bucketsize[] = { 500, 1000, 2500, 5000 };
    // int bucketsize[] = {5000, 10000, 50000, 100000, 500000, 1000000,
10000000};

    // int bucketsize[] = { 10, 100, 250, 350, 500, 750, 900, 1000, 2500, 5000};
// best 1000 - 100000
    // int bucketsize[] = {5000, 10000, 25000, 50000, 100000, 20000, 250000,
500000};

    int bl = 10;
    int bucketsize[] = { 100, 1000, 2500, 3500, 5000, 7500, 9000, 10000, 25000,
50000 }; //  best 100000 - 10000000

    // for(int length = 1000; length < 101000; length += 1000) { // 1
    // for(int length = 10000000; length < 1000000010; length += 10) {  // 3
    for(int length = 100000; length < 10100000; length += 100000) { // 2

        int *perm = new int[length];
        int *n_perm = new int[length];
        int *in = new int[length];
        int *mid = new int[length];

        int *out_p = new int[length];
        int *out_b = new int[length];
        int *out_t = new int[length];
        int **table = new int *[length];

        for(int it = 0; it < length; ++it)
            table[it] = new int[3];

        utils::random_permutation(&perm[0], &perm[length]);
        utils::populate(&in[0], &in[length]);



        papils::papi_start();

        utils::permute(&in[0], &out_p[0], &perm[0], length);

        papils::papi_stop();
        papils::papi_collect(&results_naive[0]);



        long long best_cost_bucket = -1;
        long long best_cost_table = -1;
        long long curr_cost = -1;

        for(int id = 0; id < bl; ++id) {

            int bs = bucketsize[id];
            if(bs > length)
                break;

            int n_buckets = (int)(floor(length / bs + 1));

            int *buckets = new int[n_buckets];

            buckets[0] = 0;
            for(int it = 1; it < n_buckets; ++it)
                buckets[it] = buckets[it - 1] + bs;



            papils::papi_start();

            utils::bucket(&in[0], &out_b[0], &perm[0], &buckets[0], &mid[0],
&n_perm[0], length, bs);

            papils::papi_stop();
            papils::papi_collect(&results_bucket[0]);



            curr_cost = compute_cost(0, results_bucket);

            if(best_cost_bucket == -1) {
                results_b = results_bucket;
                best_cost_bucket = curr_cost;
            } else {
                if(curr_cost < best_cost_bucket) {
                    results_b = results_bucket;
                    best_cost_bucket = curr_cost;
                }
            }

            int abs[] = { bs, bs / 2 };



            papils::papi_start();

            utils::multi_layer_preprocessing(&perm[0], &table[0], &abs[0],
length, 0, 2);

            papils::papi_stop();
            papils::papi_collect(&results_preprocess[0]);



            for(int it = 0; it < length; ++it)
                mid[it] = in[it];

            papils::papi_start();

            utils::multi_layer_table(&mid[0], &out_t[0], &table[0], length, 0,
2);

            papils::papi_stop();
            papils::papi_collect(&results_table[0]);



            curr_cost = compute_cost(0, results_table);

            if(best_cost_table == -1) {
                results_p = results_preprocess;
                results_t = results_table;
                best_cost_table = curr_cost;
            } else {
                if(curr_cost < best_cost_table) {
                    results_p = results_preprocess;
                    results_t = results_table;
                    best_cost_table = curr_cost;
                }
            }
        }

        utils::verify(&out_p[0], &out_b[0], length);
        utils::verify(&out_p[0], &out_t[0], length);

        std::string input;

        input += std::to_string(length) + "   ";
        // + std::to_string(results_naive[0]) + "   " +
std::to_string(results_bucket[0]) + "   ";
        for(int i = 0; i < 4; ++i)
            input += std::to_string(results_naive[i]) + "   ";
        for(int i = 0; i < 4; ++i)
            input += std::to_string(results_b[i]) + "   ";
        for(int i = 0; i < 4; ++i)
            input += std::to_string(results_p[i]) + "   ";
        for(int i = 0; i < 4; ++i)
            input += std::to_string(results_t[i]) + "   ";
        input += "\n";

        out << input;

        std::cout << length << " Completed " << std::endl;
    }

    out.close();

    return 0;
}
*/

/*
if(best_cost == -1) {
               time = elapsed_bucket;
               results = results_bucket;
               best_cost = 0;
               best_bucket = bs;
           } else {
               if(results_bucket[0] < results_naive[0] && results_bucket[1] <
results_naive[1] &&
                   results_bucket[2] < results_naive[2]) {
                   if(results_bucket[0] < results[0] && results_bucket[1] <
results[1] &&
                       results_bucket[2] < results[2]) {
                       time = elapsed_bucket;
                       results = results_bucket;
                       best_bucket = bs;
                   }
               }
           }*/
