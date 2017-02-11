//
//  best_bucket.cpp
//  project
//
//  Created by Matteo Dusefante on 19/02/16.
//  Copyright © 2016 Matteo Dusefante. All rights reserved.
//

#include "benchmark.h"

int main() {

    papils::papi_multiplex_init();
    std::ofstream out;
    out.open("output.txt", std::ios_base::app);

    int cores = std::thread::hardware_concurrency();

    int delta = 100;
    int end = 1000;
    int length = 10000000;
    int events = EVENTS;
    int layers = 1;
    int bucket1 = 10000;
    int bucket2 = 1;

    std::string input = "";

    long *out_a = new long[length]();
    long *out_p = new long[length];
    long *nb = new long[layers];

    utils::data *results_direct, *results_table, *results_parallel;

    long **table_perm = new long *[layers + 1];
    long **table_in = new long *[layers + 2];

    for(int ly = 0; ly < layers + 1; ++ly)
        table_perm[ly] = new long[length]();

    for(int ly = 0; ly < layers + 2; ++ly)
        table_in[ly] = new long[length]();

    std::cout << "Start..." << std::endl;

    input += "#start=" + std::to_string(delta) + " delta=" + std::to_string(delta) + " end=" + std::to_string(end) +
        " length=" + std::to_string(length) + " layers=" + std::to_string(layers) + "\n";

    results_direct = new utils::data();

    utils::random_permutation_table(&table_perm[0], length, layers);
    utils::populate_table(&table_in[0], length, layers);

#ifdef DEBUG
    algorithms::permute(&table_in[0][0], &out_p[0], &table_perm[0][0], length);
#endif

    // results_direct = benchmark::eigen(&table_in[0], &table_perm[0], &out_p[0], length, events);

    benchmark::omp_direct_algorithm(
        &table_in[0], &table_perm[0], &out_p[0], &out_a[0], length, events, cores, results_direct);

    for(int i = 1; i < 10 + 1; ++i) {
        for(bucket1 = delta; bucket1 < end + delta; bucket1 += delta) {

            int dt = (end - bucket1) / 10;

            int crs = i * cores;
            // int crs = i;
            // int crs = cores;

            // for(bucket2 = bucket1 + dt; bucket2 < end + dt; bucket2 += dt) {

            results_table = new utils::data();
            results_parallel = new utils::data();

            nb[0] = bucket1;
            nb[1] = bucket2; // bucket1 + i * dt; // vbs[0] * 10;

            /*
            for(int ly = 1; ly < layers; ++ly)
                vbs[ly] = vbs[ly - 1] / 10;
                 */

            // results_bucket = benchmark::omp_bucket_algorithm(
            //    &table_in[0], &table_perm[0], &vbs[0], &out_p[0], &out_a[0], length, layers, events);

            benchmark::omp_table_algorithm(
                &table_in[0], &table_perm[0], &nb[0], &out_p[0], length, layers, events, crs, results_table);

            benchmark::parallel_algorithm(&table_in[0], &table_perm[0], &nb[0], &out_p[0], &out_a[0], length, layers,
                events, crs, results_parallel);

            // input += std::to_string(nb[0]) + "   ";
            // input += std::to_string(crs) + "   ";
            // input += std::to_string(nb[0]) + "   " + std::to_string(nb[1]) + "   ";
            input += std::to_string(crs) + "   " + std::to_string(nb[0]) + "   ";

            input += std::to_string(results_direct->time) + "   ";
            for(int i = 0; i < events - 1; ++i)
                input += std::to_string(results_direct->counters[i]) + "   ";
            // for(int i = 0; i < events; ++i)
            //    input += std::to_string(results_bucket[i]) + "   ";
            input += std::to_string(results_table->time) + "   ";
            for(int i = 0; i < events - 1; ++i)
                input += std::to_string(results_table->counters[i]) + "   ";
            input += std::to_string(results_parallel->time) + "   ";
            for(int i = 0; i < events - 1; ++i)
                input += std::to_string(results_parallel->counters[i]) + "   ";

            input += "\n";

            out << input;

            input = "";

            // std::cout << length << " bucket " << bucket << " Completed\n";
            std::cout << length << " bucket1 " << nb[0] << " bucket2 " << nb[1] << " Completed\n";
            // std::cout << length << " bucket1 " << nb[0] << " Completed\n";

            delete results_table;
            delete results_parallel;
        }
    }

    delete results_direct;

    for(int ly = 0; ly < layers + 1; ++ly)
        delete[] table_perm[ly];

    for(int ly = 0; ly < layers + 2; ++ly)
        delete[] table_in[ly];

    delete[] out_a;
    delete[] out_p;
    delete[] table_in;
    delete[] table_perm;

    out.close();

    return 0;
}