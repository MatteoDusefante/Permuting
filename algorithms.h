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
         // printf("--- %ld  %ld  %ld  %ld\n", c, table_perm[ly][it],
         // bucket_size[ly + 1], bucket);
      }

      /*        for(it = 0; it < length; ++it)
                  printf(" %ld ", table_perm[ly][it]);
              printf("\n");

              for(it = 0; it < buckets[ly] * cores; ++it)
                  printf(" %ld ", table_pt[ly][it]);
              printf("\n\n");*/

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

      // table_pt[ly][0] = 0;

      /*        printf("\n\n");

              for(it = 0; it < buckets[ly] * cores; ++it)
                  printf(" %ld ", table_pt[ly][it]);
              printf("\n");*/

      /*

      T *tp = new T[buckets[ly] * cores];
      std::copy(&table_pt[ly][0], &table_pt[ly][0] + buckets[ly] * cores,
      &tp[0]);
      size_t rgh, lft, bs, deltaf, index;
      for(size_t bk = 0; bk < (size_t)nbuckets[ly]; ++bk) {

          if(bk == (size_t)nbuckets[ly] - 1)
              bs = length - bucket_size[ly] * (nbuckets[ly] - 1);

          deltaf = bs - (delta * (cores - 1));

          for(size_t c = 0; c < cores; ++c) {
              lft = bucket_size[ly] * bk + c * delta;
              if(c != cores - 1)
                  rgh = lft + delta;
              else
                  rgh = lft + deltaf;
              for(size_t it = lft; it < rgh; ++it) {
                  bucket = floor((table_perm[ly][it]) / bucket_size[ly + 1]);
                  // index = tp[c]++;
                  index = tp[bucket * cores + c]++;
                  table_perm[ly + 1][index] = table_perm[ly][it];
              }
          }
      }
      delete[] tp;
       */

      /*
      size_t crs = cores, lp, rp, rounds, temp = 0;
      for(size_t b = 0; b < buckets[ly]; ++b) {
          table_pt[ly][b] += temp;
          for(rounds = 1; rounds < cores; rounds *= 2) {
              crs = ceil((double)cores / (rounds * 2));
              printf("P  %ld, D  %ld\n", crs, rounds);
              omp_set_dynamic(0);
#pragma omp parallel for private(lp, rp, c) shared(table_pt) num_threads(crs)
              for(c = 1; c < crs + 1; ++c) {
                  lp = rounds * 2 * c - rounds - 1;
                  rp = std::min(lp + rounds, cores - 1);
                  printf("%ld,%ld   ", lp, rp);
                  if(rp != lp)
                      temp += table_pt[ly][buckets[ly] * rp + b] +=
table_pt[ly][buckets[ly] * lp + b];
              }
              printf("\n");
          }
          printf("-------------------------- %ld\n", b);
          crs = 1;
          for(rounds = cores; rounds > 0; rounds = ((double) ceil(rounds / 2)))
{
              printf("P  %ld, D  %ld\n", crs, rounds);
              omp_set_dynamic(0);
#pragma omp parallel for private(lp, rp, c) shared(table_pt) num_threads(crs)
              for(c = 1; c < crs + 1; ++c) {
                  lp = rounds * c - ((double) ceil(rounds / 2)) - 1;
                  rp = std::min(lp + rounds, cores - 2);
                  printf("%ld,%ld   ", lp, rp);
                  if((rp - lp) >= 0)
                      table_pt[ly][buckets[ly] * rp + b] +=
table_pt[ly][buckets[ly] * lp + b];
              }
              crs *= 2;
              printf("\n");
          }
          printf("%ld\n", b);
      }*/
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

/******************************************************************************************/

/*
template <typename T>
inline void parallel_bucket_preprocessing(T **table_perm,
    T **table_pt,
    T *nbuckets,
    T *bucket_size,
    T *pointers_length,
    size_t cores,
    size_t length,
    size_t layers) {

    T index, bucket;
    size_t lft, rgh, delta, deltaf, bs;

    printf("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");

    for(size_t it = 0; it < length; ++it)
        printf("%ld \t", table_perm[0][it]);
    printf("\n");

    for(size_t ly = 0; ly < layers; ++ly) {
        bs = bucket_size[ly];
        T *buckets = new T[nbuckets[ly + 1]];
        buckets[0] = 0;
        for(size_t it = 1; it < (size_t)nbuckets[ly + 1]; ++it)
            buckets[it] = buckets[it - 1] + bucket_size[ly + 1];
        for(size_t bk = 0; bk < (size_t)nbuckets[ly]; ++bk) {

            if(bk == (size_t)nbuckets[ly] - 1)
                bs = length - bucket_size[ly] * (nbuckets[ly] - 1);

            // delta = ceil(bs / cores);
            delta = floor(bs / cores) + (((bs / cores) * cores) != bs);
            deltaf = bs - (delta * (cores - 1));
            printf("LENGTH %ld   BS %ld   DELTA %ld   DELTAF %ld   NB  %ld\n",
length, bs, delta, deltaf,
nbuckets[ly]);

            for(size_t c = 0; c < cores; ++c) {
                lft = bucket_size[ly] * bk + c * delta;
                if(c != cores - 1)
                    rgh = lft + delta;
                else
                    rgh = lft + deltaf;
                printf("L %ld   R %ld   CORE %ld\n", lft, rgh, c);
                for(size_t it = lft; it < rgh; ++it) {
                    bucket = floor((table_perm[ly][it]) / bucket_size[ly + 1]);
                    // index = buckets[bucket]++;
                    table_pt[ly][bucket * cores + c + 1]++;
                    // table_perm[ly + 1][index] = table_perm[ly][it];
                    printf("BUCKET %ld   INDEX %ld   TP %ld   BS %ld   BCC %ld =
%ld * %ld + %ld + 1\n", bucket, it,
                        table_perm[ly + 1][it], bucket_size[ly + 1], bucket *
cores + c + 1, bucket, cores, c);
                }
            }
        }
        delete[] buckets;

        for(size_t it = 1; it < (size_t)pointers_length[ly]; ++it)
            table_pt[ly][it] = table_pt[ly][it] + table_pt[ly][it - 1];

        bs = bucket_size[ly];
        T *tp = new T[pointers_length[ly]];
        std::copy(&table_pt[ly][0], &table_pt[ly][0] + pointers_length[ly],
&tp[0]);
        for(size_t bk = 0; bk < (size_t)nbuckets[ly]; ++bk) {

            if(bk == (size_t)nbuckets[ly] - 1)
                bs = length - bucket_size[ly] * (nbuckets[ly] - 1);

            delta = floor(bs / cores) + (((bs / cores) * cores) != bs);
            deltaf = bs - (delta * (cores - 1));

            for(size_t c = 0; c < cores; ++c) {
                lft = bucket_size[ly] * bk + c * delta;
                if(c != cores - 1)
                    rgh = lft + delta;
                else
                    rgh = lft + deltaf;
                for(size_t it = lft; it < rgh; ++it) {
                    bucket = floor((table_perm[ly][it]) / bucket_size[ly + 1]);
                    // index = tp[c]++;
                    index = tp[bucket * cores + c]++;
                    table_perm[ly + 1][index] = table_perm[ly][it];
                }
            }
        }
        delete[] tp;
    }
    for(size_t ly = 0; ly < layers; ++ly) {
        for(size_t it = 0; it < (size_t)pointers_length[ly]; ++it)
            printf("%ld \t", table_pt[ly][it]);
        printf("\n");
    }

    printf("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
}
*/
/*
template <typename T>
inline void parallel_bucket_preprocessing(T **table_perm,
    T **table_pt,
    T *nbuckets,
    T *bucket_size,
    T *pointers_length,
    size_t cores,
    size_t length,
    size_t layers) {

    T index, bucket;
    size_t lft, rgh, delta, deltaf, bs;


    for(size_t ly = 0; ly < layers; ++ly) {
        bs = bucket_size[ly];
        T *buckets = new T[nbuckets[ly + 1]];
        buckets[0] = 0;
        for(size_t it = 1; it < (size_t)nbuckets[ly + 1]; ++it)
            buckets[it] = buckets[it - 1] + bucket_size[ly + 1];
        for(size_t bk = 0; bk < (size_t)nbuckets[ly]; ++bk) {
            if(bk == (size_t)nbuckets[ly] - 1)
                bs = length - bucket_size[ly] * (nbuckets[ly] - 1);

            delta = floor(bs / cores);
            deltaf = bs - (delta * (cores - 1));
            for(size_t c = 0; c < cores; ++c) {
                lft = bucket_size[ly] * bk + c * delta;
                if(c != cores - 1)
                    rgh = lft + delta;
                else
                    rgh = lft + deltaf;
                for(size_t it = lft; it < rgh; ++it) {
                    bucket = std::min(nbuckets[ly + 1] -
1,(T)floor((table_perm[ly][it]) / bucket_size[ly + 1]));
                    index = buckets[bucket]++;
                    table_pt[ly][bucket * cores + c + 1]++;
                    table_perm[ly + 1][index] = table_perm[ly][it];
                }
            }
        }

        delete[] buckets;

        for(size_t it = 1; it < (size_t)pointers_length[ly]; ++it)
            table_pt[ly][it] = table_pt[ly][it] + table_pt[ly][it - 1];
    }

}
*/

template <typename T>
inline void parallel_bucket_preprocessing1(T **table_perm, T **table_pt,
                                           T *nbuckets, T *bucket_size,
                                           T *pointers_length, size_t cores,
                                           size_t length, size_t layers) {

   T index, bucket;
   size_t lft, rgh, delta, deltaf, bs;

   for (size_t ly = 0; ly < layers; ++ly) {
      bs = bucket_size[ly];
      T *buckets = new T[nbuckets[ly + 1]];
      buckets[0] = 0;
      for (size_t it = 1; it < (size_t)nbuckets[ly + 1]; ++it)
         buckets[it] = buckets[it - 1] + bucket_size[ly + 1];
      for (size_t bk = 0; bk < (size_t)nbuckets[ly]; ++bk) {
         if (bk == (size_t)nbuckets[ly] - 1)
            bs = length - bucket_size[ly] * (nbuckets[ly] - 1);
         delta = floor(bs / cores) + (((bs / cores) * cores) != bs);
         deltaf = bs - (delta * (cores - 1));
         for (size_t c = 0; c < cores; ++c) {
            lft = bucket_size[ly] * bk + c * delta;
            if (c != cores - 1)
               rgh = lft + delta;
            else
               rgh = lft + deltaf;
            for (size_t it = lft; it < rgh; ++it) {
               bucket = floor((table_perm[ly][it]) / bucket_size[ly + 1]);
               table_pt[ly][bucket * cores + c + 1]++;
            }
         }
      }
      delete[] buckets;

      for (size_t it = 1; it < (size_t)pointers_length[ly]; ++it)
         table_pt[ly][it] = table_pt[ly][it] + table_pt[ly][it - 1];

      bs = bucket_size[ly];
      T *tp = new T[pointers_length[ly]];
      std::copy(&table_pt[ly][0], &table_pt[ly][0] + pointers_length[ly],
                &tp[0]);
      for (size_t bk = 0; bk < (size_t)nbuckets[ly]; ++bk) {

         if (bk == (size_t)nbuckets[ly] - 1)
            bs = length - bucket_size[ly] * (nbuckets[ly] - 1);

         delta = floor(bs / cores) + (((bs / cores) * cores) != bs);
         deltaf = bs - (delta * (cores - 1));

         for (size_t c = 0; c < cores; ++c) {
            lft = bucket_size[ly] * bk + c * delta;
            if (c != cores - 1)
               rgh = lft + delta;
            else
               rgh = lft + deltaf;
            for (size_t it = lft; it < rgh; ++it) {
               bucket = floor((table_perm[ly][it]) / bucket_size[ly + 1]);
               index = tp[bucket * cores + c]++;
               table_perm[ly + 1][index] = table_perm[ly][it];
            }
         }
      }
      delete[] tp;
   }
}

/******************************************************************************************/

template <typename T>
inline void parallel_bucket1(T **table_in, T **table_perm, T **table_pt, T *out,
                             T *nbuckets, T *bucket_size, size_t cores,
                             size_t length, size_t layers) {

   size_t c, it, index, lft, rgh, bs, delta, deltaf;

   T bucket;

   omp_set_dynamic(0);
   for (size_t ly = 0; ly < layers; ++ly) {
      bs = bucket_size[ly];
      for (size_t bk = 0; bk < (size_t)nbuckets[ly]; ++bk) {
         if (bk == (size_t)nbuckets[ly] - 1)
            bs = length - bucket_size[ly] * (nbuckets[ly] - 1);
         delta = floor(bs / cores) + (((bs / cores) * cores) != bs);
         deltaf = bs - (delta * (cores - 1));
#pragma omp parallel for private(c, it, index, lft, rgh,                       \
                                 bucket) shared(table_pt) num_threads(cores)
         for (c = 0; c < cores; ++c) {
            lft = bucket_size[ly] * bk + c * delta;
            if (c != cores - 1)
               rgh = lft + delta;
            else
               rgh = lft + deltaf;
            for (it = lft; it < rgh; ++it) {
               bucket = floor((table_perm[ly][it]) / bucket_size[ly + 1]);
               index = table_pt[ly][bucket * cores + c]++;
               table_in[ly + 1][index] = table_in[ly][it];
            }
         }
      }
   }

#pragma omp parallel for private(it) num_threads(cores)
   for (size_t it = 0; it < length; ++it)
      out[table_perm[layers][it]] = table_in[layers][it];
   return;
}

/*
template <typename T>
inline void parallel_bucket(T **table_in,
    T **table_perm,
    T **table_pt,
    T *out,
    T *nbuckets,
    T *bucket_size,
    size_t cores,
    size_t length,
    size_t layers) {

    size_t c, it, index, lft, rgh, bs, delta, deltaf;

    T bucket;

    omp_set_dynamic(0);
    for(size_t ly = 0; ly < layers; ++ly) {
        bs = bucket_size[ly];
        for(size_t bk = 0; bk < (size_t)nbuckets[ly]; ++bk) {
            if(bk == (size_t)nbuckets[ly] - 1)
                bs = length - bucket_size[ly] * (nbuckets[ly] - 1);
            delta = floor(bs / cores);
            deltaf = bs - (delta * (cores - 1));
#pragma omp parallel for private(c, it, index, lft, rgh, bucket)
shared(table_pt) num_threads(cores)
            for(c = 0; c < cores; ++c) {
                lft = bucket_size[ly] * bk + c * delta;
                if(c != cores - 1)
                    rgh = lft + delta;
                else
                    rgh = lft + deltaf;
                for(it = lft; it < rgh; ++it) {
                    bucket = std::min(nbuckets[ly + 1] -
1,(T)floor((table_perm[ly][it]) / bucket_size[ly + 1]));
                    //bucket = (T)floor((table_perm[ly][it]) / bucket_size[ly +
1]);
                    index = table_pt[ly][bucket * cores + c]++;
                    table_in[ly + 1][index] = table_in[ly][it];
                    table_perm[ly + 1][index] = table_perm[ly][it];
                }
            }
        }
    }
#pragma omp parallel for private(it) num_threads(cores)
    for(size_t it = 0; it < length; ++it)
        out[table_perm[layers][it]] = table_in[layers][it];
    return;
}
*/

/*
template <typename T>
inline void parallel_bucket(T **table_in,
    T **table_perm,
    T **table_pt,
    T *out,
    T *nbuckets,
    T *bucket_size,
    size_t cores,
    size_t length,
    size_t layers) {

    size_t c, it, index, lft, rgh, bs, delta, deltaf;

    T bucket;

    omp_set_dynamic(0);
    for(size_t ly = 0; ly < layers; ++ly) {
        for(size_t it = 0; it < length; ++it)
            printf("%ld \t", table_perm[ly][it]);
        printf("\n");
        for(size_t it = 0; it < length; ++it)
            printf("%ld \t", table_in[ly][it]);
        printf("\n");
        bs = bucket_size[ly];
        for(size_t bk = 0; bk < (size_t)nbuckets[ly]; ++bk) {
            if(bk == (size_t)nbuckets[ly] - 1)
                bs = length - bucket_size[ly] * (nbuckets[ly] - 1);
            // delta = (T)floor(bucket_size[ly] / cores);
            // delta = ceil(bs / cores);
            delta = floor(bs / cores) + (((bs / cores) * cores) != bs);
            deltaf = bs - (delta * (cores - 1));
            printf("LENGTH %ld   BS %ld   DELTA %ld   DELTAF %ld   NB  %ld\n",
length, bs, delta, deltaf,
nbuckets[ly]);
#pragma omp parallel for private(c, it, index, lft, rgh, bucket)
shared(table_pt) num_threads(cores)
            for(c = 0; c < cores; ++c) {
                lft = bucket_size[ly] * bk + c * delta;
                if(c != cores - 1)
                    rgh = lft + delta;
                else
                    rgh = lft + deltaf;
                printf("L %ld   R %ld   CORE %ld\n", lft, rgh, c);
                for(it = lft; it < rgh; ++it) {
                    // bucket = std::min(nbuckets[ly + 1] - 1,
(T)floor((table_perm[ly][it]) / bucket_size[ly +
1]));
                    bucket = floor((table_perm[ly][it]) / bucket_size[ly + 1]);
                    // printf("BUCKET  %ld  = %ld / %ld \n", bucket,
table_perm[ly][it], bucket_size[ly + 1]);
                    index = table_pt[ly][bucket * cores + c]++;
                    // for(size_t ly = 0; ly < layers; ++ly) {
                    //    for(size_t it = 0; it < cores * bucket_size[ly]; ++it)
                    //        printf("%ld \t", table_pt[ly][it]);
                    //    printf("\n");
                    //}

                    table_in[ly + 1][index] = table_in[ly][it];
                    //table_perm[ly + 1][index] = table_perm[ly][it];
                    printf("BUCKET %ld   INDEX %ld   TP %ld   BS %ld   BCC %ld =
%ld * %ld + %ld + 1\n", bucket,
index,
                        table_perm[ly + 1][index], bucket_size[ly + 1], bucket *
cores + c + 1, bucket, cores, c);
                }
            }
#pragma omp barrier
        }

        printf("=====\n");
    }

    for(size_t it = 0; it < length; ++it)
        printf("%ld \t", table_perm[layers][it]);
    printf("\n\n");
    for(size_t it = 0; it < length; ++it)
        printf("%ld \t", table_in[layers][it]);
    printf("\n\n");
#pragma omp parallel for private(it) num_threads(cores)
    for(size_t it = 0; it < length; ++it)
        out[table_perm[layers][it]] = table_in[layers][it];
    return;
}
*/
/*
template <typename T>
inline void parallel_bucket(T **table_in,
    T **table_perm,
    T **table_pt,
    T *out,
    T *bucket_size,
    size_t cores,
    size_t length,
    size_t layers) {

    size_t c, it, lft, rgh, delta, deltaf, index;

    T bucket;

    delta = floor(length / cores) + (((length / cores) * cores) != length);
    deltaf = length - (delta * (cores - 1));

    omp_set_dynamic(0);
    for(size_t ly = 0; ly < layers; ++ly) {
#pragma omp parallel for private(c, it, index, lft, rgh, bucket)
shared(table_pt) num_threads(cores)
schedule(static)
        for(c = 0; c < cores; ++c) {
            lft = c * delta;
            if(c != cores - 1)
                rgh = lft + delta;
            else
                rgh = lft + deltaf;
#pragma omp parallel for private(index, bucket) num_threads(cores)
            for(it = lft; it < rgh; ++it) {
                bucket = floor(table_perm[ly][it] / bucket_size[ly + 1]);
                index = table_pt[ly][bucket * cores + c]++;
                table_in[ly + 1][index] = table_in[ly][it ];
                // table_perm[ly + 1][index] = table_perm[ly][it];
            }
        }
    }
#pragma omp parallel for private(it) num_threads(cores)
    for(size_t it = 0; it < length; ++it)
        out[table_perm[layers][it]] = table_in[layers][it];
    return;
     *
     *
     * template <typename T>
inline void parallel_bucket_preprocessing(T **table_perm,
    T **table_pt,
    T *bucket_size,
    T *pointers_length,
    T *buckets,
    size_t cores,
    size_t length,
    size_t layers,
    size_t delta) {

    T index, bucket;
    size_t lft, rgh, deltaf;

    size_t c, it;

    omp_set_dynamic(0);
    omp_set_num_threads(cores);

    for(size_t ly = 0; ly < layers; ++ly) {
// for(c = 0; c < cores; ++c) {
*
lft = c * delta;
if(c != cores - 1)
rgh = lft + delta;
else
rgh = lft + deltaf;
#pragma omp parallel for private(c, it, bucket) shared(table_pt)
schedule(static, delta)
        for(it = 0; it < length; ++it) {
            c = omp_get_thread_num();
            // for(it = lft; it < rgh; ++it) {
            bucket = floor(table_perm[ly][it] / bucket_size[ly + 1]);
            // table_pt[ly][bucket * cores + c + 1]++;
            table_pt[ly][c * buckets[ly] + bucket + 1]++;
            // index = table_pt[ly][bucket * cores + c]++;
            // table_perm[ly + 1][index] = table_perm[ly][it];  // NO GOOD HERE!
        }
        //}

        *
        T *tp = new T[pointers_length[ly] + 1];
        std::copy(&table_pt[ly][0], &table_pt[ly][0] + pointers_length[ly] + 1,
&tp[0]);

        // for(it = 1; it < pointers_length[ly]; ++it)
        //    table_pt[ly][it] = table_pt[ly][it] + table_pt[ly][it - 1];

        for(c = 0; c < cores; ++c)
            for(auto b = 0; b < buckets[ly]; ++b)
                printf("%d ", table_pt[ly][c * buckets[ly] + b]);
        printf("\t\t");
        printf("\n");

        size_t temp = 0;
        for(int b = 0; b < buckets[ly]; ++b) {
            table_pt[ly][b] += temp;
            for(int c = 1; c < cores; ++c) {
                table_pt[ly][buckets[ly] * c + b] += table_pt[ly][buckets[ly] *
(c - 1) + b];
            }
            temp = table_pt[ly][buckets[ly] * (cores - 1) + b];
        }

        for(c = 0; c < cores; ++c)
            for(auto b = 0; b < buckets[ly]; ++b)
                printf("%d ", table_pt[ly][c * buckets[ly] + b]);
        printf("\t\t");
        printf("\n");


        *

        for(c = 0; c < cores; ++c)
            for(auto b = 0; b < buckets[ly]; ++b)
                table_pt[ly][c * buckets[ly] + b] = tp[b * cores + c];


        // std::copy(&table_pt[ly][0], &table_pt[ly][0] + pointers_length[ly],
&tp[0]);

        for(c = 0; c < cores; ++c) {
            lft = c * delta;
            if(c != cores - 1)
                rgh = lft + delta;
            else
                rgh = lft + deltaf;
            for(size_t it = lft; it < rgh; ++it) {
                bucket = floor((table_perm[ly][it]) / bucket_size[ly + 1]);
                // index = tp[c * buckets[ly] + bucket]++;
                // index = tp[bucket * cores + c]++;
                index = tp[buckets[ly] * c + bucket]++;
                //table_perm[ly + 1][index] = table_perm[ly][it];
            }
        }
        delete[] tp;*
    }
}
}*/

/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/
/******************************************************************************************/

template <typename T, typename P>
inline void permute(T *in, T *out, P *perm, size_t length) {

   for (size_t it = 0; it < length; ++it)
      out[perm[it]] = in[it];
   return;
}

/******************************************************************************************/

template <typename T>
inline void fast_bucketize(T **table_in, T **table_perm, T **buckets, T *out,
                           T *bucket_size, size_t length, size_t layers) {

   T bucket, index;
   for (size_t ly = 0; ly < layers; ++ly) {
      for (size_t it = 0; it < length; ++it) {
         bucket = (T)floor(table_perm[ly][it] / bucket_size[ly]);
         index = buckets[ly][bucket]++;
         table_in[ly + 1][index] = table_in[ly][it];
         table_perm[ly + 1][index] = table_perm[ly][it];
      }
   }

   for (size_t it = 0; it < length; ++it)
      out[table_perm[layers][it]] = table_in[layers][it];
   return;
}

/******************************************************************************************/

template <typename T>
inline void fast_multi_layer_table(T **table_out, T **table, size_t length,
                                   size_t layers) {

   for (size_t ly = 0; ly < layers + 1; ++ly)
      for (size_t it = 0; it < length; ++it)
         table_out[ly + 1][table[ly][it]] = table_out[ly][it];

   return;
}
}
