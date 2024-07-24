/*
 *  get_obsdatastats.c
 *  MuLoNeTests
 *
 *  Created by sebas on Mon Feb 24 2003.
 *
 */

#include "get_obsdatastats.h"

#include "zutil.h"
#include "log.h"
int get_obsstats(
    FILE *file_output,
    SGZip *file_output_gz,
    FILE *file_mask,
    // FILE *file_logerr,SGZip *file_logerr_gz,
    int n_samp, 
    long int n_site,
    long int *n_realsite, 
    char **names, 
    char *DNA_matr, 
    double *matrix_sizepos,
    double *matrix_segrpos, 
    char **matrix_pol, 
    long int **matrix_freq,
    long int **matrix_pos, 
    double *length_al, 
    long int *length_seg,
    // int *vint_perpop_nsam, 
    // int npops, 
    double *svratio, 
    double *missratio,
    // int include_unknown, 
    double *sum_sam, 
    double **tcga, 
    long int **matrix_sv,
    long int *nmhits, 
    // int output, 
    // char *ploidy, 
    // int outgroup_presence,
    double *nsites1_pop, 
    double *nsites1_pop_outg,
    double *nsites2_pop, 
    double *nsites2_pop_outg,
    double *nsites3_pop, 
    double *nsites3_pop_outg,
    double *anx, 
    double *bnx,
    double *anxo, 
    double *bnxo, 
    double **lengthamng, 
    double **lengthamng_outg,
    long int *mhitbp, 
    char **matrix_pol_tcga, 
    long int formatfile,
    mstatspop_args_t *args)
{
  /*int *unic = 0;*/
  /*long int mhits = 0;*/
  /*long int *mhitbp = 0;*/
  long int maxnsite = 128;
  long int maxbialsites = 256;
  long int maxnsamp = 32;
  int *nnsam = 0;
  /* int nnnsam; */
  long int *unic;

  long int transitions;
  long int transversions;

  double algsites;
  long int bial_sites = 0;

  int m_0, k_0, m_1, k_1;
  int v, w, x, y, z, mis, p, q, k, v2, y2, z2, xm;
  int a, aloutg, al1, al2, fal1, fal2, fal[5];

  long int xx;

  /*int v0,v1,y0,y1,b,c;*/
  /*char *strin;*/
  int d;
  double _sites;

  long int *mvbp;
  long int mv = 0;
  int nsamtot;
  /*long int totalmis=0;*/
  long int mis2, miso;
  /*int alleles,fao[5];*/
  int al3, fal3, al4, fal4;

  nsamtot = 0;
  for (x = 0; x < args->npops; x++)
    nsamtot += args->vint_perpop_nsam[x];

  /* Two pointers indicating the samples that are current samples and the outgroup samples */
  if ((nnsam = (int *)calloc(maxnsamp, sizeof(int))) == 0)
  {
    // fprintf(file_logerr,"Error: memory not reallocated. get_obsstats.1");
    log_fatal("Error: memory not reallocated, nnsam get_obsstats.1");
    return (0);
  }
  /* matrix of polymorphisms: only 0 and 1 */
  if ((*matrix_pol = (char *)calloc(maxnsamp * maxbialsites, sizeof(char))) == 0)
  {
    // fprintf(file_logerr,"Error: memory not reallocated. get_obsstat.3");
    log_fatal("Error: memory not reallocated, matrix_pol get_obsstat.3");
    return (0);
  }
  if ((*matrix_pol_tcga = (char *)calloc(maxnsamp * maxbialsites, sizeof(char))) == 0)
  {
    // fprintf(file_logerr,"Error: memory not reallocated. get_obsstat.3");
    log_fatal("Error: memory not reallocated, matrix_pol_tcga get_obsstat.3");
    return (0);
  }
  /* indicates the position and the frequency */
  if ((*matrix_pos = (long int *)calloc(maxbialsites, sizeof(long int))) == 0)
  {
    // fprintf(file_logerr,"Error: memory not reallocated. get_obsstat.4");
    log_fatal("Error: memory not reallocated, matrix_pos get_obsstat.4");
    return (0);
  }
  if ((*matrix_freq = (long int *)calloc(maxbialsites, sizeof(long int))) == 0)
  {
    // fprintf(file_logerr,"Error: memory not reallocated. get_obsstat.5");
    log_fatal("Error: memory not reallocated, matrix_freq get_obsstat.5");
    return (0);
  }
  if ((*matrix_sv = (long int *)calloc(maxbialsites, sizeof(long int))) == 0)
  {
    // fprintf(file_logerr,"Error: memory not reallocated. get_obsstat.7");
    log_fatal("Error: memory not reallocated, matrix_sv get_obsstat.7");
    return (0);
  }
  if ((unic = (long int *)calloc(maxnsamp, sizeof(long int))) == 0)
  {
    // fprintf(file_logerr,"Error: memory not reallocated. get_obsstat.6");
    log_fatal("Error: memory not reallocated, unic get_obsstat.6");
    return (0);
  }
  /* indicates the position of the mhits */
  /*
  if((mhitbp = (long int *) calloc (maxnsite, sizeof(long int))) == 0) {
    fprintf(file_logerr,"Error: memory not reallocated. get_obsstat.6");
    return(0);
  }
  */
  if ((mvbp = (long int *)calloc(maxnsite, sizeof(long int))) == 0)
  {
    // fprintf(file_logerr,"Error: memory not reallocated. get_obsstat.13");
    log_fatal("Error: memory not reallocated, mvbp get_obsstat.13");
    return (0);
  }

  if (n_samp > maxnsamp)
  {
    /* Reallocation in case the value be larger than specified */
    if ((*matrix_pol = (char *)realloc(*matrix_pol, (n_samp * maxbialsites) * sizeof(char))) == 0)
    {
      // fprintf(file_logerr,"Error: memory not reallocated. get_obsstat.7b");
      log_fatal("Error: memory not reallocated, matrix_pol get_obsstat.7b");
      return (0);
    }
    if ((*matrix_pol_tcga = (char *)realloc(*matrix_pol_tcga, (n_samp * maxbialsites) * sizeof(char))) == 0)
    {
      // fprintf(file_logerr,"Error: memory not reallocated. get_obsstat.7b");
      log_fatal("Error: memory not reallocated, matrix_pol_tcga get_obsstat.7b");
      return (0);
    }
    if ((nnsam = (int *)realloc(nnsam, n_samp * sizeof(int))) == 0)
    {
      // fprintf(file_logerr,"Error: memory not reallocated. get_obsstats.9");
      log_fatal("Error: memory not reallocated, nnsam get_obsstats.9");
      return (0);
    }
    if ((unic = (long int *)realloc(unic, n_samp * sizeof(long int))) == 0)
    {
      // fprintf(file_logerr,"Error: memory not reallocated. get_obsstats.11");
      log_fatal("Error: memory not reallocated, unic get_obsstats.11");
      return (0);
    }
    maxnsamp = n_samp;
  }
  if (n_site > maxnsite)
  {
    /*
    if((mhitbp = (long int *) realloc (mhitbp,n_site*sizeof(long int))) == 0) {
      fprintf(file_logerr,"Error: memory not reallocated. get_obsstat.8");
      return(0);
    }
    */
    if ((mvbp = (long int *)realloc(mvbp, n_site * sizeof(long int))) == 0)
    {
      // fprintf(file_logerr,"Error: memory not reallocated. get_obsstat.12");
      log_fatal("Error: memory not reallocated, mvbp get_obsstat.12");
      return (0);
    }
    maxnsite = n_site;
  }

  algsites = 0.;
  for (xx = 0; xx < n_site; xx++)
  {
    algsites += (double)matrix_sizepos[xx]; /* number of effective aligned positions, excluding gaps and mhits */
  }
  /*algsites = (double)n_site; */
  /* nnnsam = 0; */

  /* calculate number of samples in outgroup and in the current sample */
  /* if no outgroup, then all sequences are samples */
  if (n_samp < 2)
  {
    // fprintf(file_logerr," n_samples: %d .",n_samp);
    log_info("n_samples: %d .", n_samp);
    // fprintf(file_logerr," NOT ENOUGH SAMPLES.");
    log_error("NOT ENOUGH SAMPLES.");

    return (0);
  }
  for (x = 0; x < n_samp; x++)
    nnsam[x] = x;
  /* nnnsam = n_samp; */

  /*init unics*/
  for (y = 0; y < n_samp; y++)
  {
    unic[y] = 0;
  }
  /* find positions that are biallelic, excluding the rest. */
  /* IMPORTANT: multiple hits are eliminated from analysis */
  _sites = (double)0;
  transitions = 0;
  transversions = 0;
  *n_realsite = (long int)0;
  mis2 = 0;

  for (xx = 0; xx < n_site; xx++)
  {
    /*eliminate those positions that are not included in analysis*/
    if (matrix_sizepos[xx] == (double)0)
    {
      /*algsites -= (double)1-matrix_sizepos[xx];*/
      continue;
    }
    /*matrix making biallelic (and missing) variants*/
    *n_realsite += (long int)1;

    k_0 = m_1 = k_1 = y = z = v = mis = miso = 0;
    m_0 = *(DNA_matr + (((long long)n_site * (unsigned long)0) + (unsigned long)xx));
    do
    {
      w = *(DNA_matr + (((long long)n_site * (unsigned long)y) + (unsigned long)xx));
      if (w > 48 + 4)
      { /* the position y (also y == 0) is not a, g, t or c */
        mis += 1;
        if (y >= n_samp - args->vint_perpop_nsam[args->npops - 1])
          miso += 1;
        if (args->include_unknown == 1)
        { /*==*/
          if (m_0 > (48 + 4))
          { /*==*/
            if (y + 1 < nsamtot)
            {
              m_0 = *(DNA_matr + (((long long)n_site * (unsigned long)(y + 1)) + (unsigned long)xx));
            }
          }
        }
        else
        {
          z = 1;
          v = 0;
          algsites -= matrix_sizepos[xx];
          *n_realsite -= (long int)1; /*erase*/
        }
      }
      else
      {                        /* the position y is a, g, t or c */
        tcga[y][w - 49] += 1.; /*count for tcga and transition and transversions*/
        if (w != m_0 && w != k_0)
        { /* if position 0 and y are different, and also different of the k_0 */
          if (k_0 == 0)
          { /* if k_0 is 0 (initialized), the k_0 will be w */
            k_0 = w;
            k_1++; /* counting the numbers of w different from m_0 and not > 4 and not mhits */
          }
          else
          { /* if k_0 has already a value */
            if ((w <= 48 + 4) && (m_0 <= 48 + 4))
              v = 1; /* if w and m_0 are different and are also different from k_0 = mhit */
          }
        }
        else
        { /* in case w=m_0 or w=k_0 */
          if (w == m_0)
            m_1++; /* if w equal to m_0 count m_1 */
          else
            k_1++; /* if w is different from m_0 */
        }
      }
      y++;
    } while (z == 0 && y < nsamtot /*n_samp*/);

    /* mhit position */
    if (z == 0 /*valid position*/ && v == 1 /* 0 is no mhit*/ && miso < args->vint_perpop_nsam[args->npops - 1] /*outgroup present*/)
    {
      mhitbp[*nmhits] = xx + 1;
      nmhits[0] = nmhits[0] + 1;
      /*keeping mhits by taking all variants as new biallelic independent sites */
      /*(add missing where were other variants) (in case include_unknown == 1)*/
      if (args->include_unknown == 1)
      {
        aloutg = -1;
        fal[0] = fal[1] = fal[2] = fal[3] = fal[4] = 0;
        /*look at the outgroup: keep the ancestral, if more than one allele reject position, if only N also reject*/
        if (args->outgroup_presence)
        {
          for (y = nsamtot - args->vint_perpop_nsam[args->npops - 1]; y < nsamtot; y++)
          {
            a = (int)*(DNA_matr + (((long long)n_site * (unsigned long)y) + (unsigned long)xx));
            switch (a)
            {
            case 48 + 1:
              fal[0] += 1;
              aloutg = 0;
              break;
            case 48 + 2:
              fal[1] += 1;
              aloutg = 1;
              break;
            case 48 + 3:
              fal[2] += 1;
              aloutg = 2;
              break;
            case 48 + 4:
              fal[3] += 1;
              aloutg = 3;
              break;
            case 48 + 5: /*Ns*/
              fal[4] += 1;
              break;
            case 48 + 6: /*gaps=Ns*/
              fal[4] += 1;
              break;
            }
          }
        }
        if (args->outgroup_presence == 1 && (aloutg == -1 /*all Ns*/ || fal[aloutg] + fal[4] < args->vint_perpop_nsam[args->npops - 1] /*outg polym*/))
        {
          /*reject*/
          algsites -= matrix_sizepos[xx];
          *n_realsite -= (long int)1;
          /*add Ns at the ouygroup*/
          for (y = nsamtot - args->vint_perpop_nsam[args->npops - 1]; y < nsamtot; y++)
            *(DNA_matr + (((long long)n_site * (unsigned long)y) + (unsigned long)xx)) = 48 + 5;
          z = 1;
        }
        else
        {
          /*look at the rest of samples: transform all variants in new biallelic positions*/
          /*keep only the two more frequent alleles, the rest convert to Ns*/
          fal[0] = fal[1] = fal[2] = fal[3] = fal[4] = 0;
          for (y = 0; y < nsamtot - args->vint_perpop_nsam[args->npops - 1]; y++)
          {
            a = (int)*(DNA_matr + (((long long)n_site * (unsigned long)y) + (unsigned long)xx));
            switch (a)
            {
            case 48 + 1:
              fal[0] += 1;
              break;
            case 48 + 2:
              fal[1] += 1;
              break;
            case 48 + 3:
              fal[2] += 1;
              break;
            case 48 + 4:
              fal[3] += 1;
              break;
            case 48 + 5: /*Ns*/
              fal[4] += 1;
              break;
            case 48 + 6: /*gaps=Ns*/
              fal[4] += 1;
              break;
            }
          }
          /*sort the two best (or coincident with the outgroup)*/
          al1 = al2 = al3 = al4 = -1;
          fal1 = fal2 = fal3 = fal4 = 0;
          if (aloutg == -1)
          {
            al1 = 0;
            fal1 = fal[0];
            for (y = 1; y < 4; y++)
            {
              if (fal[y] > fal1)
              {
                al1 = y;
                fal1 = fal[y];
              }
            }
          }
          else
          {
            al1 = aloutg;
            fal1 = fal[aloutg];
          }
          for (y = 0; y < 4; y++)
          {
            if (y != al1)
            {
              if (fal[y] > fal2)
              {
                al2 = y;
                fal2 = fal[y];
              }
            }
          }
          for (y = 0; y < 4; y++)
          {
            if (y != al1 && y != al2)
            {
              if (fal[y] > fal3)
              {
                al3 = y;
                fal3 = fal[y];
              }
            }
          }
          for (y = 0; y < 4; y++)
          {
            if (y != al1 && y != al2 && y != al3)
            {
              if (fal[y] > fal4)
              {
                al4 = y;
                fal4 = fal[y];
              }
            }
          }
          _sites += matrix_sizepos[xx];
          /*algsites -= (double)1 - matrix_sizepos[xx];*/
          /*CALCULATE the length of each sample considering the outgroup*/
          d = 0;
          for (y2 = n_samp - args->vint_perpop_nsam[args->npops - 1]; y2 < n_samp; y2++)
          {
            if (*(DNA_matr + (((long long)n_site * (unsigned long)y2) + (unsigned long)xx)) != 48 + 5 &&
                *(DNA_matr + (((long long)n_site * (unsigned long)y2) + (unsigned long)xx)) != 48 + 6) /*6 added*/
            {
              d = 1;
            }
          }
          for (y = 0; y < nsamtot; y++)
          {
            if (*(DNA_matr + (((long long)n_site * (unsigned long)y) + (unsigned long)xx)) != 48 + 5 &&
                *(DNA_matr + (((long long)n_site * (unsigned long)y) + (unsigned long)xx)) != 48 + 6) /*6 added*/
            {
              /*the length of each sample considering the outgroup*/
              if (d == 1)
                sum_sam[y] += (double)matrix_sizepos[xx];
            }
          }
          /*Loop adding the number of bialsites (2 or 3) at this position*/
          for (y = 0; y < (int)(fal2 > 0) /*+(int)(fal3>0)+(int)(fal4>0)*/; y++)
          {
            if (matrix_segrpos[xx] /*to eliminate biallelic syn/nsyn not desired*/)
            {
              k_0 = *(DNA_matr + (((long long)n_site * (unsigned long)(nsamtot - 1) + (unsigned long)xx))) /*al1+48+1*/;
              if (y == 0)
                m_0 = al2 + 48 + 1;
              if (y == 1)
                m_0 = al3 + 48 + 1;
              if (y == 2)
                m_0 = al4 + 48 + 1;
              if ((k_0 == '1' && m_0 == '2') || (k_0 == '2' && m_0 == '1') ||
                  (k_0 == '4' && m_0 == '3') || (k_0 == '3' && m_0 == '4'))
              {
                transitions += 1;
                matrix_sv[0][bial_sites] = 1;
              }
              else
              {
                transversions += 1;
                matrix_sv[0][bial_sites] = 2;
              }
              d = 0;
              for (y2 = 0; y2 < nsamtot /*n_samp*/; y2++)
              { /* non-outgroup: 0 indicates higher frequency */
                if (*(DNA_matr + (((long long)n_site * (unsigned long)y2) + (unsigned long)xx)) == k_0)
                {
                  matrix_pol[0][((bial_sites * (nsamtot /*n_samp*/)) + y2)] = '0';
                  matrix_pol_tcga[0][((bial_sites * (nsamtot /*n_samp*/)) + y2)] = k_0;
                }
                else
                {
                  if (*(DNA_matr + (((long long)n_site * (unsigned long)y2) + (unsigned long)xx)) == m_0)
                  {
                    matrix_pol[0][((bial_sites * (nsamtot /*n_samp*/)) + y2)] = '1';
                    matrix_pol_tcga[0][((bial_sites * (nsamtot /*n_samp*/)) + y2)] = m_0;
                    d += 1;
                  }
                  else
                  {
                    matrix_pol[0][((bial_sites * (nsamtot /*n_samp*/)) + y2)] = '-';
                    matrix_pol_tcga[0][((bial_sites * (nsamtot /*n_samp*/)) + y2)] = '-';
                  }
                }
              }

              matrix_pos[0][bial_sites] = -(long int)(xx + (long int)1);
              /*positions with missing have a negative position (to exclude in haplotype calculations)*/
              mvbp[mv] = xx + 1;
              mv++;
              if (mv >= maxnsite)
              {
                if ((mvbp = (long int *)realloc(mvbp, (mv + 128) * sizeof(long int))) == 0)
                {
                  // fprintf(file_logerr,"Error: memory not reallocated. get_obsstat.12");
                  log_fatal("Error: memory not reallocated, mvbp get_obsstat.12");
                  return (0);
                }
                maxnsite = mv;
              }

              matrix_freq[0][bial_sites] = d;
              /*calculate unique variants*/
              for (y2 = 0; y2 < nsamtot /*n_samp*/; y2++)
              {
                if ((d == 1 && matrix_pol[0][((bial_sites * (nsamtot /*n_samp*/)) + y2)] == '1') ||
                    (d == nsamtot /*n_samp*/ - 1 && matrix_pol[0][((bial_sites * (nsamtot /*n_samp*/)) + y2)] == '0'))
                  unic[y2] += 1;
              }
              /* one more position */
              bial_sites++;
              /* reallocations */
              if (bial_sites == maxbialsites)
              {
                maxbialsites += 128;
                if ((*matrix_pol = realloc(*matrix_pol, (maxnsamp) * (maxbialsites) * sizeof(char))) == 0)
                {
                  // fprintf(file_logerr,"Error: memory not reallocated. get_obsstat.11");
                  log_fatal("Error: memory not reallocated, matrix_pol get_obsstat.11");
                  return (0);
                }
                if ((*matrix_pol_tcga = realloc(*matrix_pol_tcga, (maxnsamp) * (maxbialsites) * sizeof(char))) == 0)
                {
                  // fprintf(file_logerr,"Error: memory not reallocated. get_obsstat.11");
                  log_fatal("Error: memory not reallocated, matrix_pol_tcga get_obsstat.11");
                  return (0);
                }
                if ((*matrix_pos = realloc(*matrix_pos, (maxbialsites) * sizeof(long int))) == 0)
                {
                  // fprintf(file_logerr,"Error: memory not reallocated. get_obsstat.12");
                  log_fatal("Error: memory not reallocated, matrix_pos get_obsstat.12");
                  return (0);
                }
                if ((*matrix_freq = realloc(*matrix_freq, (maxbialsites) * sizeof(long int))) == 0)
                {
                  // fprintf(file_logerr,"Error: memory not reallocated. get_obsstat.13");
                  log_fatal("Error: memory not reallocated, matrix_freq get_obsstat.13");
                  return (0);
                }
                if ((*matrix_sv = realloc(*matrix_sv, (maxbialsites) * sizeof(long int))) == 0)
                {
                  // fprintf(file_logerr,"Error: memory not reallocated. get_obsstat.14");
                  log_fatal("Error: memory not reallocated, matrix_sv get_obsstat.14");
                  return (0);
                }
              }
            }
          }
        }
      }
      else
      { /*mhits excluded when missing data is not allowed*/
        *n_realsite -= (long int)1;
        algsites -= matrix_sizepos[xx];
        z = 1;
      }
    }
    if (z == 0 /*valid position*/ && v == 0 /*no mhit*/)
    {
      _sites += matrix_sizepos[xx];
      /*algsites -= (double)1 - matrix_sizepos[xx];*/
      /*CALCULATE the length of each sample considering the outgroup*/
      d = 0;
      for (y2 = n_samp - args->vint_perpop_nsam[args->npops - 1]; y2 < n_samp; y2++)
      {
        if (*(DNA_matr + (((long long)n_site * (unsigned long)y2) + (unsigned long)xx)) != 48 + 5 &&
            *(DNA_matr + (((long long)n_site * (unsigned long)y2) + (unsigned long)xx)) != 48 + 6) /*6 added*/
        {
          d = 1;
        }
      }
      for (y = 0; y < nsamtot; y++)
      {
        if (*(DNA_matr + (((long long)n_site * (unsigned long)y) + (unsigned long)xx)) != 48 + 5 &&
            *(DNA_matr + (((long long)n_site * (unsigned long)y) + (unsigned long)xx)) != 48 + 6) /*6 added*/
        {
          /*the length of each sample considering the outgroup*/
          if (d == 1)
            sum_sam[y] += (double)matrix_sizepos[xx];
        }
      }
    }
    /* do the matrix of biallelic positions: */
    if (z == 0 && v == 0 && miso < args->vint_perpop_nsam[args->npops - 1] && m_1 + mis < nsamtot /*n_samp*/ && matrix_segrpos[xx] /*to eliminate biallelic syn/nsyn not desired*/)
    {

      /*count for transitions and transversions*/
      /*
      for(y=0;y<nsamtot;y++) {
        w = *(DNA_matr+(((long long)n_site*(unsigned long)y)+(unsigned long)xx));
        tcga[y][w-49] += 1.;
      }
      */
      if (k_0 != 0 && k_0 != m_0 && v == 0 && z == 0)
      {
        if ((k_0 == '1' && m_0 == '2') || (k_0 == '2' && m_0 == '1') ||
            (k_0 == '4' && m_0 == '3') || (k_0 == '3' && m_0 == '4'))
        {
          transitions += 1;
          matrix_sv[0][bial_sites] = 1;
        }
        else
        {
          transversions += 1;
          matrix_sv[0][bial_sites] = 2;
        }
      }

      /* k_0 will be the large frequency: m_1 and k_1 are lost variables*/
      if (m_1 > k_1)
      {
        m_1 = k_0;
        k_0 = m_0;
        m_0 = m_1;
      }
      d = 0;
      for (y = 0; y < nsamtot /*n_samp*/; y++)
      { /* non-outgroup: 0 indicates higher frequency */
        if (*(DNA_matr + (((long long)n_site * (unsigned long)y) + (unsigned long)xx)) == k_0)
        {
          matrix_pol[0][((bial_sites * (nsamtot /*n_samp*/)) + y)] = '0';
          matrix_pol_tcga[0][((bial_sites * (nsamtot /*n_samp*/)) + y)] = k_0;
        }
        else
        {
          if (*(DNA_matr + (((long long)n_site * (unsigned long)y) + (unsigned long)xx)) == 48 + 5 ||
              *(DNA_matr + (((long long)n_site * (unsigned long)y) + (unsigned long)xx)) == 48 + 6)
          { /*6 addded*/
            matrix_pol[0][((bial_sites * (nsamtot /*n_samp*/)) + y)] = '-';
            matrix_pol_tcga[0][((bial_sites * (nsamtot /*n_samp*/)) + y)] = '-';
          }
          else
          {
            matrix_pol[0][((bial_sites * (nsamtot /*n_samp*/)) + y)] = '1';
            matrix_pol_tcga[0][((bial_sites * (nsamtot /*n_samp*/)) + y)] = m_0;
            d += 1;
          }
        }
      }
      if (mis == 0)
        matrix_pos[0][bial_sites] = (long int)(xx + (long int)1);
      else
      {
        matrix_pos[0][bial_sites] = -(long int)(xx + (long int)1); /*positions with missing have a negative position (to exclude in haplotype calculations)*/
        mvbp[mv] = xx + 1;
        mv++;
        if (mv >= maxnsite)
        {
          if ((mvbp = (long int *)realloc(mvbp, (mv + 128) * sizeof(long int))) == 0)
          {
            // fprintf(file_logerr,"Error: memory not reallocated. get_obsstat.12");
            log_fatal("Error: memory not reallocated, mvbp get_obsstat.12");
            return (0);
          }
          maxnsite = mv;
        }
      }
      matrix_freq[0][bial_sites] = d;
      /*calculate unics*/
      for (y = 0; y < nsamtot /*n_samp*/; y++)
      {
        if ((d == 1 && matrix_pol[0][((bial_sites * (nsamtot /*n_samp*/)) + y)] == '1') ||
            (d == nsamtot /*n_samp*/ - 1 && matrix_pol[0][((bial_sites * (nsamtot /*n_samp*/)) + y)] == '0'))
          unic[y] += 1;
      }
      /* one more position */
      bial_sites++;
      /* reallocations */
      if (bial_sites == maxbialsites)
      {
        /*
        if(maxbialsites == 32767) {
          fprintf(file_logerr,"\n Sorry, it is only accepted a maximum of 32767 biallelic sites per loci. It has been cut at position %ld.",xx+1);
          fprintf(file_logerr,"\n Sorry, it is only accepted a maximum of 32767 biallelic sites per loci. It has been cut at position %ld.",xx+1);
          return(0);
        }
        else if(maxbialsites > 32767 - 128) maxbialsites = 32767;
          else */
        maxbialsites += 128;
        if ((*matrix_pol = realloc(*matrix_pol, (maxnsamp) * (maxbialsites) * sizeof(char))) == 0)
        {
          // fprintf(file_logerr,"Error: memory not reallocated. get_obsstat.11");
          log_fatal("Error: memory not reallocated, matrix_pol get_obsstat.11");
          return (0);
        }
        if ((*matrix_pol_tcga = realloc(*matrix_pol_tcga, (maxnsamp) * (maxbialsites) * sizeof(char))) == 0)
        {
          // fprintf(file_logerr,"Error: memory not reallocated. get_obsstat.11");
          log_fatal("Error: memory not reallocated, matrix_pol_tcga get_obsstat.11");
          return (0);
        }
        if ((*matrix_pos = realloc(*matrix_pos, (maxbialsites) * sizeof(long int))) == 0)
        {
          // fprintf(file_logerr,"Error: memory not reallocated. get_obsstat.12");
          log_fatal("Error: memory not reallocated, matrix_pos get_obsstat.12");
          return (0);
        }
        if ((*matrix_freq = realloc(*matrix_freq, (maxbialsites) * sizeof(long int))) == 0)
        {
          // fprintf(file_logerr,"Error: memory not reallocated. get_obsstat.13");
          log_fatal("Error: memory not reallocated, matrix_freq get_obsstat.13");
          return (0);
        }
        if ((*matrix_sv = realloc(*matrix_sv, (maxbialsites) * sizeof(long int))) == 0)
        {
          // fprintf(file_logerr,"Error: memory not reallocated. get_obsstat.14");
          log_fatal("Error: memory not reallocated, matrix_sv get_obsstat.14");
          return (0);
        }
      }
    }
    if (args->include_unknown == 1)
    { /*redundant conditional? (unknown==1 vs mis>0 or miso>0)*/
      if (z == 0 && (mis == nsamtot || miso == args->vint_perpop_nsam[args->npops - 1]))
      {
        *n_realsite -= (long int)1;
        _sites -= matrix_sizepos[xx];
        algsites -= matrix_sizepos[xx];
      }
      else
      {
        /*totalmis += mis;*/
        if (miso < args->vint_perpop_nsam[args->npops - 1])
          mis2 += (mis - miso);
        /*printf("%ld ",xx);*/
      }
    }
  }
  /*printf("\n",xx);*/
  /*
  if(_sites != algsites) {
    printf("BUG in get_obsatastats.c: _sites != algsites");
  }
  */
  /**nmhits = mhits;*/
  if (file_output)
  {
    if (args->output == 0 || args->output == 10)
    {
      /*
      fprintf(file_logerr,"\n\nWARNING: Positions with missing values and multiple hits (including codon positions with more than 2 variants) are filtered using THE COMPLETE FILE included.");
      fprintf(file_logerr,"\n\nREADING INPUT FILE DATA:\n\n Number of Total Samples (if diploid, doubled): %d\n Valid sites: %.2f\n Multiple hits: %ld\n Polymorphic (biallelic) positions with missing samples: %ld\n Ratio missing/positions: %f\n Ratio trans/transv: %.3f.",nnnsam-!outgroup_presence,_sites,mhits,mv,(double)totalmis/(double)(nnnsam*_sites),(double)transitions/(double)transversions);
      */
      if (*nmhits && args->output == 10)
      {
        fprintf(file_output, "\n Position(s) of the multiple hits ");
        if (args->include_unknown)
          fprintf(file_output, " (minor frequency not outgroup excluded in the analysis):");
        else
          fprintf(file_output, " (excluded from the analysis):");
        for (x = 0; x < *nmhits; x++)
        {
          if (formatfile > 0)
            fprintf(file_output, " %ld", mhitbp[x] + formatfile - 1);
          else
            fprintf(file_output, " %ld", mhitbp[x]);
        }
      }
      if (mv && args->output == 10)
      {
        fprintf(file_output, "\n Position(s) of the Variable (converted biallelic if mhits) base with missing samples:");
        for (x = 0; x < mv; x++)
        {
          if (formatfile > 0)
            fprintf(file_output, " %ld", mvbp[x] + formatfile - 1);
          else
            fprintf(file_output, " %ld", mvbp[x]);
        }
      }
      fprintf(file_output, "\n\n Names from sample selected by the user (the first %d, if diploid is doubled):", (nsamtot - !args->outgroup_presence) / (int)atoi(args->ploidy));
      z = 0;
      for (y = 0; y < args->npops - !args->outgroup_presence; y++)
      {
        if (y < args->npops - 1 || args->npops == 1)
          fprintf(file_output, "\n\n Population[%d]:\n", y);
        else
          fprintf(file_output, "\n\n Outgroup: \n");
        for (x = 0; x < args->vint_perpop_nsam[y] / atoi(args->ploidy); x++)
        {
          fprintf(file_output, "\n %s", names[z]);
          z += atoi(args->ploidy);
        }
      }
      /*
      printf(".\n\nCalculating statistics...\n");
      fflush(stdout);
      */
      fprintf(file_output, "\n\nCalculating statistics...\n");
      fflush(file_output);
    }
  }
  else
  {
    if (args->output == 0 || args->output == 10)
    {
      /*printf(" n_samples: %d, valid sites: %.2f. mhits: %d, transitio/transversion: %.3f.",nnnsam,_sites,*nmhits,(double)transition/(double)transversion);*/
      /*
      printf("\n\nWARNING: Positions with missing values and multiple hits (including codon positions with more than 2 variants) are filtered using THE COMPLETE FILE included.");
      printf("\n\nREADING INPUT FILE DATA:\n\n Number of Total Samples (if diploid, doubled): %d\n Valid sites: %.2f\n Multiple hits: %ld\n Polymorphic (biallelic) positions with missing samples: %ld\n Ratio missing/positions: %f\n Ratio trans/transv: %.3f.",nnnsam-!outgroup_presence,_sites,*nmhits,mv,(double)totalmis/(double)(nnnsam*_sites),(double)transitions/(double)transversions);
      */
      if (*nmhits && args->output == 10)
      {
        fprintf(file_output, "\n Position(s) of the multiple hits ");
        if (args->include_unknown)
          fprintf(file_output, " (minor frequency not outgroup excluded in the analysis):");
        else
          fprintf(file_output, " (excluded from the analysis):");
        for (x = 0; x < *nmhits; x++)
        {
          if (formatfile > 0)
            fprintf(file_output, " %ld", mhitbp[x] + formatfile - 1);
          else
            fprintf(file_output, " %ld", mhitbp[x]);
        }
      }
      if (mv)
      {
        fprintf(file_output, "\n Position(s) of the Variable (biallelic) base with missing samples:");
        for (x = 0; x < mv; x++)
        {
          fprintf(file_output, " %ld", mvbp[x]);
        }
      }
      fprintf(file_output, "\n\n Names from samples selected by the user (the first %d, if diploid is doubled):", (nsamtot - !args->outgroup_presence) / atoi(args->ploidy));
      z = 0;
      for (y = 0; y < args->npops - !args->outgroup_presence; y++)
      {
        if (y < args->npops - 1 || args->npops == 1)
          fprintf(file_output, "\n\n Population[%d]:\n", y);
        else
          fprintf(file_output, "\n\n Outgroup: \n");
        for (x = 0; x < args->vint_perpop_nsam[y] / atoi(args->ploidy); x++)
        {
          fprintf(file_output, "\n %s", names[z]);
          z++;
        }
      }
      printf("\n");
    }
  }
  if (_sites == (double)0)
  {
    for (y = 0; y < args->npops; y++)
    {
      nsites1_pop[y] = 0.;
      nsites1_pop_outg[y] = 0.;
      nsites2_pop[y] = 0.;
      nsites2_pop_outg[y] = 0.;
      nsites3_pop[y] = 0.;
      nsites3_pop_outg[y] = 0.;
      anx[y] = bnx[y] = 0.0;
      anxo[y] = bnxo[y] = 0.0;
    }
    /*
    if(file_output)
      fprintf(file_logerr,"Not valid sites available ");
    else
      printf("Not valid sites available ");
    return 0;*/
  }
  else
  { /*print mask, in case considering missing values, in a new file named ..._mask.txt*/
    if (/*include_unknown && */ file_mask)
    {
      xm = 0;
      for (xx = 0; xx < n_site; xx++)
      {
        /*
        w = (char)0;
        for(y=0;y<nsamtot;y++) {
          w = *(DNA_matr+(((long long)n_site*(unsigned long)y)+(unsigned long)xx));
          if(w > (48+5)) {
            matrix_sizepos[xx] = (double)0;
            fprintf(file_mask,"%.3f ",matrix_sizepos[xx]);
            break;
          }
        }
        if(w <= (48+5)) */
        while (mhitbp[xm] < xx + 1 && xm < *nmhits)
        {
          xm++;
        }
        if (args->include_unknown == 0 && mhitbp[xm] == xx + 1)
          fprintf(file_mask, "%.3f ", (float)0);
        else
          fprintf(file_mask, "%.3f ", matrix_sizepos[xx]); /*print the value of each position*/
      }
      fprintf(file_mask, "\n");
      for (x = 0; x < nsamtot; x++)
      {
        for (xx = 0; xx < n_site; xx++)
        {
          /*if(matrix_sizepos[xx] > (double)0) {
              for(y=0;y<nsamtot;y++) {
                w = *(DNA_matr+(((long long)n_site*(unsigned long)y)+(unsigned long)xx));
                if(w > (48+5))
                  break;
              }
              if(y == nsamtot) {
          */
          w = *(DNA_matr + (((long long)n_site * (unsigned long)x) + (unsigned long)xx));
          if (w >= 48 + 5)
            fprintf(file_mask, "0"); /*missing*/
          else
            fprintf(file_mask, "1"); /*normal*/
                                     /*		}
                                     }*/
        }
        fprintf(file_mask, "\n");
      }
    }
    if (args->include_unknown)
    {
      /*calculate nsites1_pop, nsites2_pop and nsites1_pop_outg, nsites2_pop_outg*/
      for (y = 0; y < args->npops; y++)
      {
        nsites1_pop[y] = 0.;
        nsites1_pop_outg[y] = 0.;
        nsites2_pop[y] = 0.;
        nsites2_pop_outg[y] = 0.;
        nsites3_pop[y] = 0.;
        nsites3_pop_outg[y] = 0.;
      }
      z = 0;
      for (y = 0; y < args->npops; y++)
      {
        anx[y] = bnx[y] = 0.0;
        anxo[y] = bnxo[y] = 0.0;
        xm = 0;
        for (xx = 0; xx < n_site; xx++)
        {
          v = 0;
          for (x = z; x < z + args->vint_perpop_nsam[y]; x++)
          {
            w = *(DNA_matr + (((long long)n_site * (unsigned long)x) + (unsigned long)xx));
            if (w >= 48 + 5)
              v += 1;
          }
          while (xm < *nmhits && mhitbp[xm] < xx + 1)
          {
            xm++;
          }
          if ((matrix_sizepos[xx] > 0 /* && mhitbp[xm] != xx+1*/) /* && ((include_unknown == 1) || (include_unknown==0 && v == 0))*/)
          {
            if (v < args->vint_perpop_nsam[y])
              nsites1_pop[y] += (matrix_sizepos[xx]);
            if (v < args->vint_perpop_nsam[y] - 1)
              nsites2_pop[y] += (matrix_sizepos[xx]);
            if (v < args->vint_perpop_nsam[y] - 2)
              nsites3_pop[y] += (matrix_sizepos[xx]);
            for (k = 1; k < args->vint_perpop_nsam[y] - v; k++)
            {
              anx[y] += 1.0 / ((double)k);
              bnx[y] += 1.0 / ((double)k * (double)k);
            }
            p = 0;
            if (args->outgroup_presence == 1)
            {
              for (x = nsamtot - 1; x >= nsamtot - args->vint_perpop_nsam[args->npops - 1]; x--)
              {
                q = *(DNA_matr + (((long long)n_site * (unsigned long)x) + (unsigned long)xx));
                if (q >= 48 + 5)
                  p += 1;
              }
              if (p < args->vint_perpop_nsam[args->npops - 1] && v < args->vint_perpop_nsam[y] /* && mhitbp[xm] != xx+1*/) /*{*/
                nsites1_pop_outg[y] += (matrix_sizepos[xx]);                             /*printf("%ld ",xx);}*/
              if (p < args->vint_perpop_nsam[args->npops - 1] && v < args->vint_perpop_nsam[y] - 1 /* && mhitbp[xm] != xx+1*/)
                nsites2_pop_outg[y] += (matrix_sizepos[xx]);
              if (p < args->vint_perpop_nsam[args->npops - 1] && v < args->vint_perpop_nsam[y] - 2 /* && mhitbp[xm] != xx+1*/)
                nsites3_pop_outg[y] += (matrix_sizepos[xx]);
              if (p < args->vint_perpop_nsam[args->npops - 1])
              {
                for (k = 1; k < args->vint_perpop_nsam[y] - v; k++)
                {
                  anxo[y] += 1.0 / ((double)k);
                  bnxo[y] += 1.0 / ((double)k * (double)k);
                }
              }
            }
          }
        }
        z += args->vint_perpop_nsam[y];
        if (nsites2_pop[y] > 0)
        {
          anx[y] /= (double)nsites2_pop[y];
          bnx[y] /= (double)nsites2_pop[y];
        }
        if (nsites2_pop_outg[y] > 0)
        {
          anxo[y] /= (double)nsites2_pop_outg[y];
          bnxo[y] /= (double)nsites2_pop_outg[y];
        }
      }
      /*calculate lengthamng*/
      for (y = 0; y < args->npops - 1; y++)
      {
        for (y2 = y + 1; y2 < args->npops; y2++)
        {
          lengthamng[y][y2] = 0.;
          lengthamng_outg[y][y2] = 0.;
        }
      }
      xm = 0;
      for (xx = 0; xx < n_site; xx++)
      {
        z = 0;
        for (y = 0; y < args->npops - 1; y++)
        {
          v = 0;
          for (x = z; x < z + args->vint_perpop_nsam[y]; x++)
          {
            w = *(DNA_matr + (((long long)n_site * (unsigned long)x) + (unsigned long)xx));
            if (w >= 48 + 5)
              v += 1;
          }
          z2 = z;
          x = y;
          while (x < y + 1)
          {
            z2 += args->vint_perpop_nsam[x];
            x++;
          }
          for (y2 = y + 1; y2 < args->npops; y2++)
          {
            v2 = 0;
            for (x = z2; x < z2 + args->vint_perpop_nsam[y2]; x++)
            {
              w = *(DNA_matr + (((long long)n_site * (unsigned long)x) + (unsigned long)xx));
              if (w >= 48 + 5)
                v2 += 1;
            }
            while (mhitbp[xm] < xx + 1 && xm < *nmhits)
            {
              xm++;
            }
            p = 0;
            if ((/*include_unknown == 1 && */ matrix_sizepos[xx] > 0) /* || (include_unknown==0 && v == 0 && matrix_sizepos[xx] > 0)*/)
            {
              if (args->outgroup_presence == 1)
              {
                for (x = nsamtot - 1; x >= nsamtot - args->vint_perpop_nsam[args->npops - 1]; x--)
                {
                  q = *(DNA_matr + (((long long)n_site * (unsigned long)x) + (unsigned long)xx));
                  if (q >= 48 + 5)
                    p += 1;
                }
                if (p < args->vint_perpop_nsam[args->npops - 1] && v < args->vint_perpop_nsam[y] && v2 < args->vint_perpop_nsam[y2] /* && mhitbp[xm] != xx+1*/)
                  lengthamng_outg[y][y2] += (matrix_sizepos[xx]);
              }
              if (v < args->vint_perpop_nsam[y] && v2 < args->vint_perpop_nsam[y2] /* && mhitbp[xm] != xx+1*/)
                lengthamng[y][y2] += (matrix_sizepos[xx]);
            }
            z2 += args->vint_perpop_nsam[y2];
          }
          z += args->vint_perpop_nsam[y];
        }
      }
    }
    else
    {
      for (y = 0; y < args->npops; y++)
      {
        for (k = 1; k < args->vint_perpop_nsam[y]; k++)
        {
          anx[y] += 1.0 / ((double)k);
          bnx[y] += 1.0 / ((double)k * (double)k);
        }
        if (args->vint_perpop_nsam[y] > 0)
          nsites1_pop[y] = algsites;
        if (args->vint_perpop_nsam[y] > 1)
          nsites2_pop[y] = algsites;
        if (args->vint_perpop_nsam[y] > 2)
          nsites3_pop[y] = algsites;
        for (y2 = y + 1; y2 < args->npops; y2++)
        {
          lengthamng[y][y2] = algsites;
          lengthamng_outg[y][y2] = algsites;
        }

        if (args->outgroup_presence == 1)
        {
          for (k = 1; k < args->vint_perpop_nsam[y]; k++)
          {
            anxo[y] += 1.0 / ((double)k);
            bnxo[y] += 1.0 / ((double)k * (double)k);
          }
          if (args->vint_perpop_nsam[y] > 0)
            nsites1_pop_outg[y] = algsites;
          if (args->vint_perpop_nsam[y] > 1)
            nsites2_pop_outg[y] = algsites;
          if (args->vint_perpop_nsam[y] > 2)
            nsites3_pop_outg[y] = algsites;
        }
      }
    }

    /**/
  }

  fflush(file_output);
  *length_al = algsites;
  *length_seg = bial_sites;

  if (transversions)
    *svratio = (double)transitions / (double)transversions;
  else
    *svratio = (double)-10000;
  /**missratio = (double)totalmis/(double)(nnnsam*_sites);*/
  if (((double)(n_samp - args->vint_perpop_nsam[args->npops - 1]) * n_realsite[0]))
    *missratio = (double)mis2 / ((double)(n_samp - args->vint_perpop_nsam[args->npops - 1]) * n_realsite[0]);
  else
    *missratio = 1.0;

  free(unic);
  /*free(mhitbp);*/
  free(nnsam);
  free(mvbp);

  return 1;
}
