/***************************************************************************

  gendynGssEALF.c

  Author: Takashi U. Ito
  e-mail: tuito@post.j-parc.jp

***************************************************************************/

/***************************************************************************
 *   Copyright (C) 2024 by Takashi U. Ito                             *
 *   tuito@post.j-parc.jp                                                  *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_statistics.h>
#define DIM 3
#define gmmm 0.08516103701  //  rad/G
#define ndiv 500
#define tmax 10.0
//#define ndiv 5000
//#define tmax 100.0
#define N_TH 32  // number of threads
#define ntry_th 312500  // ntry per thread; ntry = N_TH*ntry_th = 1E7

///////////// Global variables //////////////
// time interval
double dt;
// sigma in G (static component)
double sigma_sta;
// sigma in G (dynamic component)
double sigma;
//field flip rate (<< 1.0/dt)
double nu;
//longitudinal field in G
double lf;
// initial muon spin polarization
gsl_vector * pol0;
// external field in Gauss
gsl_vector * Bext;
// random seeds
gsl_rng *r_global[N_TH];
gsl_rng *r_global2[N_TH];
// matrix for containing results
gsl_matrix *p[N_TH];
///////////////////////////////////////////


int multivariate_gaussian(double sigma, gsl_vector* result, int id)
{
  gsl_vector * mu = gsl_vector_calloc(DIM);
  gsl_matrix * Sigma = gsl_matrix_calloc(DIM, DIM);
  gsl_matrix * L = gsl_matrix_calloc(DIM, DIM);

  gsl_vector_set_all(mu, 0.0);
  gsl_matrix_set_identity(L);
  gsl_matrix_scale(L,sigma);

  gsl_ran_multivariate_gaussian(r_global[id], mu, L, result);

  gsl_vector_free(mu);
  gsl_matrix_free(Sigma);
  gsl_matrix_free(L);
  
  return 0;
  
}

// ####### Rodrigues' rotation formula ########
int rot(gsl_vector * b, double dt, gsl_matrix* result){
  gsl_vector * n = gsl_vector_calloc(DIM);
  gsl_vector_memcpy(n,b);
  gsl_vector_scale(n,1.0/gsl_blas_dnrm2(b));
  double theta=gmmm*gsl_blas_dnrm2(b)*dt;
  gsl_matrix_set(result,0,0, cos(theta)+pow(gsl_vector_get(n,0),2.0)*(1.0-cos(theta)));
  gsl_matrix_set(result,0,1, gsl_vector_get(n,0)*gsl_vector_get(n,1)*(1.0-cos(theta)) - gsl_vector_get(n,2)*sin(theta));
  gsl_matrix_set(result,0,2, gsl_vector_get(n,2)*gsl_vector_get(n,0)*(1.0-cos(theta)) + gsl_vector_get(n,1)*sin(theta));
  gsl_matrix_set(result,1,0, gsl_vector_get(n,0)*gsl_vector_get(n,1)*(1.0-cos(theta)) + gsl_vector_get(n,2)*sin(theta));
  gsl_matrix_set(result,1,1, cos(theta)+pow(gsl_vector_get(n,1),2.0)*(1.0-cos(theta)));
  gsl_matrix_set(result,1,2, gsl_vector_get(n,1)*gsl_vector_get(n,2)*(1.0-cos(theta)) - gsl_vector_get(n,0)*sin(theta));
  gsl_matrix_set(result,2,0, gsl_vector_get(n,2)*gsl_vector_get(n,0)*(1.0-cos(theta)) - gsl_vector_get(n,1)*sin(theta));
  gsl_matrix_set(result,2,1, gsl_vector_get(n,1)*gsl_vector_get(n,2)*(1.0-cos(theta)) + gsl_vector_get(n,0)*sin(theta));
  gsl_matrix_set(result,2,2, cos(theta)+pow(gsl_vector_get(n,2),2.0)*(1.0-cos(theta)));
  gsl_vector_free(n);

  return 0;
}
// ###################


void *calc_pol(void *arg)
{
  int id = *(int*)arg; //thread id
  int i=0, j=0, k=0;
  gsl_vector * pol = gsl_vector_calloc(DIM);
  gsl_vector * polold = gsl_vector_calloc(DIM);
  gsl_vector * tmpv = gsl_vector_calloc(DIM);
  gsl_vector * tmpv2 = gsl_vector_calloc(DIM);
  gsl_vector * Bint = gsl_vector_calloc(DIM);
  gsl_vector * Btot = gsl_vector_calloc(DIM);
  gsl_vector * Bint_sta = gsl_vector_calloc(DIM);
  gsl_matrix * R = gsl_matrix_calloc(DIM, DIM);
  
  for(k=0;k<ntry_th;k++){
    if((k%10000)==0) printf("%d: #%d/%d\n", id,k,ntry_th);
    multivariate_gaussian(sigma_sta,Bint_sta,id);
      gsl_vector_memcpy(pol,pol0);
      multivariate_gaussian(sigma,Bint,id);
      gsl_vector_memcpy(Btot,Bint);
      gsl_vector_add(Btot,Bint_sta);
      gsl_vector_add(Btot,Bext);
      for(i=0;i<ndiv;i++){
	if( gsl_rng_uniform(r_global2[id]) < (1.0-exp(-nu*dt)) ){ //Poisson's distribution (n!=0)
	  multivariate_gaussian(sigma,Bint,id);
	  gsl_vector_memcpy(Btot,Bint);
	  gsl_vector_add(Btot,Bint_sta);
	  gsl_vector_add(Btot,Bext);
	}
	gsl_vector_memcpy(polold,pol);
	rot(Btot,dt,R);
	gsl_blas_dgemv(CblasNoTrans,1.0,R,polold,0.0,pol);
	gsl_matrix_get_row(tmpv,p[id],i+1);
	gsl_vector_memcpy(tmpv2,pol);
	gsl_vector_scale(tmpv2,1.0/((double)ntry_th));
	gsl_vector_add(tmpv,tmpv2);
	gsl_matrix_set_row(p[id],i+1,tmpv);
      }
  }
}


int main(int argc, char *argv[]){

  if(argc !=5){
    printf("usage: $this Delta(us-1) Q nu/Delta LF/(Delta/gamma)\n");
    exit(1);
  }

  // time interval
  dt = (double)tmax / (double)ndiv;
  // sigma in G (static component)
  sigma_sta = (double)atof(argv[1])*pow(1.0-(double)atof(argv[2]),0.5)/gmmm;
  // sigma in G (dynamic component)
  sigma = (double)atof(argv[1])*pow((double)atof(argv[2]),0.5)/gmmm;
  //field flip rate (< 1.0/dt)
  nu = (double)atof(argv[3])*(double)atof(argv[1]);
  //longitudinal field in G
  lf = (double)atof(argv[4])*(double)atof(argv[1])/gmmm;


  //////////// initialization ////////////////
  // initial muon spin polarization
  pol0 = gsl_vector_calloc(DIM);
  gsl_vector_set_all(pol0, 0.0);
  gsl_vector_set( pol0, 2, 1.0);  //[0,0,1]

  
  // external field in Gauss
  Bext = gsl_vector_calloc(DIM);
  gsl_vector_set(Bext, 2, lf); //[0.0, 0.0, 0.0*B]

  // random seeds
  int i,j;
  for(i=0;i<N_TH;i++){
    r_global[i] = gsl_rng_alloc(gsl_rng_default);
    r_global2[i] = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(r_global[i], time(NULL)+i);
    gsl_rng_set(r_global2[i], time(NULL)+N_TH+1+i);    
  }

  // matrix for containing results
  for(i=0;i<N_TH;i++){  
    p[i] = gsl_matrix_calloc(ndiv+1, DIM);
    gsl_matrix_set_all(p[i], 0.0);
    gsl_matrix_set_row(p[i], 0, pol0);
  }
  //////////// end of initialization ////////////////

  //////////// create threads ////////////////  
  pthread_t thread[N_TH];
  int arg[N_TH];
  for(i=0; i<N_TH; i++){
    arg[i]=i;
    pthread_create(&thread[i], NULL, &calc_pol, (void*)&arg[i]);
  }

  for(i=0; i<N_TH; i++){
    pthread_join(thread[i], NULL);
  }

  
  //////////// average p[i]s ////////////////
  gsl_matrix *p_ave = gsl_matrix_calloc(ndiv+1, DIM);
  gsl_matrix_set_all(p_ave, 0.0);
  for(i=0; i<N_TH; i++){
    gsl_matrix_scale(p[i],1.0/(double)N_TH);
    gsl_matrix_add(p_ave,p[i]);
  }
  
  //
  char fname[128];
  sprintf(fname, "%s%s%s%s%s%s%s%s", argv[1], "-" ,argv[2], "-", argv[3], "-", argv[4], ".dat");
  FILE *fp;
  fp = fopen(fname,"w");
  printf("%s%s%s","Writing results on ", fname, ".\n");
  
  gsl_vector * tmpv = gsl_vector_calloc(DIM);
  for (i = 0; i < ndiv+1; i++) {
    gsl_matrix_get_row(tmpv,p_ave,i);
    fprintf(fp,"%lf %lf %lf %lf\n", i*dt*(double)atof(argv[1]), gsl_vector_get(tmpv,0),gsl_vector_get(tmpv,1),gsl_vector_get(tmpv,2)); //delta*t, px, py, pz

  }
  
  fclose(fp);

  return 0;
}




