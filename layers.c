/*************************************************************************
 layers.c
 ---------------------------------------------------------------------
                        Copyright (c) 2020-2022
                Andreia P. Guerreiro <andreia.guerreiro@tecnico.ulisboa.pt>
             
 This program is free software (software libre); you can redistribute
 it and/or modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2 of the
 License, or (at your option) any later version.
 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.
 You should have received a copy of the GNU General Public License
 along with this program; if not, you can obtain a copy of the GNU
 General Public License at:
                 http://www.gnu.org/copyleft/gpl.html
 or by writing to:
           Free Software Foundation, Inc., 59 Temple Place,
                 Suite 330, Boston, MA 02111-1307 USA
 ----------------------------------------------------------------------
 Reference:

    [1] Guerreiro, Andreia P.; Manquinho, Vasco; Figueira, Jose Rui. "Exact
        hypervolume subset selection through incremental computations". Computers
        & Operations Research 136 (2021): 105471. http://dx.doi.org/10.1016/j.cor.
        2021.105471.

*************************************************************************/


#include <ilcplex/cplex.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "timer.h"
#include "../HVC/hvc.h"
#include "../HVC/hv-plus.h"
#include "../EAF3D/eaf.h"
#include "../gHSS-master/gHSS.h"


#include <stdarg.h>


static void
   free_and_null     (char **ptr);


int _printSummary = 1;
int _startAlg = 0; //default


void setVerboseMode(int ps){
    _printSummary = ps;
}

void setCPLEXStartAlg(int sa){
    _startAlg = sa;
}
   
/* ------------------------------------------------------------------- */
static void sortPoints(double * pts, int n, int d, int * ix);
/* ------------------------------------------------------------------- */

  


/* -------------------------------------------------------------------
 * Auxiliar functions
 --------------------------------------------------------------------- */

static double * realloc_adouble(double * a, int oldsz, int n){
    if(oldsz >= n)
        return a;
    if(oldsz == 0)
        return (double *) malloc(n * sizeof(double));
    return realloc(a, n * sizeof(double));
}

static int * realloc_aint(int * a, int oldsz, int n){
    if(oldsz >= n)
        return a;
    if(oldsz == 0)
        return (int *) malloc(n * sizeof(int));
    return realloc(a, n * sizeof(int));
}

static char * achar(char * a, int oldsz, int n){
    if(oldsz >= n)
        return a;
    if(oldsz == 0)
        return (char *) malloc(n * sizeof(char));
    return realloc(a, n * sizeof(char));
}


static double * new_realloc_adouble(int n){
    double * a = (double *) malloc(n * sizeof(double));  
    return a;
}


static double * new_ones(int n){
    double * ones = new_realloc_adouble(n);
    for(int i = 0; i < n; i++)
        ones[i] = 1;
    return ones;
}

    
static void arange(int * a, int n){
    for(int i = 0; i < n; i++) a[i] = i;
}
    

static void fill_ones(double * ones, int n){
    for(int i = 0; i < n; i++)
        ones[i] = 1;
}

static void refill_int(int * a, int oldsz, int n, int v){
    for(int i = oldsz; i < n; i++)
        a[i] = v;
}

static void refill_char(char * a, int oldsz, int n, char c){
    for(int i = oldsz; i < n; i++)
        a[i] = c;
}


static double countDigits(int n){
    if(n == 0) return 1;
    return floor(log10(n) + 1);
}

    
    
static char ** varNames(int lvl, int n){
    char ** names = (char **) malloc(n * sizeof(double));
    int lvlsz = countDigits(lvl);
    int nsz = countDigits(n);
    for(int i = 0; i < n; i++){
        names[i] = (char *) malloc((3+lvlsz+nsz)*sizeof(char));
        sprintf(names[i], "x%d_%d", lvl, i);
    }
    return names;
}

static char * binaryCtype(char * ctypes, int oldsz, int n){
    if(oldsz >= n) return ctypes;
    int b = 0;
    if(oldsz == 0)
        ctypes = (char *) malloc(n * sizeof(char));
    else{
        ctypes = realloc(ctypes, n * sizeof(char));
        b = oldsz;
    }
    for(int i = b; i < n; i++)
        ctypes[i] = 'B';
    return ctypes;
}


int allones(double * a, int n){
    int sum = 0;
    for(int i = 0; i < n; i++)
        if(round(a[i]) <= 0) return 0; 
    return 1;
    
}


int countones(double * a, int n){
    int sum = 0;
    for(int i = 0; i < n; i++)
        sum += round(a[i]);
    return sum;
    
}
   
/* ------------------------------------------------------------------- */

   
double * hvc2d(double * pts, int n, double * ref, double * contribs){
    
    //sort
    int * ix = (int *) malloc(n * sizeof(int));
    sortPoints(pts, n, 2, ix);
    
    if(n == 1){
        contribs[0] = (ref[1]-pts[1]) * (ref[0] - pts[0]);
    }else{
        contribs[ix[0]] = (ref[1]-pts[2*ix[0]+1]) * (pts[2*ix[1]] - pts[2*ix[0]]);
        for(int i = 1; i < n-1; i++){
            contribs[ix[i]] = (pts[2*ix[i-1]+1] - pts[2*ix[i]+1]) * (pts[2*ix[i+1]] - pts[2*ix[i]]);
        }
        contribs[ix[n-1]] = (pts[2*ix[n-2]+1] - pts[2*ix[n-1]+1]) * (ref[0] - pts[2*ix[n-1]]);
    }
    free(ix);
}
   

static double hvsubset(double * pts, int n, int k, int d, double * ref, double * x){
    double * subset = (double *) malloc(d*k * sizeof(double));
    double hv;
    int i,j, c = 0;
    for(i = 0; i < n; i++){
        if(x[i] > 0.4){
            for(j = 0; j < d; j++){
                subset[d*c+j] = pts[d*i+j];
            }
            c++;
        }
    }
    hv = hvplus(subset, d, k, ref, 1);
    free(subset);
    return hv;
}

   
double * contributions(double * pts1, double * pts2, int n1, int n2, int d, double * ref, double * contribs){
    double * pts = (double *) malloc((n1+n2)*d*sizeof(double));
    int n = n1 + n2;
    
    for(int i = 0; i < n1*d; i++) pts[i] = pts1[i];
    for(int i = n1*d, j = 0; j < n2*d; i++, j++) pts[i] = pts2[j];
    
    contribs = realloc_adouble(contribs, n1, n);
    
    
    if(d == 2){
        hvc2d(pts, n, ref, contribs);
    }else if(d <= 4){
        hvc(pts, d, n, ref, contribs, 1);
    }else{
        fprintf(stderr, "Not implemented yet!\n");
        exit(-1);
    }
    free(pts);
    
    int factor = 1000000;
    for(int i = 0; i < n; i++) contribs[i] *= factor;
    
    return contribs;
}

double hvOfSelected(double * pts, int n, int d, double * ref, int * selected, int k){
    double * A = (double *) malloc(k*d * sizeof(double));
    double hv = 0;
    for(int i = 0, c = 0; i < n; i++){
        if(selected[i]){
            for(int j = 0; j < d; j++)
                A[c*d+j] = pts[i*d+j];
            c++;
        }
    }
    
    hv = hvplus(A, d, k, ref, 1);
    free(A);
    return hv;
}


/* ------------------------------------------------------------------- */


/* -------------------------------------------------------------------
 * Temporary functions -- should be moved to appropriate files or deleted
 --------------------------------------------------------------------- */


static int compare_points(const void * p1, const void * p2){
    double x1 = **(double **)p1;
    double x2 = **(double **)p2;
    if(x1 < x2) return -1;
    if(x2 < x1) return 1;
    return 0;
}


static void sortPoints(double * pts, int n, int d, int * ix){
    double ** ps = (double **) malloc(n * sizeof(double *));
    for(int i = 0; i < n; i++) 
        ps[i] = &pts[d*i];
    
    qsort((void *) ps, n, sizeof(double *), compare_points);
    
    for(int i = 0; i < n; i++){
        ix[i] = (ps[i]-pts)/d;
    }

    
    free(ps);
}


static int leq(double * a, double * b, int d){
    for(int i = 0; i < d; i++){
        if(a[i] > b[i]) return 0;
    }
    return 1;
}



double * ecdf(double * pts, int n, int d, int l, int * m){
    int sz = n+1;
    double * lvlpoints = (double *) malloc(d*sz * sizeof(double));
    
    if(d == 3){
        *m = 0;
        int malloced = sz;
        int * sets = (int *) malloc(n * sizeof(int));
        int reqlvls[1] = {l+1};
        
        for(int i = 0; i < n; i++){
            sets[i] = 1;
        }
        
        lvlpoints = eaf3dNew(pts, d, sets, n, reqlvls, 1, m, lvlpoints);
        
        free(sets);

    }else{
        printf("ECDF not implemented yet for d>3!\n");
        exit(1);
    }
        
    return lvlpoints;
}



static int getdominators(double * pts, int dvaroffset, int n, double * inters, int m, int d, int * dominators, int * dcnt){
    int c = 0;
    for(int i = 0; i < m; i++){
        dcnt[i] = 0;
        for(int j = 0; j < n; j++){
            
            if(leq(&pts[j*d], &inters[i*d], d)){
                dominators[c] = dvaroffset + j;
                dcnt[i]++;
                c++;
            }
        }
    }
    return c;
}



/* ------------------------------------------------------------------- */


/* -------------------------------------------------------------------
 * Layer's computation
 --------------------------------------------------------------------- */

typedef struct layers_t {
    double ** layers;
    double ** layersc;
    int *     lysizes;
    int *     clysizes;
    
    int t;
    int d;
    double * ref;
    
    int ** domrs;
    int ** domrc;
    int * domrssize;
    
    double * ecdftimes;
    double * contribstimes;
    double * domrstimes;
    
    double ecdftime;
    double contribstime;
    double domrstime;

} layers_t;


static layers_t * layers_alloc(int n, int d){
    layers_t * lyt = (layers_t *) malloc(sizeof(layers_t));
    
    lyt->layers   = (double **) malloc((n+1) * sizeof(double *));
    lyt->layersc  = (double **) malloc(n * sizeof(double *));
    lyt->lysizes  = (int *) malloc((n+1) * sizeof(int));
    lyt->clysizes = (int *) malloc((n+1) * sizeof(int));
    
    lyt->ref = (double *) malloc(d * sizeof(double));
    
    
    lyt->domrs = (int **) malloc(n * sizeof(int *));
    lyt->domrc = (int **) malloc(n * sizeof(int *));
    lyt->domrssize = (int *) malloc(n * sizeof(int));
    
    lyt->ecdftimes     = (double *) malloc(n * sizeof(double));
    lyt->contribstimes = (double *) malloc(n * sizeof(double));
    lyt->domrstimes    = (double *) malloc(n * sizeof(double));
    
    return lyt;
}


static layers_t * layers_init(double * pts, int n, int d, double * ref){
    layers_t * lyt = layers_alloc(n, d);
    int i;
    
    lyt->layers[0] = (double *) malloc(n*d * sizeof(double));
    for(i = 0; i < n*d; i++)
        lyt->layers[0][i] = pts[i];
    
    lyt->lysizes[0] = n;
    lyt->clysizes[0] = 0;
    lyt->clysizes[1] = n;
    
    lyt->t = 0;
    lyt->d = d;
    
    for(i = 0; i < d; i++)
        lyt->ref[i] = ref[i];
    
    lyt->domrs[0] = NULL;
    lyt->domrc[0] = NULL;
    lyt->domrssize[0] = 0;
    
    lyt->ecdftimes[0] = 0;
    lyt->domrstimes[0] = 0;
    
    lyt->ecdftime = 0;
    lyt->contribstime = 0;
    lyt->domrstime = 0;
    
    return lyt;
}


static void layers_update(layers_t * lyt, int t){
    
    if(lyt->t > t)
        return;
    
    int l;
    int n = lyt->lysizes[0];
    int d = lyt->d;
    int ni = 0;
    double * ref = lyt->ref;
    
    double tmptime;

    
    
    for(l = lyt->t; l < t; l++){
        
        //--
        tmptime = Timer_elapsed_virtual();
        
        if(l < n-1){
            lyt->layers[l+1] = ecdf(lyt->layers[0], n, d, l+1, &ni);
            lyt->lysizes[l+1] = ni;
            lyt->clysizes[l+2] = lyt->clysizes[l+1] + ni;
                
            lyt->ecdftimes[l+1] = Timer_elapsed_virtual() - tmptime;
            lyt->ecdftime += lyt->ecdftimes[l+1];
            
        }else{
            lyt->lysizes[l+1] = 0;
            lyt->layers[l+1] = NULL;
        }
        
        //--
        tmptime = Timer_elapsed_virtual();
        ni = lyt->lysizes[l];
        lyt->layersc[l] = (double *) malloc(ni * sizeof(double));
        lyt->layersc[l] = contributions(lyt->layers[l], lyt->layers[l+1], lyt->lysizes[l], lyt->lysizes[l+1], d, ref, lyt->layersc[l]);
        
        lyt->contribstimes[l] = Timer_elapsed_virtual() - tmptime;
        lyt->contribstime += lyt->contribstimes[l];
        //--
    
        
        //--
        tmptime = Timer_elapsed_virtual();
        
        if(l > 0){
            lyt->domrs[l] = (int *) malloc((ni+1) * n * sizeof(int)); 
            lyt->domrc[l] = (int *) malloc(ni * sizeof(int));
            if(l >= 0){
                lyt->domrssize[l] = getdominators(lyt->layers[0], lyt->clysizes[0], n, lyt->layers[l], ni, d, lyt->domrs[l], lyt->domrc[l]);
            }else{
                lyt->domrssize[l] = getdominators(lyt->layers[l-1], lyt->clysizes[l-1], lyt->lysizes[l-1], lyt->layers[l], ni, d, lyt->domrs[l], lyt->domrc[l]);
            
            lyt->domrstimes[l] = Timer_elapsed_virtual() - tmptime;
            lyt->domrstime += lyt->domrstimes[l];
            //--
            }
        }
        
    }
    lyt->t = t;
    
}



static double layers_hv_level(layers_t * lyt, int l, double * x0, double * xl){
    int ni = lyt->lysizes[l];
    int i, j;
    int covered, domrc, c, nc;
    double hv = 0;
    
    if(l > lyt->t){
        exit(1);
    }
    
    c = 0;
    nc = 0;
    for(i = 0; i < ni; i++){
        covered = 0;
        domrc = lyt->domrc[l][i];
        for(j = 0; j < domrc && !covered; j++){
            covered = x0[lyt->domrs[l][c + j]];
        }
        c += domrc;
        
        if(covered){
            xl[i] = 1;
            hv += lyt->layersc[l][i];
        }else{
            xl[i] = 0;
        }
        nc += covered;
    }
        
    return hv;
}


static double layers_uncoveredhv_level(layers_t * lyt, int l, double * x0, double * xl){
    int ni = lyt->lysizes[l];
    int i, j;
    int covered, domrc, c, nc;
    double hv = 0;
    
    if(l > lyt->t){
        exit(1);
    }
    
    c = 0;
    nc = 0;
    for(i = 0; i < ni; i++){
        covered = 0;
        domrc = lyt->domrc[l][i];
        for(j = 0; j < domrc && !covered; j++){
            covered = x0[lyt->domrs[l][c + j]];
        }
        c += domrc;
        
        if(covered){
            xl[i] = 1;
        }else{
            xl[i] = 0;
            hv += lyt->layersc[l][i];
        }
        nc += covered;
    }
        
    return hv;
}


static int layers_solution_stoplevel(layers_t * lyt, int k, int * selected, int ** xg){
    int n = lyt->lysizes[0];
    int i, j, l, ni;
    int covered, domrc, c, stop, nc;
    
    for(i = 0; i < n; i++) xg[0][i] = 0;
    for(i = 0; i < k; i++) xg[0][selected[i]] = 1;
        
    if(n == k)
        return 0;
    
    stop = 0;
    for(l = 1; l < n-k+1 && !stop; l++){
        layers_update(lyt, l+1);
        ni = lyt->lysizes[l];
        xg[l] = (int *) malloc(ni * sizeof(int));
        
        c = 0;
        nc = 0;
        for(i = 0; i < ni; i++){
            covered = 0;
            domrc = lyt->domrc[l][i];
            for(j = 0; j < domrc && !covered; j++){
                if(lyt->domrs[l][c + j] >= n)
                    covered = xg[l-1][lyt->domrs[l][c + j] - lyt->clysizes[l-1]];
                else 
                    covered = xg[0][lyt->domrs[l][c + j]];
            }
            c += domrc;
            xg[l][i] = covered;
            
            nc += covered;
        }
        
        if(nc == ni){
            stop = 1;
            break;
        }
    }
    return l;
}




static void layers_free(layers_t * lyt){
    
    int i;
    for(i = 0; i < lyt->t; i++){
        free(lyt->layers[i]);
        free(lyt->layersc[i]);
        free(lyt->domrs[i]);
        free(lyt->domrc[i]);
    }
    free(lyt->layers[lyt->t]);

    free(lyt->layers);
    free(lyt->layersc);
    free(lyt->lysizes);
    free(lyt->clysizes);
    
    free(lyt->ref);
    
    free(lyt->domrs);
    free(lyt->domrc);
    free(lyt->domrssize);
    
    free(lyt->ecdftimes);
    free(lyt->contribstimes);
    free(lyt->domrstimes);
    
    free(lyt);
}


/* ------------------------------------------------------------------- */


void myinitCPLEX(CPXENVptr env, CPXLPptr lp){
    int           status;
   /* Turn on output to the screen */
   CPXsetintparam(env, CPX_PARAM_THREADS, 1);
   CPXchgobjsen(env, lp, CPX_MAX);
   

   CPXsetdblparam(env, CPXPARAM_MIP_Tolerances_MIPGap, 0);   
   CPXsetintparam(env, CPX_PARAM_STARTALG, _startAlg);
   
}


static int setupConstraints(int * dominators, int * rmatcnt, int n, int * rmatbeg, int * rmatind, double * rmatval, int varoffset){
        int sz = rmatcnt[0] + 1; 
    int sum = 0;
    for(int i = 0; i < n-1; i++){
        rmatbeg[i] = sum;
        sum += sz;
        sz = rmatcnt[i+1] + 1;
    }
    rmatbeg[n-1] = sum;
    int cdom = 0;
    int c = 0;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < rmatcnt[i]; j++){
            rmatind[c] = dominators[cdom++];
            rmatval[c] = -1;
            c++;
        }
        rmatind[c] = varoffset + i;
        rmatval[c] = 1;
        c++;
    }

    return c;
}



double layersWarmStart(layers_t * lyt, int k, int ** xg, int startlevel, int * selected, double lb){
    
    Timer_start ();
    
    int n = lyt->lysizes[0];
    int d = lyt->d;
    double * ref = lyt->ref;
    int objsen    = CPX_MAX;
    
    CPXENVptr     env = NULL;
    CPXLPptr      lp = NULL;
    int           status;
    
    double * ones = new_ones(n);
    char ** colnames;
    char * ctypes;
    
    char * senses = achar(NULL, 0, n);
    refill_char(senses, 0, n, 'L');

    char *rowname[1];
    rowname[0] = "card";
    int Rt = 1, Ct = n;
    double * rhs = realloc_adouble(NULL, 0, Rt);
    int * rmatcnt = NULL;
    int * rmatbeg = realloc_aint(NULL, 0, Rt);
    int * rmatind = realloc_aint(NULL, 0, Ct);
    int * tmp = NULL;
    double * rmatval = realloc_adouble(NULL, 0, Ct);
    
    double * x = NULL;
    double * xsol = realloc_adouble(NULL, 0, n);
    double objv = -1;
    double hvRemaining = 0;
    double lbsol = 0;
    
    double tmptime, cpxmodeltime = 0, cpxsolvetime = 0;
    int nskip = 0;
  
    int ni = n;
    int finish = 0;
    
    int l;
    int ncoefs;
    int useWarmStart = 1;
      
    //--
    tmptime = Timer_elapsed_virtual();
    /* Initialize the CPLEX environment */
    env = CPXopenCPLEX (&status);
    lp = CPXcreateprob (env, &status, "layers_prog");
    myinitCPLEX(env, lp);
    cpxmodeltime += Timer_elapsed_virtual() - tmptime;
    //--
    
    layers_update(lyt, startlevel);
    
    for(l = 0; l < n && !finish; l++){
        
        layers_update(lyt, l+1);

        ni = lyt->lysizes[l];
        
        tmptime = Timer_elapsed_virtual();
        //--
        ctypes = binaryCtype(ctypes, ((l > 0) ? lyt->lysizes[l-1] : 0), ni);
        colnames = varNames(l, ni);
        status = CPXnewcols(env, lp, ni, lyt->layersc[l], NULL, NULL, ctypes, colnames);
        //--
        cpxmodeltime += Timer_elapsed_virtual() - tmptime;            
        //--
        tmptime = Timer_elapsed_virtual();
        
        if(l == 0){
            rmatbeg[0] = 0; arange(rmatind, n); fill_ones(rmatval, n);
            rhs[0] = k;
            CPXaddrows(env, lp, 0, 1, n, rhs, senses, rmatbeg, rmatind, ones, colnames, rowname);
                        
        }else{
            
            //--
            rmatbeg = realloc_aint(rmatbeg, Rt, ni);
            rmatind = realloc_aint(rmatind, Ct, (n+1)*ni);
            rmatval = realloc_adouble(rmatval, Ct, (n+1)*ni);
            senses = achar(senses, Rt, ni);
            refill_char(senses, Rt, ni, 'L');
            
            Rt = (ni > Rt) ? ni : Rt;
            Ct = ((n+1)*ni > Rt) ? (n+1)*ni : Rt;

            tmp = lyt->domrs[l];
            rmatcnt = lyt->domrc[l];
            ncoefs = lyt->domrssize[l];
            
            ncoefs = setupConstraints(tmp, rmatcnt, ni, rmatbeg, rmatind, rmatval, lyt->clysizes[l]);
            CPXaddrows(env, lp, 0, ni, ncoefs, NULL, senses, rmatbeg, rmatind, rmatval, NULL, NULL);
            
        }
        cpxmodeltime += Timer_elapsed_virtual() - tmptime;
        //--
        
        
        x = realloc_adouble(x, lyt->clysizes[l], lyt->clysizes[l+1]);
            
        if(l >= startlevel){
            
            if(xg != NULL && useWarmStart){
                rmatbeg[0] = 0;
                arange(rmatind, lyt->clysizes[l+1]);
                
                
                for(int i = 0, h=0; i <= l; i++){
                    for(int j = 0; j < lyt->lysizes[i]; j++, h++)
                        x[h] = (double) xg[i][j];
                }
                
                int effortlevel = 3;
                char * startname = "mystart";
                CPXaddmipstarts(env, lp, 1, lyt->clysizes[l+1], rmatbeg, rmatind, x, &effortlevel, &startname);
                
                int nmipstarts = CPXgetnummipstarts(env, lp);
                useWarmStart = 0;
            }else{
                
                int nmipstarts = CPXgetnummipstarts(env, lp);
                char * startname = "mystart";
                int j; 
                if(!CPXgetmipstartindex(env, lp, startname, &j)){
                    CPXdelmipstarts(env, lp, j, j);
                }
            }
            
            
            
            
            hvRemaining = hvplus(lyt->layers[l+1], d, lyt->lysizes[l+1], ref, 1);
            if(objv > 0){
                objv += layers_hv_level(lyt, l, xsol, &x[lyt->clysizes[l]]);
            }

            tmptime = Timer_elapsed_virtual();
            //--
        
            status = CPXmipopt (env, lp);
            //--
            cpxsolvetime += Timer_elapsed_virtual() - tmptime;
                            
            CPXgetx(env, lp, x, lyt->clysizes[l], lyt->clysizes[l+1]-1);
                            
            CPXgetobjval(env, lp, &objv);
            
            CPXgetx(env, lp, xsol, 0, n-1);
            
            finish = allones(x, ni);
            
            lbsol = hvsubset(lyt->layers[0], n, k, d, ref, xsol);
            
            if(lbsol >= lb){
                lb = lbsol;
            }
        }
        
        free(colnames);

   }
    

    if ( lp != NULL ) {
        status = CPXfreeprob (env, &lp);
    }
    if ( env != NULL ) {
        status = CPXcloseCPLEX (&env);
    }
    
    for(int i = 0; i < n; i++){
        selected[i] = (int) round(xsol[i]);
    }
    
    free(ones);
    free(ctypes);
    
    free(rmatind);
    free(rmatval);
    free(rmatbeg);
    free(senses);
    free(rhs);
    
    free(x);
    free(xsol);
    
    return hvOfSelected(lyt->layers[0], n, d, ref, selected, k);
}



double layersAlg(double * pts, int n, int d, int k, double * ref, int * selected, int startlevel){
    
    Timer_start ();

    layers_t * lyt = layers_init(pts, n, d, ref);
    int ** xg = NULL;
    double hv;
    
    double tmptime;
    double greedytime = 0, gsoltime = 0;
    
    hv = layersWarmStart(lyt, k, xg, startlevel, selected, 0);
    
    layers_free(lyt);
    
    return hv;
}


double layersGreedyStart(double * pts, int n, int d, int k, double * ref, int * selected){
    
    Timer_start ();

    layers_t * lyt = layers_init(pts, n, d, ref);
    int gstoplevel, gistoplevel;
    
    int ** xg = (int **) malloc(n * sizeof(int *));
    xg[0] = (int *) malloc(n * sizeof(int));
    
    int ** xgi = (int **) malloc(n * sizeof(int *));
    xgi[0] = (int *) malloc(n * sizeof(int));
    
    double hv, ghv, gihv;
    
    double tmptime;
    double greedytime = 0, gsoltime = 0;
    
    //--
    tmptime = Timer_elapsed_virtual();
    double * contribs = (double *) malloc(n * sizeof(double));
    int * greedySelected = (int *) malloc(n * sizeof(int));
    ghv = gHSSD(pts, d, n, k, ref, contribs, selected, 1);
    greedytime += Timer_elapsed_virtual() - tmptime;
    //--
    
    //--
    tmptime = Timer_elapsed_virtual();
    gstoplevel = layers_solution_stoplevel(lyt, k, selected, xg);
    gsoltime += Timer_elapsed_virtual() - tmptime;
    //--    
    
    //--
    tmptime = Timer_elapsed_virtual();
    gihv = greedyhss(pts, d, n, k, ref, contribs, selected);
    greedytime += Timer_elapsed_virtual() - tmptime;
    //--
    
    //--
    tmptime = Timer_elapsed_virtual();
    gistoplevel = layers_solution_stoplevel(lyt, k, selected, xgi);
    gsoltime += Timer_elapsed_virtual() - tmptime;
        
    //
    
    if(ghv > gihv){
        hv = layersWarmStart(lyt, k, xg, gstoplevel, selected, ghv);
    }else{
        hv = layersWarmStart(lyt, k, xgi, gistoplevel, selected, gihv);
    }
        
    free(greedySelected);
    free(contribs);

    for(int i = 0; i <= gstoplevel; i++) free(xg[i]);
    free(xg);
    
    
    for(int i = 0; i <= gistoplevel; i++) free(xgi[i]);
    free(xgi);

    
    lyt->contribstime += greedytime;
    
    layers_free(lyt);
        
    return hv;
}


static void
free_and_null (char **ptr)
{
   if ( *ptr != NULL ) {
      free (*ptr);
      *ptr = NULL;
   }
}

