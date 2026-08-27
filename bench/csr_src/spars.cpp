/*
   This code is a modified version of an algorithm
   forming part of the software program Finite
   Element Method Magnetics (FEMM), authored by
   David Meeker. The original software code is
   subject to the Aladdin Free Public Licence
   version 8, November 18, 1999. For more information
   on FEMM see www.femm.info. This modified version
   is not endorsed in any way by the original
   authors of FEMM.

   This software has been modified to use the C++
   standard template libraries and remove all Microsoft (TM)
   MFC dependent code to allow easier reuse across
   multiple operating system platforms.

   Date Modified: 2011 - 11 - 10
   By: Richard Crozier
   Contact: richard.crozier@yahoo.co.uk
*/

#include "femmcomplex.h"
#include "spars.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <utility>

using std::swap;

#define KLUDGE
#define BENCH_TIMING


CEntry::CEntry()
{
    next=NULL;
    x=0;
    c=0;
}

CBigLinProb::CBigLinProb()
{
    n=0;
    // Best guess for relaxation parameter
    Lambda = 1.5;
}

CBigLinProb::~CBigLinProb()
{
    if (n==0) return;

    int i;
    CEntry *uo,*ui;

    free(b);
    free(P);
    free(R);
    free(V);
    free(U);
    free(Z);

    for(i=0; i<n; i++)
    {
        ui=M[i];
        do
        {
            uo=ui;
            ui=uo->next;
            delete uo;
        }
        while(ui!=NULL);
    }

    free(M);
    free(Q);
    n = 0;
}

int CBigLinProb::Create(int d, int bw)
{
    int i;

    bdw=bw;
    b=(double *)calloc(d,sizeof(double));
    V=(double *)calloc(d,sizeof(double));
    P=(double *)calloc(d,sizeof(double));
    R=(double *)calloc(d,sizeof(double));
    U=(double *)calloc(d,sizeof(double));
    Z=(double *)calloc(d,sizeof(double));

    M=(CEntry **)calloc(d,sizeof(CEntry *));
    n=d;

    for(i=0; i<d; i++)
    {
        M[i] = new CEntry;
        M[i]->c = i;
    }
    Q = (int *)  calloc(d,sizeof(int));

    return 1;
}

void CBigLinProb::Put(double v, int p, int q)
{
    CEntry *e,*l = NULL;

    if (q<p)
        swap(p,q);

    e = M[p];

    while ((e->c < q) && (e->next != NULL))
    {
        l = e;
        e = e->next;
    }

    if (e->c == q)
    {
        e->x = v;
        return;
    }

    CEntry *m = new CEntry;

    if ((e->next == NULL) && (q > e->c))
    {
        e->next = m;
        m->c = q;
        m->x = v;
    }
    else
    {
        l->next = m;
        m->next = e;
        m->c = q;
        m->x = v;
    }
    return;
}

double CBigLinProb::Get(int p, int q)
{
    if (q < p)
    {
        swap(p,q);
    }

    CEntry *e = M[p];
    while ((e->c < q) && (e->next != NULL))
    {
        e = e->next;
    }

    if (e->c == q) return e->x;

    return 0;
}

void CBigLinProb::AddTo(double v, int p, int q)
{
	Put(Get(p,q)+v,p,q);
}

void CBigLinProb::FlattenMatrix()
{
    int i,k;
    CEntry *e;

    csrDiag.resize(n);
    csrRowStart.resize(n+1);

    csrRowStart[0]=0;
    for(i=0; i<n; i++)
    {
        // the first entry of each row is always the diagonal
        csrDiag[i]=M[i]->x;
        k=0;
        for(e=M[i]->next; e!=NULL; e=e->next) k++;
        csrRowStart[i+1]=csrRowStart[i]+k;
    }

    csrCol.resize(csrRowStart[n]);
    csrVal.resize(csrRowStart[n]);
    for(i=0,k=0; i<n; i++)
    {
        for(e=M[i]->next; e!=NULL; e=e->next,k++)
        {
            csrCol[k]=e->c;
            csrVal[k]=e->x;
        }
    }
}

void CBigLinProb::MultA(double *X, double *Y)
{
    int i,k;

    if ((int)csrDiag.size() != n) FlattenMatrix();

    for(i=0; i<n; i++) Y[i]=csrDiag[i]*X[i];

    for(i=0; i<n; i++)
    {
        const double xi=X[i];
        double yi=0;
        for(k=csrRowStart[i]; k<csrRowStart[i+1]; k++)
        {
            yi+=csrVal[k]*X[csrCol[k]];
            Y[csrCol[k]]+=csrVal[k]*xi;
        }
        Y[i]+=yi;
    }
}

double CBigLinProb::Dot(double *X, double *Y)
{
    int i;
    double z;

    for(i=0,z=0; i<n; i++) z+=X[i]*Y[i];

    return z;
}

void CBigLinProb::MultPC(const double *X, double *Y)
{
    // Jacobi preconditioner:
    //	int i;
    // for(i=0;i<n;i++) Y[i]=X[i]/M[i]->x;

    // SSOR preconditioner:
    int i,k;
    double c;

    if ((int)csrDiag.size() != n) FlattenMatrix();

    c= Lambda*(2.-Lambda);
    for(i=0; i<n; i++) Y[i]=X[i]*c;

    // invert Lower Triangle;
    for(i=0; i<n; i++)
    {
        Y[i]/= csrDiag[i];
        const double yl = Y[i] * Lambda;
        for(k=csrRowStart[i]; k<csrRowStart[i+1]; k++)
        {
            Y[csrCol[k]] -= csrVal[k] * yl;
        }
    }

    for(i=0; i<n; i++) Y[i]*=csrDiag[i];

    // invert Upper Triangle
    for(i=n-1; i>=0; i--)
    {
        double yi = Y[i];
        for(k=csrRowStart[i]; k<csrRowStart[i+1]; k++)
        {
            yi -= csrVal[k] * Y[csrCol[k]] * Lambda;
        }
        Y[i] = yi / csrDiag[i];
    }
}

bool CBigLinProb::PCGSolve(int flag)
{
    int i;
    double res,res_o,res_new;
    double er,del,rho,pAp;

#ifdef BENCH_TIMING
    clock_t bench_t0=clock();
    int bench_iters=0;
#endif

    // copy the assembled matrix into flat arrays for fast traversal;
    // must be redone on every call because the nonlinear solvers
    // modify the matrix between calls.
    FlattenMatrix();

    // quick check for most obvious sign of singularity;
    for(i=0; i<n; i++) if(csrDiag[i]==0)
        {
            fprintf(stderr,"singular flag tripped at %i of %i\n", i,n);
            return 0;
        }

    // initialize progress bar;
//	TheView->SetDlgItemText(IDC_FRAME1,"Conjugate Gradient Solver");
//	TheView->m_prog1.SetPos(0);
    printf("Conjugate Gradient Solver\n");

    // residual with V=0
    MultPC(b,Z);
    res_o=Dot(Z,b);
    if(res_o==0) return true;

    // if flag is false, initialize V with zeros;
    if (flag==0) for(i=0; i<n; i++) V[i]=0;

    // form residual;
    MultA(V,R);
    for(i=0; i<n; i++) R[i]=b[i]-R[i];

    // form initial search direction;
    MultPC(R,Z);
    for(i=0; i<n; i++) P[i]=Z[i];
    res=Dot(Z,R);

    // do iteration;
    do
    {
        // step i)
        MultA(P,U);
        pAp=Dot(P,U);
        del=res/pAp;

        for(i=0; i<n; i++)
        {
            // step ii)
            V[i]+=(del*P[i]);

            // step iii)
            R[i]-=(del*U[i]);
        }

        // step iv)
#ifdef BENCH_TIMING
        bench_iters++;
#endif
        MultPC(R,Z);
        res_new=Dot(Z,R);
        rho=res_new/res;
        res=res_new;

        // step v)
        for(i=0; i<n; i++) P[i]=Z[i]+(rho*P[i]);

        // have we converged yet?
        er=sqrt(res/res_o);
//        prg2=(int) (20.*log10(er)/(log10(Precision)));
//        if(prg2>prg1)
//        {
//            prg1=prg2;
//            prg2=(prg1*5);
//            if(prg2>100) prg2=100;
//			TheView->m_prog1.SetPos(prg2);
//			TheView->InvalidateRect(NULL, FALSE);
//			TheView->UpdateWindow();
//        }

    }
    while(er>Precision);

#ifdef BENCH_TIMING
    fprintf(stderr,"BENCH pcg_iters=%i pcg_seconds=%.3f\n",
            bench_iters,(double)(clock()-bench_t0)/CLOCKS_PER_SEC);
#endif

    return true;
}

void CBigLinProb::SetValue(int i, double x)
{
    int k,fst,lst;
    double z;

    if(bdw==0)
    {
        fst=0;
        lst=n;
    }
    else
    {
        fst=i-bdw;
        if (fst<0) fst=0;
        lst=i+bdw;
        if (lst>n) lst=n;
    }

    for(k=fst; k<lst; k++)
    {
        z=Get(k,i);
        if(z!=0)
        {
            b[k]=b[k]-(z*x);
            if(i!=k) Put(0.,k,i);
        }
    }
    b[i]=Get(i,i)*x;
}

void CBigLinProb::Wipe()
{
    int i;
    CEntry *e;

    for(i=0; i<n; i++)
    {
        b[i]=0.;
        e=M[i];
        do
        {
            e->x=0;
            e=e->next;
        }
        while(e!=NULL);
    }
}

void CBigLinProb::AntiPeriodicity(int i, int j)
{
    int k,fst,lst;
    double v1,v2,c;

#ifdef KLUDGE
    int tmpbdw=bdw;
    bdw=0;
#endif

    if (j<i)
        swap(j,i);

    if(bdw==0)
    {
        fst=0;
        lst=n;
    }
    else
    {
        fst=i-bdw;
        if (fst<0) fst=0;
        lst=j+bdw;
        if (lst>n) lst=n;
    }

    for(k=fst; k<lst; k++)
    {
        if((k!=i) && (k!=j))
        {
            v1=Get(k,i);
            v2=Get(k,j);
            if ((v1!=0) || (v2!=0))
            {
                c=(v1-v2)/2.;
                Put(c,k,i);
                Put(-c,k,j);
            }
        }
        if((k==i+bdw) && (k<j-bdw) && (bdw!=0)) k=j-bdw;
    }

    c=0.5*(Get(i,i)+Get(j,j));
    Put(c,i,i);
    Put(c,j,j);

    c=0.5*(b[i]-b[j]);
    b[i]=c;
    b[j]=-c;

#ifdef KLUDGE
    bdw=tmpbdw;
#endif
}

void CBigLinProb::Periodicity(int i, int j)
{
    int k,fst,lst;
    double v1,v2,c;

#ifdef KLUDGE
    int tmpbdw=bdw;
    bdw=0;
#endif

    if (j<i)
        swap(j,i);

    if(bdw==0)
    {
        fst=0;
        lst=n;
    }
    else
    {
        fst=i-bdw;
        if (fst<0) fst=0;
        lst=j+bdw;
        if (lst>n) lst=n;
    }

    for(k=fst; k<lst; k++)
    {
        if((k!=i) && (k!=j))
        {
            v1=Get(k,i);
            v2=Get(k,j);
            if ((v1!=0) || (v2!=0))
            {
                c=(v1+v2)/2.;
                Put(c,k,i);
                Put(c,k,j);
            }
        }
        if((k==i+bdw) && (k<j-bdw) && (bdw!=0)) k=j-bdw;
    }

    c=(Get(i,i)+Get(j,j))/2.;
    Put(c,i,i);
    Put(c,j,j);

    c=0.5*(b[i]+b[j]);
    b[i]=c;
    b[j]=c;

#ifdef KLUDGE
    bdw=tmpbdw;
#endif
}


// a diagnostic routine to check whether that the bandwidth of the
// constructed matrix is actually consistent with a priori bandwidth.
void CBigLinProb::ComputeBandwidth()
{
    CEntry *e;
    int k,bw,maxbw;

    for(maxbw=0,k=0; k<n; k++)
    {
        e=M[k];
        while(e->next != NULL) e=e->next;
        bw=e->c - k;
        if (bw>maxbw) maxbw=bw;
    }

//	MsgBox("Assumed Bandwidth = %i\nActual Bandwidth = %i",bdw,maxbw);

    printf("Assumed Bandwidth = %i\nActual Bandwidth = %i", bdw, maxbw);
}
