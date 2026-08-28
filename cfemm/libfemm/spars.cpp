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

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <utility>

using std::swap;


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

    // keep the column adjacency current if it has been built
    if (!colRows.empty()) colRows[q].push_back(std::make_pair(p,m));
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
    // single row walk; the old implementation was Put(Get(p,q)+v,p,q),
    // which searched the row twice
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
        e->x += v;
        return;
    }

    CEntry *m = new CEntry;
    m->c = q;
    m->x = v;

    if ((e->next == NULL) && (q > e->c))
    {
        e->next = m;
    }
    else
    {
        l->next = m;
        m->next = e;
    }

    // keep the column adjacency current if it has been built
    if (!colRows.empty()) colRows[q].push_back(std::make_pair(p,m));
}

void CBigLinProb::SyncColumnAdjacency()
{
    // once built, colRows is maintained incrementally by the insertion
    // paths (Put/AddTo), so it never goes stale
    if (!colRows.empty()) return;

    colRows.resize(n);
    for(int p=0; p<n; p++)
    {
        for(CEntry *e=M[p]->next; e!=NULL; e=e->next)
        {
            colRows[e->c].push_back(std::make_pair(p,e));
        }
    }
}

void CBigLinProb::CollectColumnRows(int i, std::vector<int> &scratch)
{
    // rows above the diagonal holding an entry in column i...
    for(size_t k=0; k<colRows[i].size(); k++)
    {
        scratch.push_back(colRows[i][k].first);
    }
    // ...and, via symmetry, the columns of row i's own entries
    for(CEntry *e=M[i]->next; e!=NULL; e=e->next)
    {
        scratch.push_back(e->c);
    }
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

    return true;
}

void CBigLinProb::SetValue(int i, double x)
{
    // visit just the structural entries of column i (rows above the
    // diagonal) and row i (columns right of the diagonal) instead of
    // scanning every row within the bandwidth
    SyncColumnAdjacency();

    for(size_t k=0; k<colRows[i].size(); k++)
    {
        CEntry *e=colRows[i][k].second;
        if (e->x != 0)
        {
            b[colRows[i][k].first] -= e->x * x;
            e->x = 0.;
        }
    }
    for(CEntry *e=M[i]->next; e!=NULL; e=e->next)
    {
        if (e->x != 0)
        {
            b[e->c] -= e->x * x;
            e->x = 0.;
        }
    }
    b[i]=M[i]->x * x;
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
    double v1,v2,c;

    if (j<i)
        swap(j,i);

    // visit just the rows holding a structural entry in column i or
    // column j.  Older versions scanned every row of the matrix here
    // (a banded scan was defeated by a KLUDGE that forced bdw=0,
    // because earlier (Anti)Periodicity calls create entries outside
    // the a-priori bandwidth which a banded scan would then miss).
    SyncColumnAdjacency();

    std::vector<int> ks;
    CollectColumnRows(i,ks);
    CollectColumnRows(j,ks);
    std::sort(ks.begin(),ks.end());
    ks.erase(std::unique(ks.begin(),ks.end()),ks.end());

    for(size_t idx=0; idx<ks.size(); idx++)
    {
        int k=ks[idx];
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
    }

    c=0.5*(Get(i,i)+Get(j,j));
    Put(c,i,i);
    Put(c,j,j);

    c=0.5*(b[i]-b[j]);
    b[i]=c;
    b[j]=-c;
}

void CBigLinProb::Periodicity(int i, int j)
{
    double v1,v2,c;

    if (j<i)
        swap(j,i);

    // see the comment in AntiPeriodicity
    SyncColumnAdjacency();

    std::vector<int> ks;
    CollectColumnRows(i,ks);
    CollectColumnRows(j,ks);
    std::sort(ks.begin(),ks.end());
    ks.erase(std::unique(ks.begin(),ks.end()),ks.end());

    for(size_t idx=0; idx<ks.size(); idx++)
    {
        int k=ks[idx];
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
    }

    c=(Get(i,i)+Get(j,j))/2.;
    Put(c,i,i);
    Put(c,j,j);

    c=0.5*(b[i]+b[j]);
    b[i]=c;
    b[j]=c;
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
