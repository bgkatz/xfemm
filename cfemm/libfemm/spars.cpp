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

#ifdef XFEMM_HAVE_EIGEN
#include <Eigen/Sparse>

#endif

using std::swap;

#ifdef XFEMM_HAVE_EIGEN
// State for the experimental Eigen LDLT direct solver.  The linked-list
// matrix is symmetric with only the upper triangle stored; the row-major
// (diagonal + strictly-upper) CSR arrays are exactly the column-major
// lower triangle of the same matrix, which is what SimplicialLDLT<...,
// Eigen::Lower> consumes.
class CBigLinProbDirect
{
public:
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>, Eigen::Lower> ldlt;
    Eigen::SparseMatrix<double> A;
    std::vector<int> rowStart;
    std::vector<int> cols;
    std::vector<double> vals;
    bool analyzed = false;
    bool factorized = false;
};

// Runtime selection of the solver.  The sparse direct solver (mode 1) is
// the default; set XFEMM_DIRECT=0 in the environment to force the legacy
// PCG solver, or XFEMM_DIRECT=2 for the experimental variant that reuses
// a stale factorization as a CG preconditioner across Newton iterations.
static int directSolveMode()
{
    static int mode = -1;
    if (mode < 0)
    {
        const char *e = getenv("XFEMM_DIRECT");
        mode = (e != nullptr) ? atoi(e) : 1;
        if (mode < 0 || mode > 2) mode = 1;
    }
    return mode;
}
#endif


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

#ifdef XFEMM_HAVE_EIGEN
    delete directSolver;
    directSolver = nullptr;
#endif
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

bool CBigLinProb::SolveDirect(int flag)
{
#ifndef XFEMM_HAVE_EIGEN
    (void)flag;
    return false;
#else
    if (directSolver == nullptr) directSolver = new CBigLinProbDirect;
    CBigLinProbDirect *ds = directSolver;

    printf("Sparse Direct Solver\n");

    // assemble the (diagonal + strictly-upper) CSR arrays into a single
    // sorted structure and copy into an Eigen sparse matrix
    {
        int nnz = csrRowStart[n] + n;
        ds->rowStart.resize(n+1);
        ds->cols.resize(nnz);
        ds->vals.resize(nnz);
        int k = 0;
        for(int i=0; i<n; i++)
        {
            ds->rowStart[i] = k;
            ds->cols[k] = i;
            ds->vals[k] = csrDiag[i];
            k++;
            for(int p=csrRowStart[i]; p<csrRowStart[i+1]; p++)
            {
                ds->cols[k] = csrCol[p];
                ds->vals[k] = csrVal[p];
                k++;
            }
        }
        ds->rowStart[n] = k;

        Eigen::Map<const Eigen::SparseMatrix<double> >
                Amap(n, n, k, ds->rowStart.data(), ds->cols.data(), ds->vals.data());
        ds->A = Amap;
    }

    Eigen::Map<Eigen::VectorXd> vb(b,n), vV(V,n), vR(R,n), vZ(Z,n);

    // Solve A z = b with the current factorization, verify the solution
    // against an explicit residual, and polish it with iterative
    // refinement if needed.  SimplicialLDLT does no pivoting and has no
    // convergence criterion of its own: on a near-indefinite matrix
    // (nonlinear Newton Jacobians can go there via B-H spline overshoot)
    // it reports Success but loses many digits, where PCG would simply
    // have kept iterating.  Refinement recovers those digits at the cost
    // of one mat-vec and one triangular solve per pass.  The candidate
    // lives in Z until it passes the residual test, so a hopeless
    // factorization leaves V (the Newton warm start) intact for the PCG
    // fallback.  Returns false if refinement stalls.
    auto solveWithRefinement = [&]() -> bool
    {
        vZ = ds->ldlt.solve(vb);

        double relres = 0;
        for (int pass=0; pass<6; pass++)
        {
            MultA(Z,U);
            double num=0, den=0;
            for(int j=0; j<n; j++)
            {
                R[j] = b[j]-U[j];
                num += R[j]*R[j];
                den += b[j]*b[j];
            }
            relres = (den==0) ? 0 : sqrt(num/den);
            if (!(relres==relres)) break;       // NaN: hopeless
            if (relres <= Precision) break;
            vZ += ds->ldlt.solve(vR);
        }

        if (relres==relres && relres <= Precision)
        {
            vV = vZ;
            return true;
        }
        fprintf(stderr,"direct solve residual %.1e exceeds Precision "
                       "after refinement; falling back to PCG\n", relres);
        return false;
    };

    if (directSolveMode() == 1 || !ds->factorized)
    {
        // factorize the current matrix and solve directly
        if (!ds->analyzed)
        {
            ds->ldlt.analyzePattern(ds->A);
            ds->analyzed = true;
        }
        ds->ldlt.factorize(ds->A);
        if (ds->ldlt.info() != Eigen::Success)
        {
            fprintf(stderr,"LDLT factorization failed; falling back to PCG\n");
            return false;
        }
        ds->factorized = true;
        return solveWithRefinement();
    }

    // mode 2 with an existing factorization: conjugate gradient on the
    // *current* matrix, preconditioned with the stale factorization.
    // The matrix drifts slowly between Newton iterations, so this
    // usually converges in a handful of iterations.
    int i,iters=0;
    double res,res_o,res_new,er=0,del,rho,pAp;
    const int maxStaleIters = 60;

    vZ = ds->ldlt.solve(vb);
    res_o = Dot(Z,b);
    if (res_o==0) return true;

    if (flag==0) for(i=0; i<n; i++) V[i]=0;

    MultA(V,R);
    for(i=0; i<n; i++) R[i]=b[i]-R[i];

    vZ = ds->ldlt.solve(vR);
    for(i=0; i<n; i++) P[i]=Z[i];
    res = Dot(Z,R);

    do
    {
        MultA(P,U);
        pAp = Dot(P,U);
        del = res/pAp;
        for(i=0; i<n; i++)
        {
            V[i] += del*P[i];
            R[i] -= del*U[i];
        }
        vZ = ds->ldlt.solve(vR);
        res_new = Dot(Z,R);
        rho = res_new/res;
        res = res_new;
        for(i=0; i<n; i++) P[i]=Z[i]+(rho*P[i]);
        er = sqrt(fabs(res/res_o));
        iters++;
    }
    while(er>Precision && iters<maxStaleIters);

    if (er>Precision)
    {
        // stale preconditioner has drifted too far; refactorize the
        // current matrix and finish directly
        ds->ldlt.factorize(ds->A);
        if (ds->ldlt.info() != Eigen::Success)
        {
            fprintf(stderr,"LDLT refactorization failed; falling back to PCG\n");
            return false;
        }
        return solveWithRefinement();
    }

    return true;
#endif
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

#ifdef XFEMM_HAVE_EIGEN
    // experimental direct-solve path, opt-in via XFEMM_DIRECT env var;
    // falls through to PCG if the factorization fails
    if (directSolveMode() != 0)
    {
        if (SolveDirect(flag)) return true;
    }
#endif

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
