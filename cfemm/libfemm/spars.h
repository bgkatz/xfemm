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

#ifndef SPARS_H
#define SPARS_H

#include <vector>

class CEntry
{
public:

    double x;				// value stored in the entry
    int c;					// column that the entry lives in
    CEntry *next;			// pointer to next entry in row;
    CEntry();

private:
};

// opaque state for the optional Eigen-based direct solver (see spars.cpp)
class CBigLinProbDirect;


class CBigLinProb
{
public:

    // data members

    double *V;				// solution
    double *P;				// search direction;
    double *R;				// residual;
    double *U;				// A * P;
    double *Z;
    double *b;				// RHS of linear equation
    CEntry **M;				// pointer to list of matrix entries;
    int n;					// dimensions of the matrix;
    int bdw;				// Optional matrix bandwidth parameter;
    double Precision;		// error tolerance for solution
    double Lambda;			// relaxation factor;

    int *Q; ///< Used by esolver and hsolver.

    // member functions

    // constructor
    CBigLinProb();
    // destructor
    ~CBigLinProb();
    virtual int Create(int d, int bw);	// initialize the problem
    void Put(double v, int p, int q);
    // use to create/set entries in the matrix
    double Get(int p, int q);
    bool PCGSolve(int flag);	// flag==true if guess for V present;
    void MultPC(const double *X, double *Y);
    void AddTo(double v, int p, int q);
    void MultA(double *X, double *Y);
    void SetValue(int i, double x);
    void Periodicity(int i, int j);
    void AntiPeriodicity(int i, int j);
    void Wipe();
    double Dot(double *X, double *Y);
    void ComputeBandwidth();

    // Pointer to the (p,q) entry, created (with value 0) if absent.
    // Entries are never deleted or moved, so the pointer stays valid for
    // the life of the matrix; assembly loops that revisit the same
    // entries every Newton iteration can add through it directly instead
    // of repeating the row walk AddTo does.
    CEntry *Entry(int p, int q);

    // Cache of the iteration-invariant part of a nonlinear system.
    // SaveLinearPart() records the current value of every entry and of
    // b; RestoreLinearPart() puts the matrix back into exactly that
    // state (entries created since are zeroed).  A Newton loop assembles
    // the linear elements once, saves, and then per iteration restores
    // and adds only the nonlinear elements -- see FSolver::Static2D.
    void SaveLinearPart();
    void RestoreLinearPart();

//		CFknDlg *TheView;

private:

    // Flat (CSR-style) copy of the matrix used by the iterative solver
    // kernels.  The linked-list representation in M is convenient for
    // incremental assembly but is slow to traverse; PCGSolve copies it
    // into these contiguous arrays before iterating.  The diagonal is
    // stored separately in csrDiag; csrCol/csrVal hold the strictly
    // upper-triangular entries of row i in
    // [csrRowStart[i], csrRowStart[i+1]).
    std::vector<int> csrRowStart;
    std::vector<int> csrCol;
    std::vector<double> csrVal;
    std::vector<double> csrDiag;

    // copy the linked-list matrix M into the CSR arrays above
    void FlattenMatrix();

    // Column adjacency: for each column q, the rows p<q holding an entry
    // (p,q), as (row, entry pointer) pairs.  Together with row q's own
    // list this gives all structural neighbors of node q, letting
    // SetValue and (Anti)Periodicity visit only the O(degree) entries of
    // a column instead of scanning O(n) rows.  Built lazily by
    // SyncColumnAdjacency; once built it is kept current by the
    // insertion paths (Put/AddTo), and entries are never deleted or
    // moved, so the pointers stay valid for the life of the matrix.
    std::vector<std::vector<std::pair<int,CEntry*> > > colRows;

    // build colRows from the current matrix structure if not built yet
    void SyncColumnAdjacency();

    // SaveLinearPart snapshot: parallel arrays of entry pointer / value,
    // plus a copy of b
    std::vector<CEntry*> linEntries;
    std::vector<double> linVals;
    std::vector<double> linB;

    // collect the rows of all structural entries in column i (both the
    // p<i side from colRows and the p>i side from row i's list) into
    // scratch, excluding row i itself
    void CollectColumnRows(int i, std::vector<int> &scratch);

    // Experimental sparse-direct (Eigen LDLT) solve path.  Only active
    // when the library is built with XFEMM_HAVE_EIGEN and the
    // XFEMM_DIRECT environment variable is set:
    //   XFEMM_DIRECT=1  factorize and solve directly on every call
    //   XFEMM_DIRECT=2  factorize on the first call; later calls run CG
    //                   preconditioned with the (stale) factorization and
    //                   refactorize only if that converges slowly
    // Returns false (falling back to PCG) if unavailable or if the
    // factorization fails.
    bool SolveDirect(int flag);
    CBigLinProbDirect *directSolver = nullptr;
};

#endif
