################################################
# MEW Deleted text not relevant to source
## Original README in SuiteSparse-5.1.0/CSparse
################################################



CSparse: a Concise Sparse Matrix package.
Copyright (c) 2006-2017, Timothy A. Davis.
http://www.suitesparse.com

Refer to "Direct Methods for Sparse Linear Systems," Timothy A. Davis,
SIAM, Philadelphia, 2006.  No detailed user guide is included in this
package; the user guide is the book itself.

The algorithms contained in CSparse have been chosen with five goals in mind:
(1) they must embody much of the theory behind sparse matrix algorithms,
(2) they must be either asymptotically optimal in their run time and memory
    usage or be fast in practice,
(3) they must be concise so as to be easily understood and short enough to
    print in the book,
(4) they must cover a wide spectrum of matrix operations, and
(5) they must be accurate and robust.
The focus is on direct methods; iterative methods and solvers for
eigenvalue problems are beyond the scope of this package.

No detailed user guide is included in this package; the user guide is the
book itself.  Some indication of how to call the CSparse C routines is given
by the M-files in the MATLAB/CSparse directory.

Complex matrices are not supported, except for methods that operate only
on the nonzero pattern of a matrix.  A complex version of CSparse appears
as a separate package, CXSparse ("Concise Extended Sparse matrix package").

The performance of the sparse factorization methods in CSparse will not be
competitive with UMFPACK or CHOLMOD, but the codes are much more concise and
easy to understand (see the above goals).  Other methods are competitive.


--------------------------------------------------------------------------------
Contents:
--------------------------------------------------------------------------------

README.txt      this file
../cs.c         primary CSparse source files (C only, no MATLAB) (moved from ./Source)

--------------------------------------------------------------------------------
./Doc:          license and change log
--------------------------------------------------------------------------------

ChangeLog       changes in CSparse since first release
License.txt     license

--------------------------------------------------------------------------------
../cs.c:       Primary source code for CSparse (moved from ./Source)
--------------------------------------------------------------------------------

cs_add.c        add sparse matrices
cs_amd.c        approximate minimum degree
cs_chol.c       sparse Cholesky
cs_cholsol.c    x=A\b using sparse Cholesky
cs_compress.c   convert a triplet form to compressed-column form
cs_counts.c     column counts for Cholesky and QR
cs_cumsum.c     cumulative sum
cs_dfs.c        depth-first-search
cs_dmperm.c     Dulmage-Mendelsohn permutation
cs_droptol.c    drop small entries from a sparse matrix
cs_dropzeros.c  drop zeros from a sparse matrix
cs_dupl.c       remove (and sum) duplicates
cs_entry.c      add an entry to a triplet matrix
cs_ereach.c     nonzero pattern of Cholesky L(k,:) from etree and triu(A(:,k))
cs_etree.c      find elimination tree
cs_fkeep.c      drop entries from a sparse matrix
cs_gaxpy.c      sparse matrix times dense matrix
cs.h            include file for CSparse
cs_happly.c     apply Householder reflection
cs_house.c      compute Householder reflection
cs_ipvec.c      x(p)=b
cs_leaf.c       determine if j is a leaf of the skeleton matrix and find lca
cs_load.c       load a sparse matrix from a file
cs_lsolve.c     x=L\b
cs_ltsolve.c    x=L'\b
cs_lu.c         sparse LU factorization
cs_lusol.c      x=A\b using sparse LU factorization
cs_malloc.c     memory manager
cs_maxtrans.c   maximum transveral (permutation for zero-free diagonal)
cs_multiply.c   sparse matrix multiply
cs_norm.c       sparse matrix norm
cs_permute.c    permute a sparse matrix
cs_pinv.c       invert a permutation vector
cs_post.c       postorder an elimination tree
cs_print.c      print a sparse matrix
cs_pvec.c       x=b(p)
cs_qr.c         sparse QR
cs_qrsol.c      solve a least-squares problem
cs_randperm.c   random permutation
cs_reach.c      find nonzero pattern of x=L\b for sparse L and b
cs_scatter.c    scatter a sparse vector
cs_scc.c        strongly-connected components
cs_schol.c      symbolic Cholesky
cs_spsolve.c    x=G\b where G, x, and b are sparse, and G upper/lower triangular
cs_sqr.c        symbolic QR (also can be used for LU)
cs_symperm.c    symmetric permutation of a sparse matrix
cs_tdfs.c       depth-first-search of a tree
cs_transpose.c  transpose a sparse matrix
cs_updown.c     sparse rank-1 Cholesky update/downate
cs_usolve.c     x=U\b
cs_util.c       various utilities (allocate/free matrices, workspace, etc)
cs_utsolve.c    x=U'\b

