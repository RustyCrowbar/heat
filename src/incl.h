#pragma once


    //
    //  Purpose:
    //
    //    ADJUST_BACKWARD_EULER adjusts the system for the backward Euler term.
    //
    //  Discussion:
    //
    //    The input linear system
    //
    //      A * U = F
    //
    //    is appropriate for the equation
    //
    //      -Uxx - Uyy = RHS
    //
    //    We need to modify the matrix A and the right hand side F to
    //    account for the approximation of the time derivative in
    //
    //      Ut - Uxx - Uyy = RHS
    //
    //    by the backward Euler approximation:
    //
    //      Ut approximately equal to ( U - Uold ) / dT
    //
    //  Parameters:
    //
    //    Input, int NODE_NUM, the number of nodes.
    //
    //    Input, double NODE_XY[2*NODE_NUM], the coordinates of nodes.
    //
    //    Input, int NNODES, the number of nodes used to form one element.
    //
    //    Input, int ELEMENT_NUM, the number of elements.
    //
    //    Input, int ELEMENT_NODE[NNODES*ELEMENT_NUM];
    //    ELEMENT_NODE(I,J) is the global index of local node I in element J.
    //
    //    Input, int QUAD_NUM, the number of quadrature points used in assembly.
    //
    //    Input, double WQ[QUAD_NUM], quadrature weights.
    //
    //    Input, double XQ[QUAD_NUM*ELEMENT_NUM],
    //    YQ[QUAD_NUM*ELEMENT_NUM], the
    //    coordinates of the quadrature points in each element.
    //
    //    Input, double ELEMENT_AREA[ELEMENT_NUM], the area of elements.
    //
    //    Input, int IB, the half-bandwidth of the matrix.
    //
    //    Input, double TIME, the current time.
    //
    //    Input, double TIME_STEP_SIZE, the size of the time step.
    //
    //    Input, double U_OLD[NODE_NUM], the finite element
    //    coefficients for the solution at the previous time.
    //
    //    Input/output, double A[(3*IB+1)*NODE_NUM], the NODE_NUM
    //    by NODE_NUM coefficient matrix, stored in a compressed format.
    //
    //    Input/output, double F[NODE_NUM], the right hand side.
    //
void adjust_backward_euler(int node_num, double node_xy[], int nnodes,
	int element_num, int element_node[], int quad_num, double wq[],
	double xq[], double yq[], double element_area[], int ib, double time,
	double time_step_size, double u_old[], double a[], double f[]);

    //
    //  Purpose:
    //
    //    ADJUST_BOUNDARY modifies the linear system for boundary conditions.
    //
    //  Parameters:
    //
    //    Input, int NODE_NUM, the number of nodes.
    //
    //    Input, double NODE_XY[2*NODE_NUM], the coordinates of nodes.
    //
    //    Input, int NODE_BOUNDARY[NODE_NUM], is
    //    0, if a node is an interior node;
    //    1, if a node is a Dirichlet boundary node.
    //
    //    Input, int IB, the half-bandwidth of the matrix.
    //
    //    Input, double TIME, the current time.
    //
    //    Input/output, double A[(3*IB+1)*NODE_NUM], the NODE_NUM by NODE_NUM
    //    coefficient matrix, stored in a compressed format.
    //    On output, A has been adjusted for boundary conditions.
    //
    //    Input/output, double F[NODE_NUM], the right hand side.
    //    On output, F has been adjusted for boundary conditions.
    //
void adjust_boundary(int node_num, double node_xy[], int node_boundary[],
	int ib, double time, double a[], double f[]);

    //
    //  Purpose:
    //
    //    AREA_SET sets the area of each element.
    //
    //  Discussion:
    //
    //    The areas of the elements are needed in order to adjust
    //    the integral estimates produced by the quadrature formulas.
    //
    //  Parameters:
    //
    //    Input, int NODE_NUM, the number of nodes.
    //
    //    Input, double NODE_XY[2*NODE_NUM], the
    //    coordinates of the nodes.
    //
    //    Input, int NNODES, the number of local nodes per element.
    //
    //    Input, int ELEMENT_NUM, the number of elements.
    //
    //    Input, int ELEMENT_NODE[NNODES*ELEMENT_NUM];
    //    ELEMENT_NODE(I,J) is the global index of local node I in element J.
    //
    //    Output, double ELEMENT_AREA[ELEMENT_NUM], the area of elements.
    //

void area_set(int node_num, double node_xy[], int nnodes, int element_num,
	int element_node[], double element_area[]);

    //
    //  Purpose:
    //
    //    ASSEMBLE assembles the matrix and right-hand side using piecewise quadratics.
    //
    //  Discussion:
    //
    //    The matrix is known to be banded.  A special matrix storage format
    //    is used to reduce the space required.  Details of this format are
    //    discussed in the routine DGB_FA.
    //
    //  Parameters:
    //
    //    Input, int NODE_NUM, the number of nodes.
    //
    //    Input, double NODE_XY[2*NODE_NUM], the X and Y coordinates of nodes.
    //
    //    Input, int NNODES, the number of nodes used to form one element.
    //
    //    Input, int ELEMENT_NUM, the number of elements.
    //
    //    Input, int ELEMENT_NODE[NNODES*ELEMENT_NUM]; ELEMENT_NODE(I,J) is the global
    //    index of local node I in element J.
    //
    //    Input, int QUAD_NUM, the number of quadrature points used in assembly.
    //
    //    Input, double WQ[QUAD_NUM], quadrature weights.
    //
    //    Input, double XQ[QUAD_NUM*ELEMENT_NUM], YQ[QUAD_NUM*ELEMENT_NUM], the X and Y
    //    coordinates of the quadrature points in each element.
    //
    //    Input, double ELEMENT_AREA[ELEMENT_NUM], the area of each element.
    //
    //    Input, int IB, the half-bandwidth of the matrix.
    //
    //    Input, double TIME, the current time.
    //
    //    Output, double A[(3*IB+1)*NODE_NUM], the NODE_NUM by NODE_NUM
    //    coefficient matrix, stored in a compressed format.
    //
    //    Output, double F[NODE_NUM], the right hand side.
    //
    //  Local parameters:
    //
    //    Local, double BI, DBIDX, DBIDY, the value of some basis function
    //    and its first derivatives at a quadrature point.
    //
    //    Local, double BJ, DBJDX, DBJDY, the value of another basis
    //    function and its first derivatives at a quadrature point.
    //
void assemble(int node_num, double node_xy[], int nnodes, int element_num,
	int element_node[], int quad_num, double wq[], double xq[], double yq[],
	double element_area[], int ib, double time, double a[], double f[]);

    //
    //  Purpose:
    //
    //    BANDWIDTH determines the bandwidth of the coefficient matrix.
    //
    //  Parameters:
    //
    //    Input, int NNODES, the number of local nodes per element.
    //
    //    Input, int ELEMENT_NUM, the number of elements.
    //
    //    Input, int ELEMENT_NODE[NNODES*ELEMENT_NUM]; ELEMENT_NODE(I,J) is the global
    //    index of local node I in element J.
    //
    //    Input, int NODE_NUM, the number of nodes.
    //
    //    Output, int BANDWIDTH, the half bandwidth of the matrix.
    //
int bandwidth(int nnodes, int element_num, int element_node[], int node_num);

    //
    //  Purpose:
    //
    //    COMPARE compares the exact and computed solution at the nodes.
    //
    //  Discussion:
    //
    //    This is a rough comparison, done only at the nodes.  Such a pointwise
    //    comparison is easy, because the value of the finite element
    //    solution is exactly the value of the finite element coefficient
    //    associated with that node.
    //
    //  Parameters:
    //
    //    Input, int NODE_NUM, the number of nodes.
    //
    //    Input, double NODE_XY[2*NODE_NUM], the nodes.
    //
    //    Input, double TIME, the current time.
    //
    //    Input, double F[NUNK], the solution vector of the finite
    //    element system.
    //
void compare(int node_num, double node_xy[], double time, double u[]);

    //
    //  Purpose:
    //
    //    DGB_FA performs a LINPACK-style PLU factorization of an DGB matrix.
    //
    //  Discussion:
    //
    //    The DGB storage format is used for an M by N banded matrix, with lower bandwidth ML
    //    and upper bandwidth MU.  Storage includes room for ML extra superdiagonals,
    //    which may be required to store nonzero entries generated during Gaussian
    //    elimination.
    //
    //    The original M by N matrix is "collapsed" downward, so that diagonals
    //    become rows of the storage array, while columns are preserved.  The
    //    collapsed array is logically 2*ML+MU+1 by N.
    //
    //    The two dimensional array can be further reduced to a one dimensional
    //    array, stored by columns.
    //
    //  Reference:
    //
    //    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
    //    LINPACK User's Guide,
    //    SIAM, 1979
    //
    //  Parameters:
    //
    //    Input, int N, the order of the matrix.
    //    N must be positive.
    //
    //    Input, int ML, MU, the lower and upper bandwidths.
    //    ML and MU must be nonnegative, and no greater than N-1.
    //
    //    Input/output, double A[(2*ML+MU+1)*N], the matrix in band storage.
    //    On output, A has been overwritten by the LU factors.
    //
    //    Output, int PIVOT[N], the pivot vector.
    //
    //    Output, int SGB_FA, singularity flag.
    //    0, no singularity detected.
    //    nonzero, the factorization failed on the INFO-th step.
    //
int dgb_fa(int n, int ml, int mu, double a[], int pivot[]);

    //
    //  Purpose:
    //
    //    DGB_PRINT_SOME prints some of a DGB matrix.
    //
    //  Discussion:
    //
    //    The DGB storage format is used for an M by N banded matrix, with lower bandwidth ML
    //    and upper bandwidth MU.  Storage includes room for ML extra superdiagonals,
    //    which may be required to store nonzero entries generated during Gaussian
    //    elimination.
    //
    //    The original M by N matrix is "collapsed" downward, so that diagonals
    //    become rows of the storage array, while columns are preserved.  The
    //    collapsed array is logically 2*ML+MU+1 by N.
    //
    //    The two dimensional array can be further reduced to a one dimensional
    //    array, stored by columns.
    //
    //  Parameters:
    //
    //    Input, int M, the number of rows of the matrix.
    //    M must be positive.
    //
    //    Input, int N, the number of columns of the matrix.
    //    N must be positive.
    //
    //    Input, int ML, MU, the lower and upper bandwidths.
    //    ML and MU must be nonnegative, and no greater than min(M,N)-1..
    //
    //    Input, double A[(2*ML+MU+1)*N], the SGB matrix.
    //
    //    Input, int ILO, JLO, IHI, JHI, designate the first row and
    //    column, and the last row and column to be printed.
    //
    //    Input, char *TITLE, a title to print.
void dgb_print_some(int m, int n, int ml, int mu, double a[],
	int ilo, int jlo, int ihi, int jhi, char *title);

    //
    //  Purpose:
    //
    //    DGB_SL solves a system factored by DGB_FA.
    //
    //  Discussion:
    //
    //    The DGB storage format is used for an M by N banded matrix, with lower bandwidth ML
    //    and upper bandwidth MU.  Storage includes room for ML extra superdiagonals,
    //    which may be required to store nonzero entries generated during Gaussian
    //    elimination.
    //
    //    The original M by N matrix is "collapsed" downward, so that diagonals
    //    become rows of the storage array, while columns are preserved.  The
    //    collapsed array is logically 2*ML+MU+1 by N.
    //
    //    The two dimensional array can be further reduced to a one dimensional
    //    array, stored by columns.
    //
    //  Reference:
    //
    //    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
    //    LINPACK User's Guide,
    //    SIAM, 1979
    //
    //  Parameters:
    //
    //    Input, int N, the order of the matrix.
    //    N must be positive.
    //
    //    Input, int ML, MU, the lower and upper bandwidths.
    //    ML and MU must be nonnegative, and no greater than N-1.
    //
    //    Input, double A[(2*ML+MU+1)*N], the LU factors from DGB_FA.
    //
    //    Input, int PIVOT[N], the pivot vector from DGB_FA.
    //
    //    Input, double B[N], the right hand side vector.
    //
    //    Input, int JOB.
    //    0, solve A * x = b.
    //    nonzero, solve A' * x = b.
    //
    //    Output, double DGB_SL[N], the solution.
    //
double *dgb_sl(int n, int ml, int mu, double a[], int pivot[], double b[],
	int job);

    //
    //  Purpose:
    //
    //    ELEMENT_WRITE writes the elements to a file.
    //
    //  Parameters:
    //
    //    Input, int NNODES, the number of nodes used to form one element.
    //
    //    Input, int ELEMENT_NUM, the number of elements.
    //
    //    Input, int ELEMENT_NODE[NNODES*ELEMENT_NUM]; ELEMENT_NODE(I,J) is the global
    //    index of local node I in element J.
    //
    //    Input, char *OUTPUT_FILENAME, the name of the file
    //    in which the data should be stored.
    //
void element_write(int nnodes, int element_num, int element_node[],
	char *triangulation_txt_file_name);

    //
    //  Purpose:
    //
    //    ERRORS calculates the error in the L2 and H1-seminorm.
    //
    //  Discussion:
    //
    //    This routine uses a 13 point quadrature rule in each element,
    //    in order to estimate the values of
    //
    //      EL2 = Sqrt ( Integral ( U(x,y) - Uh(x,y) )**2 dx dy )
    //
    //      EH1 = Sqrt ( Integral ( Ux(x,y) - Uhx(x,y) )**2 +
    //                            ( Uy(x,y) - Uhy(x,y) )**2 dx dy )
    //
    //    Here U is the exact solution, and Ux and Uy its spatial derivatives,
    //    as evaluated by a user-supplied routine.
    //
    //    Uh, Uhx and Uhy are the computed solution and its spatial derivatives,
    //    as specified by the computed finite element solution.
    //
    //  Parameters:
    //
    //    Input, double ELEMENT_AREA[ELEMENT_NUM], the area of each element.
    //
    //    Input, int ELEMENT_NODE[NNODES*ELEMENT_NUM]; ELEMENT_NODE(I,J) is the global
    //    index of local node I in element J.
    //
    //    Input, double NODE_XY[2*NODE_NUM], the X and Y coordinates of nodes.
    //
    //    Input, double U[NUNK], the coefficients of the solution.
    //
    //    Input, int ELEMENT_NUM, the number of elements.
    //
    //    Input, int NNODES, the number of nodes used to form one element.
    //
    //    Input, int NODE_NUM, the number of nodes.
    //
    //    Input, real ( kind = 8 ) TIME, the current time.
    //
    //    Output, double precision *EL2, the L2 error.
    //
    //    Output, double precision *EH1, the H1 seminorm error.
    //
    //  Local Parameters:
    //
    //    Local, double AR, the weight for a given quadrature point
    //    in a given element.
    //
    //    Local, double BI, DBIDX, DBIDY, a basis function and its first
    //    derivatives evaluated at a particular quadrature point.
    //
    //    Local, double EH1, the H1 seminorm error.
    //
    //    Local, double EL2, the L2 error.
    //
    //    Local, int NQE, the number of points in the quadrature rule.
    //    This is actually fixed at 13.
    //
    //    Local, double UEX, UEXX, UEXY, the exact solution and its first
    //    derivatives evaluated at a particular quadrature point.
    //
    //    Local, double UH, UHX, UHY, the computed solution and its first
    //    derivatives evaluated at a particular quadrature point.
    //
    //    Local, double WQE[NQE], stores the quadrature weights.
    //
    //    Local, double X, Y, the coordinates of a particular
    //    quadrature point.
    //
    //    Local, double XQE[NQE], YQE[NQE], stores the location
    //    of quadrature points in a given element.
    //
void errors(double element_area[], int element_node[], double node_xy[],
	double u[], int element_num, int nnodes, int node_num, double time,
	double *el2, double *eh1);

    //
    //  Purpose:
    //
    //    EXACT_U calculates the exact solution and its first derivatives.
    //
    //  Discussion:
    //
    //    It is assumed that the user knows the exact solution and its
    //    derivatives.  This, of course, is NOT true for a real computation.
    //    But for this code, we are interested in studying the convergence
    //    behavior of the approximations, and so we really need to assume
    //    we know the correct solution.
    //
    //    As a convenience, this single routine is used for several purposes:\
    //
    //    * it supplies the initial value function H(X,Y,T);
    //    * it supplies the boundary value function G(X,Y,T);
    //    * it is used by the COMPARE routine to make a node-wise comparison
    //      of the exact and approximate solutions.
    //    * it is used by the ERRORS routine to estimate the integrals of
    //      the L2 and H1 errors of approximation.
    //
    //    DUDX and DUDY are only needed for the ERRORS calculation.
    //
    //  Parameters:
    //
    //    Input, int NODE_NUM, the number of nodes at which
    //    a value is desired.
    //
    //    Input, double NODE_XY(2,NODE_NUM), the coordinates of
    //    the points where a value is desired.
    //
    //    Input, double TIME, the current time.
    //
    //    Output, double U[NODE_NUM], the exact solution.
    //
    //    Output, double DUDX[NODE_NUM], DUDY[NODE_NUM],
    //    the X and Y derivatives of the exact solution.
    //
void exact_u(int node_num, double node_xy[], double time,
	double u_exact[], double dudx_exact[], double dudy_exact[]);

    //
    //  Purpose:
    //
    //    FILE_NAME_INC increments a partially numeric file name.
    //
    //  Discussion:
    //
    //    It is assumed that the digits in the name, whether scattered or
    //    connected, represent a number that is to be increased by 1 on
    //    each call.  If this number is all 9's on input, the output number
    //    is all 0's.  Non-numeric letters of the name are unaffected.
    //
    //    If the input string contains no digits, a blank string is returned.
    //
    //    If a blank string is input, then an error condition results.
    //
    //  Example:
    //
    //      Input            Output
    //      -----            ------
    //      "a7to11.txt"     "a7to12.txt"  (typical case.  Last digit incremented)
    //      "a7to99.txt"     "a8to00.txt"  (last digit incremented, with carry.)
    //      "a9to99.txt"     "a0to00.txt"  (wrap around)
    //      "cat.txt"        " "           (no digits to increment)
    //      " "              STOP!         (error)
    //
    //  Parameters:
    //
    //    Input/output, character *FILE_NAME, (a pointer to) the character string
    //    to be incremented.
    //
void file_name_inc(char *file_name);

    //
    //  Purpose:
    //
    //    GRID_T6 produces a grid of pairs of 6 node triangles.
    //
    //  Example:
    //
    //    Input:
    //
    //      NX = 4, NY = 3
    //
    //    Output:
    //
    //      ELEMENT_NODE =
    //         1,  3, 15,  2,  9,  8;
    //        17, 15,  3, 16,  9, 10;
    //         3,  5, 17,  4, 11, 10;
    //        19, 17,  5, 18, 11, 12;
    //         5,  7, 19,  6, 13, 12;
    //        21, 19,  7, 20, 13, 14;
    //        15, 17, 29, 16, 23, 22;
    //        31, 29, 17, 30, 23, 24;
    //        17, 19, 31, 18, 25, 24;
    //        33, 31, 19, 32, 25, 26;
    //        19, 21, 33, 20, 27, 26;
    //        35, 33, 21, 34, 27, 28.
    //
    //  Diagram:
    //
    //   29-30-31-32-33-34-35
    //    |\ 8  |\10  |\12  |
    //    | \   | \   | \   |
    //   22 23 24 25 26 27 28
    //    |   \ |   \ |   \ |
    //    |  7 \|  9 \| 11 \|
    //   15-16-17-18-19-20-21
    //    |\ 2  |\ 4  |\ 6  |
    //    | \   | \   | \   |
    //    8  9 10 11 12 13 14
    //    |   \ |   \ |   \ |
    //    |  1 \|  3 \|  5 \|
    //    1--2--3--4--5--6--7
    //
    //  Parameters:
    //
    //    Input, int NX, NY, controls the number of elements along the
    //    X and Y directions.  The number of elements will be
    //    2 * ( NX - 1 ) * ( NY - 1 ).
    //
    //    Input, int NNODES, the number of local nodes per element.
    //
    //    Input, int ELEMENT_NUM, the number of elements.
    //
    //    Output, int ELEMENT_NODE[NNODES*ELEMENT_NUM];
    //    ELEMENT_NODE(I,J) is the index of the I-th node of the J-th element.
    //
void grid_t6(int nx, int ny, int nnodes, int element_num, int element_node[]);

    //
    //  Purpose:
    //
    //    I4_MAX returns the maximum of two I4's.
    //
    //  Parameters:
    //
    //    Input, int I1, I2, are two ints to be compared.
    //
    //    Output, int I4_MAX, the larger of I1 and I2.
    //
int i4_max(int i1, int i2);

    //
    //  Purpose:
    //
    //    I4_MIN returns the smaller of two I4's.
    //
    //  Parameters:
    //
    //    Input, int I1, I2, two ints to be compared.
    //
    //    Output, int I4_MIN, the smaller of I1 and I2.
    //
int i4_min (int i1, int i2);

    //
    //  Purpose:
    //
    //    I4VEC_PRINT_SOME prints "some" of an I4VEC.
    //
    //  Discussion:
    //
    //    The user specifies MAX_PRINT, the maximum number of lines to print.
    //
    //    If N, the size of the vector, is no more than MAX_PRINT, then
    //    the entire vector is printed, one entry per line.
    //
    //    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
    //    followed by a line of periods suggesting an omission,
    //    and the last entry.
    //
    //  Parameters:
    //
    //    Input, int N, the number of entries of the vector.
    //
    //    Input, int A[N], the vector to be printed.
    //
    //    Input, int MAX_PRINT, the maximum number of lines to print.
    //
    //    Input, char *TITLE, an optional title.
    //
void i4vec_print_some(int n, int a[], int max_print, char *title);

    //
    //  Purpose:
    //
    //    NODE_BOUNDARY_SET assigns an unknown value index at each node.
    //
    //  Discussion:
    //
    //    Every node is assigned a value which indicates whether it is
    //    an interior node, or a boundary node.
    //
    //  Example:
    //
    //    On a simple 5 by 5 grid, where the nodes are numbered starting
    //    at the lower left, and increasing in X first, we would have the
    //    following values of NODE_BOUNDARY:
    //
    //       1  1  1  1  1
    //       1  0  0  0  1
    //       1  0  0  0  1
    //       1  0  0  0  1
    //       1  1  1  1  1
    //
    //  Parameters:
    //
    //    Input, int NX, NY, the number of elements in the X and Y directions.
    //
    //    Input, int NODE_NUM, the number of nodes.
    //
    //    Output, int NODE_BOUNDARY[NODE_NUM], is
    //    0, if a node is an interior node;
    //    1, if a node is a Dirichlet boundary node.
    //
int *node_boundary_set(int nx, int ny, int node_num);

    //
    //  Purpose:
    //
    //    NODES_PLOT plots a pointset.
    //
    //  Parameters:
    //
    //    Input, char *FILE_NAME, the name of the file to create.
    //
    //    Input, int NODE_NUM, the number of nodes.
    //
    //    Input, double NODE_XY[2*NODE_NUM], the nodes.
    //
    //    Input, bool NODE_LABEL, is TRUE if the nodes are to be labeled.
    //
void nodes_plot(char *file_name, int node_num, double node_xy[],
	bool node_label);

    //
    //  Purpose:
    //
    //    NODES_WRITE writes the nodes to a file.
    //
    //  Parameters:
    //
    //    Input, int NODE_NUM, the number of nodes.
    //
    //    Input, double NODE_XY[2*NODE_NUM], the X and Y coordinates of nodes.
    //
    //    Input, char *OUTPUT_FILENAME, the name of the file
    //    in which the data should be stored.
    //
void nodes_write(int node_num, double node_xy[], char *output_filename);

    //
    //  Purpose:
    //
    //    QBF evaluates the quadratic basis functions.
    //
    //  Discussion:
    //
    //    This routine assumes that the "midpoint" nodes are, in fact,
    //    exactly the average of the two extreme nodes.  This is NOT true
    //    for a general quadratic triangular element.
    //
    //    Assuming this property of the midpoint nodes makes it easy to
    //    determine the values of (R,S) in the reference element that
    //    correspond to (X,Y) in the physical element.
    //
    //    Once we know the (R,S) coordinates, it's easy to evaluate the
    //    basis functions and derivatives.
    //
    //  The physical element T6:
    //
    //    In this picture, we don't mean to suggest that the bottom of
    //    the physical triangle is horizontal.  However, we do assume that
    //    each of the sides is a straight line, and that the intermediate
    //    points are exactly halfway on each side.
    //
    //    |
    //    |
    //    |        3
    //    |       / \
    //    |      /   \
    //    Y     6     5
    //    |    /       \
    //    |   /         \
    //    |  1-----4-----2
    //    |
    //    +--------X-------->
    //
    //  Reference element T6:
    //
    //    In this picture of the reference element, we really do assume
    //    that one side is vertical, one horizontal, of length 1.
    //
    //    |
    //    |
    //    1  3
    //    |  |\
    //    |  | \
    //    S  6  5
    //    |  |   \
    //    |  |    \
    //    0  1--4--2
    //    |
    //    +--0--R--1-------->
    //
    //  Parameters:
    //
    //    Input, double X, Y, the (global) coordinates of the point
    //    at which the basis function is to be evaluated.
    //
    //    Input, int ELEMENT, the index of the element which contains the point.
    //
    //    Input, int INODE, the local index, between 0 and 5, that
    //    specifies which basis function is to be evaluated.
    //
    //    Input, double NODE_XY[2*NODE_NUM], the nodes.
    //
    //    Input, int ELEMENT_NODE[NNODES*ELEMENT_NUM];
    //    ELEMENT_NODE(I,J) is the global index of local node I in element J.
    //
    //    Input, int ELEMENT_NUM, the number of elements.
    //
    //    Input, int NNODES, the number of nodes used to form one element.
    //
    //    Input, int NODE_NUM, the number of nodes.
    //
    //    Output, double *B, *DBDX, *DBDY, the value of the basis function
    //    and its X and Y derivatives at (X,Y).
    //
void qbf(double x, double y, int element, int inode, double node_xy[],
	int element_node[], int element_num, int nnodes,
	int node_num, double *bb, double *bx, double *by);

    //
    //  Purpose:
    //
    //    QUAD_A sets the quadrature rule for assembly.
    //
    //  Parameters:
    //
    //    Input, double NODE_XY[2*NODE_NUM], the nodes.
    //
    //    Input, int ELEMENT_NODE[NNODES*ELEMENT_NUM];
    //    ELEMENT_NODE(I,J) is the global index of local node I in element J.
    //
    //    Input, int ELEMENT_NUM, the number of elements.
    //
    //    Input, int NODE_NUM, the number of nodes.
    //
    //    Input, int NNODES, the number of nodes used to form one element.
    //
    //    Output, double WQ[3], quadrature weights.
    //
    //    Output, double XQ[3*ELEMENT_NUM], YQ[3*ELEMENT_NUM], the
    //    coordinates of the quadrature points in each element.
    //
void quad_a(double node_xy[], int element_node[], int element_num,
	int node_num, int nnodes, double wq[], double xq[], double yq[]);

    //
    //  Purpose:
    //
    //    QUAD_E sets a quadrature rule for the error calculation.
    //
    //  Parameters:
    //
    //    Input, double NODE_XY[2*NODE_NUM], the X and Y coordinates of nodes.
    //
    //    Input, int ELEMENT_NODE[NNODES*ELEMENT_NUM]; ELEMENT_NODE(I,J) is the global
    //    index of local node I in element J.
    //
    //    Input, int ELEMENT, the index of the element for which the quadrature
    //    points are to be computed.
    //
    //    Input, int ELEMENT_NUM, the number of elements.
    //
    //    Input, int NNODES, the number of nodes used to form one element.
    //
    //    Input, int NODE_NUM, the number of nodes.
    //
    //    Input, int NQE, the number of points in the quadrature rule.
    //    This is actually fixed at 13.
    //
    //    Output, double WQE[NQE], the quadrature weights.
    //
    //    Output, double XQE[NQE], YQE[NQE], the X and Y coordinates
    //    of the quadrature points.
    //
void quad_e(double node_xy[], int element_node[],
	int element, int element_num, int nnodes, int node_num,
	int nqe, double wqe[], double xqe[], double yqe[]);

    //
    //  Purpose:
    //
    //    R8_HUGE returns a "huge" R8.
    //
    //  Parameters:
    //
    //    Output, double R8_HUGE, a "huge" value.
    //
double r8_huge(void);

    //
    //  Purpose:
    //
    //    R8_MAX returns the maximum of two R8's.
    //
    //  Parameters:
    //
    //    Input, double X, Y, the quantities to compare.
    //
    //    Output, double R8_MAX, the maximum of X and Y.
    //
double r8_max(double x, double y);

    //
    //  Purpose:
    //
    //    R8_MIN returns the minimum of two R8's.
    //
    //  Parameters:
    //
    //    Input double X, Y, the quantities to compare.
    //
    //    Output, double R8_MIN, the minimum of X and Y.
    //
double r8_min(double x, double y);

    //
    //  Purpose:
    //
    //    R8_NINT returns the nearest integer to an R8.
    //
    //  Example:
    //
    //        X        R8_NINT
    //
    //      1.3         1
    //      1.4         1
    //      1.5         1 or 2
    //      1.6         2
    //      0.0         0
    //     -0.7        -1
    //     -1.1        -1
    //     -1.6        -2
    //
    //  Parameters:
    //
    //    Input, double X, the value.
    //
    //    Output, int R8_NINT, the nearest integer to X.
    //
int r8_nint(double x);

    //
    //  Purpose:
    //
    //    R8VEC_PRINT_SOME prints "some" of an R8VEC.
    //
    //  Discussion:
    //
    //    An R8VEC is a vector of R8 values.
    //
    //  Parameters:
    //
    //    Input, int N, the number of entries of the vector.
    //
    //    Input, double A[N], the vector to be printed.
    //
    //    Input, integer I_LO, I_HI, the first and last indices to print.
    //    The routine expects 1 <= I_LO <= I_HI <= N.
    //
    //    Input, char *TITLE, an optional title.
    //
void r8vec_print_some(int n, double a[], int i_lo, int i_hi, char *title);

    //
    //  Purpose:
    //
    //    RHS gives the right-hand side of the differential equation.
    //
    //  Discussion:
    //
    //    The function specified here depends on the problem being
    //    solved.  This is one of the routines that a user will
    //    normally want to change.
    //
    //  Parameters:
    //
    //    Input, double X, Y, the coordinates of a point
    //    in the region, at which the right hand side of the
    //    differential equation is to be evaluated.
    //
    //    Input, double TIME, the current time.
    //
    //    Output, double RHS, the value of the right
    //    hand side of the differential equation at (X,Y).
    //
double rhs(double x, double y, double time);

    //
    //  Purpose:
    //
    //    S_LEN_TRIM returns the length of a string to the last nonblank.
    //
    //  Parameters:
    //
    //    Input, char *S, a pointer to a string.
    //
    //    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
    //    If S_LEN_TRIM is 0, then the string is entirely blank.
    //
int s_len_trim(char *s);

    //
    //  Purpose:
    //
    //    SOLUTION_WRITE writes the solution to a file.
    //
    //  Parameters:
    //
    //    Input, int NODE_NUM, the number of nodes.
    //
    //    Input, double U[NODE_NUM], the coefficients of the solution.
    //
    //    Input, char *U_FILE_NAME, the name of the file
    //    in which the data should be stored.
    //
void solution_write(int node_num, double u[], char *u_file_name);

    //
    //  Purpose:
    //
    //    TIMESTAMP prints the current YMDHMS date as a time stamp.
    //
    //  Example:
    //
    //    May 31 2001 09:45:54 AM
    //
    //  Parameters:
    //
    //    None
    //
void timestamp(void);

    //
    //  Purpose:
    //
    //    TRIANGULATION_ORDER6_PLOT plots a 6-node triangulation of a pointset.
    //
    //  Discussion:
    //
    //    The triangulation is most usually a Delaunay triangulation,
    //    but this is not necessary.
    //
    //    This routine has been specialized to deal correctly ONLY with
    //    a mesh of 6 node elements, with the property that starting
    //    from local node 1 and traversing the edges of the element will
    //    result in encountering local nodes 1, 4, 2, 5, 3, 6 in that order.
    //
    //  Parameters:
    //
    //    Input, char *FILE_NAME, the name of the file to create.
    //
    //    Input, int NODE_NUM, the number of nodes.
    //
    //    Input, double precision NODE_XY[2*NODE_NUM], the nodes.
    //
    //    Input, int TRI_NUM, the number of triangles.
    //
    //    Input, int TRIANGLE_NODE[6*TRI_NUM], lists, for each triangle,
    //    the indices of the points that form the vertices and midsides
    //    of the triangle.
    //
    //    Input, int NODE_SHOW:
    //    0, do not show nodes;
    //    1, show nodes;
    //    2, show nodes and label them.
    //
    //    Input, int TRIANGLE_SHOW:
    //    0, do not show triangles;
    //    1, show triangles;
    //    2, show triangles and label them.
    //
void triangulation_order6_plot(char *file_name, int node_num,
	double node_xy[], int tri_num, int triangle_node[], int node_show,
	int triangle_show);

    //
    //  Purpose:
    //
    //    XY_SET sets the XY coordinates of the nodes.
    //
    //  Discussion:
    //
    //    The nodes are laid out in an evenly spaced grid, in the unit square.
    //
    //    The first node is at the origin.  More nodes are created to the
    //    right until the value of X = 1 is reached, at which point
    //    the next layer is generated starting back at X = 0, and an
    //    increased value of Y.
    //
    //  Parameters:
    //
    //    Input, int NX, NY, the number of elements in the X and
    //    Y direction.
    //
    //    Input, int NODE_NUM, the number of nodes.
    //
    //    Input, double XL, XR, YB, YT, the X coordinates of
    //    the left and right sides of the rectangle, and the Y coordinates
    //    of the bottom and top of the rectangle.
    //
    //    Output, double NODE_XY[2*NODE_NUM], the nodes.
    //
void xy_set(int nx, int ny, int node_num, double xl, double xr, double yb,
	double yt, double node_xy[]);
