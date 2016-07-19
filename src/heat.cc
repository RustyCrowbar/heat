#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "incl.hh"

# define NNODES 6
# define QUAD_NUM 3
# define NX 25
# define NY 25

# define ELEMENT_NUM ( NX - 1 ) * ( NY - 1 ) * 2
# define NODE_NUM ( 2 * NX - 1 ) * ( 2 * NY - 1 )

#define TIME_END	2.0
#define NB_ITERATIONS	80

void adjust_backward_euler(/*int node_num, */double node_xy[], int nnodes, int element_num, int element_node[], int quad_num, double wq[], double xq[], double yq[], double element_area[], int ib, /*double time, */double time_step_size, double u_old[], double a[], double f[])
{
    int basis;
    double bi;
    double bj;
    double dbidx;
    double dbidy;
    double dbjdx;
    double dbjdy;
    int element;
    //int i;
    //int ip;
    //int ipp;
    int j;
    int node;
    int quad;
    int test;
    double w;
    double x;
    double y;

    for(element = 0; element < element_num; element++)
    {
        for(quad = 0; quad < quad_num; quad++)
        {
            x = xq[quad+element*quad_num];
            y = yq[quad+element*quad_num];
            w = element_area[element] * wq[quad];

            for (test = 0; test < nnodes; test++)
            {
                node = element_node[test+element*nnodes];
                qbf(x, y, element, test, node_xy, element_node,	/*element_num, */nnodes, /*node_num, */&bi, &dbidx, &dbidy);
                //
                //  Carry the U_OLD term to the right hand side.
                //
                f[node] += w * bi * u_old[node] / time_step_size;
                //
                //  Modify the diagonal entries of A.
                //
                for (basis = 0; basis < nnodes; basis++)
                {
                    j = element_node[basis+element*nnodes];
                    qbf(x, y, element, basis, node_xy, element_node,/* element_num, */nnodes, /*node_num, */&bj, &dbjdx, &dbjdy);
                    a[node-j+2*ib+j*(3*ib+1)] = a[node-j+2*ib+j*(3*ib+1)] + w * bi * bj / time_step_size;
                }
            }
        }
    }
    return;
}

void adjust_boundary(int node_num, double node_xy[], int node_boundary[], int ib, double time, double a[], double f[])
{
    double *dudx_exact;
    double *dudy_exact;
    //int i;
    int j;
    int jhi;
    int jlo;
    int node;
    double *u_exact;
    //
    //  Get the exact solution at every node.
    //
    u_exact = new double[node_num];
    dudx_exact = new double[node_num];
    dudy_exact = new double[node_num];
    exact_u(node_num, node_xy, time, u_exact, dudx_exact, dudy_exact);

    for (node = 0; node < node_num; node++)
    {
        if (node_boundary[node] != 0)
        {
            jlo = std::max(node - ib, 0);
            jhi = std::min(node + ib, node_num - 1);
            for (j = jlo; j <= jhi; j++)
                a[node-j+2*ib+j*(3*ib+1)] = 0.0;
            a[node-node+2*ib+node*(3*ib+1)] = 1.0;
            f[node] = u_exact[node];
        }
    }

    delete [] u_exact;
    delete [] dudx_exact;
    delete [] dudy_exact;
    return;
}

void area_set(/*int node_num, */double node_xy[], int nnodes, int element_num, int element_node[], double element_area[])
{
    int element;
    int i1;
    int i2;
    int i3;
    double x1;
    double x2;
    double x3;
    double y1;
    double y2;
    double y3;

    for (element = 0; element < element_num; element++)
    {
        i1 = element_node[0+element*nnodes];
        x1 = node_xy[0+i1*2];
        y1 = node_xy[1+i1*2];

        i2 = element_node[1+element*nnodes];
        x2 = node_xy[0+i2*2];
        y2 = node_xy[1+i2*2];

        i3 = element_node[2+element*nnodes];
        x3 = node_xy[0+i3*2];
        y3 = node_xy[1+i3*2];

        element_area[element] = 0.5 * fabs(y1 * (x2 - x3) + y2 * (x3 - x1) + y3 * (x1 - x2));
    }
    return;
}

void assemble(int node_num, double node_xy[], int nnodes, int element_num, int element_node[], int quad_num, double wq[], double xq[], double yq[], double element_area[], int ib, double time, double a[], double f[])
{
    double aij;
    int basis;
    double bi;
    double bj;
    double dbidx;
    double dbidy;
    double dbjdx;
    double dbjdy;
    int element;
    int i;
    int j;
    int node;
    int quad;
    int test;
    double w;
    double x;
    double y;
    //
    //  Initialize the arrays to zero.
    //
    for (i = 0; i < node_num; i++)
        f[i] = 0.0;

    for (j = 0; j < node_num; j++)
        for (i = 0; i < 3*ib + 1; i++)
            a[i+j*(3*ib+1)] = 0.0;
    //
    //  The actual values of A and F are determined by summing up
    //  contributions from all the elements.
    //
    for (element = 0; element < element_num; element++)
    {
        for (quad = 0; quad < quad_num; quad++)
        {
            x = xq[quad+element*quad_num];
            y = yq[quad+element*quad_num];
            w = element_area[element] * wq[quad];

            for (test = 0; test < nnodes; test++)
            {
                node = element_node[test+element*nnodes];
                qbf(x, y, element, test, node_xy, element_node,	/*element_num, */nnodes, /*node_num, */&bi, &dbidx, &dbidy);
                f[node] = f[node] + w * rhs(x, y, time) * bi;
                //
                //  We are about to compute a contribution associated with the
                //  I-th test function and the J-th basis function, and add this
                //  to the entry A(I,J).
                //
                //  Because of the compressed storage of the matrix, the element
                //  will actually be stored in A(I-J+2*IB+1,J).
                //
                //  An extra complication: we are storing the array as a vector.
                //
                //  Therefore, we ACTUALLY store the entry in A[I-J+2*IB+1-1 + J * (3*IB+1)];
                //
                for(basis = 0; basis < nnodes; basis++)
                {
                    j = element_node[basis + element * nnodes];
                    qbf(x, y, element, basis, node_xy, element_node, /*element_num, */nnodes, /*node_num, */&bj, &dbjdx, &dbjdy);
                    aij = dbidx * dbjdx + dbidy * dbjdy;
                    a[node-j+2*ib+j*(3*ib+1)] += w * aij;
                }
            }
        }
    }
    return;
}

int bandwidth(int nnodes, int element_num, int element_node[]/*, int node_num*/)
{
    int element;
    int i;
    int iln;
    //int in;
    int j;
    int jln;
    //int jn;
    int nhba;

    nhba = 0;

    for (element = 0; element < element_num; element++)
        for (iln = 0; iln < nnodes; iln++)
        {
            i = element_node[iln+element*nnodes];
            for (jln = 0; jln < nnodes; jln++)
            {
                j = element_node[jln+element*nnodes];
                nhba = std::max(nhba, j - i);
            }
        }
    return nhba;
}

void compare(int node_num, double node_xy[], double time, double u[])
{
    double *dudx_exact;
    double *dudy_exact;
    int node;
    double *u_exact;

    u_exact = new double[node_num];
    dudx_exact = new double[node_num];
    dudy_exact = new double[node_num];

    exact_u (node_num, node_xy, time, u_exact, dudx_exact, dudy_exact);

    std::cout << std::endl << "COMPARE:" << std::endl
	<< "  Compare computed and exact solutions at the nodes." << std::endl << std::endl
	<< "         X           Y          U           U" << std::endl
	<< "                              exact       computed" << std::endl;

    for (node = 0; node < node_num; node++)
        std::cout << std::setw(12) << node_xy[0+node*2] << "  " << std::setw(12) << node_xy[1+node*2] << "  " << std::setw(12) << u_exact[node] << "  " << std::setw(12) << u[node] << "\n";

    delete [] u_exact;
    delete [] dudx_exact;
    delete [] dudy_exact;
    return;
}

int dgb_fa(int n, int ml, int mu, double a[], int pivot[])
{
    int col = 2 * ml + mu + 1;
    int i;
    int i0;
    int j;
    int j0;
    int j1;
    int ju;
    int jz;
    int k;
    int l;
    int lm;
    int m;
    int mm;
    double t;

    m = ml + mu + 1;
    //
    //  Zero out the initial fill-in columns.
    //
    j0 = mu + 2;
    j1 = std::min(n, m) - 1;

    for (jz = j0; jz <= j1; jz++)
    {
        i0 = m + 1 - jz;
        for (i = i0; i <= ml; i++)
            a[i-1+(jz-1)*col] = 0.0;
    }

    jz = j1;
    ju = 0;

    for (k = 1; k <= n-1; k++)
    {
        //
        //  Zero out the next fill-in column.
        //
        jz = jz + 1;
        if (jz <= n)
            for (i = 1; i <= ml; i++)
                a[i-1+(jz-1)*col] = 0.0;
        //
        //  Find L = pivot index.
        //
        lm = std::min(ml, n-k);
        l = m;

        for (j = m+1; j <= m + lm; j++)
            if (fabs(a[l-1+(k-1)*col]) < fabs(a[j-1+(k-1)*col]))
                l = j;

        pivot[k-1] = l + k - m;
        //
        //  Zero pivot implies this column already triangularized.
        //
        if (a[l-1+(k-1)*col] == 0.0)
        {
            std::cout << "\n" << "DGB_FA - Fatal error!\n  Zero pivot on step " << k << "\n";
            return k;
        }
        //
        //  Interchange if necessary.
        //
        t = a[l-1+(k-1)*col];
        a[l-1+(k-1)*col] = a[m-1+(k-1)*col];
        a[m-1+(k-1)*col] = t;
        //
        //  Compute multipliers.
        //
        for (i = m+1; i <= m+lm; i++)
            a[i-1+(k-1)*col] = - a[i-1+(k-1)*col] / a[m-1+(k-1)*col];
        //
        //  Row elimination with column indexing.
        //
        ju = std::max(ju, mu + pivot[k-1]);
        ju = std::min(ju, n);
        mm = m;

        for (j = k+1; j <= ju; j++)
        {
            l = l - 1;
            mm = mm - 1;
            if (l != mm)
            {
                t = a[l-1+(j-1)*col];
                a[l-1+(j-1)*col] = a[mm-1+(j-1)*col];
                a[mm-1+(j-1)*col] = t;
            }
            for (i = 1; i <= lm; i++)
                a[mm+i-1+(j-1)*col] = a[mm+i-1+(j-1)*col] + a[mm-1+(j-1)*col] * a[m+i-1+(k-1)*col];
        }
    }

    pivot[n-1] = n;
    if (a[m-1+(n-1)*col] == 0.0)
    {
        std::cout << "\nDGB_FA - Fatal error!\n Zero pivot on step " << n << "\n";
        return n;
    }
    return 0;
}
void dgb_print_some(int m, int n, int ml, int mu, double a[], int ilo, int jlo, int ihi, int jhi, std::string title)
{
# define INCX 5

    int col = 2 * ml + mu + 1;
    int i;
    int i2hi;
    int i2lo;
    int j;
    int j2hi;
    int j2lo;

    if (0 < s_len_trim(title))
        std::cout << "\n" << title << "\n";
    //
    //  Print the columns of the matrix, in strips of 5.
    //
    for (j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX)
    {
        j2hi = j2lo + INCX - 1;
        j2hi = std::min(j2hi, n);
        j2hi = std::min(j2hi, jhi);
        std::cout << std::endl << "  Col: ";
        for (j = j2lo; j <= j2hi; j++)
            std::cout << std::setw(7) << j << "       ";
        std::cout << std::endl
		<< "  Row" << std::endl << "  ---" << std::endl;
        //
        //  Determine the range of the rows in this strip.
        //
        i2lo = std::max({ilo, 1, j2lo - mu});
        i2hi = std::min({ihi, m, j2hi + ml});

        for (i = i2lo; i <= i2hi; i++)
        {
            //
            //  Print out (up to) 5 entries in row I, that lie in the current strip.
            //
            std::cout << std::setw(6) << i << "  ";
            for ( j = j2lo; j <= j2hi; j++ )
            {
                if ( ml < i-j || mu < j-i )
                    std::cout << "            ";
                else
                    std::cout << std::setw(10) << a[i-j+ml+mu+(j-1)*col] << "  ";
            }
            std::cout << std::endl;
        }
    }
    std::cout << std::endl;
    return;
# undef INCX
}

double *dgb_sl(int n, int ml, int mu, double a[], int pivot[], double b[], int job)
{
    int col = 2 * ml + mu + 1;
    int i;
    int k;
    int l;
    int la;
    int lb;
    int lm;
    int m;
    double t;
    double *x;

    x = new double[n];

    for ( i = 0; i < n; i++ )
        x[i] = b[i];
    //
    m = mu + ml + 1;
    //
    //  Solve A * x = b.
    //
    if (job == 0)
    {
        //
        //  Solve L * Y = B.
        //
        if (1 <= ml)
        {
            for (k = 1; k <= n-1; k++)
            {
                lm = std::min(ml, n-k);
                l = pivot[k-1];

                if (l != k)
                {
                    t      = x[l-1];
                    x[l-1] = x[k-1];
                    x[k-1] = t;
                }
                for (i = 1; i <= lm; i++)
                    x[k+i-1] = x[k+i-1] + x[k-1] * a[m+i-1+(k-1)*col];
            }
        }
        //
        //  Solve U * X = Y.
        //
        for (k = n; 1 <= k; k--)
        {
            x[k-1] = x[k-1] / a[m-1+(k-1)*col];
            lm = std::min ( k, m ) - 1;
            la = m - lm;
            lb = k - lm;
            for ( i = 0; i <= lm-1; i++ )
                x[lb+i-1] = x[lb+i-1] - x[k-1] * a[la+i-1+(k-1)*col];
        }
    }
    //
    //  Solve A' * X = B.
    //
    else
    {
        //
        //  Solve U' * Y = B.
        //
        for ( k = 1; k <= n; k++ )
        {
            lm = std::min(k, m) - 1;
            la = m - lm;
            lb = k - lm;
            for (i = 0; i <= lm-1; i++)
                x[k-1] = x[k-1] - x[lb+i-1] * a[la+i-1+(k-1)*col];
            x[k-1] = x[k-1] / a[m-1+(k-1)*col];
        }
        //
        //  Solve L' * X = Y.
        //
        if (1 <= ml)
        {
            for (k = n-1; 1 <= k; k--)
            {
                lm = std::min ( ml, n-k );
                for (i = 1; i <= lm; i++)
                    x[k-1] = x[k-1] + x[k+i-1] * a[m+i-1+(k-1)*col];
                l = pivot[k-1];

                if (l != k)
                {
                    t      = x[l-1];
                    x[l-1] = x[k-1];
                    x[k-1] = t;
                }
            }
        }
    }
    return x;
}

void element_write(int nnodes, int element_num, int element_node[], std::string output_filename)
{
    int element;
    int i;
    std::ofstream output;

    output.open(output_filename.c_str());

    if (!output)
    {
        std::cout << std::endl;
        std::cout << "ELEMENT_WRITE - Warning!" << std::endl;
        std::cout << "  Could not write the node file." << std::endl;
        return;
    }

    for(element = 0; element < element_num; element++)
    {
        for (i = 0; i < nnodes; i++)
            output << std::setw(8) << element_node[i+element*nnodes] << "  ";
        output << std::endl;
    }
    output.close();
    return;
}

void errors(double element_area[], int element_node[], double node_xy[], double u[], int element_num, int nnodes, /*int node_num, */double time, double *el2, double *eh1)
{
# define NQE 13

    double ar;
    double bi;
    double dbidx;
    double dbidy;
    double dudx_exact[1];
    double dudxh;
    double dudy_exact[1];
    double dudyh;
    int element;
    int i;
    int in1;
    //int ip;
    int quad;
    double u_exact[1];
    double uh;
    double wqe[NQE];
    double x;
    //double x1;
    double xqe[NQE];
    double xy[2];
    double y;
    //double y1;
    double yqe[NQE];

    *el2 = 0.0;
    *eh1 = 0.0;
    //
    //  For each element, retrieve the nodes, area, quadrature weights,
    //  and quadrature points.
    //
    for(element = 0; element < element_num; element++)
    {
        quad_e(node_xy, element_node, element, /*element_num, */nnodes, /*node_num, NQE, */wqe, xqe, yqe);
        //
        //  For each quadrature point, evaluate the computed solution and its X and
        //  Y derivatives.
        //
        for(quad = 0; quad < NQE; quad++)
        {
            ar = element_area[element] * wqe[quad];
            x = xqe[quad];
            y = yqe[quad];

            uh = 0.0;
            dudxh = 0.0;
            dudyh = 0.0;

            for(in1 = 0; in1 < nnodes; in1++)
            {
                i = element_node[in1+element*nnodes];

                qbf(x, y, element, in1, node_xy, element_node, /*element_num, */nnodes, /*node_num, */&bi, &dbidx, &dbidy);

                uh    = uh    + bi    * u[i];
                dudxh = dudxh + dbidx * u[i];
                dudyh = dudyh + dbidy * u[i];
            }
            //
            //  Evaluate the exact solution and its X and Y derivatives.
            //
            xy[0] = x;
            xy[1] = y;

            exact_u(1, xy, time, u_exact, dudx_exact, dudy_exact);
            //
            //  Add the weighted value at this quadrature point to the quadrature sum.
            //
            *el2 = *el2 + ar * pow (( uh - u_exact[0]), 2 );

            *eh1 = *eh1 + ar * (pow(dudxh - dudx_exact[0], 2) + pow(dudyh - dudy_exact[0], 2));
        }
    }

    *el2 = sqrt(*el2);
    *eh1 = sqrt(*eh1);

    std::cout << std::setw(14) << time
        << std::setw(14) << *el2
        << std::setw(14) << *eh1 << std::endl;

    return;
# undef NQE
}

void exact_u(int node_num, double node_xy[], double time, double u[], double dudx[], double dudy[])
{
# define PI 3.141592653589793

    int node;
    double x;
    double y;

    for(node = 0; node < node_num; node++)
    {
        x = node_xy[0+node*2];
        y = node_xy[1+node*2];

        u[node]    =      sin ( PI * x ) * sin ( PI * y ) * exp ( - time );
        dudx[node] = PI * cos ( PI * x ) * sin ( PI * y ) * exp ( - time );
        dudy[node] = PI * sin ( PI * x ) * cos ( PI * y ) * exp ( - time );
    }

    return;
# undef PI
}

void file_name_inc(std::string file_name)
{
    char c;
    int change;
    int i;
    int lens;

    lens = s_len_trim ( file_name );

    if ( lens <= 0 )
    {
        std::cout << "\n";
        std::cout << "FILE_NAME_INC - Fatal error!\n";
        std::cout << "  Input file name is blank.\n";
        exit ( 1 );
    }

    change = 0;

    for ( i = lens-1; 0 <= i; i-- )
    {
        c = file_name[i];

        if ( '0' <= c && c <= '9' )
        {
            change = change + 1;
            if ( c == '9' )
            {
                c = '0';
                file_name[i] = c;
            }
            else
            {
                c = c + 1;
                file_name[i] = c;
                return;
            }
        }
    }

    if (!change)
        file_name = "";

    return;
}

void grid_t6(int nx, int ny, int nnodes, /*int element_num, */int element_node[])
{
    int c;
    int e;
    int element;
    int i;
    int j;
    int n;
    int ne;
    int nw;
    int s;
    int se;
    int sw;
    int w;

    element = 0;

    for ( j = 1; j <= ny - 1; j++ )
    {
        for ( i = 1; i <= nx - 1; i++ )
        {
            sw = ( j - 1 ) * 2 * ( 2 * nx - 1 ) + 2 * i - 2;
            w  = sw + 1;
            nw = sw + 2;

            s  = sw + 2 * nx - 1;
            c  = s + 1;
            n  = s + 2;

            se = s  + 2 * nx - 1;
            e  = se + 1;
            ne = se + 2;

            element_node[0+element*nnodes] = sw;
            element_node[1+element*nnodes] = se;
            element_node[2+element*nnodes] = nw;
            element_node[3+element*nnodes] = s;
            element_node[4+element*nnodes] = c;
            element_node[5+element*nnodes] = w;
            element = element + 1;

            element_node[0+element*nnodes] = ne;
            element_node[1+element*nnodes] = nw;
            element_node[2+element*nnodes] = se;
            element_node[3+element*nnodes] = n;
            element_node[4+element*nnodes] = c;
            element_node[5+element*nnodes] = e;
            element = element + 1;
        }
    }

    return;
}

void i4vec_print_some(int n, int a[], int max_print, std::string title)
{
    int i;

    if ( max_print <= 0 )
        return;

    if ( n <= 0 )
        return;

    if ( 0 < s_len_trim ( title ) )
    {
        std::cout << "\n";
        std::cout << title << "\n";
        std::cout << "\n";
    }

    if ( n <= max_print )
    {
        for ( i = 0; i < n; i++ )
        {
            std::cout << std::setw(6)  << i + 1 << "  "
                << std::setw(10) << a[i] << "\n";
        }
    }
    else if ( 3 <= max_print )
    {
        for ( i = 0; i < max_print-2; i++ )
        {
            std::cout << std::setw(6)  << i + 1 << "  "
                << std::setw(10) << a[i]  << "\n";
        }
        std::cout << "......  ..............\n";
        i = n - 1;
        std::cout << std::setw(6)  << i + 1 << "  "
            << std::setw(10) << a[i]  << "\n";
    }
    else
    {
        for ( i = 0; i < max_print-1; i++ )
        {
            std::cout << std::setw(6)  << i + 1 << "  "
                << std::setw(10) << a[i]  << "\n";
        }
        i = max_print - 1;
        std::cout << std::setw(6)  << i + 1 << "  "
            << std::setw(10) << a[i]  << "...more entries...\n";
    }

    return;
}

int *node_boundary_set(int nx, int ny, int node_num)
{
    int i;
    int j;
    int node;
    int *node_boundary;

    node_boundary = new int[node_num];

    node = 0;

    for ( j = 1; j <= 2 * ny - 1; j++ )
    {
        for ( i = 1; i <= 2 * nx - 1; i++ )
        {
            if (j == 1 ||
                j == 2 * ny - 1 ||
                i == 1 ||
                i == 2 * nx - 1 )
            {
                node_boundary[node] = 1;
            }
            else
            {
                node_boundary[node] = 0;
            }

            ++node;
        }
    }
    return node_boundary;
}

void nodes_plot(std::string file_name, int node_num, double node_xy[]/*, bool node_label*/)
{
    int circle_size;
    int delta;
    std::ofstream file_unit;
    //int i;
    int node;
    double x_max;
    double x_min;
    int x_ps;
    int x_ps_max = 576;
    int x_ps_max_clip = 594;
    int x_ps_min = 36;
    int x_ps_min_clip = 18;
    double x_scale;
    double y_max;
    double y_min;
    int y_ps;
    int y_ps_max = 666;
    int y_ps_max_clip = 684;
    int y_ps_min = 126;
    int y_ps_min_clip = 108;
    double y_scale;
    //
    //  We need to do some figuring here, so that we can determine
    //  the range of the data, and hence the height and width
    //  of the piece of paper.
    //
    x_max = -r8_huge ( );
    for ( node = 0; node < node_num; node++ )
    {
        if ( x_max < node_xy[0+node*2] )
        {
            x_max = node_xy[0+node*2];
        }
    }
    x_min = r8_huge ( );
    for ( node = 0; node < node_num; node++ )
    {
        if ( node_xy[0+node*2] < x_min )
        {
            x_min = node_xy[0+node*2];
        }
    }
    x_scale = x_max - x_min;

    x_max = x_max + 0.05 * x_scale;
    x_min = x_min - 0.05 * x_scale;
    x_scale = x_max - x_min;

    y_max = -r8_huge ( );
    for ( node = 0; node < node_num; node++ )
    {
        if ( y_max < node_xy[1+node*2] )
        {
            y_max = node_xy[1+node*2];
        }
    }
    y_min = r8_huge ( );
    for ( node = 0; node < node_num; node++ )
    {
        if ( node_xy[1+node*2] < y_min )
        {
            y_min = node_xy[1+node*2];
        }
    }
    y_scale = y_max - y_min;

    y_max = y_max + 0.05 * y_scale;
    y_min = y_min - 0.05 * y_scale;
    y_scale = y_max - y_min;

    if ( x_scale < y_scale )
    {
        delta = (int)std::round ( ( double ) ( x_ps_max - x_ps_min )
                * ( y_scale - x_scale ) / ( 2.0 * y_scale ) );

        x_ps_max = x_ps_max - delta;
        x_ps_min = x_ps_min + delta;

        x_ps_max_clip = x_ps_max_clip - delta;
        x_ps_min_clip = x_ps_min_clip + delta;

        x_scale = y_scale;
    }
    else if ( y_scale < x_scale )
    {
        delta = (int)std::round ( ( double ) ( y_ps_max - y_ps_min )
                * ( x_scale - y_scale ) / ( 2.0 * x_scale ) );

        y_ps_max = y_ps_max - delta;
        y_ps_min = y_ps_min + delta;

        y_ps_max_clip = y_ps_max_clip - delta;
        y_ps_min_clip = y_ps_min_clip + delta;

        y_scale = x_scale;
    }

    file_unit.open(file_name.c_str());

    if ( !file_unit )
    {
        std::cout << "\n";
        std::cout << "POINTS_PLOT - Fatal error!\n";
        std::cout << "  Could not open the output EPS file.\n";
        exit ( 1 );
    }

    file_unit << "%!PS-Adobe-3.0 EPSF-3.0\n";
    file_unit << "%%Creator: nodes_plot.C\n";
    file_unit << "%%Title: " << file_name << "\n";
    file_unit << "%%Pages: 1\n";
    file_unit << "%%BoundingBox:  "
        << x_ps_min << "  "
        << y_ps_min << "  "
        << x_ps_max << "  "
        << y_ps_max << "\n";
    file_unit << "%%Document-Fonts: Times-Roman\n";
    file_unit << "%%LanguageLevel: 1\n";
    file_unit << "%%EndComments\n";
    file_unit << "%%BeginProlog\n";
    file_unit << "/inch {72 mul} def\n";
    file_unit << "%%EndProlog\n";
    file_unit << "%%Page:      1     1\n";
    file_unit << "save\n";
    file_unit << "%\n";
    file_unit << "% Set the RGB line color to very light gray.\n";
    file_unit << "%\n";
    file_unit << " 0.9000 0.9000 0.9000 setrgbcolor\n";
    file_unit << "%\n";
    file_unit << "% Draw a gray border around the page.\n";
    file_unit << "%\n";
    file_unit << "newpath\n";
    file_unit << x_ps_min << "  "
        << y_ps_min << "  moveto\n";
    file_unit << x_ps_max << "  "
        << y_ps_min << "  lineto\n";
    file_unit << x_ps_max << "  "
        << y_ps_max << "  lineto\n";
    file_unit << x_ps_min << "  "
        << y_ps_max << "  lineto\n";
    file_unit << x_ps_min << "  "
        << y_ps_min << "  lineto\n";
    file_unit << "stroke\n";
    file_unit << "%\n";
    file_unit << "% Set RGB line color to black.\n";
    file_unit << "%\n";
    file_unit << " 0.0000 0.0000 0.0000 setrgbcolor\n";
    file_unit << "%\n";
    file_unit << "%  Set the font and its size:\n";
    file_unit << "%\n";
    file_unit << "/Times-Roman findfont\n";
    file_unit << "0.50 inch scalefont\n";
    file_unit << "setfont\n";
    file_unit << "%\n";
    file_unit << "%  Print a title:\n";
    file_unit << "%\n";
    file_unit << "%  210  702 moveto\n";
    file_unit << "%(Pointset) show\n";
    file_unit << "%\n";
    file_unit << "% Define a clipping polygon\n";
    file_unit << "%\n";
    file_unit << "newpath\n";
    file_unit << x_ps_min_clip << "  "
        << y_ps_min_clip << "  moveto\n";
    file_unit << x_ps_max_clip << "  "
        << y_ps_min_clip << "  lineto\n";
    file_unit << x_ps_max_clip << "  "
        << y_ps_max_clip << "  lineto\n";
    file_unit << x_ps_min_clip << "  "
        << y_ps_max_clip << "  lineto\n";
    file_unit << x_ps_min_clip << "  "
        << y_ps_min_clip << "  lineto\n";
    file_unit << "clip newpath\n";
    //
    //  Draw the nodes.
    //
    if ( node_num <= 200 )
    {
        circle_size = 5;
    }
    else if ( node_num <= 500 )
    {
        circle_size = 4;
    }
    else if ( node_num <= 1000 )
    {
        circle_size = 3;
    }
    else if ( node_num <= 5000 )
    {
        circle_size = 2;
    }
    else
    {
        circle_size = 1;
    }

    file_unit << "%\n";
    file_unit << "%  Draw filled dots at each node:\n";
    file_unit << "%\n";
    file_unit << "%  Set the color to blue:\n";
    file_unit << "%\n";
    file_unit << "0.000  0.150  0.750  setrgbcolor\n";
    file_unit << "%\n";

    for ( node = 0; node < node_num; node++ )
    {
        x_ps = ( int ) (
                ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )
                  + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) )
                / ( x_max                     - x_min ) );

        y_ps = ( int ) (
                ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )
                  + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) )
                / ( y_max                     - y_min ) );

        file_unit << "newpath  "
            << x_ps << "  "
            << y_ps << "  "
            << circle_size << " 0 360 arc closepath fill\n";
    }
    //
    //  Label the nodes.
    //
    file_unit << "%\n";
    file_unit << "%  Label the nodes:\n";
    file_unit << "%\n";
    file_unit << "%  Set the color to darker blue:\n";
    file_unit << "%\n";
    file_unit << "0.000  0.250  0.850  setrgbcolor\n";
    file_unit << "/Times-Roman findfont\n";
    file_unit << "0.20 inch scalefont\n";
    file_unit << "setfont\n";

    file_unit << "%\n";

    for ( node = 0; node < node_num; node++ )
    {
        x_ps = ( int ) (
                ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )
                  + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) )
                / ( x_max                     - x_min ) );

        y_ps = ( int ) (
                ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )
                  + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) )
                / ( y_max                     - y_min ) );

        file_unit << "newpath  "
            << x_ps     << "  "
            << y_ps + 5 << "  moveto ("
            << node     << ") show\n";
    }

    file_unit << "%\n";
    file_unit << "restore showpage\n";
    file_unit << "%\n";
    file_unit << "% End of page\n";
    file_unit << "%\n";
    file_unit << "%%Trailer\n";
    file_unit << "%%EOF\n";

    file_unit.close ( );

    return;
}

void nodes_write(int node_num, double node_xy[], std::string output_filename)
{
    int node;
    std::ofstream output;
    double x;
    double y;

    output.open(output_filename.c_str());

    if (!output)
    {
        std::cout << "\n";
        std::cout << "NODES_WRITE - Warning!\n";
        std::cout << "  Could not write the node file.\n";
        return;
    }

    for ( node = 0; node < node_num; node++ )
    {
        x = node_xy[0+node*2];
        y = node_xy[1+node*2];

        output << std::setw(8)  << x << "  "
            << std::setw(8)  << y << "\n";
    }

    output.close ( );

    return;
}

void qbf(double x, double y, int element, int inode, double node_xy[], int element_node[], /*int element_num, */int nnodes, /*int node_num, */double *b, double *dbdx, double *dbdy)
{
    double dbdr;
    double dbds;
    double det;
    double drdx;
    double drdy;
    double dsdx;
    double dsdy;
    int i;
    double r;
    double s;
    double xn[6];
    double yn[6];

    for ( i = 0; i < 6; i++ )
    {
        xn[i] = node_xy[0+element_node[i+element*nnodes]*2];
        yn[i] = node_xy[1+element_node[i+element*nnodes]*2];
    }
    //
    //  Determine the (R,S) coordinates corresponding to (X,Y).
    //
    //  What is happening here is that we are solving the linear system:
    //
    //    ( X2-X1  X3-X1 ) * ( R ) = ( X - X1 )
    //    ( Y2-Y1  Y3-Y1 )   ( S )   ( Y - Y1 )
    //
    //  by computing the inverse of the coefficient matrix and multiplying
    //  it by the right hand side to get R and S.
    //
    //  The values of dRdX, dRdY, dSdX and dSdY are easily from the formulas
    //  for R and S.
    //
    det =   ( xn[1] - xn[0] ) * ( yn[2] - yn[0] )
        - ( xn[2] - xn[0] ) * ( yn[1] - yn[0] );

    r = ( ( yn[2] - yn[0] ) * ( x     - xn[0] )
            + ( xn[0] - xn[2] ) * ( y     - yn[0] ) ) / det;

    drdx = ( yn[2] - yn[0] ) / det;
    drdy = ( xn[0] - xn[2] ) / det;

    s = ( ( yn[0] - yn[1] ) * ( x     - xn[0] )
            + ( xn[1] - xn[0] ) * ( y     - yn[0] ) ) / det;

    dsdx = ( yn[0] - yn[1] ) / det;
    dsdy = ( xn[1] - xn[0] ) / det;
    //
    //  The basis functions can now be evaluated in terms of the
    //  reference coordinates R and S.  It's also easy to determine
    //  the values of the derivatives with respect to R and S.
    //
    if ( inode == 0 )
    {
        *b   =   2.0 *     ( 1.0 - r - s ) * ( 0.5 - r - s );
        dbdr = - 3.0 + 4.0 * r + 4.0 * s;
        dbds = - 3.0 + 4.0 * r + 4.0 * s;
    }
    else if ( inode == 1 )
    {
        *b   =   2.0 * r * ( r - 0.5 );
        dbdr = - 1.0 + 4.0 * r;
        dbds =   0.0;
    }
    else if ( inode == 2 )
    {
        *b   =   2.0 * s * ( s - 0.5 );
        dbdr =   0.0;
        dbds = - 1.0               + 4.0 * s;
    }
    else if ( inode == 3 )
    {
        *b   =   4.0 * r * ( 1.0 - r - s );
        dbdr =   4.0 - 8.0 * r - 4.0 * s;
        dbds =           - 4.0 * r;
    }
    else if ( inode == 4 )
    {
        *b   =   4.0 * r * s;
        dbdr =                           4.0 * s;
        dbds =             4.0 * r;
    }
    else if ( inode == 5 )
    {
        *b   =   4.0 * s * ( 1.0 - r - s );
        dbdr =                         - 4.0 * s;
        dbds =   4.0 - 4.0 * r - 8.0 * s;
    }
    else
    {
        std::cout << "\n";
        std::cout << "QBF - Fatal error!\n";
        std::cout << "  Request for local basis function INODE = " << inode << "\n";
        exit ( 1 );
    }
    //
    //  We need to convert the derivative information from (R(X,Y),S(X,Y))
    //  to (X,Y) using the chain rule.
    //
    *dbdx = dbdr * drdx + dbds * dsdx;
    *dbdy = dbdr * drdy + dbds * dsdy;

    return;
}

void quad_a(double node_xy[], int element_node[], int element_num, /*int node_num, */int nnodes, double wq[], double xq[], double yq[])
{
    int element;
    int ip1;
    int ip2;
    int ip3;
    double x1;
    double x2;
    double x3;
    double y1;
    double y2;
    double y3;

    wq[0] = 1.0 / 3.0;
    wq[1] = wq[0];
    wq[2] = wq[0];

    for ( element = 0; element < element_num; element++ )
    {
        ip1 = element_node[0+element*nnodes];
        ip2 = element_node[1+element*nnodes];
        ip3 = element_node[2+element*nnodes];

        x1 = node_xy[0+ip1*2];
        x2 = node_xy[0+ip2*2];
        x3 = node_xy[0+ip3*2];

        y1 = node_xy[1+ip1*2];
        y2 = node_xy[1+ip2*2];
        y3 = node_xy[1+ip3*2];

        xq[0+element*3] = 0.5 * ( x1 + x2 );
        xq[1+element*3] = 0.5 * ( x2 + x3 );
        xq[2+element*3] = 0.5 * ( x1 + x3 );

        yq[0+element*3] = 0.5 * ( y1 + y2 );
        yq[1+element*3] = 0.5 * ( y2 + y3 );
        yq[2+element*3] = 0.5 * ( y1 + y3 );
    }
    return;
}

void quad_e(double node_xy[], int element_node[], int element, /*int element_num, */int nnodes, /*int node_num, int nqe, */double wqe[], double xqe[], double yqe[])
{
    int i;
    int ii;
    int iii;
    int ip1;
    int ip2;
    int ip3;
    double x1;
    double x2;
    double x3;
    double y1;
    double y2;
    double y3;
    double z1;
    double z2;
    double z3;
    double z4;
    double z5;
    double z6;
    double z7;

    for ( i = 0; i < 3; i++ )
    {
        wqe[i] = 0.175615257433204;
        ii = i + 3;
        wqe[ii] = 0.053347235608839;
        ii = i + 6;
        iii = ii + 3;
        wqe[ii] = 0.077113760890257;
        wqe[iii] = wqe[ii];
    }

    wqe[12] = -0.14957004446767;

    z1 = 0.479308067841923;
    z2 = 0.260345966079038;
    z3 = 0.869739794195568;
    z4 = 0.065130102902216;
    z5 = 0.638444188569809;
    z6 = 0.312865496004875;
    z7 = 0.048690315425316;

    ip1 = element_node[0+element*nnodes];
    ip2 = element_node[1+element*nnodes];
    ip3 = element_node[2+element*nnodes];

    x1 = node_xy[0+ip1*2];
    x2 = node_xy[0+ip2*2];
    x3 = node_xy[0+ip3*2];

    y1 = node_xy[1+ip1*2];
    y2 = node_xy[1+ip2*2];
    y3 = node_xy[1+ip3*2];

    xqe[ 0] = z1 * x1 + z2 * x2 + z2 * x3;
    yqe[ 0] = z1 * y1 + z2 * y2 + z2 * y3;
    xqe[ 1] = z2 * x1 + z1 * x2 + z2 * x3;
    yqe[ 1] = z2 * y1 + z1 * y2 + z2 * y3;
    xqe[ 2] = z2 * x1 + z2 * x2 + z1 * x3;
    yqe[ 2] = z2 * y1 + z2 * y2 + z1 * y3;
    xqe[ 3] = z3 * x1 + z4 * x2 + z4 * x3;
    yqe[ 3] = z3 * y1 + z4 * y2 + z4 * y3;
    xqe[ 4] = z4 * x1 + z3 * x2 + z4 * x3;
    yqe[ 4] = z4 * y1 + z3 * y2 + z4 * y3;
    xqe[ 5] = z4 * x1 + z4 * x2 + z3 * x3;
    yqe[ 5] = z4 * y1 + z4 * y2 + z3 * y3;
    xqe[ 6] = z5 * x1 + z6 * x2 + z7 * x3;
    yqe[ 6] = z5 * y1 + z6 * y2 + z7 * y3;
    xqe[ 7] = z5 * x1 + z7 * x2 + z6 * x3;
    yqe[ 7] = z5 * y1 + z7 * y2 + z6 * y3;
    xqe[ 8] = z6 * x1 + z5 * x2 + z7 * x3;
    yqe[ 8] = z6 * y1 + z5 * y2 + z7 * y3;
    xqe[ 9] = z6 * x1 + z7 * x2 + z5 * x3;
    yqe[ 9] = z6 * y1 + z7 * y2 + z5 * y3;
    xqe[10] = z7 * x1 + z5 * x2 + z6 * x3;
    yqe[10] = z7 * y1 + z5 * y2 + z6 * y3;
    xqe[11] = z7 * x1 + z6 * x2 + z5 * x3;
    yqe[11] = z7 * y1 + z6 * y2 + z5 * y3;
    xqe[12] = ( x1 + x2 + x3 ) / 3.0;
    yqe[12] = ( y1 + y2 + y3 ) / 3.0;

    return;
}

double r8_huge()
{
    return (double)HUGE_VAL;
}

void r8vec_print_some(int n, double a[], int i_lo, int i_hi, std::string title)
{
    int i;

    if ( 0 < s_len_trim ( title ) )
    {
        std::cout << "\n";
        std::cout << title << "\n";
    }

    std::cout << "\n";
    for ( i = std::max ( 1, i_lo ); i <= std::min ( n, i_hi ); i++ )
    {
        std::cout << "  " << std::setw(8)  << i       << "  "
            << "  " << std::setw(14) << a[i-1]  << "\n";
    }

    return;
}

double rhs(double x, double y, double time)
{
# define PI 3.141592653589793

    double ut;
    double uxx;
    double uyy;
    double value;

    ut =            - sin ( PI * x ) * sin ( PI * y ) * exp ( - time );
    uxx = - PI * PI * sin ( PI * x ) * sin ( PI * y ) * exp ( - time );
    uyy = - PI * PI * sin ( PI * x ) * sin ( PI * y ) * exp ( - time );

    value = ut - uxx - uyy;

    return value;
# undef PI
}

int s_len_trim(std::string s)
{
    return s.find(' ');
    /*int n;
    char* t;

    n = s.size();
    t = s + strlen (s) - 1;

    while (0 < n)
    {
        if (*t != ' ')
            return n;
        t--;
        n--;
    }

    return n;*/
}

void solution_write(int node_num, double u[], std::string u_file_name)
{
    int node;
    std::ofstream u_file;

    u_file.open(u_file_name.c_str());

    if ( !u_file )
    {
        std::cout << "\n";
        std::cout << "SOLUTION_WRITE - Warning!\n";
        std::cout << "  Could not write the solution file \""
            << u_file_name << "\".\n";
        return;
    }

    const int nb_nodes_x = NX * 2 - 1;
    for ( node = 0; node < node_num; node++ )
    {
	if ((node % nb_nodes_x) == 0 && (node > 0))
		u_file << std::endl << u[node] << " ";
	else
		u_file << u[node] << " ";
    }

    u_file.close ( );

    return;
}

void timestamp()
{
# define TIME_SIZE 40

    char time_buffer[TIME_SIZE];
    struct tm *tm;
    size_t len;
    time_t now;

    now = time(NULL); //Original code
    //time(NULL); //We have a segfault, here.
    tm = localtime(&now);
    len = strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);
    std::cout << time_buffer << "\n";
    //std::cout << "DELETE_ME" << endl;
    return;
# undef TIME_SIZE
}

void triangulation_order6_plot(std::string file_name, int node_num, double node_xy[], int tri_num, int triangle_node[], int node_show, int triangle_show)
{
    double ave_x;
    double ave_y;
    int circle_size;
    int delta;
    //int e;
    std::ofstream file_unit;
    int i;
    int ip1;
    int node;
    int order[6] = { 0, 3, 1, 4, 2, 5 };
    int triangle;
    double x_max;
    double x_min;
    int x_ps;
    int x_ps_max = 576;
    int x_ps_max_clip = 594;
    int x_ps_min = 36;
    int x_ps_min_clip = 18;
    double x_scale;
    double y_max;
    double y_min;
    int y_ps;
    int y_ps_max = 666;
    int y_ps_max_clip = 684;
    int y_ps_min = 126;
    int y_ps_min_clip = 108;
    double y_scale;
    //
    //  We need to do some figuring here, so that we can determine
    //  the range of the data, and hence the height and width
    //  of the piece of paper.
    //
    x_max = -r8_huge ( );
    for ( node = 0; node < node_num; node++ )
        if ( x_max < node_xy[0+node*2] )
            x_max = node_xy[0+node*2];
    x_min = r8_huge ( );
    for ( node = 0; node < node_num; node++ )
        if ( node_xy[0+node*2] < x_min )
            x_min = node_xy[0+node*2];
    x_scale = x_max - x_min;

    x_max = x_max + 0.05 * x_scale;
    x_min = x_min - 0.05 * x_scale;
    x_scale = x_max - x_min;

    y_max = -r8_huge ( );
    for ( node = 0; node < node_num; node++ )
        if ( y_max < node_xy[1+node*2] )
            y_max = node_xy[1+node*2];
    y_min = r8_huge ( );
    for ( node = 0; node < node_num; node++ )
        if ( node_xy[1+node*2] < y_min )
            y_min = node_xy[1+node*2];
    y_scale = y_max - y_min;

    y_max = y_max + 0.05 * y_scale;
    y_min = y_min - 0.05 * y_scale;
    y_scale = y_max - y_min;

    if ( x_scale < y_scale )
    {
        delta = (int)std::round ( ( double ) ( x_ps_max - x_ps_min )
                * ( y_scale - x_scale ) / ( 2.0 * y_scale ) );

        x_ps_max = x_ps_max - delta;
        x_ps_min = x_ps_min + delta;

        x_ps_max_clip = x_ps_max_clip - delta;
        x_ps_min_clip = x_ps_min_clip + delta;

        x_scale = y_scale;
    }
    else if ( y_scale < x_scale )
    {
        delta = (int)std::round ( ( double ) ( y_ps_max - y_ps_min )
                * ( x_scale - y_scale ) / ( 2.0 * x_scale ) );

        y_ps_max = y_ps_max - delta;
        y_ps_min = y_ps_min + delta;

        y_ps_max_clip = y_ps_max_clip - delta;
        y_ps_min_clip = y_ps_min_clip + delta;

        y_scale = x_scale;
    }

    file_unit.open(file_name.c_str());

    if ( !file_unit )
    {
        std::cout << "\n";
        std::cout << "TRIANGULATION_ORDER6_PLOT - Fatal error!\n";
        std::cout << "  Could not open the output EPS file.\n";
        exit ( 1 );
    }

    file_unit << "%!PS-Adobe-3.0 EPSF-3.0\n";
    file_unit << "%%Creator: triangulation_order6_plot.C\n";
    file_unit << "%%Title: " << file_name << "\n";
    file_unit << "%%Pages: 1\n";
    file_unit << "%%BoundingBox:  "
        << x_ps_min << "  "
        << y_ps_min << "  "
        << x_ps_max << "  "
        << y_ps_max << "\n";
    file_unit << "%%Document-Fonts: Times-Roman\n";
    file_unit << "%%LanguageLevel: 1\n";
    file_unit << "%%EndComments\n";
    file_unit << "%%BeginProlog\n";
    file_unit << "/inch {72 mul} def\n";
    file_unit << "%%EndProlog\n";
    file_unit << "%%Page:      1     1\n";
    file_unit << "save\n";
    file_unit << "%\n";
    file_unit << "% Set the RGB line color to very light gray.\n";
    file_unit << "%\n";
    file_unit << " 0.9000 0.9000 0.9000 setrgbcolor\n";
    file_unit << "%\n";
    file_unit << "% Draw a gray border around the page.\n";
    file_unit << "%\n";
    file_unit << "newpath\n";
    file_unit << x_ps_min << "  "
        << y_ps_min << "  moveto\n";
    file_unit << x_ps_max << "  "
        << y_ps_min << "  lineto\n";
    file_unit << x_ps_max << "  "
        << y_ps_max << "  lineto\n";
    file_unit << x_ps_min << "  "
        << y_ps_max << "  lineto\n";
    file_unit << x_ps_min << "  "
        << y_ps_min << "  lineto\n";
    file_unit << "stroke\n";
    file_unit << "%\n";
    file_unit << "% Set RGB line color to black.\n";
    file_unit << "%\n";
    file_unit << " 0.0000 0.0000 0.0000 setrgbcolor\n";
    file_unit << "%\n";
    file_unit << "%  Set the font and its size:\n";
    file_unit << "%\n";
    file_unit << "/Times-Roman findfont\n";
    file_unit << "0.50 inch scalefont\n";
    file_unit << "setfont\n";
    file_unit << "%\n";
    file_unit << "%  Print a title:\n";
    file_unit << "%\n";
    file_unit << "%  210  702 moveto\n";
    file_unit << "%(Pointset) show\n";
    file_unit << "%\n";
    file_unit << "% Define a clipping polygon\n";
    file_unit << "%\n";
    file_unit << "newpath\n";
    file_unit << x_ps_min_clip << "  "
        << y_ps_min_clip << "  moveto\n";
    file_unit << x_ps_max_clip << "  "
        << y_ps_min_clip << "  lineto\n";
    file_unit << x_ps_max_clip << "  "
        << y_ps_max_clip << "  lineto\n";
    file_unit << x_ps_min_clip << "  "
        << y_ps_max_clip << "  lineto\n";
    file_unit << x_ps_min_clip << "  "
        << y_ps_min_clip << "  lineto\n";
    file_unit << "clip newpath\n";
    //
    //  Draw the nodes.
    //
    if ( node_num <= 200 )
        circle_size = 5;
    else if ( node_num <= 500 )
        circle_size = 4;
    else if ( node_num <= 1000 )
        circle_size = 3;
    else if ( node_num <= 5000 )
        circle_size = 2;
    else
        circle_size = 1;
    if ( 1 <= node_show )
    {
        file_unit << "%\n";
        file_unit << "%  Draw filled dots at each node:\n";
        file_unit << "%\n";
        file_unit << "%  Set the color to blue:\n";
        file_unit << "%\n";
        file_unit << "0.000  0.150  0.750  setrgbcolor\n";
        file_unit << "%\n";

        for ( node = 0; node < node_num; node++ )
        {
            x_ps = ( int ) (
                    ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )
                      + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) )
                    / ( x_max                     - x_min ) );

            y_ps = ( int ) (
                    ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )
                      + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) )
                    / ( y_max                     - y_min ) );

            file_unit << "newpath  "
                << x_ps << "  "
                << y_ps << "  "
                << circle_size << " 0 360 arc closepath fill\n";
        }
    }
    //
    //  Label the nodes.
    //
    if ( 2 <= node_show )
    {
        file_unit << "%\n";
        file_unit << "%  Label the nodes:\n";
        file_unit << "%\n";
        file_unit << "%  Set the color to darker blue:\n";
        file_unit << "%\n";
        file_unit << "0.000  0.250  0.850  setrgbcolor\n";
        file_unit << "/Times-Roman findfont\n";
        file_unit << "0.20 inch scalefont\n";
        file_unit << "setfont\n";

        file_unit << "%\n";

        for ( node = 0; node < node_num; node++ )
        {
            x_ps = ( int ) (
                    ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )
                      + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) )
                    / ( x_max                     - x_min ) );

            y_ps = ( int ) (
                    ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )
                      + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) )
                    / ( y_max                     - y_min ) );

            file_unit << "newpath  "
                << x_ps     << "  "
                << y_ps + 5 << "  moveto ("
                << node     << ") show\n";
        }
    }
    //
    //  Draw the triangles.
    //
    if ( 1 <= triangle_show )
    {
        file_unit << "%\n";
        file_unit << "%  Set the RGB color to red.\n";
        file_unit << "%\n";
        file_unit << "0.900  0.200  0.100 setrgbcolor\n";
        file_unit << "%\n";
        file_unit << "%  Draw the triangles.\n";
        file_unit << "%\n";

        for ( triangle = 0; triangle < tri_num; triangle++ )
        {
            node = triangle_node[order[0]+triangle*6];

            x_ps = ( int ) (
                    ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )
                      + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) )
                    / ( x_max                     - x_min ) );

            y_ps = ( int ) (
                    ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )
                      + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) )
                    / ( y_max                     - y_min ) );

            file_unit << "newpath  " << x_ps << "  " << y_ps << "  moveto\n";

            for ( i = 1; i <= 6; i++ )
            {
                ip1 = ( i % 6 );
                node = triangle_node[order[ip1]+triangle*6];

                x_ps = ( int ) (
                        ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )
                          + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) )
                        / ( x_max                     - x_min ) );

                y_ps = ( int ) (
                        ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )
                          + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) )
                        / ( y_max                     - y_min ) );

                file_unit << x_ps << "  " << y_ps << "  lineto\n";
            }
            file_unit << "stroke\n";
        }
    }
    //
    //  Label the triangles.
    //
    if ( 2 <= triangle_show )
    {
        file_unit << "%\n";
        file_unit << "%  Label the triangles.\n";
        file_unit << "%\n";
        file_unit << "%  Set the RGB color to darker red.\n";
        file_unit << "%\n";
        file_unit << "0.950  0.250  0.150 setrgbcolor\n";
        file_unit << "/Times-Roman findfont\n";
        file_unit << "0.20 inch scalefont\n";
        file_unit << "setfont\n";
        file_unit << "%\n";

        for ( triangle = 0; triangle < tri_num; triangle++ )
        {
            ave_x = 0.0;
            ave_y = 0.0;

            for ( i = 0; i < 6; i++ )
            {
                node = triangle_node[i+triangle*6];
                ave_x = ave_x + node_xy[0+node*2];
                ave_y = ave_y + node_xy[1+node*2];
            }

            ave_x = ave_x / 6.0;
            ave_y = ave_y / 6.0;

            x_ps = ( int ) (
                    ( ( x_max - ave_x         ) * ( double ) ( x_ps_min )
                      + (       + ave_x - x_min ) * ( double ) ( x_ps_max ) )
                    / ( x_max         - x_min ) );

            y_ps = ( int ) (
                    ( ( y_max - ave_y         ) * ( double ) ( y_ps_min )
                      + (         ave_y - y_min ) * ( double ) ( y_ps_max ) )
                    / ( y_max         - y_min ) );

            file_unit << std::setw(4) << x_ps << "  "
                << std::setw(4) << y_ps << "  "
                << "moveto (" << triangle << ") show\n";
        }
    }

    file_unit << "%\n";
    file_unit << "restore showpage\n";
    file_unit << "%\n";
    file_unit << "% End of page\n";
    file_unit << "%\n";
    file_unit << "%%Trailer\n";
    file_unit << "%%EOF\n";

    file_unit.close();

    return;
}

void xy_set(int nx, int ny, /*int node_num, */double xl, double xr, double yb, double yt, double node_xy[])
{
    int i;
    int j;

    for(j = 0; j < 2*ny - 1; j++)
        for(i = 0; i < 2*nx - 1; i++)
        {
            node_xy[0+(i+j*(2*nx-1))*2] =
                (double(2 * nx - i - 2) * xl
                  + double (i) * xr)
                / double(2 * nx - 2);

            node_xy[1+(i+j*(2*nx-1))*2] =
                (double(2 * ny - j - 2) * yb
                  + double(j) * yt)
                / double(2 * ny - 2);
        }
    return;
}
