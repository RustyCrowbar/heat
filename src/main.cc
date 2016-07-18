#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "incl.hh"

//  Purpose:
//
//    MAIN is the main program for FEM2D_HEAT_RECTANGLE.
//
//  Discussion:
//
//    FEM2D_HEAT_RECTANGLE solves
//
//      dUdT - Laplacian U(X,Y,T) = F(X,Y,T)
//
//    in a rectangular region in the plane.
//
//    Along the boundary of the region, Dirichlet conditions
//    are imposed:
//
//      U(X,Y,T) = G(X,Y,T)
//
//    At the initial time T_INIT, the value of U is given at all points
//    in the region:
//
//      U(X,Y,T) = H(X,Y,T)
//
//    The code uses continuous piecewise quadratic basis functions on
//    triangles determined by a uniform grid of NX by NY points.
//
//    The backward Euler approximation is used for the time derivatives.
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Local parameters:
//
//    Local, double A[(3*IB+1)*NODE_NUM], the coefficient matrix.
//
//    Local, double EH1, the H1 seminorm error.
//
//    Local, double EL2, the L2 error.
//
//    Local, double ELEMENT_AREA[ELEMENT_NUM], the area of each element.
//
//    Local, int ELEMENT_NODE[ELEMENT_NUM*NNODES]; ELEMENT_NODE(I,J) is the
//    global node index of the local node J in element I.
//
//    Local, int ELEMENT_NUM, the number of elements.
//
//    Local, double F[NODE_NUM], the right hand side.
//
//    Local, int IB, the half-bandwidth of the matrix.
//
//    Local, integer NODE_BOUNDARY[NODE_NUM], is
//    0, if a node is an interior node;
//    1, if a node is a Dirichlet boundary node.
//
//    Local, int NNODES, the number of nodes used to form one element.
//
//    Local, double NODE_XY[2*NODE_NUM], the X and Y coordinates of nodes.
//
//    Local, int NX, the number of points in the X direction.
//
//    Local, int NY, the number of points in the Y direction.
//
//    Local, int QUAD_NUM, the number of quadrature points used for assembly.
//
//    Local, double TIME, the current time.
//
//    Local, double TIME_FINAL, the final time.
//
//    Local, double TIME_INIT, the initial time.
//
//    Local, double TIME_OLD, the time at the previous time step.
//
//    Local, int TIME_STEP_NUM, the number of time steps to take.
//
//    Local, double TIME_STEP_SIZE, the size of the time steps.
//
//    Local, double U[NODE_NUM], the finite element coefficients
//    defining the solution at the current time.
//
//    Local, double WQ[QUAD_NUM], quadrature weights.
//
//    Local, double XL, XR, YB, YT, the X coordinates of
//    the left and right sides of the rectangle, and the Y coordinates
//    of the bottom and top of the rectangle.
//
//    Local, double XQ[QUAD_NUM*ELEMENT_NUM], YQ[QUAD_NUM*ELEMENT_NUM], the X and Y
//    coordinates of the quadrature points in each element.
//

int main(void)
{
# define NNODES 6
# define QUAD_NUM 3
# define NX 5
# define NY 5

# define ELEMENT_NUM ( NX - 1 ) * ( NY - 1 ) * 2
# define NODE_NUM ( 2 * NX - 1 ) * ( 2 * NY - 1 )

    double *a;
    //double dt;
    double *dudx_exact;
    double *dudy_exact;
    double eh1;
    double el2;
    //int element;
    double element_area[ELEMENT_NUM];
    //bool *element_mask;
    int element_node[NNODES*ELEMENT_NUM];
    double *f;
    //int i;
    int ib;
    int ierr;
    int job;
    //int local;
    int node;
    int *node_boundary;
    std::string node_eps_file_name = "rectangle_nodes.eps";
    std::string node_txt_file_name = "rectangle_nodes.txt";
    bool node_label;
    int node_show;
    double node_xy[2*NODE_NUM];
    int *pivot;
    double time;
    std::string time_file_name = "rectangle_time.txt";
    double time_final;
    double time_init;
    double time_old;
    int time_step;
    int time_step_num;
    double time_step_size;
    std::ofstream time_unit;
    int triangle_show;
    std::string triangulation_eps_file_name = "rectangle_elements.eps";
    std::string triangulation_txt_file_name = "rectangle_elements.txt";
    double *u;
    double *u_exact;
    std::string u_file_name = "rectangle_u0000.txt";
    double *u_old;
    double wq[QUAD_NUM];
    double xl = 0.0;
    double xq[QUAD_NUM*ELEMENT_NUM];
    double xr = 1.0;
    double yb = 0.0;
    double yq[QUAD_NUM*ELEMENT_NUM];
    double yt = 1.0;

    timestamp();

    std::cout << "\n";
    std::cout << "FEM2D_HEAT_RECTANGLE\n";
    std::cout << "  C++ version\n";
    std::cout << "\n";
    std::cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
    std::cout << "\n";
    std::cout << "  Solution of the time-dependent heat equation\n";
    std::cout << "  on a unit box in 2 dimensions.\n";
    std::cout << "\n";
    std::cout << "  Ut - Uxx - Uyy = F(x,y,t) in the box\n";
    std::cout << "        U(x,y,t) = G(x,y,t) for (x,y) on the boundary.\n";
    std::cout << "        U(x,y,t) = H(x,y,t) for t = T_INIT.\n";
    std::cout << "\n";
    std::cout << "  The finite element method is used, with piecewise\n";
    std::cout << "  quadratic basis functions on 6 node triangular\n";
    std::cout << "  elements.\n";
    std::cout << "\n";
    std::cout << "  The corner nodes of the triangles are generated by an\n";
    std::cout << "  underlying grid whose dimensions are\n";
    std::cout << "\n";
    std::cout << "  NX =                 " << NX << "\n";
    std::cout << "  NY =                 " << NY << "\n";
    std::cout << "\n";
    std::cout << "  Number of nodes    = " << NODE_NUM << "\n";
    std::cout << "  Number of elements = " << ELEMENT_NUM << "\n";
    //
    //  Set the coordinates of the nodes.
    //
    xy_set(NX, NY, /*NODE_NUM, */xl, xr, yb, yt, node_xy);
    //
    //  Organize the nodes into a grid of 6-node triangles.
    //
    grid_t6(NX, NY, NNODES, /*ELEMENT_NUM, */element_node);
    //
    //  Set the quadrature rule for assembly.
    //
    quad_a(node_xy, element_node, ELEMENT_NUM, /*NODE_NUM, */NNODES, wq, xq, yq);
    //
    //  Determine the areas of the elements.
    //
    area_set(/*NODE_NUM, */node_xy, NNODES, ELEMENT_NUM, element_node, element_area);
    //
    //  Determine which nodes are boundary nodes and which have a
    //  finite element unknown.  Then set the boundary values.
    //
    node_boundary = node_boundary_set(NX, NY, NODE_NUM);

#ifndef NDEBUG
    i4vec_print_some(NODE_NUM, node_boundary, 10, "  NODE_BOUNDARY:");
#endif
    //
    //  Determine the bandwidth of the coefficient matrix.
    //
    ib = bandwidth(NNODES, ELEMENT_NUM, element_node/*, NODE_NUM*/);

    std::cout << "\n";
    std::cout << "  The matrix half bandwidth is " << ib << "\n";
    std::cout << "  The matrix row size is       " << 3 * ib + 1 << "\n";
    //
    //  Make an EPS picture of the nodes.
    //
    if (NX <= 10 && NY <= 10)
    {
        node_label = true;
        nodes_plot(node_eps_file_name, NODE_NUM, node_xy/*, node_label*/);

        std::cout << "\n";
        std::cout << "FEM2D_HEAT_RECTANGLE:\n";
        std::cout << "  Wrote an EPS file\n";
        std::cout << "    \"" << node_eps_file_name << "\".\n";
        std::cout << "  containing a picture of the nodes.\n";
    }
    //
    //  Write the nodes to an ASCII file that can be read into MATLAB.
    //
    nodes_write (NODE_NUM, node_xy, node_txt_file_name);

    std::cout << "\n";
    std::cout << "FEM2D_HEAT_RECTANGLE:\n";
    std::cout << "  Wrote an ASCII node file\n";
    std::cout << "    " << node_txt_file_name << "\n";
    std::cout << "  of the form\n";
    std::cout << "    X(I), Y(I)\n";
    std::cout << "  which can be used for plotting.\n";
    //
    //  Make a picture of the elements.
    //
    if (NX <= 10 && NY <= 10)
    {
        node_show = 2;
        triangle_show = 2;

        triangulation_order6_plot(triangulation_eps_file_name, NODE_NUM, node_xy, ELEMENT_NUM, element_node, node_show, triangle_show);
        std::cout << "\n";
        std::cout << "FEM2D_HEAT_RECTANGLE:\n";
        std::cout << "  Wrote an EPS file\n";
        std::cout << "    \"" << triangulation_eps_file_name << "\".\n";
        std::cout << "  containing a picture of the elements.\n";
    }
    //
    //  Write the elements to a file that can be read into MATLAB.
    //
    element_write(NNODES, ELEMENT_NUM, element_node, triangulation_txt_file_name);

    std::cout << "\n";
    std::cout << "FEM2D_HEAT_RECTANGLE:\n";
    std::cout << "  Wrote an ASCII element file\n";
    std::cout << "    \"" << triangulation_txt_file_name << "\".\n";
    std::cout << "  of the form\n";
    std::cout << "    Node(1) Node(2) Node(3) Node(4) Node(5) Node(6)\n";
    std::cout << "  which can be used for plotting.\n";
    //
    //  Set time stepping quantities.
    //
    time_init = 0.0;
    time_final = 0.5;
    time_step_num = 10;
    time_step_size = (time_final - time_init) / (double)(time_step_num);
    //
    //  Allocate space.
    //
    a = new double[(3*ib+1)*NODE_NUM];
    dudx_exact = new double[NODE_NUM];
    dudy_exact = new double[NODE_NUM];
    f = new double[NODE_NUM];
    pivot = new int[NODE_NUM];
    u = new double[NODE_NUM];
    u_exact = new double[NODE_NUM];
    u_old = new double[NODE_NUM];
    //
    //  Set the value of U at the initial time.
    //
    time = time_init;
    exact_u(NODE_NUM, node_xy, time, u_exact, dudx_exact, dudy_exact);

    for(node = 0; node < NODE_NUM; node++)
        u[node] = u_exact[node];

    time_unit.open(time_file_name.c_str());

    if (!time_unit)
    {
        std::cout << "\n";
        std::cout << "FEM2D_HEAT_RECTANGLE- Warning!\n";
        std::cout << "  Could not write the time file \"" << time_file_name << "\".\n";
        exit(1);
    }

    time_unit << "  " << std::setw(14) << time << "\n";

    solution_write(NODE_NUM, u, u_file_name);
    //
    //  Time looping.
    //
    std::cout << "\n";
    std::cout << "     Time        L2 Error       H1 Error\n";
    std::cout << "\n";

    for(time_step = 1; time_step <= time_step_num; time_step++ )
    {
        time_old = time;
        for ( node = 0; node < NODE_NUM; node++ )
            u_old[node] = u[node];
        delete [] u;

        time = ((double)(time_step_num - time_step) * time_init	+ (double)(time_step) * time_final) / (double)(time_step_num);
        //
        //  Assemble the coefficient matrix A and the right-hand side F of the
        //  finite element equations.
        //
        assemble(NODE_NUM, node_xy, NNODES, ELEMENT_NUM, element_node, QUAD_NUM, wq, xq, yq, element_area, ib, time, a, f);

#ifndef NDEBUG
        dgb_print_some(NODE_NUM, NODE_NUM, ib, ib, a, 10, 1, 12, 25, "  Initial block of coefficient matrix A:");
        r8vec_print_some(NODE_NUM, f, 1, 10, "  Part of the right hand side F:");
#endif
        //
        //  Modify the coefficient matrix and right hand side to account for the dU/dt
        //  term, which we are treating using the backward Euler formula.
        //
        adjust_backward_euler(/*NODE_NUM, */node_xy, NNODES, ELEMENT_NUM, element_node, QUAD_NUM, wq, xq, yq, element_area, ib, /*time, */time_step_size, u_old, a, f);

#ifndef NDEBUG
        dgb_print_some(NODE_NUM, NODE_NUM, ib, ib, a, 10, 1, 12, 25, "  A after DT adjustment:");
        r8vec_print_some ( NODE_NUM, f, 1, 10, "  F after DT adjustment:");
#endif
        //
        //  Modify the coefficient matrix and right hand side to account for
        //  boundary conditions.
        //
        adjust_boundary(NODE_NUM, node_xy, node_boundary, ib, time, a, f);

#ifndef NDEBUG
        dgb_print_some(NODE_NUM, NODE_NUM, ib, ib, a, 10, 1, 12, 25, "  A after BC adjustment:");
        r8vec_print_some(NODE_NUM, f, 1, 10, "  F after BC adjustment:");
#endif
        //
        //  Solve the linear system using a banded solver.
        //
        ierr = dgb_fa(NODE_NUM, ib, ib, a, pivot);

        if (ierr != 0)
        {
            std::cout << "\n";
            std::cout << "FEM2D_HEAT_RECTANGLE - Error!\n";
            std::cout << "  DGB_FA returned an error condition.\n";
            std::cout << "\n";
            std::cout << "  The linear system was not factored, and the\n";
            std::cout << "  algorithm cannot proceed.\n";
            exit (1);
        }

        job = 0;
        u = dgb_sl (NODE_NUM, ib, ib, a, pivot, f, job);

#ifndef NDEBUG
        r8vec_print_some(NODE_NUM, u, 1, 10, "  Part of the solution vector:");
#endif
        //
        //  Calculate error using 13 point quadrature rule.
        //
        errors(element_area, element_node, node_xy, u, ELEMENT_NUM, NNODES, /*NODE_NUM, */time, &el2, &eh1);
        //
        //  Compare the exact and computed solutions just at the nodes.
        //
#ifndef NDEBUG
        compare(NODE_NUM, node_xy, time, u);
#endif
        //
        //  Increment the file name, and write the new solution.
        //
        time_unit << std::setw(14) << time << "\n";
        file_name_inc(u_file_name);
        solution_write(NODE_NUM, u, u_file_name);
    }
    //
    //  Deallocate memory.
    //
    delete [] a;
    delete [] dudx_exact;
    delete [] dudy_exact;
    delete [] f;
    delete [] node_boundary;
    delete [] pivot;
    delete [] u;
    delete [] u_exact;
    delete [] u_old;

    time_unit.close();
    //
    //  Terminate.
    //
    std::cout << "\n";
    std::cout << "FEM2D_HEAT_RECTANGLE:\n";
    std::cout << "  Normal end of execution.\n";

    std::cout << "\n";
    timestamp();

    return 0;
# undef ELEMENT_NUM
# undef NNODES
# undef NODE_NUM
# undef NX
# undef NY
# undef QUAD_NUM
}
