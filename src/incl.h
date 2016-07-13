#pragma once

int main(void);

void adjust_backward_euler(int node_num, double node_xy[], int nnodes, int element_num, int element_node[], int quad_num, double wq[], double xq[], double yq[], double element_area[], int ib, double time,double time_step_size, double u_old[], double a[], double f[]);

void adjust_boundary(int node_num, double node_xy[], int node_boundary[], int ib, double time, double a[], double f[]);

void area_set(int node_num, double node_xy[], int nnodes, int element_num, int element_node[], double element_area[]);

void assemble(int node_num, double node_xy[], int nnodes, int element_num, int element_node[], int quad_num, double wq[], double xq[], double yq[], double element_area[], int ib, double time, double a[], double f[]);

int bandwidth(int nnodes, int element_num, int element_node[], int node_num);

void compare(int node_num, double node_xy[], double time, double u[]);

int dgb_fa(int n, int ml, int mu, double a[], int pivot[]);

void dgb_print_some(int m, int n, int ml, int mu, double a[], int ilo, int jlo, int ihi, int jhi, char *title);

double *dgb_sl(int n, int ml, int mu, double a[], int pivot[], double b[], int job);

void element_write(int nnodes, int element_num, int element_node[],char *triangulation_txt_file_name);

void errors(double element_area[], int element_node[], double node_xy[], double u[], int element_num, int nnodes, int node_num, double time, double *el2, double *eh1);

void exact_u(int node_num, double node_xy[], double time, double u_exact[], double dudx_exact[], double dudy_exact[]);

void file_name_inc(char *file_name);

void grid_t6(int nx, int ny, int nnodes, int element_num, int element_node[]);

int i4_max(int i1, int i2);

int i4_min (int i1, int i2);

void i4vec_print_some(int n, int a[], int max_print, char *title);

int *node_boundary_set(int nx, int ny, int node_num);

void nodes_plot(char *file_name, int node_num, double node_xy[], bool node_label);

void nodes_write(int node_num, double node_xy[], char *output_filename);

void qbf(double x, double y, int element, int inode, double node_xy[], int element_node[], int element_num, int nnodes, int node_num, double *bb, double *bx, double *by);

void quad_a(double node_xy[], int element_node[], int element_num, int node_num, int nnodes, double wq[], double xq[], double yq[]);

void quad_e(double node_xy[], int element_node[], int element, int element_num, int nnodes, int node_num, int nqe, double wqe[], double xqe[], double yqe[]);

double r8_huge(void);

double r8_max(double x, double y);

double r8_min(double x, double y);

int r8_nint(double x);

void r8vec_print_some(int n, double a[], int i_lo, int i_hi, char *title);

double rhs(double x, double y, double time);

int s_len_trim(char *s);

void solution_write(int node_num, double u[], char *u_file_name);

void timestamp(void);

void triangulation_order6_plot(char *file_name, int node_num, double node_xy[], int tri_num, int triangle_node[], int node_show, int triangle_show);

void xy_set(int nx, int ny, int node_num, double xl, double xr, double yb, double yt, double node_xy[]);
