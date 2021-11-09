#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
//#include <cs.h>

#define RECTANGLE 3
#define TRAPEZOID_ROW 2
#define TRAPEZOID_COLUMN 1

const int NET_COLUMN = 6;       //net_line 的列数

const double PERCENTAGE_DIV_0_7 = 0.7;
const double PERCENTAGE_INCREASING = 1.1;
double dielectric_constant = 3.9;
int cnt_x = 0;
int cnt_y = 0;

//函数声明

// TODO
// return the min value
double compare_min(double a1, double a2, double a3, double a4);
// TODO
// return the max value
double compare_max(double a1, double a2, double a3, double a4);

// TODO
double** storage_net(double(*net_input)[16]);

// TODO
double** minimum_spacing(double net_input[][16], double** net_outline);

// all coordinates sorting
void quick_sort(double* arr, int low, int high);

// TODO
// ceil
int ceil_numbers_x(double a);

// extend related
double* extend_net_coordinate_y(double** block_net_p, int* count_epitaxial_del_cloumn_y);
double* extend_net_coordinate_x(double** block_net_p, int* count_epitaxial_del_cloumn_x);

// mesh related
double* mesh_divide_x_p(double* ext_coor_x_p, double** net, double** block, int* cnt_clmn_x, double** ext_block_coordinate, double** ext_net_coordinate, double** net_div_min);
double* mesh_divide_y_p(double* ext_coor_y_p, double** net, double** block, int* cnt_clmn_y, double** ext_block_coordinate, double** ext_net_coordinate, double** net_div_min);

// TODO
double *assembly_matrix(double *inter_x_p, double *inter_y_p, double **ext_block, double **net_outline, double **ext_net_coordinate);

// begin after solver
double** calculation_capacitance(double* inter_x_p, double* inter_y_p, double** ext_net_coordinate, double** net_div_min, double* vector_x);


// TODO
double boundary[4] = { -10, 0, 10, 9.9 };

double min_spac_x = 0;
double min_spac_y = 0;

// TODO
// may cause error when they are not global var
int cnt_epitaxial_del_cloumn_x = 0;
int cnt_epitaxial_del_cloumn_y = 0;


// mesh related
double DIV_NUMS = 8;
double EXTEND_MESH_x = 8;
double EXTEND_MESH_y = 8;

// parse related
double numbs_trapezoid = 0;
double net_input_row =  9;
double net_input[9][16] = { {-0.017, 0.072, 0.017, 0.139, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                           {-0.114, 0.072, -0.078, 0.139, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                           {0.078, 0.072, 0.114, 0.139, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                           {-0.155, 0.105, -0.133, 0.128, -0.156, 0.128, -0.132, 0.152, -0.157, 0.152, -0.131, 0.176, -0.158, 0.176, -0.13, 0.2},
                           {0.133, 0.105, 0.155, 0.128, 0.132, 0.128, 0.156, 0.152, 0.131, 0.152, 0.157, 0.176, 0.13, 0.176, 0.158, 0.2},
                           {-0.059, 0.105, -0.037, 0.128, -0.06, 0.128, -0.036, 0.152, -0.061, 0.152, -0.035, 0.176, -0.062, 0.176, -0.034, 0.2},
                           {0.037, 0.105, 0.059, 0.128, 0.036, 0.128, 0.060, 0.152, 0.035, 0.152, 0.061, 0.176, 0.034, 0.176, 0.062, 0.200},
                           {-0.061, 0.247, -0.024, 0.261, -0.062, 0.261, -0.023, 0.276, -0.063, 0.276, -0.022, 0.290, -0.064, 0.290, -0.021, 0.304},
                           {0.024, 0.247, 0.061, 0.261, 0.023, 0.261, 0.062, 0.276, 0.022, 0.276, 0.063, 0.290, 0.021, 0.290, 0.064, 0.304} };

int main(void)
{

    //开始时间
    long int start = clock();

    // TODO
    // opt first
    double** storage_net_p = storage_net(net_input);
    double** block_net_p = minimum_spacing(net_input, storage_net_p);

    double* extend_net_coordinate_x_p = extend_net_coordinate_x(block_net_p, &cnt_epitaxial_del_cloumn_x);
    double* extend_net_coordinate_y_p = extend_net_coordinate_y(block_net_p, &cnt_epitaxial_del_cloumn_y);

    //计算外扩的net在网格上的坐标
    double** ext_net_coordinate = (double**)malloc(net_input_row * sizeof(double));
    if (ext_net_coordinate == NULL)
    {
        printf("Failed to allocate memory\n");
    }
    else
    {
        for (int i = 0; i < net_input_row; ++i)
        {
            ext_net_coordinate[i] = (double*)malloc(4 * sizeof(double));
            if (ext_net_coordinate[i] == NULL)
            {
                printf("Failed to allocate memory\n");
            }
        }
    }

    //计算外扩的block在网格上的坐标
    double** ext_block_coordinate = (double**)malloc((net_input_row + 3 * numbs_trapezoid) * sizeof(double));
    if (ext_block_coordinate == NULL)
    {
        printf("Failed to allocate memory\n");
    }
    else
    {
        for (int i = 0; i < net_input_row + 3 * numbs_trapezoid; ++i)
        {
            ext_block_coordinate[i] = (double*)malloc(4 * sizeof(double));
            if (ext_block_coordinate[i] == NULL)
            {
                printf("Failed to allocate memory\n");
            }
        }
    }

    double** net_div_min = (double**)malloc(net_input_row * sizeof(double));
    if (net_div_min == NULL)
    {
        printf("net_div_min failed to allocate memory\n");
    }
    else
    {
        for (int i = 0; i < net_input_row; ++i)
        {
            net_div_min[i] = (double*)malloc(4 * sizeof(double));
            if (net_div_min[i] == NULL)
            {
                printf(" net_div_min Failed to allocate memory\n");
            }
        }
    }

    double* inter_x_p = mesh_divide_x_p(extend_net_coordinate_x_p, storage_net_p, block_net_p, &cnt_epitaxial_del_cloumn_x, ext_block_coordinate, ext_net_coordinate, net_div_min);
    double* inter_y_p = mesh_divide_y_p(extend_net_coordinate_y_p, storage_net_p, block_net_p, &cnt_epitaxial_del_cloumn_y, ext_block_coordinate, ext_net_coordinate, net_div_min);

    free(extend_net_coordinate_x_p);
    free(extend_net_coordinate_y_p);

    double * vector_x = assembly_matrix(inter_x_p, inter_y_p, ext_block_coordinate, storage_net_p, ext_net_coordinate);

    double** capacitance_p = calculation_capacitance(inter_x_p, inter_y_p, ext_net_coordinate, net_div_min, vector_x);

   
    free(storage_net_p);
    free(block_net_p);
    free(ext_net_coordinate);
    free(ext_block_coordinate);

    long int finish = clock(); //结束时间
    printf("\ntime: %fs\n", (double) (finish - start) / CLOCKS_PER_SEC);

    return 0;
}

double** calculation_capacitance(double* inter_x, double* inter_y, double** ext_net_coordinate, double** net_div_min, double* vector_x)
{

    int length_z_range = 6;
    int length_gradient_range = 3;

    int max_row = 0;
    int max_clumn = 0;
    for (int u = 0; u < net_input_row; u++)
    {
        if (max_row < round(ext_net_coordinate[u][3]) - round(ext_net_coordinate[u][1]) + 2 * length_z_range - 1)
            max_row = round(ext_net_coordinate[u][3]) - round(ext_net_coordinate[u][1]) + 2 * length_z_range - 1;
        if (max_clumn < round(ext_net_coordinate[u][2]) - round(ext_net_coordinate[u][0]) + 2 * length_z_range - 1)
            max_clumn = round(ext_net_coordinate[u][2]) - round(ext_net_coordinate[u][0]) + 2 * length_z_range - 1;
    }

    double** electric_potential = (double**)malloc((max_row + 1) * sizeof(double));
    if (electric_potential == NULL)
    {
        printf("electric_potential failed to allocate memory\n");
    }
    else
    {
        for (int i = 0; i < (max_row + 1); ++i)
        {
            electric_potential[i] = (double*)malloc((max_clumn + 1) * sizeof(double));
            if (electric_potential[i] == NULL)
            {
                printf("electric_potential[i] failed to allocate memory\n");
            }
        }
    }

    double* gradient_x_upward = (double*)malloc(max_row * sizeof(double));
    if (gradient_x_upward == NULL)
    {
        printf("gradient_x_upward failed to allocate memory\n");
    }

    double* gradient_x_down = (double*)malloc(max_row * sizeof(double));
    if (gradient_x_down == NULL)
    {
        printf("gradient_x_down failed to allocate memory\n");
    }

    double* gradient_y_left = (double*)malloc(max_clumn * sizeof(double));
    if (gradient_y_left == NULL)
    {
        printf("gradient_y_left failed to allocate memory\n");
    }

    double* gradient_y_rigt = (double*)malloc(max_clumn * sizeof(double));
    if (gradient_y_rigt == NULL)
    {
        printf("gradient_y_rigt failed to allocate memory\n");
    }

    double** capacitance = (double**)malloc(3 * sizeof(double));
    if (capacitance == NULL)
    {
        printf("capacitance failed to allocate memory\n");
    }
    else
    {
        for (int i = 0; i < 3; ++i)
        {
            capacitance[i] = (double*)malloc(net_input_row * sizeof(double));
            if (capacitance[i] == NULL)
            {
                printf("capacitance[i] failed to allocate memory\n");
            }
        }
    }

    int temp_j = 0;
    int temp_x_1 = 0;
    int temp_x_2 = 0;
    int temp_y_1 = 0;
    int temp_y_2 = 0;
    int p = 0;
    int q = 0;
    int s = 0;
    int t = 0;
    for (int u = 0; u < net_input_row; u++)
    {

        //将包括net的范围的电势按位置坐标存进块状矩阵z中
        for (int j = 0; j < round(ext_net_coordinate[u][3]) - round(ext_net_coordinate[u][1]) + 2 * length_z_range - 1; j++)
        {
            for (int i = 0; i < round(ext_net_coordinate[u][2]) - round(ext_net_coordinate[u][0]) + 2 * length_z_range - 1; i++)
            {
                electric_potential[j][i] = vector_x[(int)((j + round(ext_net_coordinate[u][1]) - length_z_range) * cnt_x + i + round(ext_net_coordinate[u][0]) - length_z_range + 1)];
            }
        }

        //求梯度
        temp_j = round(ext_net_coordinate[u][3] - round(ext_net_coordinate[u][1])) + 2 * length_z_range - 1 - length_gradient_range;
        temp_x_1 = round(ext_net_coordinate[u][2] - round(ext_net_coordinate[u][0])) + 2 * length_z_range - 1 - length_gradient_range - 1;
        temp_x_2 = round(ext_net_coordinate[u][2] - round(ext_net_coordinate[u][0])) + 2 * length_z_range - length_gradient_range;
        for (int j = length_gradient_range; j < temp_j; j++)
        {
            gradient_x_upward[j] = (electric_potential[j][length_gradient_range + 1] - electric_potential[j][length_gradient_range - 1]) / (2 * net_div_min[u][0]);
            gradient_x_down[j] = (electric_potential[j][temp_x_1] - electric_potential[j][temp_x_2]) / (2 * net_div_min[u][2]);
        }

        temp_j = round(ext_net_coordinate[u][2] - round(ext_net_coordinate[u][0])) + 2 * length_z_range - 1 - length_gradient_range;
        temp_y_1 = round(ext_net_coordinate[u][3] - round(ext_net_coordinate[u][1])) + 2 * length_z_range - 1 - length_gradient_range - 1;
        temp_y_2 = round(ext_net_coordinate[u][3] - round(ext_net_coordinate[u][1])) + 2 * length_z_range - length_gradient_range;

        for (int j = length_gradient_range; j < temp_j; j++)
        {
            gradient_y_left[j] = (electric_potential[length_gradient_range + 1][j] - electric_potential[length_gradient_range - 1][j]) / (2 * net_div_min[u][1]);
            gradient_y_rigt[j] = (electric_potential[temp_y_1][j] - electric_potential[temp_y_2][j]) / (2 * net_div_min[u][3]);
        }

        for (int i = length_gradient_range; i < temp_x_1; i++)
        {
            if (gradient_y_left[i] < 0)
            {
                gradient_y_left[i] = -gradient_y_left[i];
            }
            if (gradient_y_left[i + 1] < 0)
            {
                gradient_y_left[i + 1] = -gradient_y_left[i + 1];
            }
            p = p + inter_x[(int)(i + round(ext_net_coordinate[u][0]) - 7)] * (gradient_y_left[i] - gradient_y_left[i + 1]) / 2;
        }

        for (int i = length_gradient_range; i < temp_x_1; i++)
        {
            if (gradient_y_rigt[i] < 0)
            {
                gradient_y_rigt[i] = -gradient_y_rigt[i];
            }
            if (gradient_y_rigt[i + 1] < 0)
            {
                gradient_y_rigt[i + 1] = -gradient_y_rigt[i + 1];
            }
            q = q + inter_x[(int)(i + round(ext_net_coordinate[u][0]) - 7)] * (gradient_y_rigt[i] - gradient_y_rigt[i + 1]) / 2;
        }

        for (int i = length_gradient_range; i < temp_y_1; i++)
        {
            if (gradient_x_upward[i] < 0)
            {
                gradient_x_upward[i] = -gradient_x_upward[i];
            }
            if (gradient_x_upward[i + 1] < 0)
            {
                gradient_x_upward[i + 1] = -gradient_x_upward[i + 1];
            }
            s = s + inter_y[(int)(i + round(ext_net_coordinate[u][1]) - 7)] * (gradient_x_upward[i] - gradient_x_upward[i + 1]) / 2;
        }

        for (int i = length_gradient_range; i < temp_y_1; i++)
        {
            if (gradient_x_down[i] < 0)
            {
                gradient_x_down[i] = -gradient_x_down[i];
            }
            if (gradient_x_down[i + 1] < 0)
            {
                gradient_x_down[i + 1] = -gradient_x_down[i + 1];
            }
            t = t + inter_y[(int)(i + round(ext_net_coordinate[u][1]) - 7)] * (gradient_x_down[i] - gradient_x_down[i + 1]) / 2;
        }

        capacitance[0][u] = dielectric_constant * 8.854 * (q + p + t + s) / 1000;
        capacitance[2][u] = capacitance[0][u] / capacitance[1][u];
    }
    free(electric_potential);
    free(gradient_x_upward);
    free(gradient_x_down);
    free(gradient_y_left);
    free(gradient_y_rigt);

    free(vector_x);

    return capacitance;
}

double *assembly_matrix(double *inter_x_p, double *inter_y_p, double **ext_block, double **net_outline, double **ext_net_coordinate)
{

    long int cnt_b_numbs = 0;
    long int cnt_A_3_clmn_1 = 0;
    double **matrix_A_3 = (double **)malloc((cnt_x + 1) * (cnt_y + 1) * sizeof(double));
    if (matrix_A_3 == NULL)
    {
        printf("matrix_A Failed to allocate memory\n");
    }
    else
    {
        for (int i = 0; i < (cnt_x + 1) * (cnt_y + 1); ++i)
        {
            matrix_A_3[i] = (double *)malloc(3 * sizeof(double));
            if (matrix_A_3[i] == NULL)
            {
                printf("matrix_A[i] Failed to allocate memory\n");
            }
        }
    }

    double **matrix_A_2 = (double **)malloc(((cnt_x + 1) * (cnt_y + 1) - cnt_x - 1) * sizeof(double));
    if (matrix_A_2 == NULL)
    {
        printf("matrix_A Failed to allocate memory\n");
    }
    else
    {
        for (int i = 0; i < ((cnt_x + 1) * (cnt_y + 1) - cnt_x - 1); ++i)
        {
            matrix_A_2[i] = (double *)malloc(2 * sizeof(double));
            if (matrix_A_2[i] == NULL)
            {
                printf("matrix_A[i] Failed to allocate memory\n");
            }
        }
    }

    int flag = 0;

    for (int j = 0; j < cnt_y + 1; j++)
    {
        for (int i = 0; i < cnt_x + 1; i++)
        {

            for (int k = 0; k < net_input_row + 3 * numbs_trapezoid; k++)
            {
                if ((i <= ext_block[k][2] + 0.01 * min_spac_x) && (i >= ext_block[k][0] + 0.01 * min_spac_x) && (j <= ext_block[k][3] + 0.01 * min_spac_y) && (j >= ext_block[k][1] + 0.01 * min_spac_y))
                {

                    matrix_A_3[cnt_A_3_clmn_1][1] = 1;
                    matrix_A_3[cnt_A_3_clmn_1][2] = 0;
                    matrix_A_3[cnt_A_3_clmn_1][0] = 0;

                    matrix_A_2[cnt_A_3_clmn_1][1] = 0;
                    matrix_A_2[cnt_A_3_clmn_1 - (cnt_x + 1)][0] = 0;

                    //break;
                }
            }

            if ((i == 0) || (j == 0) || (i == cnt_x) || (j == cnt_y))
            {

                matrix_A_3[cnt_A_3_clmn_1][1] = 1;
                matrix_A_3[cnt_A_3_clmn_1][2] = 0;
                matrix_A_3[cnt_A_3_clmn_1][0] = 0;

                if (cnt_A_3_clmn_1 < (cnt_y + 1) * (cnt_x + 1) - (cnt_x + 1))
                {

                    matrix_A_2[cnt_A_3_clmn_1][1] = 0;
                }

                if (cnt_A_3_clmn_1 > (cnt_x))
                {

                    matrix_A_2[cnt_A_3_clmn_1 - (cnt_x + 1)][0] = 0;
                }
            }

            else
            {
                //对角线+1
                matrix_A_3[cnt_A_3_clmn_1][2] = (inter_y_p[j] + inter_y_p[j - 1]) / inter_x_p[i];

                //对角线
                matrix_A_3[cnt_A_3_clmn_1][1] = -(inter_y_p[j] + inter_y_p[j - 1]) / inter_x_p[i] - (inter_y_p[j] + inter_y_p[j - 1]) / inter_x_p[i - 1] - (inter_x_p[i] + inter_x_p[i - 1]) / inter_y_p[j] - (inter_x_p[i] + inter_x_p[i - 1]) / inter_y_p[j - 1];

                //对角线-1
                matrix_A_3[cnt_A_3_clmn_1][0] = (inter_y_p[j] + inter_y_p[j - 1]) / inter_x_p[i - 1];

                if (cnt_A_3_clmn_1 < (cnt_y + 1) * (cnt_x + 1) - (cnt_x + 1))
                {

                    //对角线+n
                    matrix_A_2[cnt_A_3_clmn_1][1] = (inter_x_p[i] + inter_x_p[i - 1]) / inter_y_p[j];
                }

                if (cnt_A_3_clmn_1 > (cnt_x))
                {

                    //对角线-n
                    matrix_A_2[cnt_A_3_clmn_1 - (cnt_x + 1)][0] = (inter_x_p[i] + inter_x_p[i - 1]) / inter_y_p[j - 1];
                }
            }

            cnt_A_3_clmn_1++;
        }
    }

    //如果是梯形会多出1000多个

    int cnt_b_numbs_2 = (round(ext_net_coordinate[0][2]) - round(ext_net_coordinate[0][0])) * (round(ext_net_coordinate[0][3]) - round(ext_net_coordinate[0][1]));

    //这个地方需要修改
    double *coordinate_b = (double *)malloc(cnt_b_numbs_2 * sizeof(double));
    if (coordinate_b == NULL)
    {
        printf("matrix_A Failed to allocate memory\n");
    }

    cnt_b_numbs = 0;

    if ((net_outline[0][4] == 2) || (net_outline[0][5] == 1))
    {
        for (int u = 0; u < 4; u++)
        {
            for (int i = round(ext_block[u][0]); i < round(ext_block[u][2]); i++)
            {
                for (int j = round(ext_block[u][1]); j < round(ext_block[u][3]); j++)
                {
                    coordinate_b[cnt_b_numbs] = j * cnt_x + i;
                    cnt_b_numbs++;
                }
            }
        }
    }
    else
    {
        for (int i = round(ext_block[0][0]); i < round(ext_block[0][2]); i++)
        {
            for (int j = round(ext_block[0][1]); j < round(ext_block[0][3]); j++)
            {
                coordinate_b[cnt_b_numbs] = j * cnt_x + i;
                cnt_b_numbs++;
            }
        }
    }

    double *vector_x = (double *) malloc((cnt_x + 1) * (cnt_y + 1) * sizeof(double));
    if (vector_x == NULL)
    {
        printf("matrix_A Failed to allocate memory\n");
    }
    ///////////////////////开始解矩阵///////////////

    //////////////////////////////////////////////

    free(matrix_A_3);
    free(matrix_A_2);
    free(coordinate_b);
    return vector_x;
}

double* mesh_divide_y_p(double* ext_coor_y_p, double** net, double** block, int* cnt_clmn_y, double** ext_block_coordinate, double** ext_net_coordinate, double** net_div_min)
{
    //计算上边的0.7部分
    double long_cnt_high = 5;
    double long_high = pow(1.1, 5.0);
    double long_mul_high = pow(1.1, 5.0);
    double div_numbs_high = (boundary[3] - ext_coor_y_p[*cnt_clmn_y - 1]) * PERCENTAGE_DIV_0_7 / min_spac_y;

    while (((long_high - 1) / (1.1 - 1)) < div_numbs_high)
    {
        long_high = long_high * long_mul_high;
        long_cnt_high = long_cnt_high + 5;
    }
    //就按上边0.3部分
    int long_cnt_high_03 = round((div_numbs_high * 3 / 7) / pow(1.1, long_cnt_high - 5));

    double long_cnt_low = 5;
    double long_low = pow(1.1, 5);
    double long_mul_low = pow(1.1, 5);
    double div_numbs_low = (ext_coor_y_p[0] - boundary[1]) * PERCENTAGE_DIV_0_7 / min_spac_y;
    double* dividing_line_y = (double*)malloc(*cnt_clmn_y * sizeof(double));
    if (dividing_line_y == NULL)
    {
        printf("dividing_line_y Failed to allocate memory\n");
    }
    while ((long_low - 1) / (1.1 - 1) < div_numbs_low)
    {
        long_low = long_low * long_mul_low;
        long_cnt_low = long_cnt_low + 5;
    }

    int long_cnt_low_03 = round((div_numbs_low * 3 / 7) / pow(1.1, long_cnt_low - 5));
    int long_cnt_y = long_cnt_high + long_cnt_low + long_cnt_high_03 + long_cnt_low_03;
    //在这些线段里进行细密的划分，划分间隔与之前定的最小划分间隔相近，uniform_filigree_divide_numbs_x是每个的均匀细划分区间隔数，uniform_filigree_divide_interval_x是每个均�?细划分区的间�?
    int* unif_div_numbs_y = (int*)malloc(*cnt_clmn_y / 2 * sizeof(int)); //uniform_filigree_divide_numbs_y
    if (unif_div_numbs_y == NULL)
    {
        printf("unif_div_numbs_y Failed to allocate memory\n");
    }

    double* unif_div_inter_y = (double*)malloc(*cnt_clmn_y / 2 * sizeof(double)); //uniform_filigree_divide_interval_y
    if (unif_div_inter_y == NULL)
    {
        printf("unif_div_inter_y Failed to allocate memory\n");
    }

    for (int i = 0; i < *cnt_clmn_y - 1; i = i + 2)
    {
        unif_div_numbs_y[i / 2] = round((ext_coor_y_p[i + 1] - ext_coor_y_p[i]) / min_spac_y);
        unif_div_inter_y[i / 2] = (ext_coor_y_p[i + 1] - ext_coor_y_p[i]) / round((ext_coor_y_p[i + 1] - ext_coor_y_p[i]) / min_spac_y);
        long_cnt_y = long_cnt_y + unif_div_numbs_y[i / 2];
    }

    double long_cnt_middle = 5;
    double long_middle = pow(1.1, 5);
    double long_mul_middle = pow(1.1, 5);
    double div_numbs_middle = 0;

    double max_long_cnt_middle = 0;
    for (int i = 1; i < *cnt_clmn_y - 1; i = i + 2)
    {
        long_cnt_middle = 5;
        long_middle = pow(1.1, 5);
        long_mul_middle = pow(1.1, 5);
        div_numbs_middle = ((ext_coor_y_p[i + 1] - ext_coor_y_p[i]) * PERCENTAGE_DIV_0_7 / min_spac_y) / 2;
        while ((long_middle - 1) / (1.1 - 1) < div_numbs_middle)
        {
            long_middle = long_middle * long_mul_middle;
            long_cnt_middle = long_cnt_middle + 5;
        }

        long_cnt_y = long_cnt_y + 2 * long_cnt_middle + round((2 * div_numbs_middle * 3 / 7) / pow(1.1, long_cnt_middle - 5));
        if (max_long_cnt_middle < long_cnt_middle)
        {
            max_long_cnt_middle = long_cnt_middle;
        }
    }

    //这一部分的代码是对x轴方向靠近boundary(3)的最后一块区域进行划分，先非均匀划分70%的部分，剩下的均匀划分
    double long_dived_high = 0;                                                               //已经划分了的区域长度
    double* non_unif_div_numbs_high = (double*)malloc((long_cnt_high * 2) * sizeof(double)); //uniform_filigree_divide_numbs_x
    if (non_unif_div_numbs_high == NULL)
    {
        printf("non_unif_div_numbs_high Failed to allocate memory\n");
    }
    int cnt_non_unif_inter_high = 0;

    non_unif_div_numbs_high[cnt_non_unif_inter_high] = unif_div_inter_y[ceil_numbers_x((double)*cnt_clmn_y - 2.0) / 2];
    int unif_div_numbs_high = 0;
    double long_unif_dived_high = 0;
    if (boundary[3] - (EXTEND_MESH_x / 2) * min_spac_y > ext_coor_y_p[*cnt_clmn_y - 1])
    {
        //非均匀划分
        while (long_dived_high < PERCENTAGE_DIV_0_7 * boundary[3] - ext_coor_y_p[*cnt_clmn_y - 1])
        {
            long_dived_high = long_dived_high + non_unif_div_numbs_high[cnt_non_unif_inter_high];
            non_unif_div_numbs_high[cnt_non_unif_inter_high + 1] = non_unif_div_numbs_high[cnt_non_unif_inter_high] * 1.1;

            cnt_non_unif_inter_high++;
        }
        //均匀划分
        unif_div_numbs_high = round((boundary[3] - ext_coor_y_p[*cnt_clmn_y - 1] - long_dived_high) / non_unif_div_numbs_high[cnt_non_unif_inter_high - 1]);
        long_unif_dived_high = (boundary[3] - ext_coor_y_p[*cnt_clmn_y - 1] - long_dived_high) / unif_div_numbs_high;
    }

    ////这一部分的代码是对x轴方向靠近boundary(0)的最后一块区域进行划分，先非均匀划分70%的部分，剩下的均匀划分
    double long_dived_low = 0;                                                        //已经划分了的区域长度
    double* non_unif_div_numbs_low = (double*)malloc(long_cnt_low * sizeof(double)); //uniform_filigree_divide_numbs_x
    if (non_unif_div_numbs_low == NULL)
    {
        printf("non_unif_div_numbs_low Failed to allocate memory\n");
    }
    int cnt_non_unif_inter_low = 0;
    non_unif_div_numbs_low[cnt_non_unif_inter_low] = unif_div_inter_y[0];
    int unif_div_numbs_low = 0;
    double long_unif_dived_low = 0;
    double* inter_y = (double*)malloc(long_cnt_y * sizeof(double));
    if (inter_y == NULL)
    {
        printf("inter_y Failed to allocate memory\n");
    }
    int cnt_dividing_line_y = 0;
    if (boundary[1] + (EXTEND_MESH_x / 2) * min_spac_y < ext_coor_y_p[0])
    {
        //非均匀划分
        while (long_dived_low < PERCENTAGE_DIV_0_7 * (ext_coor_y_p[0] - boundary[1]))
        {
            long_dived_low = long_dived_low + non_unif_div_numbs_low[cnt_non_unif_inter_low];

            non_unif_div_numbs_low[cnt_non_unif_inter_low + 1] = non_unif_div_numbs_low[cnt_non_unif_inter_low] * 1.1;
            cnt_non_unif_inter_low++;
        }
        //均匀划分
        unif_div_numbs_low = round((ext_coor_y_p[0] - boundary[1] - long_dived_low) / non_unif_div_numbs_low[cnt_non_unif_inter_low - 1]);
        long_unif_dived_low = (ext_coor_y_p[0] - boundary[1] - long_dived_low) / unif_div_numbs_low;

        for (int i = 0; i < unif_div_numbs_low; i++)
        {
            inter_y[cnt_y] = long_unif_dived_low;

            cnt_y++;
        }

        for (int i = 0; i < cnt_non_unif_inter_low; i++)
        {

            inter_y[cnt_y] = non_unif_div_numbs_low[cnt_non_unif_inter_low - i - 1];

            cnt_y++;
        }

        // %下面两行专门来记录细划分区边界的格点坐标
        dividing_line_y[cnt_dividing_line_y] = cnt_y;
        cnt_dividing_line_y++;
    }

    free(non_unif_div_numbs_low);
    double long_dived = 0;                                    //已经划分了的区域长度
    double* array_y = (double*)malloc(100 * sizeof(double)); //uniform_filigree_divide_numbs_x
    if (array_y == NULL)
    {
        printf("array_y Failed to allocate memory\n");
    }

    //double array_y[100];
    int cnt_temp = 0;
    int unif_div_numbs_middle;
    double long_unif_dived_middle;
    for (int i = 1; i < *cnt_clmn_y - 1; i = i + 2)
    {
        // 确定非均匀的划分区间
        long_dived = 0;
        cnt_temp = 0;
        array_y[0] = min_spac_y;

        // 两边70%区域同时进行非均匀划分
        while ((ext_coor_y_p[i + 1] - ext_coor_y_p[i]) * PERCENTAGE_DIV_0_7 > 2 * long_dived)
        {

            long_dived = long_dived + array_y[cnt_temp];
            array_y[cnt_temp + 1] = array_y[cnt_temp] * 1.1;
            cnt_temp++;
        }

        ////中间部分粗均匀划分
        unif_div_numbs_middle = round((ext_coor_y_p[i + 1] - ext_coor_y_p[i] - 2 * long_dived) / array_y[cnt_temp - 1]);
        long_unif_dived_middle = (ext_coor_y_p[i + 1] - ext_coor_y_p[i] - 2 * long_dived) / unif_div_numbs_middle;

        //填装中间细均匀部分
        for (int j = 0; j < unif_div_numbs_y[(i - 1) / 2]; j++)
        {
            inter_y[cnt_y] = unif_div_inter_y[(i - 1) / 2];
            cnt_y++;
        }

        // %下面两行专门来记录细划分区边界的格点坐标
        dividing_line_y[cnt_dividing_line_y] = cnt_y;
        cnt_dividing_line_y++;
        //填充非均匀间隔到inter_y
        for (int j = 0; j < cnt_temp; j++)
        {
            inter_y[cnt_y] = array_y[j];
            cnt_y++;
        }

        //填充粗均匀间隔
        for (int j = 0; j < unif_div_numbs_middle; j++)
        {
            inter_y[cnt_y] = long_unif_dived_middle;
            cnt_y++;
        }

        //填充非均匀划分间隔
        for (int j = 0; j < cnt_temp; j++)
        {
            inter_y[cnt_y] = array_y[cnt_temp - j - 1];
            cnt_y++;
        }

        // %下面两行专门来记录细划分区边界的格点坐标
        dividing_line_y[cnt_dividing_line_y] = cnt_y;
        cnt_dividing_line_y++;
    }

    // free(array_y);

    //最后一均匀划分填充
    //可能有问题
    for (int j = 0; j < unif_div_numbs_y[(*cnt_clmn_y - 2) / 2]; j++)
    {
        inter_y[cnt_y] = unif_div_inter_y[(*cnt_clmn_y - 1) / 2];
        cnt_y++;
    }
    free(array_y);
    free(unif_div_numbs_y);
    // %下面两行专门来记录细划分区边界的格点坐标
    dividing_line_y[cnt_dividing_line_y] = cnt_y;
    cnt_dividing_line_y++;

    //填充整个inter_y的靠近boundary(2)部分
    if (boundary[2] - (EXTEND_MESH_x / 2) * min_spac_y > ext_coor_y_p[*cnt_clmn_y - 1])
    {
        for (int i = 0; i < cnt_non_unif_inter_high; i++)
        {
            inter_y[cnt_y] = non_unif_div_numbs_high[i];
            cnt_y++;
        }

        for (int i = 0; i < unif_div_numbs_high; i++)
        {
            inter_y[cnt_y] = long_unif_dived_high;
            cnt_y++;
        }
    }

    ////
    free(non_unif_div_numbs_high);
    //计算外扩的net在网格上的坐标

    for (int i = 0; i < net_input_row; i++)
    {
        for (int j = 0; j < *cnt_clmn_y - 1; j = j + 2)
        {
            if ((net[i][1] > ext_coor_y_p[j]) && (net[i][1] < ext_coor_y_p[j + 1]))
            {
                ext_net_coordinate[i][1] = (net[i][1] - ext_coor_y_p[j]) / unif_div_inter_y[j / 2] + dividing_line_y[j];
                net_div_min[i][1] = unif_div_inter_y[j / 2];
            }
            if ((net[i][3] > ext_coor_y_p[j]) && (net[i][3] < ext_coor_y_p[j + 1]))
            {
                ext_net_coordinate[i][3] = (net[i][3] - ext_coor_y_p[j]) / unif_div_inter_y[j / 2] + dividing_line_y[j];
                net_div_min[i][3] = unif_div_inter_y[j / 2];
            }
        }
    }

    //计算外扩的block在网格上的坐标

    for (int i = 0; i < net_input_row + 3 * numbs_trapezoid; i++)
    {
        for (int j = 0; j < *cnt_clmn_y - 1; j = j + 2)
        {
            if ((block[i][1] > ext_coor_y_p[j]) && (block[i][1] < ext_coor_y_p[j + 1]))
            {

                ext_block_coordinate[i][1] = (block[i][1] - ext_coor_y_p[j]) / unif_div_inter_y[j / 2] + dividing_line_y[j];
            }
            if ((block[i][3] > ext_coor_y_p[j]) && (block[i][0] < ext_coor_y_p[j + 1]))
            {
                ext_block_coordinate[i][3] = (block[i][3] - ext_coor_y_p[j]) / unif_div_inter_y[j / 2] + dividing_line_y[j];
            }
        }
    }

    free(unif_div_inter_y);
    free(dividing_line_y);
    return inter_y;
}

double* mesh_divide_x_p(double* ext_coor_x_p, double** net, double** block, int* cnt_clmn_x, double** ext_block_coordinate, double** ext_net_coordinate, double** net_div_min)
{

    //计算右边0.7部分
    double long_cnt_right = 5;
    double long_right = pow(1.1, 5.0);
    double long_mul_right = pow(1.1, 5.0);
    double div_numbs_right = (boundary[2] - ext_coor_x_p[*cnt_clmn_x - 1]) * PERCENTAGE_DIV_0_7 / min_spac_x;

    while (((long_right - 1) / (1.1 - 1)) < div_numbs_right)
    {
        long_right = long_right * long_mul_right;
        long_cnt_right = long_cnt_right + 5;
    }
    //就按右边0.3部分
    int long_cnt_right_03 = round((div_numbs_right * 3 / 7) / pow(1.1, long_cnt_right - 5));
    //计算左边0.7部分
    double long_cnt_left = 5;
    double long_left = pow(1.1, 5);

    double long_mul_left = pow(1.1, 5);
    double div_numbs_left = (ext_coor_x_p[1] - boundary[0]) * PERCENTAGE_DIV_0_7 / min_spac_x;

    while ((long_left - 1) / (1.1 - 1) < div_numbs_left)
    {
        long_left = long_left * long_mul_left;
        long_cnt_left = long_cnt_left + 5;
    }

    int long_cnt_left_03 = ((div_numbs_left * 3 / 7) / pow(1.1, long_cnt_left - 5));
    int long_cnt_x = long_cnt_right + long_cnt_left + long_cnt_right_03 + long_cnt_left_03;

    double* unif_div_numbs_x = (double*)malloc((*cnt_clmn_x - 2) * sizeof(double)); //uniform_filigree_divide_numbs_x
    if (unif_div_numbs_x == NULL)
    {
        printf("unif_div_numbs_x Failed to allocate memory\n");
    }

    double* unif_div_inter_x = (double*)malloc((*cnt_clmn_x - 2) * sizeof(double)); //uniform_filigree_divide_interval_x
    if (unif_div_inter_x == NULL)
    {
        printf("unif_div_inter_x Failed to allocate memory\n");
    }

    for (int i = 0; i < *cnt_clmn_x - 1; i = i + 2)
    {

        unif_div_numbs_x[i / 2] = round((ext_coor_x_p[i + 1] - ext_coor_x_p[i]) / min_spac_x);
        unif_div_inter_x[i / 2] = (ext_coor_x_p[i + 1] - ext_coor_x_p[i]) / round((ext_coor_x_p[i + 1] - ext_coor_x_p[i]) / min_spac_x);

        long_cnt_x = long_cnt_x + unif_div_numbs_x[i / 2];
    }

    double long_cnt_middle = 5;
    double long_middle = pow(1.1, 5);
    double long_mul_middle = pow(1.1, 5);
    double div_numbs_middle = 0;

    double max_long_cnt_middle = 0;
    for (int i = 1; i < *cnt_clmn_x - 2; i = i + 2)
    {
        long_cnt_middle = 5;
        long_middle = pow(1.1, 5);
        long_mul_middle = pow(1.1, 5);
        div_numbs_middle = ((ext_coor_x_p[i + 1] - ext_coor_x_p[i]) * PERCENTAGE_DIV_0_7 / min_spac_x) / 2;
        while ((long_middle - 1) / (1.1 - 1) < div_numbs_middle)
        {
            long_middle = long_middle * long_mul_middle;
            long_cnt_middle = long_cnt_middle + 5;
        }

        long_cnt_x = long_cnt_x + 2 * long_cnt_middle + round((2 * div_numbs_middle * 3 / 7) / pow(1.1, long_cnt_middle - 5));
        if (max_long_cnt_middle < long_cnt_middle)
        {
            max_long_cnt_middle = long_cnt_middle;
        }
    }

    //这一部分的代码是对x轴方向靠近boundary(2)的最后一块区域进行划分，先非均匀划分70%的部分，剩下的均匀划分
    double long_dived_right = 0;                                                              //已经划分了的区域长度
    double* non_unif_div_numbs_right = (double*)malloc(long_cnt_right * 2 * sizeof(double)); //uniform_filigree_divide_numbs_x
    if (non_unif_div_numbs_right == NULL)
    {
        printf("non_unif_div_numbs_right Failed to allocate memory\n");
    }
    int cnt_non_unif_inter_right = 0;
    double tt = *cnt_clmn_x - 2.0;
    non_unif_div_numbs_right[cnt_non_unif_inter_right] = unif_div_inter_x[ceil_numbers_x(tt) / 2];
    int unif_div_numbs_right = 0;
    double long_unif_dived_right = 0;
    if (boundary[2] - (EXTEND_MESH_x / 2) * min_spac_x > ext_coor_x_p[*cnt_clmn_x - 1])
    {
        //非均匀划分
        while (long_dived_right < PERCENTAGE_DIV_0_7 * boundary[2] - ext_coor_x_p[*cnt_clmn_x - 1])
        {
            long_dived_right = long_dived_right + non_unif_div_numbs_right[cnt_non_unif_inter_right];
            non_unif_div_numbs_right[cnt_non_unif_inter_right + 1] = non_unif_div_numbs_right[cnt_non_unif_inter_right] * 1.1;

            cnt_non_unif_inter_right++;
        }
        //均匀划分
        unif_div_numbs_right = round((boundary[2] - ext_coor_x_p[*cnt_clmn_x - 1] - long_dived_right) / non_unif_div_numbs_right[cnt_non_unif_inter_right - 1]);
        long_unif_dived_right = (boundary[2] - ext_coor_x_p[*cnt_clmn_x - 1] - long_dived_right) / unif_div_numbs_right;
    }

    ////这一部分的代码是对x轴方向靠近boundary(0)的最后一块区域进行划分，先非均匀划分70%的部分，剩下的均匀划分
    double long_dived_left = 0; //已经划分了的区域长度

    double* non_unif_div_numbs_left = (double*)malloc(long_cnt_left * sizeof(double)); //uniform_filigree_divide_numbs_x
    if (non_unif_div_numbs_left == NULL)
    {
        printf("non_unif_div_numbs_left Failed to allocate memory\n");
    }
    int cnt_non_unif_inter_left = 0;

    non_unif_div_numbs_left[cnt_non_unif_inter_left] = unif_div_inter_x[0];

    int unif_div_numbs_left = 0;
    double long_unif_dived_left = 0;

    double* inter_x = (double*)malloc(long_cnt_x * sizeof(double));
    if (inter_x == NULL)
    {
        printf("inter_x Failed to allocate memory\n");
    }
    int* dividing_line_x = (int*)malloc(*cnt_clmn_x * sizeof(int));
    if (dividing_line_x == NULL)
    {
        printf("dividing_line_x Failed to allocate memory\n");
    }

    int cnt_dividing_line_x = 0;

    //int cnt_x = 0;
    if (boundary[0] + (EXTEND_MESH_x / 2) * min_spac_x < ext_coor_x_p[0])
    {
        //非均匀划分
        while (long_dived_left < PERCENTAGE_DIV_0_7 * (ext_coor_x_p[0] - boundary[0]))
        {
            long_dived_left = long_dived_left + non_unif_div_numbs_left[cnt_non_unif_inter_left];

            non_unif_div_numbs_left[cnt_non_unif_inter_left + 1] = non_unif_div_numbs_left[cnt_non_unif_inter_left] * 1.1;
            cnt_non_unif_inter_left++;
        }
        //均匀划分
        unif_div_numbs_left = round((ext_coor_x_p[0] - boundary[0] - long_dived_left) / non_unif_div_numbs_left[cnt_non_unif_inter_left - 1]);
        long_unif_dived_left = (ext_coor_x_p[0] - boundary[0] - long_dived_left) / unif_div_numbs_left;

        for (int i = 0; i < unif_div_numbs_left; i++)
        {
            inter_x[cnt_x] = long_unif_dived_left;

            cnt_x++;
        }

        for (int i = 0; i < cnt_non_unif_inter_left; i++)
        {

            inter_x[cnt_x] = non_unif_div_numbs_left[cnt_non_unif_inter_left - i - 1];

            cnt_x++;
        }

        // %下面两行专门来记录细划分区边界的格点坐标
        dividing_line_x[cnt_dividing_line_x] = cnt_x;
        cnt_dividing_line_x++;
    }
    free(non_unif_div_numbs_left);

    //下面是个循环，用来填充inter_x中间的部分
    double long_dived = 0; //已经划分了的区域长度

    double* array_x = (double*)malloc((max_long_cnt_middle + 1) * sizeof(double)); //uniform_filigree_divide_numbs_x
    if (array_x == NULL)
    {
        printf("array_x Failed to allocate memory\n");
    }

    //double array_x[50];
    int cnt_temp = 0;
    int unif_div_numbs_middle;
    double long_unif_dived_middle;

    for (int i = 1; i < *cnt_clmn_x - 1; i = i + 2)
    {
        // 确定非均匀的划分区间
        long_dived = 0;
        cnt_temp = 0;
        array_x[0] = min_spac_x;

        // 两边70%区域同时进行非均匀划分
        while ((ext_coor_x_p[i + 1] - ext_coor_x_p[i]) * PERCENTAGE_DIV_0_7 > 2 * long_dived)
        {

            long_dived = long_dived + array_x[cnt_temp];
            array_x[cnt_temp + 1] = array_x[cnt_temp] * 1.1;
            cnt_temp++;
        }

        ////中间部分粗均匀划分
        unif_div_numbs_middle = round((ext_coor_x_p[i + 1] - ext_coor_x_p[i] - 2 * long_dived) / array_x[cnt_temp - 1]);
        long_unif_dived_middle = (ext_coor_x_p[i + 1] - ext_coor_x_p[i] - 2 * long_dived) / unif_div_numbs_middle;

        //填装中间细均匀部分
        for (int j = 0; j < unif_div_numbs_x[(i - 1) / 2]; j++)
        {
            inter_x[cnt_x] = unif_div_inter_x[(i - 1) / 2];
            cnt_x++;
        }

        // %下面两行专门来记录细划分区边界的格点坐标
        dividing_line_x[cnt_dividing_line_x] = cnt_x;
        cnt_dividing_line_x++;
        //填充非均匀间隔到inter_x
        for (int j = 0; j < cnt_temp; j++)
        {
            inter_x[cnt_x] = array_x[j];
            cnt_x++;
        }

        //填充粗均匀间隔
        for (int j = 0; j < unif_div_numbs_middle; j++)
        {
            inter_x[cnt_x] = long_unif_dived_middle;
            cnt_x++;
        }

        //填充非均匀划分间隔
        for (int j = 0; j < cnt_temp; j++)
        {
            inter_x[cnt_x] = array_x[cnt_temp - j - 1];
            cnt_x++;
        }

        // %下面两行专门来记录细划分区边界的格点坐标
        dividing_line_x[cnt_dividing_line_x] = cnt_x;
        cnt_dividing_line_x++;
    }

    //最后一均匀划分填充
    //可能有问题
    for (int j = 0; j < unif_div_numbs_x[(*cnt_clmn_x - 1) / 2]; j++)
    {
        inter_x[cnt_x] = unif_div_inter_x[(*cnt_clmn_x - 1) / 2];
        cnt_x++;
    }
    free(array_x);

    free(unif_div_numbs_x);
    // %下面两行专门来记录细划分区边界的格点坐标
    dividing_line_x[cnt_dividing_line_x] = cnt_x;
    cnt_dividing_line_x++;

    //填充整个inter_x的靠近boundary(2)部分
    if (boundary[2] - (EXTEND_MESH_x / 2) * min_spac_x > ext_coor_x_p[*cnt_clmn_x - 1])
    {
        for (int i = 0; i < cnt_non_unif_inter_right; i++)
        {
            inter_x[cnt_x] = non_unif_div_numbs_right[i];
            cnt_x++;
        }

        for (int i = 0; i < unif_div_numbs_right; i++)
        {
            inter_x[cnt_x] = long_unif_dived_right;
            cnt_x++;
        }
    }

    free(non_unif_div_numbs_right);
    //free(array_x);
    //计算外扩的net在网格上的坐标
    for (int i = 0; i < net_input_row; i++)
    {
        for (int j = 0; j < *cnt_clmn_x - 1; j = j + 2)
        {
            if ((net[i][0] > ext_coor_x_p[j]) && (net[i][0] < ext_coor_x_p[j + 1]))
            {
                ext_net_coordinate[i][0] = (net[i][0] - ext_coor_x_p[j]) / unif_div_inter_x[j / 2] + dividing_line_x[j];
                net_div_min[i][0] = unif_div_inter_x[j / 2];
            }
            if ((net[i][2] > ext_coor_x_p[j]) && (net[i][0] < ext_coor_x_p[j + 1]))
            {
                ext_net_coordinate[i][2] = (net[i][2] - ext_coor_x_p[j]) / unif_div_inter_x[j / 2] + dividing_line_x[j];
                net_div_min[i][2] = unif_div_inter_x[j / 2];
            }
        }
    }

    for (int i = 0; i < net_input_row + 3 * numbs_trapezoid; i++)
    {
        for (int j = 0; j < *cnt_clmn_x - 1; j = j + 2)
        {
            if ((block[i][0] > ext_coor_x_p[j]) && (block[i][0] < ext_coor_x_p[j + 1]))
            {
                ext_block_coordinate[i][0] = (block[i][0] - ext_coor_x_p[j]) / unif_div_inter_x[j / 2] + dividing_line_x[j];
            }
            if ((block[i][2] > ext_coor_x_p[j]) && (block[i][0] < ext_coor_x_p[j + 1]))
            {
                ext_block_coordinate[i][2] = (block[i][2] - ext_coor_x_p[j]) / unif_div_inter_x[j / 2] + dividing_line_x[j];
            }
        }
    }

    free(unif_div_inter_x);
    free(dividing_line_x);
    return inter_x;
}

double* extend_net_coordinate_x(double** block_net_p, int* count_epitaxial_del_cloumn_x)
{
    int net_numbs = net_input_row + numbs_trapezoid * 3;
    int net_numbs_2 = 2 * (net_input_row + numbs_trapezoid * 3);

    int count_row = 0;
    //申请向外延申的内存
    double** epitaxial_extend_x = (double**)malloc(net_numbs_2 * sizeof(double));
    if (epitaxial_extend_x == NULL)
    {
        printf("epitaxial_extend_x Failed to allocate memory\n");
    }
    else
        for (int i = 0; i < net_numbs_2; ++i)
        {
            epitaxial_extend_x[i] = (double*)malloc(2 * sizeof(double));
            if (epitaxial_extend_x[i] == NULL)
            {
                printf("epitaxial_extend_x[i] Failed to allocate memory\n");
            }
        }
    //对所有矩形的四个边进行外扩
    for (int i = 0; i < net_numbs_2; i++)
    {

        epitaxial_extend_x[i][0] = block_net_p[count_row][0] - min_spac_x * EXTEND_MESH_x;
        epitaxial_extend_x[i][1] = block_net_p[count_row][0] + min_spac_x * EXTEND_MESH_x;

        i++;
        epitaxial_extend_x[i][0] = block_net_p[count_row][2] - min_spac_x * EXTEND_MESH_x;
        epitaxial_extend_x[i][1] = block_net_p[count_row][2] + min_spac_x * EXTEND_MESH_x;
        count_row++;
    }

    int count_column_x = 0;
    int* overlap_egments_x = NULL;
    overlap_egments_x = (int*)malloc(sizeof(int) * net_numbs_2);
    if (overlap_egments_x == NULL)
    {
        printf("overlap_egments_x failed to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i < net_numbs_2; i++)
    {
        for (int j = i + 1; j < net_numbs_2; j++)
        {

            if ((epitaxial_extend_x[i][0] >= epitaxial_extend_x[j][0]) && (epitaxial_extend_x[i][1] <= epitaxial_extend_x[j][1]))
            {
                overlap_egments_x[count_column_x] = i;
                count_column_x++;
                break;
            }
            if ((epitaxial_extend_x[i][0] <= epitaxial_extend_x[j][0]) && (epitaxial_extend_x[i][1] >= epitaxial_extend_x[j][1]))
            {
                overlap_egments_x[count_column_x] = i;
                epitaxial_extend_x[j][0] = epitaxial_extend_x[i][0];
                epitaxial_extend_x[j][1] = epitaxial_extend_x[i][1];
                count_column_x++;
                break;
            }
            if ((epitaxial_extend_x[i][0] <= epitaxial_extend_x[j][0]) && (epitaxial_extend_x[j][0] <= epitaxial_extend_x[i][1]) && (epitaxial_extend_x[i][1] <= epitaxial_extend_x[j][1]))
            {
                overlap_egments_x[count_column_x] = i;
                count_column_x++;
                epitaxial_extend_x[j][0] = epitaxial_extend_x[i][0];

                break;
            }
            if ((epitaxial_extend_x[j][0] <= epitaxial_extend_x[i][0]) && (epitaxial_extend_x[i][0] <= epitaxial_extend_x[j][1]) && (epitaxial_extend_x[j][1] <= epitaxial_extend_x[i][1]))
            {
                overlap_egments_x[count_column_x] = i;
                count_column_x++;
                epitaxial_extend_x[j][1] = epitaxial_extend_x[i][1];
                break;
            }
        }
    }

    double* epitaxial_extend_del_x = (double*)malloc(sizeof(double) * net_numbs_2);
    if (epitaxial_extend_del_x == NULL)
    {
        printf("epitaxial_exte  nd_del_x Failed to allocate memory\n");
        exit(1);
    }

    int temp = 0;
    for (int i = 0; i < net_numbs_2; i++)
    {

        for (int j = 0; j < count_column_x; j++)
        {
            if (overlap_egments_x[j] == i)
            {
                temp = 1;
            }
        }
        if (temp == 0)
        {

            epitaxial_extend_del_x[*count_epitaxial_del_cloumn_x] = epitaxial_extend_x[i][0];
            *count_epitaxial_del_cloumn_x = *count_epitaxial_del_cloumn_x + 1;
            epitaxial_extend_del_x[*count_epitaxial_del_cloumn_x] = epitaxial_extend_x[i][1];
            *count_epitaxial_del_cloumn_x = *count_epitaxial_del_cloumn_x + 1;
        }
        else
            temp = 0;
    }

    free(overlap_egments_x);
    free(epitaxial_extend_x);

    quick_sort(epitaxial_extend_del_x, 0, (*count_epitaxial_del_cloumn_x - 1));

    if (epitaxial_extend_del_x[0] < boundary[0] + min_spac_x * round(EXTEND_MESH_x / 2))
    {
        epitaxial_extend_del_x[0] = boundary[0];
    }
    else
    {
        epitaxial_extend_del_x[0] = epitaxial_extend_del_x[0] + min_spac_x * round(EXTEND_MESH_x / 2);
    }

    if (epitaxial_extend_del_x[*count_epitaxial_del_cloumn_x - 1] > boundary[2] - min_spac_x * round(EXTEND_MESH_x / 2))
    {
        epitaxial_extend_del_x[*count_epitaxial_del_cloumn_x - 1] = boundary[2];
    }
    else
    {
        epitaxial_extend_del_x[*count_epitaxial_del_cloumn_x - 1] = epitaxial_extend_del_x[*count_epitaxial_del_cloumn_x - 1] - min_spac_x * round(EXTEND_MESH_x / 2);
    }

    for (int i = 1; i < *count_epitaxial_del_cloumn_x - 1; i++)
    {
        if ((i + 1) % 2 != 0)
        {
            epitaxial_extend_del_x[i] = epitaxial_extend_del_x[i] + min_spac_x * round(EXTEND_MESH_x / 2);
        }
        else
        {
            epitaxial_extend_del_x[i] = epitaxial_extend_del_x[i] - min_spac_x * round(EXTEND_MESH_x / 2);
        }
    }

    return epitaxial_extend_del_x;
}

double* extend_net_coordinate_y(double** block_net_p, int* count_epitaxial_del_cloumn_y)
{
    int net_numbs = net_input_row + numbs_trapezoid * 3;
    int net_numbs_2 = 2 * (net_input_row + numbs_trapezoid * 3);
    int count_row = 0;

    double** epitaxial_extend_y = (double**)malloc(net_numbs_2 * sizeof(double*));
    if (epitaxial_extend_y == NULL)
    {
        printf("epitaxial_extend_y Failed to allocate memory\n");
        exit(1);
    }

    else
    {
        for (int i = 0; i < net_numbs_2; ++i)
        {
            epitaxial_extend_y[i] = (double*)malloc(2 * sizeof(double));
            if (epitaxial_extend_y[i] == NULL)
            {
                printf("epitaxial_extend_y[i] Failed to allocate memory\n");
                exit(1);
            }
        }
    }

    for (int i = 0; i < net_numbs_2; i++)
    {

        epitaxial_extend_y[i][0] = block_net_p[count_row][1] - min_spac_y * EXTEND_MESH_y;
        epitaxial_extend_y[i][1] = block_net_p[count_row][1] + min_spac_y * EXTEND_MESH_y;
        i++;
        epitaxial_extend_y[i][0] = block_net_p[count_row][3] - min_spac_y * EXTEND_MESH_y;
        epitaxial_extend_y[i][1] = block_net_p[count_row][3] + min_spac_y * EXTEND_MESH_y;
        count_row++;
    }

    int count_column_y = 0;
    int* overlap_egments_y = NULL;
    overlap_egments_y = (int*)malloc(sizeof(int) * net_numbs_2);
    if (overlap_egments_y == NULL)
    {
        printf("Failed to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i < net_numbs_2; i++)
    {
        for (int j = i + 1; j < net_numbs_2; j++)
        {

            if ((epitaxial_extend_y[i][0] >= epitaxial_extend_y[j][0]) && (epitaxial_extend_y[i][1] <= epitaxial_extend_y[j][1]))
            {
                overlap_egments_y[count_column_y] = i;
                count_column_y++;
                break;
            }
            if ((epitaxial_extend_y[i][0] <= epitaxial_extend_y[j][0]) && (epitaxial_extend_y[i][1] >= epitaxial_extend_y[j][1]))
            {
                overlap_egments_y[count_column_y] = i;
                epitaxial_extend_y[j][0] = epitaxial_extend_y[i][0];
                epitaxial_extend_y[j][1] = epitaxial_extend_y[i][1];
                count_column_y++;
                break;
            }
            if ((epitaxial_extend_y[i][0] <= epitaxial_extend_y[j][0]) && (epitaxial_extend_y[j][0] <= epitaxial_extend_y[i][1]) && (epitaxial_extend_y[i][1] <= epitaxial_extend_y[j][1]))
            {
                overlap_egments_y[count_column_y] = i;
                count_column_y++;
                epitaxial_extend_y[j][0] = epitaxial_extend_y[i][0];

                break;
            }
            if ((epitaxial_extend_y[j][0] <= epitaxial_extend_y[i][0]) && (epitaxial_extend_y[i][0] <= epitaxial_extend_y[j][1]) && (epitaxial_extend_y[j][1] <= epitaxial_extend_y[i][1]))
            {
                overlap_egments_y[count_column_y] = i;
                count_column_y++;
                epitaxial_extend_y[j][1] = epitaxial_extend_y[i][1];

                break;
            }
        }
    }

    double* epitaxial_extend_del_y = (double*)malloc(sizeof(double) * net_numbs_2);
    if (epitaxial_extend_del_y == NULL)
    {
        printf("epitaxial_extend_del_y  Failed to allocate memory\n");
        exit(1);
    }
    int temp = 0;
    for (int i = 0; i < net_numbs_2; i++)
    {

        for (int j = 0; j < count_column_y + 1; j++)
        {
            if (i == overlap_egments_y[j])
            {
                temp = 1;
            }
        }
        if (temp == 0)
        {

            epitaxial_extend_del_y[*count_epitaxial_del_cloumn_y] = epitaxial_extend_y[i][0];
            *count_epitaxial_del_cloumn_y = *count_epitaxial_del_cloumn_y + 1;
            epitaxial_extend_del_y[*count_epitaxial_del_cloumn_y] = epitaxial_extend_y[i][1];
            *count_epitaxial_del_cloumn_y = *count_epitaxial_del_cloumn_y + 1;
        }
        else
            temp = 0;
    }

    free(overlap_egments_y);
    free(epitaxial_extend_y);

    quick_sort(epitaxial_extend_del_y, 0, *count_epitaxial_del_cloumn_y - 1);

    if (epitaxial_extend_del_y[0] < boundary[1] + min_spac_y * round(EXTEND_MESH_y / 2))
    {
        epitaxial_extend_del_y[0] = boundary[1];
    }

    else
    {
        epitaxial_extend_del_y[0] = epitaxial_extend_del_y[0] + min_spac_y * round(EXTEND_MESH_y / 2);
    }

    if (epitaxial_extend_del_y[*count_epitaxial_del_cloumn_y - 1] > boundary[3] - min_spac_y * round(EXTEND_MESH_y / 2))
    {
        epitaxial_extend_del_y[*count_epitaxial_del_cloumn_y - 1] = boundary[3];
    }
    else
    {
        epitaxial_extend_del_y[*count_epitaxial_del_cloumn_y - 1] = epitaxial_extend_del_y[*count_epitaxial_del_cloumn_y - 1] - min_spac_y * round(EXTEND_MESH_y / 2);
    }

    for (int i = 1; i < *count_epitaxial_del_cloumn_y - 1; i++)
    {
        if ((i + 1) % 2 != 0)
        {

            epitaxial_extend_del_y[i] = epitaxial_extend_del_y[i] + min_spac_y * round(EXTEND_MESH_y / 2);
        }
        else
        {
            epitaxial_extend_del_y[i] = epitaxial_extend_del_y[i] - min_spac_y * round(EXTEND_MESH_y / 2);
        }
    }

    return epitaxial_extend_del_y;
}

void quick_sort(double* arr, int low, int high)
{
    if (low < high)
    {
        int i = low;
        int j = high;
        double k = arr[low];
        while (i < j)
        {
            while (i < j && arr[j] >= k)
            {
                j--;
            }

            if (i < j)
            {
                arr[i++] = arr[j];
            }

            while (i < j && arr[i] < k)
            {
                i++;
            }

            if (i < j)
            {
                arr[j--] = arr[i];
            }
        }

        arr[i] = k;

        quick_sort(arr, low, i - 1);
        quick_sort(arr, i + 1, high);
    }
}

// TODO
// re-store the TRPZ shapes
// obtain the minimum spaces between the max block and the min block
double** minimum_spacing(double net_input[][16], double** net_outline)
{
    int temp_t = net_input_row + numbs_trapezoid * 3;
    double temp_min;
    double minsize_net_x, minsize_net_y;
    //申请block_net内存
    double** block_net = (double**)malloc((net_input_row + numbs_trapezoid * 3) * sizeof(double));
    if (block_net == NULL)
    {
        printf("block_net  failed to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i < (net_input_row + numbs_trapezoid * 3); i++)
    {
        block_net[i] = (double*)malloc(4 * sizeof(double));
        if (block_net[i] == NULL)
        {
            printf("block_net  failed to allocate memory\n");
        }
    }

    //block_net就是将net_input展开，行数变成net_input_row + numbs_trapezoid * 3
    int count_block_row = 0;
    for (int i = 0; i < net_input_row; i++)
    {
        if (net_input[i][4] == net_input[i][6])
        {

            block_net[count_block_row][0] = net_input[i][0];
            block_net[count_block_row][1] = net_input[i][1];
            block_net[count_block_row][2] = net_input[i][2];
            block_net[count_block_row][3] = net_input[i][3];
            count_block_row++;
        }
        else
        {
            for (int j = 0; j < 4; j++)
            {
                block_net[count_block_row][0] = net_input[i][0 + j * 4];
                block_net[count_block_row][1] = net_input[i][1 + j * 4];
                block_net[count_block_row][2] = net_input[i][2 + j * 4];
                block_net[count_block_row][3] = net_input[i][3 + j * 4];

                count_block_row++;
            }
        }
    }

    minsize_net_x = block_net[0][2] - block_net[0][0];
    minsize_net_y = block_net[0][3] - block_net[0][1];

    //求梯形最大块的一边与最小块一边间的距离
    for (int i = 0; i < net_input_row; i++)
    {
        if (net_outline[i][5] == TRAPEZOID_COLUMN)
        {
            temp_min = net_outline[i][4] - net_outline[i][0];
            if (temp_min < minsize_net_x)
            {
                minsize_net_x = temp_min;
            }
        }
        if (net_outline[i][5] == TRAPEZOID_ROW)
        {
            temp_min = net_outline[i][4] - net_outline[i][1];
            if (temp_min < minsize_net_y)
            {
                minsize_net_y = temp_min;
            }
        }
    }

    // 所有块中的最小边长比较
    for (int i = 1; i < temp_t; i++)
    {
        temp_min = block_net[i][3] - block_net[i][1];
        if (temp_min < minsize_net_y)
        {
            minsize_net_y = temp_min;
        }
    }

    for (int i = 1; i < temp_t; i++)
    {
        temp_min = block_net[i][2] - block_net[i][0];
        if (temp_min < minsize_net_x)
        {
            minsize_net_x = temp_min;
        }
    }

    min_spac_x = minsize_net_x / DIV_NUMS;
    min_spac_y = minsize_net_y / DIV_NUMS;

    return block_net;
}

// TODO
// re-store the shapes
// these codes may be useless
double** storage_net(double(*net_input)[16])
{
    //申请net_outline的内存
    double** net_outline = (double**)malloc(net_input_row * sizeof(double));
    if (net_outline == NULL)
    {
        printf("net_outline  failed to allocate memory\n");
        exit(1);
    }
    else
    {
        for (int i = 0; i < net_input_row; ++i)
        {
            net_outline[i] = (double*)malloc(NET_COLUMN * sizeof(double));
            if (net_outline[i] == NULL)
            {
                printf("net_outline  failed to allocate memory\n");
                exit(1);
            }
        }
    }

    for (int i = 0; i < net_input_row; i++)
    {
        //如果是矩形
        if (net_input[i][4] == net_input[i][6])
        {
            net_outline[i][0] = net_input[i][0];
            net_outline[i][1] = net_input[i][1];
            net_outline[i][2] = net_input[i][2];
            net_outline[i][3] = net_input[i][3];
            net_outline[i][5] = RECTANGLE;
        }
        //如果是梯形
        else
        {

            net_outline[i][0] = compare_min(net_input[i][0], net_input[i][4], net_input[i][8], net_input[i][12]);
            net_outline[i][1] = compare_min(net_input[i][1], net_input[i][5], net_input[i][9], net_input[i][13]);
            net_outline[i][2] = compare_max(net_input[i][2], net_input[i][6], net_input[i][10], net_input[i][14]);
            net_outline[i][3] = compare_max(net_input[i][3], net_input[i][7], net_input[i][11], net_input[i][15]);
            if ((net_input[i][3] == net_input[i][5]) || (net_input[i][3] == net_input[i][9]) || (net_input[i][3] == net_input[i][13]) || (net_input[i][7] == net_input[i][1]) || (net_input[i][7] == net_input[i][9]) || (net_input[i][7] == net_input[i][13]))
            {
                //这块有什么用？？？
                net_outline[i][4] = compare_max(net_input[i][0], net_input[i][4], net_input[i][8], net_input[i][12]);
                net_outline[i][5] = TRAPEZOID_COLUMN;
            }
            else
            {
                net_outline[i][4] = compare_max(net_input[i][1], net_input[i][5], net_input[i][9], net_input[i][13]);
                net_outline[i][5] = TRAPEZOID_ROW;
            }
            //可以通过解析得到
            numbs_trapezoid = numbs_trapezoid + 1;
        }
    }

    return net_outline;
}

double compare_min(double a1, double a2, double a3, double a4)
{
    double temp;
    if (a1 < a2)
        if (a1 < a3)
            if (a1 < a4)
                temp = a1;
            else
                temp = a4;
        else if (a3 < a4)
            temp = a3;
        else
            temp = a4;
    else if (a2 < a3)
        if (a2 < a4)
            temp = a2;
        else
            temp = a4;
    else if (a3 < a4)
        temp = a3;
    else
        temp = a4;

    return temp;
}

double compare_max(double a1, double a2, double a3, double a4)
{
    double temp;
    if (a1 > a2)
        if (a1 > a3)
            if (a1 > a4)
                temp = a1;
            else
                temp = a4;
        else if (a3 > a4)
            temp = a3;
        else
            temp = a4;
    else if (a2 > a3)
        if (a2 > a4)
            temp = a2;
        else
            temp = a4;
    else if (a3 > a4)
        temp = a3;
    else
        temp = a4;

    return temp;
}

int ceil_numbers_x(double a)
{
    int y;
    int x;
    x = a;

    if (x + 0.1 < a)
    {

        y = a;
        y = a + 1;
    }
    else
    {
        y = x;
    }

    return y;
}

