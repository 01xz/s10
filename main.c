#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#define FILE_BUFFER_LENGTH 300000
#pragma warning(push)
#pragma warning(disable:6385)
#pragma warning(disable:6386)

#define RECTANGLE           3
#define TRAPEZOID_ROW       2
#define TRAPEZOID_COLUMN    1
#define DIV_NUMS            30
double EXTEND_MESH_x = 15;
double EXTEND_MESH_y = 15;
const float  PERCENTAGE_DIV_0_7 = 0.7;
#define PERCENTAGE_INCREASING 1.1
double boundary[4] = {-10,0,10,9.9};
//double dielectric_constant = 3.9;
//double number_of_net = 2;
double net_input_row = 9;
int trapezoid = 0;
int cnt_x = 0;
int cnt_y = 0;

//Declare functions
double compare_min(double a1, double a2, double a3, double a4);
double compare_max(double a1, double a2, double a3, double a4);
//int round_numbers(double x);
//int cpy = aa(double net_input, 2, &minsize_net_x_direction, &minsize_net_x_direction);
//void minimum_spacing(double net_input[][16], int net_input_row, double* minsize_net_x_direction, double* minsize_net_y_direction);
int** storage_net(double(*net_inpu)[16]);
int** minimum_spacing(double net_input[][16], double** net_outline,double* min_spacing_net_x, double* min_spacing_net_y);
int extend_net_coordinate_x(double** block_net_p, double* min_spacing_net_x,int *count_epitaxial_del_cloumn_x);
int extend_net_coordinate_y(double** block_net_p, double* min_spacing_net_y,int* count_epitaxial_del_cloumn_y);
void quick_sort(double arr[],int start, int end);
int partition(int arr[], int low, int high);
int  round_numbers(double x);
int ceil_numbers_x(double a);
void mesh_divide_x_p(double* ext_coor_x_p, double** net, double** block, double* min_spac_x, int* cnt_clmn_x, double* ext_block_coordinate, double* ext_net_coordinate);
void mesh_divide_y_p(double* extend_net_coordinate_y_p, double** net, double** block, double* min_spacing_net_y, int* count_epitaxial_del_cloumn_y, double* ext_block_coordinate, double* ext_net_coordinate);
double net_input[9][16] = { {-0.017, 0.072,  0.017, 0.139,  0,     0,      0,     0,      0,     0,      0,     0,      0,     0,      0,     0    },
                            {-0.114, 0.072, -0.078, 0.139,  0,     0,      0,     0,      0,     0,      0,     0,      0,     0,      0,     0    },
                            { 0.078, 0.072,  0.114, 0.139,  0,     0,      0,     0,      0,     0,      0,     0,      0,     0,      0,     0    },
                            {-0.155, 0.105, -0.133, 0.128, -0.156, 0.128, -0.132, 0.152, -0.157, 0.152, -0.131, 0.176, -0.158, 0.176, -0.13,  0.2  },
                            { 0.133, 0.105,  0.155, 0.128,  0.132, 0.128,  0.156, 0.152,  0.131, 0.152,  0.157, 0.176,  0.13,  0.176,  0.158, 0.2  },
                            {-0.059, 0.105, -0.037, 0.128, -0.06,  0.128, -0.036, 0.152, -0.061, 0.152, -0.035, 0.176, -0.062, 0.176, -0.034, 0.2  },
                            { 0.037, 0.105,  0.059, 0.128,  0.036, 0.128,  0.060, 0.152,  0.035, 0.152,  0.061, 0.176,  0.034, 0.176,  0.062, 0.200},
                            {-0.061, 0.247, -0.024, 0.261, -0.062, 0.261, -0.023, 0.276, -0.063, 0.276, -0.022, 0.290, -0.064, 0.290, -0.021, 0.304},
                            { 0.024, 0.247,  0.061, 0.261,  0.023, 0.261,  0.062, 0.276,  0.022, 0.276,  0.063, 0.290,  0.021, 0.290,  0.064, 0.304} };



int main(void)
{
    double* min_spacing_net_x = NULL;
    double* min_spacing_net_y = NULL;

    int * count_epitaxial_del_cloumn_x = 0;
    int * count_epitaxial_del_cloumn_y = 0;


    



   double** storage_net_p = storage_net(net_input);
   
   double** block_net_p   = minimum_spacing(net_input, storage_net_p,&min_spacing_net_x, &min_spacing_net_y);


   //计算外扩的net在网格上的坐标
   double** ext_net_coordinate = (double**)malloc(net_input_row * sizeof(double));
   if (ext_net_coordinate == NULL)
   {
       printf("Failed to allocate memory\n");

   }
   for (int i = 0; i < net_input_row; ++i)
   {
       ext_net_coordinate[i] = (double*)malloc(4 * sizeof(double));
       if (ext_net_coordinate == NULL)
       {
           printf("Failed to allocate memory\n");
       }
   }

   


   //计算外扩的block在网格上的坐标
   double** ext_block_coordinate = (double**)malloc(cnt_x * sizeof(double));
   if (ext_block_coordinate == NULL)
   {
       printf("Failed to allocate memory\n");

   }
   for (int i = 0; i < cnt_x; ++i)
   {
       ext_block_coordinate[i] = (double*)malloc(4 * sizeof(double));
       if (ext_net_coordinate == NULL)
       {
           printf("Failed to allocate memory\n");

       }

   }










   double* extend_net_coordinate_x_p = extend_net_coordinate_x(block_net_p, &min_spacing_net_x,&count_epitaxial_del_cloumn_x);
   double* extend_net_coordinate_y_p = extend_net_coordinate_y(block_net_p, &min_spacing_net_y,&count_epitaxial_del_cloumn_y);
   
   mesh_divide_x_p(extend_net_coordinate_x_p, storage_net_p,block_net_p, &min_spacing_net_x,& count_epitaxial_del_cloumn_x, ext_block_coordinate,ext_net_coordinate);
   mesh_divide_y_p(extend_net_coordinate_y_p, storage_net_p, block_net_p, &min_spacing_net_y,&count_epitaxial_del_cloumn_y,ext_block_coordinate, ext_net_coordinate);
   



 
   free(block_net_p);
   free(storage_net_p);
   free(extend_net_coordinate_x_p);
   free(extend_net_coordinate_y_p);

   printf("%d", 1);
   return 0;

}





void mesh_divide_x_p(double* ext_coor_x_p, double** net, double** block, double* min_spac_x, int* cnt_clmn_x, double **ext_block_coordinate, double **ext_net_coordinate)
{

   

    //计算右边0.7部分
    int long_cnt_right = 5;
    double long_right = pow(1.1,5.0);
    double long_mul_right = pow(1.1, 5.0);
    double div_numbs_right = (boundary[2] - ext_coor_x_p[*cnt_clmn_x - 1])* PERCENTAGE_DIV_0_7 /(*min_spac_x);
    
    while (((long_right - 1) / (1.1 - 1)) < div_numbs_right)
    {
        long_right = long_right * long_mul_right;
        long_cnt_right = long_cnt_right + 5;
    }
    //就按右边0.3部分
    int long_cnt_right_03 = round_numbers((div_numbs_right * 3 / 7) / pow(1.1, long_cnt_right - 5));
   

    //计算左边0.7部分
    int long_cnt_left = 5;
    double long_left = pow(1.1, 5);
    
    double long_mul_left = pow(1.1, 5);
    double div_numbs_left = (ext_coor_x_p[1] - boundary[0])* PERCENTAGE_DIV_0_7 / (*min_spac_x);
    
    while ((long_left - 1) / (1.1 - 1) < div_numbs_left)
    {
        long_left = long_left * long_mul_left;
        long_cnt_left = long_cnt_left + 5;


    }

    int long_cnt_left_03 = round_numbers((div_numbs_left * 3 / 7) / pow(1.1, long_cnt_left - 5));
    int long_cnt_x = long_cnt_right + long_cnt_left + long_cnt_right_03 + long_cnt_left_03;

  

    int* unif_div_numbs_x = (double*)malloc((*cnt_clmn_x - 2) * sizeof(double*));  //uniform_filigree_divide_numbs_x 
    double* unif_div_inter_x = (double*)malloc((*cnt_clmn_x - 2) * sizeof(double*));//uniform_filigree_divide_interval_x


    for (int i = 0; i < *cnt_clmn_x - 1; i = i + 2)
    {

        unif_div_numbs_x[i / 2] = round_numbers((ext_coor_x_p[i + 1] - ext_coor_x_p[i]) / (*min_spac_x));
        unif_div_inter_x[i / 2] = (ext_coor_x_p[i + 1] - ext_coor_x_p[i]) / round_numbers((ext_coor_x_p[i + 1] - ext_coor_x_p[i]) / (*min_spac_x));
       
        long_cnt_x = long_cnt_x + unif_div_numbs_x[i/2];
     
    }


    int long_cnt_middle = 5;
    double long_middle= pow(1.1, 5);
    double long_mul_middle = pow(1.1, 5);
    double div_numbs_middle = 0;

    for (int i = 1; i < *cnt_clmn_x - 2; i = i + 2)
    {
        long_cnt_middle = 5;
        long_middle = pow(1.1, 5);
        long_mul_middle = pow(1.1, 5);
        div_numbs_middle = ((ext_coor_x_p[i + 1] - ext_coor_x_p[i]) * PERCENTAGE_DIV_0_7 / (*min_spac_x)) / 2;
        while ((long_middle - 1) / (1.1 - 1) < div_numbs_middle)
        {
            long_middle = long_middle * long_mul_middle;
            long_cnt_middle = long_cnt_middle + 5;

        }

        
        long_cnt_x = long_cnt_x+2 * long_cnt_middle + round_numbers((2 * div_numbs_middle * 3 / 7) / pow(1.1, long_cnt_middle - 5));

    }



    //这一部分的代码是对x轴方向靠近boundary(2)的最后一块区域进行划分，先非均匀划分70%的部分，剩下的均匀划分
    double long_dived_right = 0;//已经划分了的区域长度
    double* non_unif_div_numbs_right = (double*)malloc(long_cnt_right * sizeof(double));  //uniform_filigree_divide_numbs_x 
    //double* non_unif_div_numbs_left = (double*)malloc(long_cnt_left * sizeof(double*));  //uniform_filigree_divide_numbs_x 
    int cnt_non_unif_inter_right = 0;
    non_unif_div_numbs_right[cnt_non_unif_inter_right] = unif_div_inter_x[ceil_numbers_x(*cnt_clmn_x - 2)/2];
    //non_unif_div_numbs_right[cnt_non_unif_inter_right] = unif_div_inter_x(ceil_numbers_x(* cnt_clmn_x - 2)/2);
    int unif_div_numbs_right = 0;
    double long_unif_dived_right = 0;
    if(boundary[2] - (EXTEND_MESH_x/2)*(*min_spac_x) > ext_coor_x_p[* cnt_clmn_x - 1])
    {
        //非均匀划分
        while(long_dived_right < PERCENTAGE_DIV_0_7*boundary[2] - ext_coor_x_p[* cnt_clmn_x - 1])
        {
            long_dived_right = long_dived_right + non_unif_div_numbs_right[cnt_non_unif_inter_right];
            non_unif_div_numbs_right[cnt_non_unif_inter_right + 1] = non_unif_div_numbs_right[cnt_non_unif_inter_right] * 1.1;
        
            cnt_non_unif_inter_right++;

        }
        //均匀划分
        unif_div_numbs_right = round_numbers((boundary[2] - ext_coor_x_p[* cnt_clmn_x - 1]-long_dived_right)/ non_unif_div_numbs_right[cnt_non_unif_inter_right - 1]);
        long_unif_dived_right = (boundary[2] - ext_coor_x_p[*cnt_clmn_x - 1] - long_dived_right) / unif_div_numbs_right;

    }
    //free(non_unif_div_numbs_right);
    ////这一部分的代码是对x轴方向靠近boundary(0)的最后一块区域进行划分，先非均匀划分70%的部分，剩下的均匀划分
    double long_dived_left = 0;//已经划分了的区域长度
    
    double* non_unif_div_numbs_left = (double*)malloc(long_cnt_left * sizeof(double));  //uniform_filigree_divide_numbs_x 
    int cnt_non_unif_inter_left = 0;
    //free(non_unif_div_numbs_left);
    
    non_unif_div_numbs_left[cnt_non_unif_inter_left] = unif_div_inter_x[0];
    
    int unif_div_numbs_left = 0;
    double long_unif_dived_left = 0;

    double* inter_x = (double*)malloc(long_cnt_x * sizeof(double*));
    double* dividing_line_x = (double*)malloc(*cnt_clmn_x * sizeof(double));
    int cnt_dividing_line_x = 0;
    
    //int cnt_x = 0;
    if(boundary[0] + (EXTEND_MESH_x/2)*(*min_spac_x) < ext_coor_x_p[0])
    {
        //非均匀划分
        while(long_dived_left < PERCENTAGE_DIV_0_7*(ext_coor_x_p[0] - boundary[0]))
        {
            long_dived_left = long_dived_left + non_unif_div_numbs_left[cnt_non_unif_inter_left];
       
            non_unif_div_numbs_left[cnt_non_unif_inter_left + 1] = non_unif_div_numbs_left[cnt_non_unif_inter_left] * 1.1;
            cnt_non_unif_inter_left++;

        }
        //均匀划分
        unif_div_numbs_left = round((ext_coor_x_p[0] - boundary[0] - long_dived_left)/ non_unif_div_numbs_left[cnt_non_unif_inter_left - 1]);
        long_unif_dived_left = (ext_coor_x_p[0] - boundary[0] - long_dived_left) / unif_div_numbs_left;





        for (int i = 0; i < unif_div_numbs_left; i++)
        {
            inter_x[cnt_x] = long_unif_dived_left;
      
            cnt_x++;

        }

        for (int i = 0; i < cnt_non_unif_inter_left; i++)
        {
            
           inter_x[cnt_x] = non_unif_div_numbs_left[cnt_non_unif_inter_left - i -1];
          
            cnt_x++;
            
        }

        // %下面两行专门来记录细划分区边界的格点坐标
        dividing_line_x[cnt_dividing_line_x] = cnt_x;
        cnt_dividing_line_x ++;
    }
    
    //下面是个循环，用来填充inter_x中间的部分
   double long_dived = 0;//已经划分了的区域长度


   double* array_x = (double*)malloc(long_cnt_x * sizeof(double*));  //uniform_filigree_divide_numbs_x 
   int cnt_temp = 0;
   int unif_div_numbs_middle;
   double long_unif_dived_middle;

   for (int i = 1; i < *cnt_clmn_x - 1; i = i + 2)
   {
       // 确定非均匀的划分区间
       long_dived = 0;
       cnt_temp = 0;
       array_x[0] = (*min_spac_x);

       // 两边70%区域同时进行非均匀划分
       while ((ext_coor_x_p[i + 1] - ext_coor_x_p[i]) * PERCENTAGE_DIV_0_7 > 2 * long_dived)
       {

           long_dived = long_dived + array_x[cnt_temp];
           array_x[cnt_temp + 1] = array_x[cnt_temp] * 1.1;
           cnt_temp++;

       }

       ////中间部分粗均匀划分
       unif_div_numbs_middle = round_numbers((ext_coor_x_p[i + 1] - ext_coor_x_p[i] - 2 * long_dived) / array_x[cnt_temp - 1]);
       long_unif_dived_middle = (ext_coor_x_p[i + 1] - ext_coor_x_p[i] - 2 * long_dived) / unif_div_numbs_middle;

       //填装中间细均匀部分
       for (int j = 0; j < unif_div_numbs_x[(i - 1)/ 2]; j++)
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
       for (int j = 0; j < cnt_temp ; j++)
       {
           inter_x[cnt_x] = array_x[cnt_temp - j -1];
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

   // %下面两行专门来记录细划分区边界的格点坐标
   dividing_line_x[cnt_dividing_line_x] = cnt_x;
   cnt_dividing_line_x++;

   //填充整个inter_x的靠近boundary(2)部分
   if (boundary[2] - (EXTEND_MESH_x / 2) * (*min_spac_x) > ext_coor_x_p[*cnt_clmn_x - 1])
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



   //计算外扩的net在网格上的坐标
   double** ext_net_coordinate = (double**)malloc(net_input_row * sizeof(double));
   if (ext_net_coordinate == NULL)
   {
       printf("Failed to allocate memory\n");
    
   }
   for (int i = 0; i < net_input_row; ++i)
   {
       ext_net_coordinate[i] = (double*)malloc(4 * sizeof(double));
       if (ext_net_coordinate == NULL)
       {
           printf("Failed to allocate memory\n");  
       }
   }
   
   int i = 0;
   int j = 0;
   for ( i = 0; i < net_input_row ; i++)
   {
       for ( j = 0; j < *cnt_clmn_x -1; j = j + 2)
       {
           if ((net[i][0] > ext_coor_x_p[j]) && (net[i][0] < ext_coor_x_p[j + 1]))
           {
               ext_net_coordinate[i][0] = (net[i][0] - ext_coor_x_p[j]) / unif_div_inter_x[j/2] + dividing_line_x[j];
               
           }
           if ((net[i][2] > ext_coor_x_p[j]) && (net[i][0] < ext_coor_x_p[j + 1]))
           {
               ext_net_coordinate[i][2] = (net[i][2] - ext_coor_x_p[j]) / unif_div_inter_x[j/2] + dividing_line_x[j];
              
           }
            
       }
   }



   //计算外扩的block在网格上的坐标
   double** ext_block_coordinate = (double**)malloc(cnt_x * sizeof(double));
   if (ext_block_coordinate == NULL)
   {
       printf("Failed to allocate memory\n");

   }
   for (int i = 0; i < cnt_x; ++i)
   {
       ext_block_coordinate[i] = (double*)malloc(4 * sizeof(double));
       if (ext_net_coordinate == NULL)
       {
           printf("Failed to allocate memory\n");

       }

   }


   for (int i = 0; i < 27; i++)
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

   

   




    
    printf("d", 1);



}




void mesh_divide_y_p(double* ext_coor_y_p, double** net, double** block, double* min_spac_y, int* cnt_clmn_y, double** ext_block_coordinate, double** ext_net_coordinate)
{
    //计算上边的0.7部分
    int long_cnt_high = 5;
    double long_high = pow(1.1, 5.0);
    double long_mul_high = pow(1.1, 5.0);
    double div_numbs_high = (boundary[3] - ext_coor_y_p[*cnt_clmn_y - 1]) * PERCENTAGE_DIV_0_7 / (*min_spac_y);

    while (((long_high - 1) / (1.1 - 1)) < div_numbs_high)
    {
        long_high = long_high * long_mul_high;
        long_cnt_high = long_cnt_high + 5;
    }
    //就按上边0.3部分
    int long_cnt_high_03 = round_numbers((div_numbs_high * 3 / 7) / pow(1.1, long_cnt_high - 5));


    int long_cnt_low = 5;
    double long_low = pow(1.1, 5);

    double long_mul_low = pow(1.1, 5);
    double div_numbs_low = (ext_coor_y_p[1] - boundary[1]) * PERCENTAGE_DIV_0_7 / (*min_spac_y);

    while ((long_low - 1) / (1.1 - 1) < div_numbs_low)
    {
        long_low = long_low * long_mul_low;
        long_cnt_low = long_cnt_low + 5;


    }

    int long_cnt_low_03 = round_numbers((div_numbs_low * 3 / 7) / pow(1.1, long_cnt_low - 5));
    int long_cnt_y = long_cnt_high + long_cnt_low + long_cnt_high_03 + long_cnt_low_03;





    //在这些线段里进行细密的划分，划分间隔与之前定的最小划分间隔相近，uniform_filigree_divide_numbs_x是每个的均匀细划分区间隔数，uniform_filigree_divide_interval_x是每个均�?细划分区的间�?

    int * unif_div_numbs_y = (double*)malloc((*cnt_clmn_y - 2) / 2 * sizeof(double*)); //uniform_filigree_divide_numbs_y
    double* unif_div_inter_y = (double*)malloc((*cnt_clmn_y - 2) / 2 * sizeof(double*));//uniform_filigree_divide_interval_y

    double* tt = (double*)malloc((*cnt_clmn_y - 2) * sizeof(double*));

    for (int i = 0; i < *cnt_clmn_y - 1; i = i + 2)
    {
        unif_div_numbs_y[i / 2] = round((ext_coor_y_p[i + 1] - ext_coor_y_p[i]) / (*min_spac_y));
        unif_div_inter_y[i / 2] = (ext_coor_y_p[i + 1] - ext_coor_y_p[i]) / round((ext_coor_y_p[i + 1] - ext_coor_y_p[i]) / (*min_spac_y));
        long_cnt_y = long_cnt_y + unif_div_numbs_y[i / 2];

    }

    int long_cnt_middle = 5;
    double long_middle = pow(1.1, 5);
    double long_mul_middle = pow(1.1, 5);
    double div_numbs_middle = 0;

    for (int i = 1; i < *cnt_clmn_y - 1; i = i + 2)
    {
        long_cnt_middle = 5;
        long_middle = pow(1.1, 5);
        long_mul_middle = pow(1.1, 5);
        div_numbs_middle = ((ext_coor_y_p[i + 1] - ext_coor_y_p[i]) * PERCENTAGE_DIV_0_7 / (*min_spac_y)) / 2;
        while ((long_middle - 1) / (1.1 - 1) < div_numbs_middle)
        {
            long_middle = long_middle * long_mul_middle;
            long_cnt_middle = long_cnt_middle + 5;

        }


        long_cnt_y = long_cnt_y + 2 * long_cnt_middle + round_numbers((2 * div_numbs_middle * 3 / 7) / pow(1.1, long_cnt_middle - 5));

    }



    //这一部分的代码是对x轴方向靠近boundary(3)的最后一块区域进行划分，先非均匀划分70%的部分，剩下的均匀划分
    double long_dived_high = 0;//已经划分了的区域长度
    double* non_unif_div_numbs_high = (double*)malloc(long_cnt_high * sizeof(double));  //uniform_filigree_divide_numbs_x 
    //double* non_unif_div_numbs_low = (double*)malloc(long_cnt_low * sizeof(double*));  //uniform_filigree_divide_numbs_x 
    int cnt_non_unif_inter_high = 0;
    non_unif_div_numbs_high[cnt_non_unif_inter_high] = unif_div_inter_y[ceil_numbers_x(*cnt_clmn_y - 2) / 2];
    //non_unif_div_numbs_high[cnt_non_unif_inter_high] = unif_div_inter_y(ceil_numbers_x(* cnt_clmn_y - 2)/2);
    int unif_div_numbs_high = 0;
    double long_unif_dived_high = 0;
    if (boundary[3] - (EXTEND_MESH_x / 2) * (*min_spac_y) > ext_coor_y_p[*cnt_clmn_y - 1])
    {
        //非均匀划分
        while (long_dived_high < PERCENTAGE_DIV_0_7 * boundary[3] - ext_coor_y_p[*cnt_clmn_y - 1])
        {
            long_dived_high = long_dived_high + non_unif_div_numbs_high[cnt_non_unif_inter_high];
            non_unif_div_numbs_high[cnt_non_unif_inter_high + 1] = non_unif_div_numbs_high[cnt_non_unif_inter_high] * 1.1;

            cnt_non_unif_inter_high++;

        }
        //均匀划分
        unif_div_numbs_high = round_numbers((boundary[2] - ext_coor_y_p[*cnt_clmn_y - 1] - long_dived_high) / non_unif_div_numbs_high[cnt_non_unif_inter_high - 1]);
        long_unif_dived_high = (boundary[3] - ext_coor_y_p[*cnt_clmn_y - 1] - long_dived_high) / unif_div_numbs_high;

    }




    ////这一部分的代码是对x轴方向靠近boundary(0)的最后一块区域进行划分，先非均匀划分70%的部分，剩下的均匀划分
    double long_dived_low = 0;//已经划分了的区域长度

    double* non_unif_div_numbs_low = (double*)malloc(long_cnt_low * sizeof(double*));  //uniform_filigree_divide_numbs_x 
    int cnt_non_unif_inter_low = 0;


    non_unif_div_numbs_low[cnt_non_unif_inter_low] = unif_div_inter_y[0];

    int unif_div_numbs_low = 0;
    double long_unif_dived_low = 0;

    double* inter_y = (double*)malloc(long_cnt_y * sizeof(double*));
    double* dividing_line_y = (double*)malloc(*cnt_clmn_y * sizeof(double*));
    int cnt_dividing_line_y = 0;

    //int cnt_y = 0;
    if (boundary[1] + (EXTEND_MESH_x / 2) * (*min_spac_y) < ext_coor_y_p[0])
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


    double long_dived = 0;//已经划分了的区域长度


    double* array_y = (double*)malloc(long_cnt_y * sizeof(double*));  //uniform_filigree_divide_numbs_x 
    int cnt_temp = 0;
    int unif_div_numbs_middle;
    double long_unif_dived_middle;

    for (int i = 1; i < *cnt_clmn_y - 1; i = i + 2)
    {
        // 确定非均匀的划分区间
        long_dived = 0;
        cnt_temp = 0;
        array_y[0] = (*min_spac_y);

        // 两边70%区域同时进行非均匀划分
        while ((ext_coor_y_p[i + 1] - ext_coor_y_p[i]) * PERCENTAGE_DIV_0_7 > 2 * long_dived)
        {

            long_dived = long_dived + array_y[cnt_temp];
            array_y[cnt_temp + 1] = array_y[cnt_temp] * 1.1;
            cnt_temp++;

        }

        ////中间部分粗均匀划分
        unif_div_numbs_middle = round_numbers((ext_coor_y_p[i + 1] - ext_coor_y_p[i] - 2 * long_dived) / array_y[cnt_temp - 1]);
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




    //最后一均匀划分填充
    //可能有问题
    for (int j = 0; j < unif_div_numbs_y[(*cnt_clmn_y - 1) / 2]; j++)
    {
        inter_y[cnt_y] = unif_div_inter_y[(*cnt_clmn_y - 1) / 2];
        cnt_y++;

    }

    // %下面两行专门来记录细划分区边界的格点坐标
    dividing_line_y[cnt_dividing_line_y] = cnt_y;
    cnt_dividing_line_y++;

    //填充整个inter_y的靠近boundary(2)部分
    if (boundary[2] - (EXTEND_MESH_x / 2) * (*min_spac_y) > ext_coor_y_p[*cnt_clmn_y - 1])
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





    //计算外扩的net在网格上的坐标
   

  
    for (int i = 0; i < net_input_row; i++)
    {
        for (int j = 0; j < *cnt_clmn_y - 1; j = j + 2)
        {
            if ((net[i][1] > ext_coor_y_p[j]) && (net[i][1] < ext_coor_y_p[j + 1]))
            {
                ext_net_coordinate[i][1] = (net[i][1] - ext_coor_y_p[j]) / unif_div_inter_y[j / 2] + dividing_line_y[j];

            }
            if ((net[i][3] > ext_coor_y_p[j]) && (net[i][3] < ext_coor_y_p[j + 1]))
            {
                ext_net_coordinate[i][3] = (net[i][3] - ext_coor_y_p[j]) / unif_div_inter_y[j / 2] + dividing_line_y[j];

            }

        }
    }



    //计算外扩的block在网格上的坐标
   

    for (int i = 0; i < 27; i++)
    {
        for (int j = 0; j < *cnt_clmn_y - 1; j = j + 2)
        {
            if ((block[i][1] > ext_coor_y_p[j]) && (block[i][1] < ext_coor_y_p[j + 1]))
            {

                ext_block_coordinate[i][1] = (block[i][1] - ext_coor_y_p[j]) / unif_div_inter_y[j / 2] + dividing_line_y[j];
            }
            if ((block[i][2] > ext_coor_y_p[j]) && (block[i][0] < ext_coor_y_p[j + 1]))
            {
                ext_block_coordinate[i][3] = (block[i][3] - ext_coor_y_p[j]) / unif_div_inter_y[j / 2] + dividing_line_y[j];

            }
        }

    }







      printf("%d", 1);



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
    else {
        y = x;

    }


    return y;

}


int round_numbers(double x)
{
    int y;
    y = (int)(x + 0.5) > (int)x ? (int)x + 1 : (int)x;

    return y;

}







int extend_net_coordinate_x(double **block_net_p , double *min_spacing_net_x,int *count_epitaxial_del_cloumn_x)
{
    int net_numbs =  net_input_row + trapezoid * 3;
    int net_numbs_2 = 2 * (net_input_row + trapezoid * 3);

    int count_row = 0;
    double**  epitaxial_extend_x = (double**)malloc(net_numbs_2 * sizeof(double*));
    if (epitaxial_extend_x == NULL)
    {
        printf("Failed to allocate memory\n");
        exit(1);
    }
    for (int i = 0; i < net_numbs_2; ++i)
    {
        epitaxial_extend_x[i] = (double*)malloc(2 * sizeof(double));
        if (epitaxial_extend_x == NULL)
        {
            printf("Failed to allocate memory\n");
            exit(1);
        }

    }  


    for (int i = 0; i < net_numbs_2; i++)
    {

        epitaxial_extend_x[i][0] = block_net_p[count_row][0] - *min_spacing_net_x * EXTEND_MESH_x;
        epitaxial_extend_x[i][1] = block_net_p[count_row][0] + *min_spacing_net_x * EXTEND_MESH_x;

        i++;
        epitaxial_extend_x[i][0] = block_net_p[count_row][2] - *min_spacing_net_x * EXTEND_MESH_x;
        epitaxial_extend_x[i][1] = block_net_p[count_row][2] + *min_spacing_net_x * EXTEND_MESH_x;
        count_row ++;
    }


    int count_column_x = 0;
    double* overlap_egments_x = NULL;
    overlap_egments_x = (double*)malloc(sizeof(double) * net_numbs_2);
    if (overlap_egments_x == NULL)
    {
        printf("Failed to allocate memory\n");
        exit(1);
    }
    if (overlap_egments_x == NULL)
    {
        printf("Failed to allocate memory\n");
        exit(1);
    }
    for (int i = 0; i < net_numbs_2; i++)
    {
        for(int j = i+1 ; j< net_numbs_2;j++)
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
            if ((epitaxial_extend_x[i][0] <= epitaxial_extend_x[j][0]) && (epitaxial_extend_x[j][0] <= epitaxial_extend_x[i][1])&& (epitaxial_extend_x[i][1] <= epitaxial_extend_x[j][1]))
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
    





    double* epitaxial_extend_del_x = (double*)malloc(sizeof(double) * net_numbs_2);//��������С�Ƕ�����
  
    int temp =0; 
    for (int i = 0; i < net_numbs_2; i++)
    {
     
        for (int j = 0; j < count_column_x + 1; j++)
        {
            if (overlap_egments_x[j] == i)
            { 
                temp = 1;
                
            }
        }
        if(temp == 0)
        {
            
            epitaxial_extend_del_x[*count_epitaxial_del_cloumn_x] = epitaxial_extend_x[i][0];
            *count_epitaxial_del_cloumn_x = *count_epitaxial_del_cloumn_x + 1;
            epitaxial_extend_del_x[*count_epitaxial_del_cloumn_x] = epitaxial_extend_x[i][1];
            *count_epitaxial_del_cloumn_x = *count_epitaxial_del_cloumn_x + 1;
        }
        else temp = 0;

    }
    


   
    free(overlap_egments_x);
    free(epitaxial_extend_x);

    quick_sort(epitaxial_extend_del_x, 0, (*count_epitaxial_del_cloumn_x - 1));
    
   



    if (epitaxial_extend_del_x[0] < boundary[0] + *min_spacing_net_x * round_numbers(EXTEND_MESH_x / 2))
    {
        epitaxial_extend_del_x[0] = boundary[0];
    }
    else
    {
        epitaxial_extend_del_x[0] = epitaxial_extend_del_x[0] + *min_spacing_net_x * round_numbers(EXTEND_MESH_x / 2);
    }

    if (epitaxial_extend_del_x[*count_epitaxial_del_cloumn_x - 1] > boundary[2] - *min_spacing_net_x * round_numbers(EXTEND_MESH_x / 2))
    {
        epitaxial_extend_del_x[*count_epitaxial_del_cloumn_x - 1] = boundary[2];
    }
    else
    {
        epitaxial_extend_del_x[*count_epitaxial_del_cloumn_x - 1] = epitaxial_extend_del_x[*count_epitaxial_del_cloumn_x - 1] - (*min_spacing_net_x) * round_numbers(EXTEND_MESH_x / 2);
    }

    



    for (int i = 1; i < *count_epitaxial_del_cloumn_x - 1; i++)
    {
        if ((i + 1) % 2 != 0)
        {
            //
            epitaxial_extend_del_x[i] = epitaxial_extend_del_x[i] + *min_spacing_net_x * round_numbers(EXTEND_MESH_x / 2);
        }
        else
        {
            
            epitaxial_extend_del_x[i] = epitaxial_extend_del_x[i] - *min_spacing_net_x * round_numbers(EXTEND_MESH_x / 2);
        }

    }



    

    return epitaxial_extend_del_x;




}



int extend_net_coordinate_y(double** block_net_p, double* min_spacing_net_y, int* count_epitaxial_del_cloumn_y)
{
    int net_numbs = net_input_row + trapezoid * 3;
    int net_numbs_2 = 2 * (net_input_row + trapezoid * 3);
    int count_row = 0;

    double** epitaxial_extend_y = (double**)malloc(net_numbs_2 * sizeof(double*));
    if (epitaxial_extend_y == NULL)
    {
        printf("�ڴ���䲻�ɹ���\n");
        exit(1);
    }

    for (int i = 0; i < net_numbs_2; ++i)
    {
        epitaxial_extend_y[i] = (double*)malloc(2 * sizeof(double));
        if (epitaxial_extend_y == NULL)
        {
            printf("�ڴ���䲻�ɹ���\n");
            exit(1);
        }

    }


    for (int i = 0; i < net_numbs_2; i++)
    {


        epitaxial_extend_y[i][0] = block_net_p[count_row][1] - *min_spacing_net_y * EXTEND_MESH_y;
        epitaxial_extend_y[i][1] = block_net_p[count_row][1] + *min_spacing_net_y * EXTEND_MESH_y;
        i++;
        epitaxial_extend_y[i][0] = block_net_p[count_row][3] - *min_spacing_net_y * EXTEND_MESH_y;
        epitaxial_extend_y[i][1] = block_net_p[count_row][3] + *min_spacing_net_y * EXTEND_MESH_y;
        count_row++;
    }

    



    int count_column_y = 0;
    double* overlap_egments_y = NULL;
    overlap_egments_y = (double*)malloc(sizeof(double) * net_numbs_2);
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

    //count_column_x = 0;

    
    int temp = 0;
    for (int i = 0; i < net_numbs_2; i++)
    {
        // printf("�ڴ���䲻�ɹ���\n");
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
        else temp = 0;
    }


   


    free(overlap_egments_y);
    free(epitaxial_extend_y);

    quick_sort(epitaxial_extend_del_y, 0, *count_epitaxial_del_cloumn_y - 1);

  





    if (epitaxial_extend_del_y[0] < boundary[1] + *min_spacing_net_y * round_numbers(EXTEND_MESH_y / 2))
    {
        epitaxial_extend_del_y[0] = boundary[1];
    }

    else
    {
        epitaxial_extend_del_y[0] = epitaxial_extend_del_y[0] + *min_spacing_net_y * round_numbers(EXTEND_MESH_y / 2);
    }

    if (epitaxial_extend_del_y[*count_epitaxial_del_cloumn_y - 1] > boundary[3] - *min_spacing_net_y * round_numbers(EXTEND_MESH_y / 2))
    {
        epitaxial_extend_del_y[*count_epitaxial_del_cloumn_y - 1] = boundary[3];
    }
    else
    {
        epitaxial_extend_del_y[*count_epitaxial_del_cloumn_y - 1] = epitaxial_extend_del_y[*count_epitaxial_del_cloumn_y - 1] - *min_spacing_net_y * round_numbers(EXTEND_MESH_y / 2);
    }

    for (int i = 1; i < *count_epitaxial_del_cloumn_y - 1; i++)
    {
        if ((i + 1) % 2 != 0)
        {

            epitaxial_extend_del_y[i] = epitaxial_extend_del_y[i] + *min_spacing_net_y * round_numbers(EXTEND_MESH_y / 2);
            
        }
        else
        {
            epitaxial_extend_del_y[i] = epitaxial_extend_del_y[i] - *min_spacing_net_y * round_numbers(EXTEND_MESH_y / 2);

        }

    }




    return epitaxial_extend_del_y;



}














void quick_sort(double *arr, int low, int high)
{
    if (low < high)
    {
        int i = low;
        int j = high;
        double k = arr[low];
        while (i < j)
        {
            while (i < j && arr[j] >= k)     // ���������ҵ�һ��С��k����
            {
                j--;
            }

            if (i < j)
            {
                arr[i++] = arr[j];
            }

            while (i < j && arr[i] < k)      // ���������ҵ�һ�����ڵ���k����
            {
                i++;
            }

            if (i < j)
            {
                arr[j--] = arr[i];
            }
        }

        arr[i] = k;

        // �ݹ����
        quick_sort(arr, low, i - 1);     // ����k���
        quick_sort(arr, i + 1, high);    // ����k�ұ�
    }
}



int** minimum_spacing(double net_input[][16], double **net_outline,double* min_spacing_net_x, double* min_spacing_net_y)
{
    double** block_net = NULL;
    double temp_min;
    double minsize_net_x, minsize_net_y;

    printf("%f\n", net_outline[3][3]);
    block_net = (double**)malloc((net_input_row + trapezoid * 3 )* sizeof(double*));
    if (block_net == NULL) {
        printf("�ڴ���䲻�ɹ���\n");
        exit(1);
    }

    for (int i = 0; i < (net_input_row + trapezoid * 3); i++)
    {
        block_net[i] = (double*)malloc(4 * sizeof(double));
        if (block_net[i] == NULL)
        {
            printf("�ڴ���䲻�ɹ���\n");

        }
       // block_net[i] = 0;
    }

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
   
    double temp_t = net_input_row + trapezoid * 3;
   
    for (int i = 1; i < temp_t; i++)
    {
       
        temp_min = block_net[i][3] - block_net[i][1];
        if (temp_min < minsize_net_y)
        {
            minsize_net_y = temp_min;

        }


    }
 
    //�����ߵ��������
    for (int i = 1; i < net_input_row + trapezoid * 3; i++)
    {
        temp_min = block_net[i][2] - block_net[i][0];
        if (temp_min < minsize_net_x)
        {
            minsize_net_x= temp_min;

        }


    }
   
   
    *min_spacing_net_x = minsize_net_x / DIV_NUMS;
    *min_spacing_net_y = minsize_net_y / DIV_NUMS;



    return block_net;
  // free(block_net);



}




//�洢net������
int ** storage_net(double (*net_input)[16])
{

    double **net_outline = NULL;  //������¼����߽�����ε����߽�
    const int NET_COLUMN = 6;
    
    int i = 0; //������ѭ������

    //���붯̬�ڴ�
    net_outline = (double**)malloc(net_input_row * sizeof(double));
    if (net_outline == NULL) {
        printf("�ڴ���䲻�ɹ���\n");
    }
    else
    {
    for (i = 0; i < net_input_row; ++i)
    {
        net_outline[i] = (double*)malloc(NET_COLUMN * sizeof(double));
        if (net_outline[i] == NULL) 
        {
            printf("�ڴ���䲻�ɹ���\n");
        }
    }
    }
   
  
    for (i = 0; i < net_input_row; i++)
    {
        if (net_input[i][4] == net_input[i][ 6])
        {
           
            net_outline[i][0] = net_input[i][0];
            net_outline[i][1] = net_input[i][1];
            net_outline[i][2] = net_input[i][2];
            net_outline[i][3] = net_input[i][3];
            net_outline[i][5] = RECTANGLE;
      
           
            

        }
        else
        {
        
            net_outline[i][0] = compare_min(net_input[i][0], net_input[i][4], net_input[i][8],  net_input[i][12]);
            net_outline[i][1] = compare_min(net_input[i][1], net_input[i][5], net_input[i][9],  net_input[i][13]);
            net_outline[i][2] = compare_max(net_input[i][2], net_input[i][6], net_input[i][10], net_input[i][12]);
            net_outline[i][3] = compare_max(net_input[i][3], net_input[i][7], net_input[i][11], net_input[i][13]);
            if ((net_input[i][3] == net_input[i][5]) || (net_input[i][3] == net_input[i][9]) || (net_input[i][3] == net_input[i][13]) || (net_input[i][7] == net_input[i][1]) || (net_input[i][7] == net_input[i][9]) || (net_input[i][7] == net_input[i][13]))
            {
                net_outline[i][4] = compare_max(net_input[i][0], net_input[i][4], net_input[i][8], net_input[i][12]);
                net_outline[i][5] = TRAPEZOID_COLUMN;
            }
            else
            { 
                net_outline[i][4] = compare_max(net_input[i][1], net_input[i][5], net_input[i][9], net_input[i][13]);
                net_outline[i][5] = TRAPEZOID_ROW;
             
              
            }
           
            trapezoid = trapezoid + 1;
        }


    }
   

    return net_outline;
}

//����Сֵ
double compare_min(double a1, double a2, double a3, double a4)
{
    double temp;
    if (a1 < a2)
        if (a1 < a3)
            if (a1 < a4)
                temp = a1;
            else temp = a4;
        else
            if (a3 < a4)
                temp = a3;
            else temp = a4;
    else
        if (a2 < a3)
            if (a2 < a4)
                temp = a2;
            else temp = a4;
        else
            if (a3 < a4)
                temp = a3;
            else
                temp = a4;

    return temp;

}



//�����ֵ
double compare_max(double a1, double a2, double a3, double a4)
{
    double temp;
    if (a1 > a2)
        if (a1 > a3)
            if (a1 > a4)
                temp = a1;
            else temp = a4;
        else
            if (a3 > a4)
                temp = a3;
            else temp = a4;
    else
        if (a2 > a3)
            if (a2 > a4)
                temp = a2;
            else temp = a4;
        else
            if (a3 > a4)
                temp = a3;
            else
                temp = a4;

    return temp;

}





















