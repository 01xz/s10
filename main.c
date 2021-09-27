#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "cholmod.h"

//float conductor_boundary(float * boundary);
int main(void) {
    // original data
    float boundary[4] = { -10.0000,0.0000,10.0000,9.9000 };
    float net0[4]     = { -0.0160,0.4800,0.0160,0.5509 };
    float net1[4]     = { -0.0800,0.4800,0.0480,0.5500 };
    float net2[4]     = { -0.0480,0.4800,0.0800,0.5500 };

    float dx          = 0.004;
    float dy          = 0.01;
    float net0_right  = 0;
    float net1_right  = 0; 
    float net2_right  = 0;
    float net0_left   = 0; 
    float net1_left   = 0;
    float net2_left   = 0;
    float net0_top    = 0;
    float net1_top    = 0; 
    float net2_top    = 0;
    float net0_bottom = 0;
    float net1_bottom = 0;
    float net2_bottom = 0;

    //强制转化成int类型后出现了截断误差？？如何解决
    int number_mesh_x    = ceil((boundary[2] - boundary[0])/0.004);
    int number_mesh_y    = ceil((boundary[3] - boundary[1])/0.01);
    int number_element_A = (number_mesh_x + 1) * (number_mesh_y + 1);

    int* pt_b = (int*)malloc(number_element_A * sizeof(int));

    float** pt_column;

    pt_column = (float**)malloc(sizeof(float*) * number_element_A);
    for (int i; i < 5; i++) {
        pt_column = (float*)malloc(sizeof(float) * 5);
    }

    float coefficient_y  = dy * dy;
    float coefficient_x  = dx * dx;
    float coefficient_xy = -2*(dx * dy);

    net0_right  = (net0[0] - boundary[0]) / dx;
    net0_left   = (net0[2] - boundary[0]) / dx;
    net1_right  = (net1[0] - boundary[0]) / dx;
    net1_left   = (net1[2] - boundary[0]) / dx;
    net2_right  = (net2[0] - boundary[0]) / dx;
    net2_left   = (net2[2] - boundary[0]) / dx;
    net0_top    = (net0[3] - boundary[1]) / dy;
    net0_bottom = (net0[2] - boundary[1]) / dy;
    net1_top    = (net1[3] - boundary[1]) / dy;
    net1_bottom = (net1[2] - boundary[1]) / dy;
    net2_top    = (net2[3] - boundary[1]) / dy;
    net2_bottom = (net2[2] - boundary[1]) / dy;

    // filling matrix A and vector b
    for (int y_i = 0; y_i <= number_mesh_y; y_i++) {
        for (int x_i = 0; x_i <= number_mesh_x; x_i++) {
            if ((x_i == 0) | (y_i == 0) | (x_i == number_mesh_x) | (y_i == number_mesh_y) | \
                ((x_i >= net0_right) && (x_i <= net0_left) && (y_i >= net0_bottom) && (y_i <= net0_top)) | \
                ((x_i >= net1_right) && (x_i <= net1_left) && (y_i >= net1_bottom) && (y_i <= net1_top)) | \
                ((x_i >= net2_right) && (x_i <= net2_left) && (y_i >= net2_bottom) && (y_i <= net2_top))) {
                if ((x_i >= net0_right) && (x_i <= net0_left) && (y_i <= net0_top) && (y_i >= net0_bottom)) {
                    pt_b[y_i * number_mesh_x + x_i] = 1;
                }
                pt_column[y_i * number_mesh_x + x_i][2]         = 1;
                pt_column[y_i * number_mesh_x + x_i + 1][3]     = 0;
                pt_column[y_i * number_mesh_x  + x_i - 1][1]    = 0;
                pt_column[(y_i + 1) * number_mesh_x  + x_i][4]  = 0;
                pt_column[(y_i - 1) * number_mesh_x + x_i][0]   = 0;
            } else {
                pt_b[y_i * number_mesh_x + x_i] = 0;
            }
            pt_column[y_i * number_mesh_x + x_i][2]       = coefficient_xy;
            pt_column[y_i * number_mesh_x + x_i + 1][3]   = coefficient_y;
            pt_column[y_i * number_mesh_x + x_i - 1][1]   = coefficient_y;
            pt_column[(y_i + 1) * number_mesh_x + x_i][4] = coefficient_x;
            pt_column[(y_i - 1) * number_mesh_x + x_i][0] = coefficient_x;
        }
    }
    return 0;
}
