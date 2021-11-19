#include "stdio.h"
#include "math.h"
#include "stdlib.h"
#define DIV 2
#define RANGE3 192.0
#define RANGE2 128.0
#define RANGE1 32.0
#define RANGE0 7.0

struct layout_field
{
	int nn;
	int ns;

	double b[4];
	double d;
	double nets[108];
	int shps[9];
};

//struct min_space
//{
//	int cnt;
//	double min_spac;
//	double *coord_del_p;
//};

typedef struct min_space
{
	int cnt;
	double min_spac;
	double *coord_del_p;
} ms;

//calculate_space
struct cal_space
{
	int flag_1;
	int flag_2;
	int numbs;
	int cnt_1;
	int cnt_2;
	int *divtype_mid_p;
	int *lencnt;
	int max_lencnt;
};

struct div_mesh
{
	int *div_int;
	int in_num; //insert_numbs
	double *ins;
};
struct nb_coord
{
	int **n_p;
	int **b_p;
};
//Declare functions
double **caluculate_coord(double nets[], int numbs_nn, int shps[]);
double compare_min(double a1, double a2, double a3, double a4);
double compare_max(double a1, double a2, double a3, double a4);
void min_apace_x(int nn, double nets[], double b[], struct min_space *min_x_p);
void min_apace_y(int nn, double nets[], double b[], struct min_space *min_y_p);
void quick_sort(double *arr, int low, int high);
void calculate_space_x(double b[], struct min_space *min_x_p, double ratio_nouni[], double ratio_uni[], struct cal_space *calspce_p);
void calculate_space_y(double b[], struct min_space *min_y_p, double ratio_nouni[], double ratio_uni[], struct cal_space *calspce_y_p);
void div_mesh_x(struct cal_space *calspce_p, double ratio_uni[], struct min_space *min_x_p, double b[], double ratio_nouni[], struct div_mesh *divmesh_x);
void div_mesh_y(struct cal_space *calspce_y_p, double ratio_uni[], struct min_space *min_y_p, double b[], double ratio_nouni[], struct div_mesh *divmesh_y);
void nb_coord(struct layout_field *field, struct min_space *min_y_p, struct min_space *min_x_p, struct div_mesh *divmesh_x, struct div_mesh *divmesh_y, double **nets_coord, struct nb_coord *nbc);

int main(void)
{
	//相邻非均匀划分之比常数
	double ratio_nouni[3] = {1.2, 1.4, 1.8};
	//当非均匀划分占比大于large_div,开始粗均匀划分，三个large_div与三个range，三个pow相对应
	double ratio_uni[3] = {0.6, 0.5, 0.4};
	struct layout_field lf = {27, 9, {-10, 0, 10, 9.9}, 3.9, {-0.017, 0.072, 0.017, 0.139, -0.114, 0.072, -0.078, 0.139, 0.078, 0.072, 0.114, 0.139, -0.155, 0.105, -0.133, 0.128, -0.156, 0.128, -0.132, 0.152, -0.157, 0.152, -0.131, 0.176, -0.158, 0.176, -0.13, 0.2, 0.133, 0.105, 0.155, 0.128, 0.132, 0.128, 0.156, 0.152, 0.131, 0.152, 0.157, 0.176, 0.13, 0.176, 0.158, 0.2, -0.059, 0.105, -0.037, 0.128, -0.06, 0.128, -0.036, 0.152, -0.061, 0.152, -0.035, 0.176, -0.062, 0.176, -0.034, 0.2, 0.037, 0.105, 0.059, 0.128, 0.036, 0.128, 0.060, 0.152, 0.035, 0.152, 0.061, 0.176, 0.034, 0.176, 0.062, 0.200, -0.061, 0.247, -0.024, 0.261, -0.062, 0.261, -0.023, 0.276, -0.063, 0.276, -0.022, 0.290, -0.064, 0.290, -0.021, 0.304, 0.024, 0.247, 0.061, 0.261, 0.023, 0.261, 0.062, 0.276, 0.022, 0.276, 0.063, 0.290, 0.021, 0.290, 0.064, 0.304}, {1, 1, 1, 4, 4, 4, 4, 4, 4}};

	struct layout_field *field = &lf;

	double **nets_coord_p = caluculate_coord(field->nets, field->ns, field->shps);
	//ms *min_x_p = malloc(sizeof(ms));
	struct min_space min_x_p;
	min_apace_x(field->nn, field->nets, field->b, &min_x_p);
	struct min_space min_y_p;
	//ms *min_y_p = malloc(sizeof(ms));
	min_apace_y(field->nn, field->nets, field->b, &min_y_p);

	struct cal_space calspce_p;
	calculate_space_x(field->b, &min_x_p, ratio_nouni, ratio_uni, &calspce_p);

	struct cal_space calspce_y_p;
	calculate_space_y(field->b, &min_y_p, ratio_nouni, ratio_uni, &calspce_y_p);
	struct div_mesh divmesh_x;
	div_mesh_x(&calspce_p, ratio_uni, &min_x_p, field->b, ratio_nouni, &divmesh_x);
	struct div_mesh divmesh_y;
	div_mesh_y(&calspce_y_p, ratio_uni, &min_y_p, field->b, ratio_nouni, &divmesh_y);

	struct nb_coord nbc;
	nb_coord(&lf, &min_y_p, &min_x_p, &divmesh_x, &divmesh_y, nets_coord_p, &nbc);
	free(divmesh_y.div_int);
	free(divmesh_x.div_int);
	free(min_y_p.coord_del_p);
	free(min_x_p.coord_del_p);
	



	free(&nets_coord_p);
	return 0;
}

// TODO
// make sure to fix the '...'
//
//void calculable_capacitor(struct layout_field *field, struct nb_coord *nbc, cholmod_triplet *T, struct div_mesh *divmesh_x, struct div_mesh *divmesh_y)
//{
//
//	int len_z = 2;
//	int len_grad = 1; //length_gradient_range = 1;
//	double q = 0;
//	double q1 = 0;
//	double p = 0;
//	double p1 = 0;
//
//	//z的大小需要比较出来
//	double max_row = 0;
//	double max_clumn = 0;
//	for (int u = 0; u < field->nn; u++)
//	{
//		if (max_row < nbc->n_p[u][3] - nbc->n_p[u][1] + 2 * len_z - 1)
//			max_row = nbc->n_p[u][3] - nbc->n_p[u][1]+ 2 * len_z - 1;
//		if (max_clumn < nbc->n_p[u][2] - nbc->n_p[u][0] + 2 * len_z - 1)
//			max_clumn = nbc->n_p[u][2] - nbc->n_p[u][0] + 2 * len_z - 1;
//	}
//
//	double **z = (double **)malloc((max_row + 1) * sizeof(double));
//	if (z == NULL)
//		printf("Failed to allocate memory\n");
//	else
//	{
//		for (int i = 0; i < (max_row + 1); ++i)
//		{
//			z[i] = (double *)malloc((max_clumn + 1) * sizeof(double));
//			if (z[i] == NULL)
//				printf("Failed to allocate memory\n");
//		}
//	}
//
//
//	double *dx1 = (double *)malloc((max_row + 1 - 2 * len_grad + 1) * sizeof(double));
//	if (dx1 == NULL)
//		printf("Failed to allocate memory\n");
//	double *dx2 = (double *)malloc((max_row + 1 - 2 * len_grad + 1) * sizeof(double));
//	if (dx2 == NULL)
//		printf("Failed to allocate memory\n");
//	double *dy1 = (double *)malloc((max_clumn + 1 - 2 * len_grad + 1) * sizeof(double));
//	if (dy1 == NULL)
//		printf("Failed to allocate memory\n");
//	double *dy2 = (double *)malloc((max_clumn + 1 - 2 * len_grad + 1) * sizeof(double));
//	if (dy2 == NULL)
//		printf("Failed to allocate memory\n");
//
//	double *C = (double *)malloc(field->ns * sizeof(double));
//	if (C == NULL)
//		printf("Failed to allocate memory\n");
//
//
//
//
//////////////////////////////////////////////////////////////////////////
//	for (int u = 0; u < field->ns; u++)
//	{
//		for (int j = 0; j < nbc->n_p[u][3] - nbc->n_p[u][1] + 2 * len_z + 1; j++)
//		{
//			for (int i = 0; i < nbc->n_p[u][2] - nbc->n_p[u][0] + 2 * len_z + 1;i++)
//			{
//				
//				z[j][i] = T->x[(j + nbc->n_p[u][1] - len_z - 2) * divmesh_x->in_num + i + nbc->n_p[u][0] - len_z - 1];
//			}
//		}
//
//		//求场强，一共四条线
//
//		
//		for (int j = len_grad; j < nbc->n_p[u][3] - nbc->n_p[u][1] + 2 * len_z + 1 - len_grad; j++)
//		{
//			
//			dx1[j] = (z[j, len_grad + 2 - 1] - z[j, len_grad - 1]) / (divmesh_x->ins[nbc->n_p[u][0] - 1 - 1] + divmesh_x->ins[nbc->n_p[u][0] - 2 - 1]);
//			dx2[j] = -(z[j, nbc->n_p[u][2] - nbc->n_p[u][0] + 2 * len_z + 1 - len_grad - 1 - 1] - z[j, nbc->n_p[u][2] - nbc->n_p[u][0] + 2 * len_z + 1 - len_grad + 1 - 1]) / (divmesh_x->ins[nbc->n_p[u][2] + 1 - 1] + divmesh_x->ins[nbc->n_p[u][2]]);
//			
//		}
//
//		for (int i = len_grad; i < nbc->n_p[u][2] - nbc->n_p[u][0] + 2 * len_z + 1 - len_grad; i++)
//		{
//			
//			dy1[i] = (z[len_grad + 2 - 1, i] - z[len_grad - 1, i]) / (divmesh_y->ins[nbc->n_p[u][1] - 1 - 1] + divmesh_y->ins[nbc->n_p[u][1] - 2 - 1]);
//			dy2[i] = -(z[nbc->n_p[u][3] - nbc->n_p[u][1] + 2 * len_z + 1 - len_grad - 1 - 1, i] - z[nbc->n_p[u][3] - nbc->n_p[u][1] + 2 * len_z + 1 - len_grad + 1 - 1, i]) / (divmesh_y->ins[nbc->n_p[u][3] - 1] + divmesh_y->ins[nbc->n_p[u][3] + 1 - 1]);
//			
//		}
//
//	
//
//		for (int v = len_grad; v < nbc->n_p[u][2] - nbc->n_p[u][0] + 2 * len_z + 1 - len_grad - 1; v++)
//		{
//			if (dy1[v] < 0)
//				dy1[v] = -dy1[v];
//			if (dy1[v + 1] < 0)
//				dy1[v + 1] = -dy1[v + 1];
//			p = p + divmesh_x->ins[v + nbc->n_p[u][0] - 1 - len_z ] * (dy1[v] + dy1[v + 1]) / 2;
//		}
//
//
//
//		for (int v = len_grad; v < nbc->n_p[u][2] - nbc->n_p[u][0] + 2 * len_z + 1 - len_grad - 1; v++)
//		{
//			if (dy2[v] < 0)
//				dy2[v] = -dy2[v];
//			if (dy2[v + 1] < 0)
//				dy2[v + 1] = -dy2[v + 1];
//			q = q + divmesh_x->ins[v + nbc->n_p[u][0] - 1 - len_z] * (dy2[v] + dy2[v + 1]) / 2;
//		}
//
//	
//		for (int v = len_grad; v < nbc->n_p[u][3] - nbc->n_p[u][1] + 2 * len_z + 1 - len_grad - 1; v++)
//		{
//			if (dx1[v] < 0)
//				dx1[v] = -dx1[v];
//			if (dx1[v + 1] < 0)
//				dx1[v + 1] = -dx1[v + 1];
//			q1 = q1 + divmesh_y->ins[v + nbc->n_p[u][1] - 1 - len_z ] * (dx1[v] + dx1[v + 1]) / 2;
//		}
//
//		for (int v = len_grad; v < nbc->n_p[u][3] - nbc->n_p[u][1] + 2 * len_z + 1 - len_grad - 1; v++)
//		{
//			if (dx2[v] < 0)
//				dx2[v] = -dx2[v];
//			if (dx2[v + 1] < 0)
//				dx2[v + 1] = -dx2[v + 1];
//			p1 = p1 + divmesh_y->ins[v + nbc->n_p[u][1] - 1 - len_z ] * (dx2[v] + dx2[v + 1]) / 2;
//		}
//
//		C[u] = field->d*8.854*(p+q+p1+q)/1000;
//	}
//
//	free(dx1);
//	free(dx2);
//	free(dy1);
//	free(dy2);
//	free(z);
//
//	free(nbc->n_p);
//}

//cholmod_sparse *cholmod_generate(struct div_mesh *divmesh_x, struct div_mesh *divmesh_y, struct layout_field *field, struct nb_coord *nbc, cholmod_common *cc)
//{
//	int flag = 1;
	int nrows, ncols, nnz;
	nrows = ncols = divmesh_x->in_num * divmesh_y->in_num;
	nnz = divmesh_x->in_num * divmesh_y->in_num * 5; //会多于非0；
	long int k = 0;
	//double u[447125];
	//double v[447125];
	//double x[447125];

	double *u = (double*)malloc(447125 * sizeof(double));
	if (u == NULL)
		printf("Failed to allocate memory\n");
	double *v = (double*)malloc(447125 * sizeof(double));
	if (v == NULL)
		printf("Failed to allocate memory\n");
	
	double *x = (double*)malloc(447125 * sizeof(double));
	if (x == NULL)
		printf("Failed to allocate memory\n");
	double c, d, e, f;

	
	for (int j = 0; j <= divmesh_y->in_num; j++)
	{
		for (int i = 0; i <= divmesh_x->in_num; i++)
		{
			if ((i == 0) || (j == 0) || (i == divmesh_x->in_num) || (j == divmesh_y->in_num))
			{
				u[k] = j * divmesh_x->in_num + i;
				v[k] = j * divmesh_x->in_num + i;
				x[k] = 1.0;
				flag = 0;
			}

			for (int q = 0; q < field->nn; q++)
			{
				if ((i <= nbc->b_p[q][2] ) && (i >= nbc->b_p[q][0]) && (j <= nbc->b_p[q][3] ) && (j >= nbc->b_p[q][1] ))
				{
					u[k] = j * divmesh_x->in_num + i;
					v[k] = j * divmesh_x->in_num + i;
					x[k] = 1.0;
					flag = 0;
				}
			}

			if (flag == 1)
			{

				c = divmesh_x->ins[i];
				d = divmesh_x->ins[i - 1];
				e = divmesh_y->ins[j];
					f = divmesh_y->ins[j - 1];

				u[k] = j * divmesh_x->in_num + i;
				v[k] = j * divmesh_x->in_num + i + 1;
				x[k] = (e + f) / c;
				k++;

				u[k] = j * divmesh_x->in_num + i;
				v[k] = j * divmesh_x->in_num + i - 1;
				x[k] = (e + f) / d;
				k++;

				u[k] = j * divmesh_y->in_num + i;
				v[k] = j * divmesh_y->in_num + i;
				x[k] = -(e + f) / c- (e + f) / d- (c + d) / e- (c + d) / f;
				k++;

				u[k] = j * divmesh_x->in_num + i;
				v[k] = (j + 1) * divmesh_x->in_num + i;
				x[k] = (c + d) / e;
				k++;

				u[k] = j * divmesh_x->in_num + i;
				v[k] = (j - 1) * divmesh_x->in_num + i;
				x[k] = (c + d) / f;
		
			}

			k++;
			flag = 1;
		}
	}

//	}
//
//	int n = 0;
//
//
//	
//
//	if (field->shps[0] == 4)
//	{
//		for (int u = 0; u < 4; u++)
//		{
//			for (int i = nbc->b_p[u][0]; i < nbc->b_p[u][2]; i++)
//			{
//				for (int j = nbc->b_p[u][1]; j < nbc->b_p[u][3]; j++))
//				{
//					bidx[n] = j * divmesh_x->in_num + i;
//					x[bidx[n]] = 1;
//					n++;
//				}
//			}
//		}
//		else
//		{
//			for (int i = nbc->b_p[0][0]; i < nbc->b_p[0][2]; i++)
//			{
//				for (int(int j = nbc->b_p[0][1]; j < nbc->b_p[0][3]; j++))
//				{
//					bidx[n] = j * divmesh_x->in_num + i;
//					x[bidx[n]] = 1;
//					n++;
//				}
//			}
//		}
//	}
//
//free(nbc->b_p);
//free(divmesh_y->ins);
//free(divmesh_x->ins);
//}
//}

void nb_coord(struct layout_field *field, struct min_space *min_y_p, struct min_space *min_x_p, struct div_mesh *divmesh_x, struct div_mesh *divmesh_y, double **nets_coord, struct nb_coord *nbc)
{
	//block_int
	nbc->b_p = (int **)malloc(field->nn * sizeof(int));
	if (nbc->b_p == NULL)
		printf("Failed to allocate memory\n");
	else
	{
		for (int i = 0; i < field->nn; ++i)
		{
			nbc->b_p[i] = (int *)malloc(4 * sizeof(int));
			if (nbc->b_p[i] == NULL)
				printf("Failed to allocate memory\n");
		}
	}

	for (int i = 0; i < field->nn; i++)
	{
		for (int j = 0; j < min_x_p->cnt; j++)
		{
			if ((field->nets[i * 4 + 0] > min_x_p->coord_del_p[j] - min_x_p->min_spac) && (field->nets[i * 4 + 0] < min_x_p->coord_del_p[j] + min_x_p->min_spac))
				nbc->b_p[i][0] = divmesh_x->div_int[j];
			if ((field->nets[i * 4 + 2] > min_x_p->coord_del_p[j] - min_x_p->min_spac) && (field->nets[i * 4 + 2] < min_x_p->coord_del_p[j] + min_x_p->min_spac))
				nbc->b_p[i][2] = divmesh_x->div_int[j];
		}
	}
	for (int i = 0; i < field->nn; i++)
	{
		for (int j = 0; j <= min_y_p->cnt; j++)
		{
			if ((field->nets[i * 4 + 1] > min_y_p->coord_del_p[j] - min_y_p->min_spac) && (field->nets[i * 4 + 1] < min_y_p->coord_del_p[j] + min_y_p->min_spac))
				nbc->b_p[i][1] = divmesh_y->div_int[j];
			if ((field->nets[i * 4 + 3] > min_y_p->coord_del_p[j] - min_y_p->min_spac) && (field->nets[i * 4 + 3] < min_y_p->coord_del_p[j] + min_y_p->min_spac))
			{
				nbc->b_p[i][3] = divmesh_y->div_int[j];
				break;
			}
		}
	}

	nbc->n_p = (int **)malloc((field->ns) * sizeof(int));
	if (nbc->n_p == NULL)
		printf("Failed to allocate memory\n");
	else
	{
		for (int i = 0; i < field->ns; ++i)
		{
			nbc->n_p[i] = (int *)malloc(4 * sizeof(int));
			if (nbc->n_p[i] == NULL)
				printf("Failed to allocate memory\n");
		}
	}

	for (int i = 0; i < field->ns; i++)
	{
		for (int j = 0; j < min_x_p->cnt; j++)
		{
			if ((nets_coord[i][0] > min_x_p->coord_del_p[j] - min_x_p->min_spac) && (nets_coord[i][0] < min_x_p->coord_del_p[j] + min_x_p->min_spac))
				nbc->n_p[i][0] = divmesh_x->div_int[j];
			if ((nets_coord[i][2] > min_x_p->coord_del_p[j] - min_x_p->min_spac) && (nets_coord[i][2] < min_x_p->coord_del_p[j] + min_x_p->min_spac))
			{
				nbc->n_p[i][2] = divmesh_x->div_int[j];
				break;
			}
		}
	}

	for (int i = 0; i < field->ns; i++)
	{
		for (int j = 0; j <= min_y_p->cnt; j++)
		{
			if ((nets_coord[i][1] > min_y_p->coord_del_p[j] - min_y_p->min_spac) && (nets_coord[i][1] < min_y_p->coord_del_p[j] + min_y_p->min_spac))
				nbc->n_p[i][1] = divmesh_y->div_int[j];
			if ((nets_coord[i][3] > (min_y_p->coord_del_p[j] - min_y_p->min_spac)) && (nets_coord[i][3] < (min_y_p->coord_del_p[j] + min_y_p->min_spac)))
			{
				nbc->n_p[i][3] = divmesh_y->div_int[j];
				break;
			}
		}
	}

}
void div_mesh_x(struct cal_space *calspce_p, double ratio_uni[], struct min_space *min_x_p, double b[], double ratio_nouni[], struct div_mesh *divmesh_x)
{
	divmesh_x->ins = (double *)malloc(calspce_p->numbs* sizeof(double));
	if (divmesh_x->ins == NULL)
		printf("Failed to create memory of divtype_mid \n");

	double sum_r = 0;
	double *insert_r = (double *)malloc((calspce_p->cnt_1+1) * sizeof(double));
	if (insert_r == NULL)
		printf("Failed to create memory of divtype_mid \n");
	int cnt_r = 0;
	insert_r[cnt_r] = min_x_p->min_spac;

	int divnums_r = 0;
	double divunif_r = 0;

	divmesh_x->in_num = 0;
	int insert_numbs = 0;
	if (calspce_p->flag_1 > 0)
	{
		while (sum_r < ratio_uni[calspce_p->flag_1 - 1] * (min_x_p->coord_del_p[0] - b[0] - RANGE0 * min_x_p->min_spac))
		{
			sum_r = sum_r + insert_r[cnt_r];
			insert_r[cnt_r + 1] = insert_r[cnt_r] * ratio_nouni[calspce_p->flag_1 - 1];
			cnt_r++;
		}
		divnums_r = (int) round((min_x_p->coord_del_p[0] - b[0] - sum_r - RANGE0 * min_x_p->min_spac) / insert_r[cnt_r - 1]);
		divunif_r = (min_x_p->coord_del_p[0] - b[0] - sum_r - RANGE0 * min_x_p->min_spac) / divnums_r;

		for (int i = 0; i < divnums_r; i++)
		{
			divmesh_x->ins[divmesh_x->in_num] = divunif_r;
			divmesh_x->in_num++;
		}
		for (int i = 0; i < cnt_r; i++)
		{
			divmesh_x->ins[divmesh_x->in_num] = insert_r[cnt_r - i - 1];
			divmesh_x->in_num++;
		}
		for (int i = 0; i < RANGE0; i++)
		{
			divmesh_x->ins[divmesh_x->in_num] = min_x_p->min_spac;
			divmesh_x->in_num++;
		}
	}
	else
	{
		divunif_r = (min_x_p->coord_del_p[0] - b[0]) / calspce_p->cnt_1;
		for (int i = 0; i < calspce_p->cnt_1; i++)
		{
			divmesh_x->ins[divmesh_x->in_num] = divunif_r;
			divmesh_x->in_num++;
		}
	}

	free(insert_r);
	//细划分区格点坐标
	divmesh_x->div_int = (int *)malloc(calspce_p->numbs * sizeof(int));
	if (divmesh_x->div_int == NULL)
		printf("Failed to create memory of divtype_mid \n");
	int cntint_x = 0;
	divmesh_x->div_int[cntint_x] = divmesh_x->in_num;
	cntint_x++;

	//
	double *x_array = (double *)malloc((calspce_p->max_lencnt+2) * sizeof(double));
	if (x_array == NULL)
		printf("Failed to create memory of divtype_mid \n");
	double sum_x = 0;
	int divnums_x = 0;
	double divunif_x = 0;
	double len_unif_x = 0;
	cnt_r = 0;
	for (int i = 0; i < min_x_p->cnt; i++)
	{
		if (calspce_p->divtype_mid_p[i] > 0)
		{
			sum_x = 0;
			cnt_r = 0;
			x_array[0] = min_x_p->min_spac;
			while ((min_x_p->coord_del_p[i + 1] - min_x_p->coord_del_p[i] - 2 * 7.0 * min_x_p->min_spac) * ratio_uni[calspce_p->divtype_mid_p[i] - 1] > 2 * sum_x)
			{
				sum_x = sum_x + x_array[cnt_r];
				x_array[cnt_r + 1] = x_array[cnt_r] * ratio_nouni[calspce_p->divtype_mid_p[i] - 1];
				cnt_r++;
			}

			divnums_x = (int)round((min_x_p->coord_del_p[i + 1] - min_x_p->coord_del_p[i] - 2 * sum_x - 2 * RANGE0 * min_x_p->min_spac) / x_array[cnt_r - 1]);
			divunif_x = (min_x_p->coord_del_p[i + 1] - min_x_p->coord_del_p[i] - 2 * sum_x - 2 * RANGE0 * min_x_p->min_spac) / divnums_x;

			//清醒的时候重新检查这里
			for (int j = 0; j < RANGE0; j++)
			{
				divmesh_x->ins[divmesh_x->in_num] = min_x_p->min_spac;
				divmesh_x->in_num++;
			}

			for (int j = 0; j < cnt_r; j++)
			{
				divmesh_x->ins[divmesh_x->in_num] = x_array[j];
				divmesh_x->in_num++;
			}
			for (int j = 0; j < divnums_x; j++)
			{
				divmesh_x->ins[divmesh_x->in_num] = divunif_x;
				divmesh_x->in_num++;
			}

			for (int j = 0; j < cnt_r; j++)
			{
				divmesh_x->ins[divmesh_x->in_num] = x_array[cnt_r - j - 1];
				divmesh_x->in_num++;
			}

			for (int j = 0; j < RANGE0; j++)
			{
				divmesh_x->ins[divmesh_x->in_num] = min_x_p->min_spac;
				divmesh_x->in_num++;
			}
			divmesh_x->div_int[cntint_x] = divmesh_x->in_num;
			cntint_x++;
		}
		else
		{
			len_unif_x = (min_x_p->coord_del_p[i + 1] - min_x_p->coord_del_p[i]) / calspce_p->lencnt[i];
			for (int k = 0; k < calspce_p->lencnt[i]; k++)
			{
				divmesh_x->ins[divmesh_x->in_num] = len_unif_x;
				divmesh_x->in_num++;
			}

			divmesh_x->div_int[cntint_x] = divmesh_x->in_num;
			cntint_x++;
		}
	}
	
	int cnt_l = 0;
	double sum_l = 0;
	double *insert_l = (double *)malloc((calspce_p->cnt_2 +2)* sizeof(double));
	if (insert_l == NULL)
		printf("Failed to create memory of divtype_mid \n");
	insert_l[cnt_l] = min_x_p->min_spac;

	int divnums_l = 0;
	double divunif_l = 0;

	if (calspce_p->flag_2 > 0)
	{
		while (sum_l < (ratio_uni[calspce_p->flag_2 - 1] * (b[2] - min_x_p->coord_del_p[min_x_p->cnt] - RANGE0 * min_x_p->min_spac)))
		{
			sum_l = sum_l + insert_l[cnt_l];
			insert_l[cnt_l + 1] = insert_l[cnt_l] * ratio_nouni[calspce_p->flag_2 - 1];
			cnt_l++;
		}
		divnums_l = (int) round((b[2] - min_x_p->coord_del_p[min_x_p->cnt] - sum_l - RANGE0 * min_x_p->min_spac) / insert_l[cnt_l - 1]);
		divunif_l = (b[2] - min_x_p->coord_del_p[min_x_p->cnt] - sum_l - RANGE0 * min_x_p->min_spac) / divnums_l;

		for (int i = 0; i < RANGE0; i++)
		{
			divmesh_x->ins[divmesh_x->in_num] = min_x_p->min_spac;
			divmesh_x->in_num++;
		}
		for (int i = 0; i < cnt_l; i++)
		{
			divmesh_x->ins[divmesh_x->in_num] = insert_l[i];
			divmesh_x->in_num++;
		}
		for (int i = 0; i < divnums_l; i++)
		{
			divmesh_x->ins[divmesh_x->in_num] = divunif_l;
			divmesh_x->in_num++;
		}
	}
	else
	{
		divunif_l = (b[2] - min_x_p->coord_del_p[min_x_p->cnt]) / calspce_p->cnt_2;
		for (int i = 0; i < calspce_p->cnt_2; i++)
		{
			divmesh_x->ins[divmesh_x->in_num] = divunif_r;
			divmesh_x->in_num++;
		}
	}
	free(x_array);
	free(insert_l);  
}
void div_mesh_y(struct cal_space *calspce_y_p, double ratio_uni[], struct min_space *min_y_p, double b[], double ratio_nouni[], struct div_mesh *divmesh_y)
{
	divmesh_y->ins = (double *)malloc(calspce_y_p->numbs * sizeof(double));
	if (divmesh_y->ins == NULL)
		printf("Failed to create memory of divtype_mid \n");

	double sum_h = 0;
	double *insert_h = (double *)malloc((calspce_y_p->cnt_2+2) * sizeof(double));
	if (insert_h == NULL)
		printf("Failed to create memory of divtype_mid \n");
	int cnt_h = 0;
	insert_h[cnt_h] = min_y_p->min_spac;

	int divnums_h = 0;
	double divunif_h = 0;

	divmesh_y->in_num = 0;
	int insert_numbs = 0;
	if (calspce_y_p->flag_2 > 0)
	{
		while (sum_h < ratio_uni[calspce_y_p->flag_2 - 1] * (min_y_p->coord_del_p[0] - b[1] - RANGE0 * min_y_p->min_spac))
		{
			sum_h = sum_h + insert_h[cnt_h];
			insert_h[cnt_h + 1] = insert_h[cnt_h] * ratio_nouni[calspce_y_p->flag_2 - 1];
			cnt_h++;
		}
		divnums_h = (int) round((min_y_p->coord_del_p[0] - b[1] - sum_h - RANGE0 * min_y_p->min_spac) / insert_h[cnt_h - 1]);
		divunif_h = (min_y_p->coord_del_p[0] - b[1] - sum_h - RANGE0 * min_y_p->min_spac) / divnums_h;

		for (int i = 0; i < divnums_h; i++)
		{
			divmesh_y->ins[divmesh_y->in_num] = divunif_h;
			divmesh_y->in_num++;
		}
		for (int i = 0; i < cnt_h; i++)
		{
			divmesh_y->ins[divmesh_y->in_num] = insert_h[cnt_h - i -1];
			divmesh_y->in_num++;
		}
		for (int i = 0; i < RANGE0; i++)
		{
			divmesh_y->ins[divmesh_y->in_num] = min_y_p->min_spac;
			divmesh_y->in_num++;
		}
	}
	else
	{
		divunif_h = (min_y_p->coord_del_p[0] - b[1]) / calspce_y_p->cnt_2;
		for (int i = 0; i < calspce_y_p->cnt_2; i++)
		{
			divmesh_y->ins[divmesh_y->in_num] = divunif_h;
			divmesh_y->in_num++;
		}
	}

	
	//细划分区格点坐标
	divmesh_y->div_int = (int *)malloc(calspce_y_p->numbs * sizeof(int));
	if (divmesh_y->ins == NULL)
		printf("Failed to create memory of divtype_mid \n");
	int cntint_y = 0;
	divmesh_y->div_int[cntint_y] = divmesh_y->in_num;
	cntint_y++;

	//
	double *y_array = (double *)malloc(calspce_y_p->max_lencnt  * sizeof(double));
	if (y_array == NULL)
		printf("Failed to create memory of divtype_mid \n");
	double sum_y = 0;
	int divnums_y = 0;
	double divunif_y = 0;
	double len_unif_y = 0;
	cnt_h = 0;
	for (int i = 0; i < min_y_p->cnt; i++)
	{
		if (calspce_y_p->divtype_mid_p[i] > 0)
		{
			sum_y = 0;
			cnt_h = 0;
			y_array[0] = min_y_p->min_spac;
			while ((min_y_p->coord_del_p[i + 1] - min_y_p->coord_del_p[i] - 2 * 7.0 * min_y_p->min_spac) * ratio_uni[calspce_y_p->divtype_mid_p[i] - 1] > 2 * sum_y)
			{
				sum_y = sum_y + y_array[cnt_h];
				y_array[cnt_h + 1] = y_array[cnt_h] * ratio_nouni[calspce_y_p->divtype_mid_p[i] - 1];
				cnt_h++;
			}

			divnums_y = (int)round((min_y_p->coord_del_p[i + 1] - min_y_p->coord_del_p[i] - 2 * sum_y - 2 * RANGE0 * min_y_p->min_spac) / y_array[cnt_h - 1]);
			divunif_y = (min_y_p->coord_del_p[i + 1] - min_y_p->coord_del_p[i] - 2 * sum_y - 2 * RANGE0 * min_y_p->min_spac) / divnums_y;

			//清醒的时候重新检查这里
			for (int j = 0; j < RANGE0; j++)
			{
				divmesh_y->ins[divmesh_y->in_num] = min_y_p->min_spac;
				divmesh_y->in_num++;
			}

			for (int j = 0; j < cnt_h; j++)
			{
				divmesh_y->ins[divmesh_y->in_num] = y_array[j];
				divmesh_y->in_num++;
			}
			for (int j = 0; j < divnums_y; j++)
			{
				divmesh_y->ins[divmesh_y->in_num] = divunif_y;
				divmesh_y->in_num++;
			}
			for (int j = 0; j < cnt_h; j++)
			{
				divmesh_y->ins[divmesh_y->in_num] = y_array[cnt_h - j - 1];
				divmesh_y->in_num++;
			}
			for (int j = 0; j < RANGE0; j++)
			{
				divmesh_y->ins[divmesh_y->in_num] = min_y_p->min_spac;
				divmesh_y->in_num++;
			}
			divmesh_y->div_int[cntint_y] = divmesh_y->in_num;
			cntint_y++;
		}
		else
		{
			len_unif_y = (min_y_p->coord_del_p[i + 1] - min_y_p->coord_del_p[i]) / calspce_y_p->lencnt[i];
			for (int k = 0; k < calspce_y_p->lencnt[i]; k++)
			{
				divmesh_y->ins[divmesh_y->in_num] = len_unif_y;
				divmesh_y->in_num++;
			}

			divmesh_y->div_int[cntint_y] = divmesh_y->in_num;
			cntint_y++;
		}
	}
	free(insert_h);
	free(y_array);
	free(calspce_y_p->lencnt);
	free(calspce_y_p->divtype_mid_p);
	int cnt_low = 0;
	double sum_low = 0;
	double *insert_low = (double *)malloc(calspce_y_p->cnt_1 * sizeof(double));
	if (insert_low == NULL)
		printf("Failed to create memory of divtype_mid \n");
	insert_low[cnt_low] = min_y_p->min_spac;

	int divnums_low = 0;
	double divunif_low = 0;

	if (calspce_y_p->flag_1 > 0)
	{
		while (sum_low < ratio_uni[calspce_y_p->flag_1 - 1] * (b[3] - min_y_p->coord_del_p[min_y_p->cnt] - RANGE0 * min_y_p->min_spac))
		{
			sum_low = sum_low + insert_low[cnt_low];
			insert_low[cnt_low + 1] = insert_low[cnt_low] * ratio_nouni[calspce_y_p->flag_1 - 1];
			cnt_low++;
		}
		divnums_low = (int)round((b[3] - min_y_p->coord_del_p[min_y_p->cnt] - sum_low - RANGE0 * min_y_p->min_spac) / insert_low[cnt_low - 1]);
		divunif_low = (b[3] - min_y_p->coord_del_p[min_y_p->cnt] - sum_low - RANGE0 * min_y_p->min_spac) / divnums_low;

		for (int i = 0; i < RANGE0; i++)
		{
			divmesh_y->ins[divmesh_y->in_num] = min_y_p->min_spac;
			divmesh_y->in_num++;
		}
		for (int i = 0; i < cnt_low; i++)
		{
			divmesh_y->ins[divmesh_y->in_num] = insert_low[i];
			divmesh_y->in_num++;
		}
		for (int i = 0; i < divnums_low; i++)
		{
			divmesh_y->ins[divmesh_y->in_num] = divunif_low;
			divmesh_y->in_num++;
		}
	}
	else
	{
		divunif_low = (b[3] - min_y_p->coord_del_p[min_y_p->cnt]) / calspce_y_p->cnt_1;
		for (int i = 0; i < calspce_y_p->cnt_1; i++)
		{
			divmesh_y->ins[divmesh_y->in_num] = divunif_low;
			divmesh_y->in_num++;
		}
	}

	free(insert_low);
}
void calculate_space_x(double b[], struct min_space *min_x_p, double ratio_nouni[], double ratio_uni[], struct cal_space *calspce_p)
{
	//判断x轴方向最后一段属于哪种划分类型
	calspce_p->flag_1 = 0;
	if (((b[2] - min_x_p->coord_del_p[min_x_p->cnt]) > ((RANGE1 - 0.01) * min_x_p->min_spac)) && (b[2] - min_x_p->coord_del_p[min_x_p->cnt]) < (RANGE2 + 0.01) * min_x_p->min_spac)
		calspce_p->flag_1 = 1;
	else if (((b[2] - min_x_p->coord_del_p[min_x_p->cnt]) > ((RANGE2 - 0.01) * min_x_p->min_spac)) && (b[2] - min_x_p->coord_del_p[min_x_p->cnt]) < (RANGE3 + 0.01) * min_x_p->min_spac)
		calspce_p->flag_1 = 2;
	else if ((b[2] - min_x_p->coord_del_p[min_x_p->cnt]) > ((RANGE3 - 0.01) * min_x_p->min_spac))
		calspce_p->flag_1 = 3;
	else if ((b[2] - min_x_p->coord_del_p[min_x_p->cnt]) < ((RANGE1 - 0.01) * min_x_p->min_spac))
		calspce_p->flag_1 = 0;

	calspce_p->numbs = 0;
	calspce_p->cnt_1 = 3;
	double len_right = 0;
	double len_mul_right = 0;
	double len_uni_right = 0;
	if (calspce_p->flag_1 > 0)
	{
		len_right = pow(ratio_nouni[calspce_p->flag_1 - 1], 3);
		len_mul_right = pow(ratio_nouni[calspce_p->flag_1 - 1], 3);
		len_uni_right = (b[2] - min_x_p->coord_del_p[min_x_p->cnt] - RANGE0 * min_x_p->min_spac) * ratio_uni[calspce_p->flag_1 - 1] / min_x_p->min_spac;
		while (((len_right - 1) / (ratio_nouni[calspce_p->flag_1 - 1] - 1)) < len_uni_right)
		{
			len_right = len_right * len_mul_right;
			calspce_p->cnt_1 = calspce_p->cnt_1 + 3;
		}
		calspce_p->numbs = calspce_p->numbs + calspce_p->cnt_1 + (int)round(len_uni_right * (1.0 - ratio_uni[calspce_p->flag_1 - 1]) / ratio_uni[calspce_p->flag_1 - 1] / pow(ratio_nouni[calspce_p->flag_1 - 1], (double)calspce_p->cnt_1 - 3.0)) + (int)RANGE0;
	}
	else
	{
		calspce_p->cnt_1 = (int)round((b[2] - min_x_p->coord_del_p[min_x_p->cnt]) / min_x_p->min_spac);
		calspce_p->numbs = calspce_p->numbs + calspce_p->cnt_1;
	}

	calspce_p->flag_2 = 0;
	if (((min_x_p->coord_del_p[0] - b[0]) > ((RANGE1 - 0.01) * min_x_p->min_spac)) && (min_x_p->coord_del_p[0] - b[0]) < (RANGE2 + 0.01) * min_x_p->min_spac)
		calspce_p->flag_2 = 1;
	else if (((min_x_p->coord_del_p[0] - b[0]) > ((RANGE2 - 0.01) * min_x_p->min_spac)) && (min_x_p->coord_del_p[0] - b[0]) < (RANGE3 + 0.01) * min_x_p->min_spac)
		calspce_p->flag_2 = 2;
	else if ((min_x_p->coord_del_p[0] - b[0]) > ((RANGE3 - 0.01) * min_x_p->min_spac))
		calspce_p->flag_2 = 3;
	else if ((min_x_p->coord_del_p[0] - b[0]) < ((RANGE1 - 0.01) * min_x_p->min_spac))
		calspce_p->flag_2 = 0;

	calspce_p->cnt_2 = 3;
	double len_left = 0;
	double len_mul_left = 0;
	double len_uni_left = 0;
	if (calspce_p->flag_2 > 0)
	{
		len_left = pow(ratio_nouni[calspce_p->flag_2 - 1], 3);
		len_mul_left = pow(ratio_nouni[calspce_p->flag_2 - 1], 3);
		len_uni_left = (min_x_p->coord_del_p[0] - b[0] - RANGE0 * min_x_p->min_spac) * ratio_uni[calspce_p->flag_2 - 1] / min_x_p->min_spac;
		while (((len_left - 1) / (ratio_nouni[calspce_p->flag_2 - 1] - 1)) < len_uni_left)
		{
			len_left = len_left * len_mul_left;
			calspce_p->cnt_2 = calspce_p->cnt_2 + 3;
		}
		calspce_p->numbs = calspce_p->numbs + calspce_p->cnt_2 + (int)round(len_uni_left * (1.0 - ratio_uni[calspce_p->flag_2 - 1]) / ratio_uni[calspce_p->flag_2 - 1] / pow(ratio_nouni[calspce_p->flag_2 - 1], calspce_p->cnt_2 - 3.0)) + (int)RANGE0;
	}
	else
	{
		calspce_p->cnt_2 = (int) round((b[0] - min_x_p->coord_del_p[min_x_p->cnt]) / min_x_p->min_spac);
		calspce_p->numbs = calspce_p->numbs + calspce_p->cnt_2;
	}

	//pow_flag_x用来储存中间的所有间隔的划分类型

	//int* divtype_mid = (int**)malloc(min_x_p->cnt * sizeof(int));
	calspce_p->divtype_mid_p = (int *)malloc(min_x_p->cnt * sizeof(int));
	if (calspce_p->divtype_mid_p == NULL)
		printf("Failed to create memory of divtype_mid \n");
	calspce_p->lencnt = (int *)malloc(min_x_p->cnt* sizeof(int));
	if (calspce_p->lencnt == NULL)
		printf("Failed to create memory of divtype_mid \n");
	double len_x = 0;
	double lenmul_x = 0;
	double long_x = 0;
	calspce_p->max_lencnt = 0;
	for (int i = 0; i < min_x_p->cnt; i++)
	{
		if (((min_x_p->coord_del_p[i + 1] - min_x_p->coord_del_p[i]) > (RANGE1 - 0.01) * min_x_p->min_spac) && ((min_x_p->coord_del_p[i + 1] - min_x_p->coord_del_p[i]) < (RANGE2 + 0.01) * min_x_p->min_spac))
			calspce_p->divtype_mid_p[i] = 1;
		else if (((min_x_p->coord_del_p[i + 1] - min_x_p->coord_del_p[i]) > (RANGE2 - 0.01) * min_x_p->min_spac) && ((min_x_p->coord_del_p[i + 1] - min_x_p->coord_del_p[i]) < (RANGE3 + 0.01) * min_x_p->min_spac))
			calspce_p->divtype_mid_p[i] = 2;
		else if ((min_x_p->coord_del_p[i + 1] - min_x_p->coord_del_p[i]) > (RANGE3 - 0.01) * min_x_p->min_spac)
			calspce_p->divtype_mid_p[i] = 3;
		else if ((min_x_p->coord_del_p[i + 1] - min_x_p->coord_del_p[i]) < (RANGE2 - 0.01) * min_x_p->min_spac)
			calspce_p->divtype_mid_p[i] = 0;

		if (  calspce_p->divtype_mid_p[i] > 0)
		{
			calspce_p->lencnt[i] = 3;
			len_x = pow(ratio_nouni[calspce_p->divtype_mid_p[i] - 1], 3);
			lenmul_x = pow(ratio_nouni[calspce_p->divtype_mid_p[i] - 1], 3);
			long_x = ((min_x_p->coord_del_p[i + 1] - min_x_p->coord_del_p[i] - 2 * RANGE0 * min_x_p->min_spac) * ratio_uni[calspce_p->divtype_mid_p[i] - 1] / min_x_p->min_spac) / 2;
			while ((len_x - 1) / (ratio_nouni[calspce_p->divtype_mid_p[i] - 1] - 1) < long_x)
			{
				len_x = len_x * lenmul_x;
				calspce_p->lencnt[i] = calspce_p->lencnt[i] + 3;
				if (calspce_p->max_lencnt < calspce_p->lencnt[i])
					calspce_p->max_lencnt = calspce_p->lencnt[i];
			}
			calspce_p->numbs = calspce_p->numbs + 2 * calspce_p->lencnt[i] + (int)round((2 * long_x * (1 - ratio_uni[calspce_p->divtype_mid_p[i] - 1]) / ratio_uni[calspce_p->divtype_mid_p[i] - 1]) / pow(ratio_nouni[calspce_p->divtype_mid_p[i] - 1], calspce_p->lencnt[i] - 3)) + 2 * (int)RANGE0;
		}
		else
		{
			calspce_p->lencnt[i] =(int) round((min_x_p->coord_del_p[i + 1] - min_x_p->coord_del_p[i]) / min_x_p->min_spac);

			calspce_p->numbs = (calspce_p->numbs + calspce_p->lencnt[i]);
		}
	}
}

void calculate_space_y(double b[], struct min_space *min_y_p, double ratio_nouni[], double ratio_uni[], struct cal_space *calspce_y_p)
{

	calspce_y_p->flag_2 = 0;
	if (((min_y_p->coord_del_p[0] - b[1]) > ((RANGE1 - 0.01) * min_y_p->min_spac)) && (min_y_p->coord_del_p[0] - b[1]) < (RANGE2 + 0.01) * min_y_p->min_spac)
		calspce_y_p->flag_2 = 1;
	else if (((min_y_p->coord_del_p[0] - b[1]) > ((RANGE2 - 0.01) * min_y_p->min_spac)) && (min_y_p->coord_del_p[0] - b[1]) < (RANGE3 + 0.01) * min_y_p->min_spac)
		calspce_y_p->flag_2 = 2;
	else if ((min_y_p->coord_del_p[0] - b[1]) > ((RANGE3 - 0.01) * min_y_p->min_spac))
		calspce_y_p->flag_2 = 3;
	else if ((min_y_p->coord_del_p[0] - b[1]) < ((RANGE1 - 0.01) * min_y_p->min_spac))
		calspce_y_p->flag_2 = 0;

	calspce_y_p->cnt_2 = 3;
	double len_low = 0;
	double len_mul_low = 0;
	double len_uni_low = 0;
	calspce_y_p->numbs = 0;
	if (calspce_y_p->flag_2 > 0)
	{
		len_low = pow(ratio_nouni[calspce_y_p->flag_2 - 1], 3);
		len_mul_low = pow(ratio_nouni[calspce_y_p->flag_2 - 1], 3);
		len_uni_low = (min_y_p->coord_del_p[0] - b[1] - RANGE0 * min_y_p->min_spac) * ratio_uni[calspce_y_p->flag_2 - 1] / min_y_p->min_spac;
		//这里明天问杨，是不是写错了，对应matlab
		//while (((len_low - 1) / (ratio_nouni[calspce_y_p->flag_2 - 1] - 1)) < len_uni_low)
		while (((len_low - 1) / (ratio_nouni[calspce_y_p->flag_2 - 1] - 1)) < len_uni_low)
		{
			len_low = len_low * len_mul_low;
			calspce_y_p->cnt_2 = calspce_y_p->cnt_2 + 3;
		}
		calspce_y_p->numbs = calspce_y_p->numbs + calspce_y_p->cnt_2 + (int)round(len_uni_low * (1 - ratio_uni[calspce_y_p->flag_2 - 1]) / ratio_uni[calspce_y_p->flag_2 - 1] / pow(ratio_nouni[calspce_y_p->flag_2 - 1], calspce_y_p->cnt_2 - 3)) + (int)RANGE0;
	}
	else
	{
		calspce_y_p->cnt_2 = (int)round((b[1] - min_y_p->coord_del_p[min_y_p->cnt]) / min_y_p->min_spac);
		calspce_y_p->numbs = calspce_y_p->numbs + calspce_y_p->cnt_2;
	}

	//判断y轴方向最后一段属于哪种划分类型
	calspce_y_p->flag_1 = 0;
	if (((b[3] - min_y_p->coord_del_p[min_y_p->cnt]) > ((RANGE1 - 0.01) * min_y_p->min_spac)) && (b[3] - min_y_p->coord_del_p[min_y_p->cnt]) < (RANGE2 + 0.01) * min_y_p->min_spac)
		calspce_y_p->flag_1 = 1;
	else if (((b[3] - min_y_p->coord_del_p[min_y_p->cnt]) > ((RANGE2 - 0.01) * min_y_p->min_spac)) && (b[3] - min_y_p->coord_del_p[min_y_p->cnt]) < (RANGE3 + 0.01) * min_y_p->min_spac)
		calspce_y_p->flag_1 = 2;
	else if ((b[3] - min_y_p->coord_del_p[min_y_p->cnt]) > ((RANGE3 - 0.01) * min_y_p->min_spac))
		calspce_y_p->flag_1 = 3;
	else if ((b[3] - min_y_p->coord_del_p[min_y_p->cnt]) < ((RANGE1 - 0.01) * min_y_p->min_spac))
		calspce_y_p->flag_1 = 0;

	//calspce_y_p->numbs = 0;
	calspce_y_p->cnt_1 = 3;
	double len_high = 0;
	double len_mul_high = 0;
	double len_uni_high = 0;
	if (calspce_y_p->flag_1 > 0)
	{
		len_high = pow(ratio_nouni[calspce_y_p->flag_1 - 1], 3);
		len_mul_high = pow(ratio_nouni[calspce_y_p->flag_1 - 1], 3);
		len_uni_high = (b[3] - min_y_p->coord_del_p[min_y_p->cnt] - RANGE0 * min_y_p->min_spac) * ratio_uni[calspce_y_p->flag_1 - 1] / min_y_p->min_spac;
		while (((len_high - 1) / (ratio_nouni[calspce_y_p->flag_1 - 1] - 1)) < len_uni_high)
		{
			len_high = len_high * len_mul_high;
			calspce_y_p->cnt_1 = calspce_y_p->cnt_1 + 3;
		}
		calspce_y_p->numbs = calspce_y_p->numbs + calspce_y_p->cnt_1 + (int)round(len_uni_high * (1.0 - ratio_uni[calspce_y_p->flag_1 - 1]) / ratio_uni[calspce_y_p->flag_1 - 1] / pow(ratio_nouni[calspce_y_p->flag_1 - 1], calspce_y_p->cnt_1 - 3)) + (int)RANGE0;
	}
	else
	{
		calspce_y_p->cnt_1 =(int) round((b[3] - min_y_p->coord_del_p[min_y_p->cnt]) / min_y_p->min_spac);
		calspce_y_p->numbs = calspce_y_p->numbs + calspce_y_p->cnt_1;
	}
	///////////////////////////////////////////////////

	//pow_flag_x用来储存中间的所有间隔的划分类型

	calspce_y_p->divtype_mid_p = (int *)malloc(min_y_p->cnt * sizeof(int));
	if (calspce_y_p->divtype_mid_p == NULL)
		printf("Failed to create memory of divtype_mid \n");

	calspce_y_p->lencnt = (int *)malloc((min_y_p->cnt + 1) * sizeof(int));
	if (calspce_y_p->lencnt == NULL)
		printf("Failed to create memory of divtype_mid \n");
	double len_y = 0;
	double lenmul_y = 0;
	double long_y = 0;
	calspce_y_p->max_lencnt = 0;
	for (int i = 0; i < min_y_p->cnt; i++)
	{
		if (((min_y_p->coord_del_p[i + 1] - min_y_p->coord_del_p[i]) > (RANGE1 - 0.01) * min_y_p->min_spac) && ((min_y_p->coord_del_p[i + 1] - min_y_p->coord_del_p[i]) < (RANGE2 + 0.01) * min_y_p->min_spac))
			calspce_y_p->divtype_mid_p[i] = 1;
		else if (((min_y_p->coord_del_p[i + 1] - min_y_p->coord_del_p[i]) > (RANGE2 - 0.01) * min_y_p->min_spac) && ((min_y_p->coord_del_p[i + 1] - min_y_p->coord_del_p[i]) < (RANGE3 + 0.01) * min_y_p->min_spac))
			calspce_y_p->divtype_mid_p[i] = 2;
		else if ((min_y_p->coord_del_p[i + 1] - min_y_p->coord_del_p[i]) > (RANGE3 - 0.01) * min_y_p->min_spac)
			calspce_y_p->divtype_mid_p[i] = 3;
		else if ((min_y_p->coord_del_p[i + 1] - min_y_p->coord_del_p[i]) < (RANGE2 - 0.01) * min_y_p->min_spac)
			calspce_y_p->divtype_mid_p[i] = 0;

		if (calspce_y_p->divtype_mid_p[i] > 0)
		{
			calspce_y_p->lencnt[i] = 3;
			len_y = pow(ratio_nouni[calspce_y_p->divtype_mid_p[i] - 1], 3);
			lenmul_y = pow(ratio_nouni[calspce_y_p->divtype_mid_p[i] - 1], 3);
			long_y = ((min_y_p->coord_del_p[i + 1] - min_y_p->coord_del_p[i] - 2 * RANGE0 * min_y_p->min_spac) * ratio_uni[calspce_y_p->divtype_mid_p[i] - 1] / min_y_p->min_spac) / 2;
			while ((len_y - 1) / (ratio_nouni[calspce_y_p->divtype_mid_p[i] - 1] - 1) < long_y)
			{
				len_y = len_y * lenmul_y;
				calspce_y_p->lencnt[i] = calspce_y_p->lencnt[i] + 3;
				if (calspce_y_p->max_lencnt < calspce_y_p->lencnt[i])
					calspce_y_p->max_lencnt = calspce_y_p->lencnt[i];
			}
			calspce_y_p->numbs = calspce_y_p->numbs + 2 * calspce_y_p->lencnt[i] + (int)round((2 * long_y * (1 - ratio_uni[calspce_y_p->divtype_mid_p[i] - 1]) / ratio_uni[calspce_y_p->divtype_mid_p[i] - 1]) / pow(ratio_nouni[calspce_y_p->divtype_mid_p[i] - 1], calspce_y_p->lencnt[i] - 3)) + 2 *(int) RANGE0;
		}
		else
		{
			calspce_y_p->lencnt[i] = (int)round((min_y_p->coord_del_p[i + 1] - min_y_p->coord_del_p[i]) / min_y_p->min_spac);

			calspce_y_p->numbs = calspce_y_p->numbs + calspce_y_p->lencnt[i];
		}
	}
}

void min_apace_x(int nn, double nets[], double b[], struct min_space *min_x_p)
{
	//struct min_space min_x_p;
	//设定根据每个独立的块来确定最小划分的初值
	double min_net_x = nets[2] - nets[0];
	//把块的边界都投影到x轴上，投并求出同一块内最小的长度
	double *coord_x = (double *)malloc(nn * 2 * sizeof(double));
	if (coord_x == NULL)
		printf("Failed to create memory of nets_coord \n");

	min_x_p->cnt = 0;
	//int cnt_x = 0;
	double temp_min = 0;
	for (int i = 0; i < nn; i++)
	{
		coord_x[min_x_p->cnt] = nets[i * 4 + 0];
		min_x_p->cnt++;
		coord_x[min_x_p->cnt] = nets[i * 4 + 2];
		min_x_p->cnt++;
		temp_min = nets[i * 4 + 2] - nets[i * 4 + 0];
		if (temp_min < min_net_x)
			min_net_x = temp_min;
	}
	//对投影的坐标排序
	quick_sort(coord_x, 0, min_x_p->cnt - 1);
	//删除重复的坐标看过法律

	min_x_p->coord_del_p = (double *)malloc(nn * 2 * sizeof(double));

	if (min_x_p->coord_del_p == NULL)
		printf("Failed to create memory of nets_coord \n");
	min_x_p->cnt = 0;
	min_x_p->coord_del_p[0] = coord_x[0];
	for (int i = 1; i < 2 * nn; i++)
	{
		if (coord_x[i] != min_x_p->coord_del_p[min_x_p->cnt])
		{
			min_x_p->coord_del_p[min_x_p->cnt + 1] = coord_x[i];
			min_x_p->cnt++;
		}
	}
	free(coord_x);
	//通过比较，获得所有相邻坐标的最小差值，用这个差值除div得到最小划分间隔
	//min adjacent coordinate
	double min_x = (min_x_p->coord_del_p[1] - min_x_p->coord_del_p[0]) / DIV;
	temp_min = 0;
	for (int i = 1; i < min_x_p->cnt - 1; i++)
	{
		temp_min = (min_x_p->coord_del_p[i + 1] - min_x_p->coord_del_p[i]) / DIV;
		if (temp_min < min_x)
			min_x = temp_min;
	}
	//之前还得得到了块内坐标的最小差值block_div_x_min和block_div_y_min，用这个最小差值/8和x_div_min，y_div_min与比较
	if (min_net_x / 8 < min_x)
		min_x = min_net_x / 8;
	if (min_x_p->coord_del_p[0] - b[0] < min_x)
		min_x = (min_x_p->coord_del_p[0] - b[0]) / DIV;
	if (b[2] - min_x_p->coord_del_p[min_x_p->cnt - 2] < min_x)
		//-2可能有问题
		min_x = (b[2] - min_x_p->coord_del_p[min_x_p->cnt - 2] - b[0]) / DIV;

	min_x_p->min_spac = min_x;
}

void min_apace_y(int nn, double nets[], double b[], struct min_space *min_y_p)
{
	//struct min_space min_y_p;
	//设定根据每个独立的块来确定最小划分的初值
	double min_net_y = nets[3] - nets[1];
	//把块的边界都投影到x轴上，投并求出同一块内最小的长度
	double *coord_y = (double *)malloc(nn * 2 * sizeof(double));
	if (coord_y == NULL)
		printf("Failed to create memory of nets_coord \n");

	min_y_p->cnt = 0;
	//int cnt_x = 0;
	double temp_min = 0;
	for (int i = 0; i < nn; i++)
	{
		coord_y[min_y_p->cnt] = nets[i * 4 + 1];
		min_y_p->cnt++;
		coord_y[min_y_p->cnt] = nets[i * 4 + 3];
		min_y_p->cnt++;
		temp_min = nets[i * 4 + 3] - nets[i * 4 + 1];
		if (temp_min < min_net_y)
			min_net_y = temp_min;
	}
	//对投影的坐标排序
	quick_sort(coord_y, 0, min_y_p->cnt - 1);

	min_y_p->coord_del_p = (double *)malloc((nn * 2+2) * sizeof(double));

	if (min_y_p->coord_del_p == NULL)
		printf("Failed to create memory of nets_coord \n");
	min_y_p->cnt = 0;
	min_y_p->coord_del_p[0] = coord_y[0];
	for (int i = 1; i < 2 * nn; i++)
	{
		if (coord_y[i] != min_y_p->coord_del_p[min_y_p->cnt])
		{
			min_y_p->coord_del_p[min_y_p->cnt + 1] = coord_y[i];
			min_y_p->cnt++;
		}
	}
	free(coord_y);
	//通过比较，获得所有相邻坐标的最小差值，用这个差值除div得到最小划分间隔
	//min adjacent coordinate
	double min_y = (min_y_p->coord_del_p[1] - min_y_p->coord_del_p[0]) / DIV;
	temp_min = 0;
	for (int i = 1; i < min_y_p->cnt - 1; i++)
	{
		temp_min = (min_y_p->coord_del_p[i + 1] - min_y_p->coord_del_p[i]) / DIV;
		if (temp_min < min_y)
			min_y = temp_min;
	}
	//之前还得得到了块内坐标的最小差值block_div_x_min和block_div_y_min，用这个最小差值/8和x_div_min，y_div_min与比较
	if (min_net_y / 8 < min_y)
		min_y = min_net_y / 8;
	if (min_y_p->coord_del_p[0] - b[1] < min_y)
		min_y = (min_y_p->coord_del_p[0] - b[1]) / DIV;
	if (b[3] - min_y_p->coord_del_p[min_y_p->cnt - 2] < min_y)
		//-2可能有问题
		min_y = (b[3] - min_y_p->coord_del_p[min_y_p->cnt - 2] - b[1]) / DIV;

	min_y_p->min_spac = min_y;
}

//Calculate the coordinate
double **caluculate_coord(double nets[], int numbs_ns, int shps[])
{

	int cnt_rect = 0;
	double **nets_coord = (double **)malloc(numbs_ns * sizeof(double));
	if (nets_coord == NULL)
		printf("Failed to create memory of nets_coord \n");
	else
	{
		for (int i = 0; i < numbs_ns; ++i)
		{
			nets_coord[i] = (double *)malloc(4 * sizeof(double));
			if (nets_coord[i] == NULL)
				printf("Failed to create memory of nets_coord[i]\n");
		}
	}

	for (int i = 0; i < numbs_ns; i++)
	{
		if (shps[i] == 1)
		{
			nets_coord[i][0] = nets[cnt_rect * 4 + 0];
			nets_coord[i][1] = nets[cnt_rect * 4 + 1];
			nets_coord[i][2] = nets[cnt_rect * 4 + 2];
			nets_coord[i][3] = nets[cnt_rect * 4 + 3];
			cnt_rect++;
		}
		else
		{
			nets_coord[i][0] = compare_min(nets[cnt_rect * 4 + 0], nets[cnt_rect * 4 + 4], nets[cnt_rect * 4 + 8], nets[cnt_rect * 4 + 12]);
			nets_coord[i][1] = nets[cnt_rect * 4 + 1];
			nets_coord[i][2] = compare_max(nets[cnt_rect * 4 + 2], nets[cnt_rect * 4 + 6], nets[cnt_rect * 4 + 10], nets[cnt_rect * 4 + 14]);
			nets_coord[i][3] = nets[cnt_rect * 4 + 15];
			cnt_rect = cnt_rect + 4;
		}
	}

	return nets_coord;
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
