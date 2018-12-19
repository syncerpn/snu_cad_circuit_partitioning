#ifndef __PARTITION_ALG_H__
#define	__PARTITION_ALG_H__

#ifdef __cplusplus
extern "C" {
#endif
	//helper
	float get_area(int cell_id);
	int cal_cut_size(int net_index, int* part_a);
	void print_partitioning_result(const char* method_name, int cut_size, int* part_mask);

	//main
    int fm_partitioning(float r, float t, int* final_part_mask); //fm
    int kl_partitioning(int* final_part_mask); //kl
    int sa_partitioning(float r, float t, float T, float alpha, int max_time, int M, int* final_part_mask); //sa

#ifdef __cplusplus
}
#endif
#endif