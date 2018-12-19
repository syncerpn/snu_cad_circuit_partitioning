/****************************************************************************
*   Title : 
*   Desc  : 
*   Author: 
*   Date  : 
*****************************************************************************/

#include <stdio.h>
#include "standard.h"
#include "parser.h"
#include "util.h"
//nghiant
#include <stdlib.h>
#include <string.h>
#include "partition_alg.h"
//nghiant_end

static void part(char *file)
{
    inputParse(file);

    /* test print, you may delete this line */
    // printCircuit();

	/*********************************
	* insert your partition code here 
	**********************************/
    randomize();
    
    //settings
    float partition_ratio = 0.5;
    float tolerance = 0.1;
    
    //SA settings
    float T = 10000;
    float alpha = 0.5; //cooling schedule factor
    int M = cellNum / 4; //number of moves at each temperature
    M = M > 0 ? M : 1;
    int max_time = M * 100; //max number of moves

    //running time
    int repeat = 10;

    int i;
    int* partition_mask = mem_calloc(sizeof(int), cellNum);
    int* tmp_partition_mask = mem_calloc(sizeof(int), cellNum);
    int cut_size = netNum + 1;
    float avg_size = 0;
    //FM Partitioning
    printf("\n");
    for (i = 0; i < repeat; ++i) {
        printf("\033[1A");
        printf("\033[K");
    	printf("Passed: %d/%d\n", i+1, repeat);
    	int tmp = fm_partitioning(partition_ratio, tolerance, tmp_partition_mask);
        avg_size += (float)tmp;
    	if (tmp < cut_size) {
    		memcpy(partition_mask, tmp_partition_mask, cellNum * sizeof(int));
    		cut_size = tmp;
    	}
    }
    print_partitioning_result("FM Algorithm", cut_size, partition_mask);
    printf("Average Cut Size: %.2f\n", avg_size/repeat);
    printf("==================\n\n");
	
	//KL Partitioning
	cut_size = netNum + 1;
    avg_size = 0;
    printf("\n");
    for (i = 0; i < repeat; ++i) {
        printf("\033[1A");
        printf("\033[K");
    	printf("Passed: %d/%d\n", i+1, repeat);
    	int tmp = kl_partitioning(tmp_partition_mask);
        avg_size += (float)tmp;
    	if (tmp < cut_size) {
    		memcpy(partition_mask, tmp_partition_mask, cellNum * sizeof(int));
    		cut_size = tmp;
    	}
    }
    print_partitioning_result("KL Algorithm", cut_size, partition_mask);
    printf("Average Cut Size: %.2f\n", avg_size/repeat);
    printf("==================\n\n");

    //SA Partitioning
    cut_size = netNum + 1;
    avg_size = 0;
    printf("\n");
    for (i = 0; i < repeat; ++i) {
        printf("\033[1A");
        printf("\033[K");
    	printf("Passed: %d/%d\n", i+1, repeat);
    	int tmp = sa_partitioning(partition_ratio, tolerance, T, alpha, max_time, M, tmp_partition_mask);
        avg_size += (float)tmp;
    	if (tmp < cut_size) {
    		memcpy(partition_mask, tmp_partition_mask, cellNum * sizeof(int));
    		cut_size = tmp;
    	}
    }
    print_partitioning_result("SA Algorithm", cut_size, partition_mask);
    printf("Average Cut Size: %.2f\n", avg_size/repeat);
    printf("==================\n\n");

    return;

    /*********************************
    * custom code ends here
    **********************************/
}

int main(int argc, char **argv)
{
    if (argc < 2)
	errorMsg("usage: part input");

    part(argv[1]);
    return 0;
}
