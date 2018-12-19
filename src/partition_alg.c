#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util.h"
#include "standard.h"
#include "partition_alg.h"

//====================  HELPER  ============================
float get_area(int cell_id) {
    return (float)(cellArray[cell_id].width * cellArray[cell_id].height);
}

int cal_cut_size(int net_index, int* part_a) {
    int i;
    NET *netptr;
    CELL *cellptr;
    netptr = netArray + net_index;
    int count_a = 0;
    int count_b = 0;
    
    for (i = 0; i < netptr->cellcount2; ++i) {
        cellptr = cellArray + netptr->cell2[i];
        part_a[cellptr->id] == 1 ? ++count_a : ++count_b;
    }

    return ((count_a > 0) && (count_b > 0));
}

void print_partitioning_result(const char* method_name, int cut_size, int* part_mask) {
    int i;
    printf("== %s ==\n", method_name);
    printf("Cut Size  : %d\n", cut_size);
    float area_a = 0;
    float area_b = 0;
    float total_size = stat.avg_cell_size * cellNum;

    for (i = 0; i < cellNum; ++i) {
        // printf("%d ", part_mask[i]);
        if (part_mask[i] == 1) area_a += get_area(i);
        else area_b += get_area(i);
    }
    printf("Area Ratio: %.3f : %.3f\n", area_a/total_size, area_b/total_size);
    printf("\n");
    return;
}

//======================  FM  ==============================
int fm_update_gain(int net_index, int* part_a, int* gain) {
    int i;
    NET *netptr;
    CELL *cellptr;
    netptr = netArray + net_index;
    int count_a = 0;
    int count_b = 0;
    
    for (i = 0; i < netptr->cellcount2; ++i) {
        cellptr = cellArray + netptr->cell2[i];
        part_a[cellptr->id] == 1 ? ++count_a : ++count_b;
    }

    if (((count_a >= 2) && (count_b == 0)) || ((count_b >= 2) && (count_a == 0))) {
        for (i = 0; i < netptr->cellcount2; ++i) {
            cellptr = cellArray + netptr->cell2[i];
            --gain[cellptr->id];
        }
    } else if ((count_a == 1) && (count_b == 1)) {
        for (i = 0; i < netptr->cellcount2; ++i) {
            cellptr = cellArray + netptr->cell2[i];
            ++gain[cellptr->id];
        }
    } else if (((count_a >= 2) && (count_b == 1)) || ((count_b >= 2) && (count_a == 1))) {
        for (i = 0; i < netptr->cellcount2; ++i) {
            cellptr = cellArray + netptr->cell2[i];
            if (((count_b == 1) && (part_a[cellptr->id] == -1)) || ((count_a == 1) && (part_a[cellptr->id] == 1))) {
                ++gain[cellptr->id];
                break;
            }
        }
    }
    return ((count_a > 0) && (count_b > 0));
}

int fm_partitioning(float r, float t, int* final_part_mask) {
    int i;

    //for random permutation at init phase
    int* index_array  = mem_calloc(sizeof(int), cellNum);

    //mask of partition: 1 and -1 are equivalent to part a and b respectively
    int* part_a       = mem_calloc(sizeof(int), cellNum);

    //area of part a
    float part_a_size = 0;

    //cell gain
    int* gain         = mem_calloc(sizeof(int), cellNum);
    for (i = 0; i < cellNum; ++i) gain[i] = 0;

    //cell gain after sorted
    int* sorted_gain  = mem_calloc(sizeof(int), cellNum);

    //lock mask
    int* locked       = mem_calloc(sizeof(int), cellNum);
    for (i = 0; i < cellNum; ++i) locked[i] = 0;

    //cut size after each iteration, cut_set[0] = init cut size
    int* cut_set      = mem_calloc(sizeof(int), cellNum + 1);
    for (i = 0; i < cellNum + 1; ++i) cut_set[i] = 0;

    //index of cell being moved after each iteration
    int* move_id      = mem_calloc(sizeof(int), cellNum);
    for (i = 0; i < cellNum; ++i) move_id[i] = -1;

    //total size of circuit
    float total_size = stat.avg_cell_size * cellNum;

    //lower and upper bounds of area constraint
    float lb_size = (r - t) * total_size;
    float ub_size = (r + t) * total_size;

    //init phase: randomize partitioning
    do {
        part_a_size = 0;
        for (i = 0; i < cellNum; ++i) index_array[i] = i;
        rand_permute(index_array, cellNum);

        for (i = 0; i < cellNum; ++i) {
            if ((part_a_size >= lb_size) && (part_a_size <= ub_size)) {
                part_a[index_array[i]] = -1;
            } else {
                float tmp_size = get_area(index_array[i]);
                tmp_size += part_a_size;
                if (tmp_size <= ub_size) {
                    part_a_size = tmp_size;
                    part_a[index_array[i]] = 1;
                } else {
                    part_a[index_array[i]] = -1;
                }
            }
        }
    } while ((part_a_size < lb_size) || (part_a_size > ub_size)); //verify the initialization
    free(index_array); //index_array is no longer needed

    //store the init state of partition
    memcpy(final_part_mask, part_a, cellNum * sizeof(int));

    //update
    int total_locked = 0;
    int min_cut_size = netNum + 1;
    int min_cut_size_id = -1;
    float best_size = part_a_size;

    do {
        for (i = 0; i < cellNum; ++i) gain[i] = 0;
        for (i = 0; i < netNum; ++i) cut_set[total_locked] += fm_update_gain(i, part_a, gain);

        //last update when total_locked == cellNum is not necessary since its cut size should be equal to init cut size
        if ((min_cut_size > cut_set[total_locked]) || ((min_cut_size == cut_set[total_locked]) && (abs(best_size/total_size - r) < abs(part_a_size/total_size - r)))) {
            min_cut_size = cut_set[total_locked];
            min_cut_size_id = total_locked - 1;
            best_size = part_a_size;
        }

        bucket_sort(gain, cellNum, stat.max_net_per_cell, sorted_gain);

        for (i = cellNum - 1; i >= 0; --i) {
            int id = sorted_gain[i];
            if (locked[id]) continue; //locked cell
            float tmp_size = part_a_size - part_a[id] * cellArray[id].width * cellArray[id].height;
            if ((tmp_size < lb_size) || (tmp_size > ub_size)) continue; //imbalance

            part_a_size = tmp_size;
            if (part_a_size/total_size < 0.4 || part_a_size/total_size > 0.6) printf("%.1f\n", part_a_size/total_size);
            locked[id] = 1;
            part_a[id] = -part_a[id];
            move_id[total_locked] = id;
            break;
        }
        if (i == -1) break; //none is chosen
        ++total_locked;
    } while (total_locked < cellNum);

    //reveal best partitioning based on cut size
    for (i = 0; i <= min_cut_size_id; ++i) final_part_mask[move_id[i]] = -final_part_mask[move_id[i]];

    //finish
    free(part_a);
    free(gain);
    free(sorted_gain);
    free(move_id);
    free(cut_set);
    free(locked);
    return min_cut_size;
}

//======================  KL  ==============================
int kl_update_gain(int cell_id, int* part_a, int* D) {
    int i, j, partner_cell_id;
    int cut = 0;
    
    D[cell_id] = 0;

    NET *netptr;
    for (i = 0; i < cellArray[cell_id].netcount2; ++i) {
    	netptr = netArray + cellArray[cell_id].net2[i];
    	for (j = 0; j < netptr->cellcount2; ++j) {
    		partner_cell_id = netptr->cell2[j];
    		if (cell_id == partner_cell_id) continue;
    		cut += (part_a[cell_id] != part_a[partner_cell_id]);
    		D[cell_id] += (part_a[cell_id] == part_a[partner_cell_id] ? -1 : 1);
    	}
    }
    return cut;
}

int kl_cal_c(int cell_i, int cell_j) {
    int i, j;
    int c_ij = 0;
    
    NET *netptr;
    for (i = 0; i < cellArray[cell_i].netcount2; ++i) {
    	netptr = netArray + cellArray[cell_i].net2[i];
    	for (j = 0; j < netptr->cellcount2; ++j) {
    		if (cell_j == netptr->cell2[j]) {
    			++c_ij;
    			break;
    		}
    	}
    }
    return c_ij;
}

int kl_find_maximal_gain(int* D, int* ref_part_a) {
	int i, j;

	int* part_a = mem_calloc(sizeof(int), cellNum);
	memcpy(part_a, ref_part_a, cellNum * sizeof(int));

	int* Ds_id = mem_calloc(sizeof(int), cellNum);
	for (i = 0; i < cellNum; ++i) Ds_id[i] = i;

	int* Ds = mem_calloc(sizeof(int), cellNum);
	memcpy(Ds, D, cellNum * sizeof(int));

	merge_sort(Ds, 0, cellNum - 1, Ds_id);

    int* locked = mem_calloc(sizeof(int), cellNum);
	for (i = 0; i < cellNum; ++i) locked[i] = 0;

    int n_gain = (cellNum + 1)/2;

    int* gain = mem_calloc(sizeof(int), n_gain);
    int* sorted_gain_id = mem_calloc(sizeof(int), n_gain);
	for (i = 0; i < n_gain; ++i) sorted_gain_id[i] = i;

    int* move_a = mem_calloc(sizeof(int), n_gain);
    int* move_b = mem_calloc(sizeof(int), n_gain);
	for (i = 0; i < n_gain; ++i) {
		move_a[i] = -1;
		move_b[i] = -1;
	}

    int g_best, g_best_compare;
    int ia_best, ib_best;
    int ia, ib, ic;
    int first;
    int no_more;

    for (i = 0; i < n_gain-1; ++i) {
    	ia_best = cellNum - 1;
    	ib_best = cellNum - 1;
    	first = 1;
    	no_more = 0;
    	while (1) {
    		if (first) {
		    	while (locked[Ds_id[ia_best]] == 1 || part_a[Ds_id[ia_best]] == -1) --ia_best;
		    	while (locked[Ds_id[ib_best]] == 1 || part_a[Ds_id[ib_best]] == 1) --ib_best;

		    	g_best = Ds[ia_best] + Ds[ib_best] - 2 * kl_cal_c(Ds_id[ia_best], Ds_id[ib_best]);

			    ia = ia_best;
			    ib = ib_best;
	    		ic = (ia > ib) ? ia : ib;

		    	first = 0;
		    }

		    ia = ia_best;
		    ib = ib_best;

		    while (ia == ia_best && ib == ib_best) {
		    	--ic;
		    	if (ic < 0) {
		    		no_more = 1;
		    		break;
		    	}
		    	if (locked[Ds_id[ic]] == 0) {
		    		if (part_a[Ds_id[ic]] == 1) ia = ic;
		    		else ib = ic;
		    	}
		    }

		    if (no_more) break;

		    g_best_compare = Ds[ia] + Ds[ib];

		    if (g_best >= g_best_compare) break; //no more

		    g_best_compare -= 2 * kl_cal_c(Ds_id[ia], Ds_id[ib]);

		    if (g_best < g_best_compare) {
		    	ia_best = ia;
		    	ib_best = ib;
		    	g_best = g_best_compare;
		    }
		}

    	//swap and lock
    	locked[Ds_id[ia_best]] = 1;
    	locked[Ds_id[ib_best]] = 1;
    	part_a[Ds_id[ia_best]] = -part_a[Ds_id[ia_best]];
    	part_a[Ds_id[ib_best]] = -part_a[Ds_id[ib_best]];

    	//update move and gain
    	move_a[i] = Ds_id[ia_best];
    	move_b[i] = Ds_id[ib_best];
    	gain[i] = g_best;

    	for (j = 0; j < cellNum; ++j) {
    		if (locked[j] == 0) {
				// kl_update_gain(j, part_a, D); //simple update, but slower
                //quick update
                int ca = (part_a[j] == part_a[move_a[i]]) ? -1 : 1;
                int cb = (part_a[j] == part_a[move_b[i]]) ? -1 : 1;
                D[j] = D[j] + 2 * ca * kl_cal_c(j,move_a[i]) + 2 * cb * kl_cal_c(j,move_b[i]);
    		}
    	}
		
		memcpy(Ds, D, cellNum * sizeof(int));
		for (j = 0; j < cellNum; ++j) Ds_id[j] = j;
		merge_sort(Ds, 0, cellNum - 1, Ds_id);
    }
    //handle the last swap in case there is an odd number of cells
    if (cellNum % 2 == 1) {
    	for (i = 0; i < cellNum; ++i) {
    		if (locked[i] == 1) continue;
			gain[n_gain - 1] = D[i];
			locked[i] = 1;
    		if (part_a[i] == 1) move_a[n_gain - 1] = i;
    		else move_b[n_gain - 1] = i;
    		part_a[i] = -part_a[i];
    		break;
    	}
    } else {
    	for (i = 0; i < cellNum; ++i) {
    		if (locked[i] == 1) continue;
    		if (part_a[i] == 1) ia = i;
    		else ib = i;
    	}
    	gain[n_gain - 1] = D[ia] + D[ib] - 2 * kl_cal_c(ia, ib);
		move_a[n_gain - 1] = ia;
		move_b[n_gain - 1] = ib;
		locked[ia] = 1;
		locked[ib] = 1;
		part_a[ia] = -part_a[ia];
		part_a[ib] = -part_a[ib];
    }

    for (i = 1; i < n_gain; ++i) gain[i] += gain[i-1];
    merge_sort(gain, 0, n_gain - 1, sorted_gain_id);

    if (gain[n_gain - 1] > 0) {
	    for (i = 0; i <= sorted_gain_id[n_gain-1]; ++i) {
	    	if (move_a[i] >= 0) ref_part_a[move_a[i]] = -ref_part_a[move_a[i]];
	    	if (move_b[i] >= 0) ref_part_a[move_b[i]] = -ref_part_a[move_b[i]];
	    }
	}

	int G_k_best = gain[n_gain - 1]; //store into intermediate to free malloc

    //finish
    free(part_a);
	free(Ds_id);
	free(Ds);
	free(locked);
	free(gain);
	free(sorted_gain_id);
	free(move_a);
	free(move_b);

    return G_k_best;
}

int kl_partitioning(int* final_part_mask) {
    int i;

    //for random permutation at init phase
    int* index_array  = mem_calloc(sizeof(int), cellNum);

    //mask of partition: 1 and -1 are equivalent to part a and b respectively
    int* part_a       = mem_calloc(sizeof(int), cellNum);

    //init phase: randomize partitioning
    for (i = 0; i < cellNum; ++i) index_array[i] = i;
    rand_permute(index_array, cellNum);

	for (i = 0; i < cellNum; ++i) {
		if (i < cellNum/2) part_a[index_array[i]] = 1;
		else part_a[index_array[i]] = -1;
	}

    free(index_array); //index_array is no longer needed

	int* D = mem_calloc(sizeof(int), cellNum);
    int G_k_best;

    do {
		for (i = 0; i < cellNum; ++i) kl_update_gain(i, part_a, D);
		G_k_best = kl_find_maximal_gain(D, part_a);
	} while (G_k_best > 0);
	
	int min_cut_size = 0;

	memcpy(final_part_mask, part_a, cellNum * sizeof(int));
    //we recalculate and return the net cut size instead of graph cut size
    for (i = 0; i < netNum; ++i) {
        min_cut_size += cal_cut_size(i, part_a);
    }

    //finish
	free(part_a);
    free(D);
    
    return min_cut_size;
}

//======================  SA  ==============================
int sa_accept(int cut_size, int cut_size_new, float T) {
    float dE = (float)(cut_size_new - cut_size);
    float r = (float)(rand() % RAND_MAX) / (RAND_MAX - 1);
    if (r < exp(-dE/T)) return 1;
    return 0;
}

int sa_partitioning(float r, float t, float T, float alpha, int max_time, int M, int* final_part_mask) {
    int i, j, k, c, ptr;

    T = T > 0 ? T : 10000; //temperature
    //check validity of alpha 0 < alpha < 1
    alpha = alpha < 1 ? (alpha > 0 ? alpha : 0.01) : 0.99;

    // printf("Initial temperature: %.3f\n", T);
    // printf("Cooling schedule alpha: %.3f\n", alpha);
    
    //for random permutation at init phase
    int* index_array  = mem_calloc(sizeof(int), cellNum);

    //mask of partition: 1 and -1 are equivalent to part a and b respectively
    int* part_a       = mem_calloc(sizeof(int), cellNum);

    //area of part a
    float part_a_size = 0;

    //total size of circuit
    float total_size = stat.avg_cell_size * cellNum;

    //lower and upper bounds of area constraint
    float lb_size = (r - t) * total_size;
    float ub_size = (r + t) * total_size;

    //init phase: randomize partitioning
    do {
        part_a_size = 0;
        for (i = 0; i < cellNum; ++i) index_array[i] = i;
        rand_permute(index_array, cellNum);

        for (i = 0; i < cellNum; ++i) {
            if ((part_a_size >= lb_size) && (part_a_size <= ub_size)) {
                part_a[index_array[i]] = -1;
            } else {
                float tmp_size = get_area(index_array[i]);
                tmp_size += part_a_size;
                if (tmp_size <= ub_size) {
                    part_a_size = tmp_size;
                    part_a[index_array[i]] = 1;
                } else {
                    part_a[index_array[i]] = -1;
                }
            }
        }
    } while ((part_a_size < lb_size) || (part_a_size > ub_size)); //verify the initialization
    free(index_array); //index_array is no longer needed

    int cut_size = 0;
    for (i = 0; i < netNum; ++i) cut_size += cal_cut_size(i, part_a);

    int cut_size_new = cut_size;
    int no_update = 0;

    float* all_size = mem_calloc(sizeof(float), cellNum);
    int* id = mem_calloc(sizeof(int), cellNum);

    for (i = 0; i < cellNum; ++i) {
        id[i] = i;
        all_size[i] = get_area(i);
    }

    merge_sort_float(all_size, 0, cellNum - 1, id);

    int time = 0;
    int single_move;

    do {
        single_move = 0;
        i = rand() % cellNum;
        float tmp_part_a_size = part_a_size - part_a[id[i]] * all_size[i];
        //if this cell moves to the other partition without breaking the area constraint, we simply approve it.
        if ((lb_size <= tmp_part_a_size) && (ub_size >= tmp_part_a_size)) single_move = 1;

        //calculate the lower and upper bounds of area of the partner cell
        float lb_b_size = lb_size - tmp_part_a_size;
        float ub_b_size = ub_size - tmp_part_a_size;

        c = 0;
        ptr = -1;

        for (k = 0; k < cellNum; ++k) {
            if (part_a[id[k]] != part_a[id[i]]) {
                //only accept the partner cell with which the swapping does not violate the area constraint
                if ((- part_a[id[k]] * all_size[k] <= ub_b_size) && (- part_a[id[k]] * all_size[k] >= lb_b_size)) {
                    if (ptr == -1) ptr = k;
                    ++c;
                } else if (ptr != -1) {
                    break;
                }
            }
        }

        c += single_move;

        if (c > 0) {
            //choose the partner cell randomly
            j = rand() % c;
            //allow a move to be made with only one cell, if it does not violate the constraint
            if (single_move && (j == c - 1)) {
                part_a[id[i]] = -part_a[id[i]];

                cut_size_new = 0;
                for (k = 0; k < netNum; ++k) {
                    cut_size_new += cal_cut_size(k, part_a);
                }

                //validate move
                if (sa_accept(cut_size, cut_size_new, T)) {
                    cut_size = cut_size_new;
                    no_update = 0;
                    part_a_size = part_a_size + part_a[id[i]] * all_size[i];
                    
                } else {
                    //if not accepted, swap back
                    part_a[id[i]] = -part_a[id[i]];
                }
            } else {
                while (j > 0) {
                    ++ptr;
                    if (part_a[id[ptr]] != part_a[id[i]]) --j;
                }
                j = ptr;
                part_a[id[i]] = -part_a[id[i]];
                part_a[id[j]] = -part_a[id[j]];

                cut_size_new = 0;
                for (k = 0; k < netNum; ++k) {
                    cut_size_new += cal_cut_size(k, part_a);
                }

                //validate move
                if (sa_accept(cut_size, cut_size_new, T)) {
                    cut_size = cut_size_new;
                    no_update = 0;
                    part_a_size = part_a_size + part_a[id[i]] * all_size[i] + part_a[id[j]] * all_size[j];

                } else {
                    //if not accepted, swap back
                    part_a[id[i]] = -part_a[id[i]];
                    part_a[id[j]] = -part_a[id[j]];
                }
            }
        }
        ++time;

        //if there is no update at this temperature, we stop the whole process
        if ((time % M) == 0) {
            if (no_update) {
                break;
            } else {
                no_update = 1;
                T *= alpha;
            }
        }
    } while (time < max_time);

    memcpy(final_part_mask, part_a, cellNum * sizeof(int));

    free(part_a);
    free(all_size);
    free(id);

    return cut_size;
}