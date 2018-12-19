/****************************************************************************
*   Title : util.c
*   Desc  : 
*   Author: Jinhwan Jeon
*   Date  : 1995.2.13, 1995.2.15
*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <stdarg.h>

void errorMsg(const char *format, ...)
{
    va_list ap;
    char buf[1024];

    va_start(ap, format);
    vsprintf(buf, format, ap);
    fprintf(stderr, "%s\n", buf);
    va_end(ap);
    exit(-1);
}

FILE *file_open(const char *name, const char *mode)
{
    FILE *fp;
    fp = fopen(name, mode);
    if (!fp) {
	errorMsg("File open error %s", name);
    }
    return fp;
}

void *mem_alloc(int len)
{
    void *ptr;
    ptr = malloc(len);
    if (!ptr) {
	errorMsg("Memory allocation error(in malloc(%d))", len);
    }
    return ptr;
}

void *mem_calloc(int len, int num)
{
    void *ptr;
    ptr = calloc(len, num);
    if (!ptr) {
	errorMsg("Memory allocation error(in calloc(%d, %d))", len, num);
    }
    return ptr;
}

int randomize()
{
    struct timeval tv;
    struct timezone tz;
    int seed;

    gettimeofday(&tv, &tz);

    seed = tv.tv_usec;

    srand(seed);
    return seed;
}

double get_ptime(double *utime, double *stime)
{
    struct rusage rusage;
    double user, system;

    getrusage(RUSAGE_SELF, &rusage);
    user = rusage.ru_utime.tv_sec * 1.0e6 + rusage.ru_utime.tv_usec;
    system = rusage.ru_stime.tv_sec * 1.0e6 + rusage.ru_stime.tv_usec;

    if (utime)
	*utime = user;
    if (stime)
	*stime = system;

    return user + system;
}

double get_utime()
{
    double utime;
    double stime;

    get_ptime(&utime, &stime);

    return utime;
}

double get_stime()
{
    double utime;
    double stime;

    get_ptime(&utime, &stime);

    return stime;
}

int sign(int val)
{
    return val >= 0 ? 1 : -1;
}

void rand_permute(int* arr, int n) {
    int i, x, y;
    int t;
    for (i = 0; i < n; ++i) {
        x = rand() % n;
        y = rand() % n;
        t = arr[x];
        arr[x] = arr[y];
        arr[y] = t;
    }
}

void bucket_sort(int* arr, int n, int p_max, int* sorted_index_array) {
	int i;
	int n_bucket = 2*p_max + 1;
	int* count = mem_calloc(sizeof(int), n_bucket);
	for (i = 0; i < n_bucket; ++i) count[i] = 0;
	for (i = 0; i < n; ++i) ++count[arr[i] + p_max];

	int* placement = mem_calloc(sizeof(int), n_bucket);	
	placement[0] = 0;
	for (i = 1; i < n_bucket; ++i) placement[i] = placement[i-1] + count[i-1];
	for (i = 0; i < n; ++i) {
		int bucket_id = arr[i] + p_max;
		int place = placement[bucket_id];
		sorted_index_array[place] = i;
		++placement[bucket_id];
	}
	free(count);
	free(placement);
}

void merge(int* arr, int l, int m, int r, int* sorted_index) {
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;
    int* L = mem_calloc(sizeof(int), n1);
    int* R = mem_calloc(sizeof(int), n2);
    int* iL = mem_calloc(sizeof(int), n1);
    int* iR = mem_calloc(sizeof(int), n2);

    for (i = 0; i < n1; ++i) {
        L[i] = arr[l + i];
        iL[i] = sorted_index[l + i];
    }
    
    for (j = 0; j < n2; ++j) {
        R[j] = arr[m + 1 + j];
        iR[j] = sorted_index[m + 1 + j];
    }

    i = 0;
    j = 0;
    k = l;
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            arr[k] = L[i];
            sorted_index[k] = iL[i];
            ++i;
        } else {
            arr[k] = R[j];
            sorted_index[k] = iR[j];
            ++j;
        }
        ++k;
    }

    while (i < n1) {
        arr[k] = L[i];
        sorted_index[k] = iL[i];
        ++i;
        ++k;
    }

    while (j < n2) {
        arr[k] = R[j];
        sorted_index[k] = iR[j];
        ++j;
        ++k;
    }
    free(L);
    free(R);
    free(iL);
    free(iR);
}

void merge_sort(int* arr, int l, int r, int* sorted_index) {
    if (l < r) {
        int m = l + (r - l)/2;
        merge_sort(arr, l, m, sorted_index);
        merge_sort(arr, m + 1, r, sorted_index);
        merge(arr, l, m, r, sorted_index);
    }
}

void merge_float(float* arr, int l, int m, int r, int* sorted_index) {
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;
    float* L = mem_calloc(sizeof(int), n1);
    float* R = mem_calloc(sizeof(int), n2);
    int* iL = mem_calloc(sizeof(int), n1);
    int* iR = mem_calloc(sizeof(int), n2);

    for (i = 0; i < n1; ++i) {
        L[i] = arr[l + i];
        iL[i] = sorted_index[l + i];
    }
    
    for (j = 0; j < n2; ++j) {
        R[j] = arr[m + 1 + j];
        iR[j] = sorted_index[m + 1 + j];
    }

    i = 0;
    j = 0;
    k = l;
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            arr[k] = L[i];
            sorted_index[k] = iL[i];
            ++i;
        } else {
            arr[k] = R[j];
            sorted_index[k] = iR[j];
            ++j;
        }
        ++k;
    }

    while (i < n1) {
        arr[k] = L[i];
        sorted_index[k] = iL[i];
        ++i;
        ++k;
    }

    while (j < n2) {
        arr[k] = R[j];
        sorted_index[k] = iR[j];
        ++j;
        ++k;
    }
    free(L);
    free(R);
    free(iL);
    free(iR);
}

void merge_sort_float(float* arr, int l, int r, int* sorted_index) {
    if (l < r) {
        int m = l + (r - l)/2;
        merge_sort_float(arr, l, m, sorted_index);
        merge_sort_float(arr, m + 1, r, sorted_index);
        merge_float(arr, l, m, r, sorted_index);
    }
}