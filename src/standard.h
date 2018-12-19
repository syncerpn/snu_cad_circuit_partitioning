/****************************************************************************
*	Title : standard.h
*	Desc  : 
*	Author: Jeon Jinhwan
*	Date  : 1995.2.13
*****************************************************************************/

#ifndef __STANDARD_H__
#define __STANDARD_H__

#include <stdio.h>

#define		MAXNAMELEN		30

/* module type */
#define		STANDARD		1
#define		PAD				2
#define		GENERAL			3
#define		FEEDTHROUGH		5

/* cell orientation */
#define		ROT0			1
#define		ROT90			2
#define		ROT180			3
#define		ROT270			4
#define		YREFLECTION		0x100
#define		XREFLECTION 	0x200

/* padside */
#define 	BOTTOM 			1
#define 	RIGHT  			2
#define 	TOP    			3
#define 	LEFT   			4
#define 	NOSIDE 			5

/* modlue structure */
typedef struct moduletype {
    char name[MAXNAMELEN];
    int type;
    int left;
    int right;
    int top;
    int bottom;
    short pincount;
    struct pintype *pin;
} MODULE;

typedef struct celltype {
    int id;                 //unique ID of each CELL
    char name[MAXNAMELEN];  //name of CELL
    int width, height;      //width and height of CELL, used to obtain size of CELL
    short netcount2;        //size of net2 array
    short *net2; 			//stores the indices of NETs connected with CELL
    						//e.g.: cellArray[10].net2[0] is 1 then the first NET connected with CELL is netArray[1]

    short module; 			//unused in partitioning
    short netcount; 		//unused in partitioning
    short *net; 			//unused in partitioning
    /* pad information */
    short padside;          //unused in partitioning
    double padpos;          //unused in partitioning
} CELL;

typedef struct pintype {
    char name[MAXNAMELEN];
    int x, y;
    int layer;
    struct pintype *equiv;
} PIN;

typedef struct nettype {
    char name[MAXNAMELEN];  //name of the net and name of the circuit (?)
    short cellcount2;       //size of cell2 array
    short *cell2; 			//stores the indices of CELLs connected with NET
    						//e.g.: netArray[10].cell2[0] is 1 then the first CELL connected with NET is cellArray[1]

    short cellcount;        //unused in partitioning
    short *cell;            //unused in partitioning

    short *pin;             //unused in partitioning
} NET;

/* statistical informations */
typedef struct stattype {
    int max_net_per_cell;
    int max_cell_per_net;
    double avg_net_per_cell;
    double avg_cell_per_net;
    int max_cell_width;
    int min_cell_width;
    int max_cell_height;
    int min_cell_height;
    double avg_cell_width;
    double avg_cell_height;
    int feed_width;
    int max_cell_size;
    int min_cell_size;
    double avg_cell_size;

    /* standard deviation */
    double sig_net_per_cell;
    double sig_cell_per_net;
    double sig_cell_width;
    double sig_cell_height;
    double sig_cell_size;
} STAT;

extern CELL *cellArray; 		//nghiant: this one is important
extern NET *netArray; 			//nghiant: this one is important
extern MODULE *moduleArray;
extern STAT stat;

extern int cellNum;             //nghiant: this one is important
extern int netNum;              //nghiant: this one is important
extern int padNum;
extern int moduleNum;

extern short *tempnets;

/* utility functions */

/*
void  errorMsg();
void* mem_alloc();
FILE* file_open();
int   randomize();
double get_utime();
double get_stime();
*/
#define		RAND()			rand()
#define		SRAND(seed)		srand(seed)

#endif
