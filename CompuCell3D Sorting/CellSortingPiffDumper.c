#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int printCompleteSortingPIF();
int printEngulfmentPIF();
int printCheckerboardPIF();

/* Lattice Dimensions */
int latticeSize = 800;

/* Gap between cells - Default: 0 */
int gap = 0;

/* Margin between lattice and tissue */
int margin = 200;

/* Cell characteristics */
int cellSize = 20;

/* Store cell type for each cell in lattice */
int lattice[20][20] = {{-1}};


int intcmp(const void *aa, const void *bb){
	const int *a = aa;
	const int *b = bb;
	return (*a < *b) ? -1 : (*a > *b);
}


int main(int argc, char* argv[]){

	if (printCompleteSortingPIF() == 0){
		// print lattice
		int i = 0, j = 0;
		int index = 0;
		int x = margin;
		int y = margin;
		
		while (i < 20){
		
			if (y + margin > latticeSize){
				printf("Error: Vertical lattice size is too small\n");
				return 1;
			}
			
			j = 0;
			x = margin;
			
			while (j < 20){
			
				if (x + margin > latticeSize){
					printf("Error: Horizontal lattice size is too small\n");
					return 1;
				}
				
				if (lattice[i][j] == 1){
					printf("%d CellU %d %d %d %d 0 0\n", index, x, x+cellSize, y, y+cellSize);
				} else if (lattice[i][j] == 2){
					printf("%d CellV %d %d %d %d 0 0\n", index, x, x+cellSize, y, y+cellSize);	
				}
				
				index++;
				x += gap + cellSize;
				j++;
			}
			
			i++;
			y += cellSize + gap;
		}
		return 0;
		
	} else {
	
		printf("Error: Mismatch between generated lattice and desired tissue composition\n");
		return 1;
	}
}


/*
	Print complete sorting PIF file 
	Random cell placement
*/
int printCompleteSortingPIF(){
	/* Tissue composition */
	int numCells = 400;
	int numCellV = 0.4*400;					// 160 red cells				
	int numCellU = numCells - numCellV;		// 240 green cells
	
	int i = 0, j = 0;
	
	// select random indexes
	srand(time(NULL));
	int flag;
	int ind = 0, k = 0;
	int randomindex = -1;
	int indexes[numCellV];
	while(ind < numCellV){
		randomindex = rand() % numCells;
		flag = 0;
		for(k=0; k<ind; k++){
			if (indexes[k] == randomindex){
				flag = 1;
			}
		}
		if (flag == 0){
			indexes[ind] = randomindex;
			ind++;
		}
	}
	qsort(indexes,numCellV,sizeof(int),intcmp);
	
	while (i < 20){
		j = 0;
		while (j < 20){
		
			lattice[i][j] = 1;
		
			for (k=0; k<numCellV; k++){
				if(indexes[k] == i+j*20){
					lattice[i][j] = 2;
				}
			}
			j++;
		}
		i++;
	}

	// check
	i = 0, j = 0;
	int UCells = 0, VCells = 0;
	while (i < 20){
		j = 0;
		while (j < 20){
			if (lattice[i][j] == 1)
				UCells++;
			if (lattice[i][j] == 2)
				VCells++;
			j++;
		}
		i++;
	}
	
	if (UCells == numCellU && VCells == numCellV)
		return 0;
	return 1;
}


/*
	Print engulfment PIF file
	CellV (Red) cells surrounding CellU (Green) cells
	CellV (Red) engulfed by CellU (Green) cells
*/
int printEngulfmentPIF(){
	/* Tissue composition */
	int numCells = 400;
	int numCellV = (20*4 - 4) + (18*4 - 4);		// 144 red cells				
	int numCellU = numCells - numCellV;			// 384 green cells
	
	int i = 0, j = 0;
	
	while (i < 20){
		j = 0;
		while (j < 20){
		
			lattice[i][j] = 1; 	// type 1: green cell
		
			if (i == 0 || i == 1 || i == 20-1 || i == 20-2)
				lattice[i][j] = 2;	// type 2: red cell
			
			if (j == 0 || j == 1 || j == 20-1 || j == 20-2)
				lattice[i][j] = 2;	// type 2: red cell	
				
			j++;
		}
		i++;
	}
	
	// check
	i = 0, j = 0;
	int UCells = 0, VCells = 0;
	while (i < 20){
		j = 0;
		while (j < 20){
			if (lattice[i][j] == 1)
				UCells++;
			if (lattice[i][j] == 2)
				VCells++;
			j++;
		}
		i++;
	}
	
	if (UCells == numCellU && VCells == numCellV)
		return 0;
	return 1;
}


/*
	Print checkerboard sorting PIF file 
	CellV (Red) cells surrounded by CellU (Green) cells 
*/
int printCheckerboardPIF(){
	/* Tissue composition */
	int numCells = 400;
	int numCellV = 10*10;					// 100 red cells				
	int numCellU = numCells - numCellV;		// 300 green cells

	int i = 0, j = 0;
	
	while (i < 20){
		j = 0;
		while (j < 20){
		
			lattice[i][j] = 1; 	// type 1: green cell
		
			if (i >= 5 && i <= 15-1 && j >= 5 && j <= 15-1)
				lattice[i][j] = 2;	// type 2: red cell
				
			j++;
		}
		i++;
	}
	
	// check
	i = 0, j = 0;
	int UCells = 0, VCells = 0;
	while (i < 20){
		j = 0;
		while (j < 20){
			if (lattice[i][j] == 1)
				UCells++;
			if (lattice[i][j] == 2)
				VCells++;
			j++;
		}
		i++;
	}
	
	if (UCells == numCellU && VCells == numCellV)
		return 0;
	return 1;
}
