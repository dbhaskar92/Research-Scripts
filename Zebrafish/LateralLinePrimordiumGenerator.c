#include<stdio.h>
#include<stdlib.h>

int printMedPrimordiumV1(int,int,int,int,int);
int printMedPrimordiumV2(int,int,int,int,int);

int main(int argc, char* argv[]){

	/* Lattice Dimensions */
	int latticeWidth = 750;
	int latticeHeight = 200; 
	
	/* Separation between "front" and "back" cells - Default: 0 */
	int split = 0;

	/* Primordium placement */
	int xstart = 20;

	/* Cell characteristics */
	int cellSize = 10;
	
	return printMedPrimordiumV2(latticeWidth, latticeHeight, xstart, split, cellSize);
}

/*
	Print small sized PLLP
	Back Configuration: 2 X 1 | 4 X 3 | 6 X 2 = 26
	Front Configuration: 6 X 3 | 4 X 4 | 2 X 1 = 36
*/
int printSmlPrimordium(int latticeWidth, int latticeHeight, int xstart, int split, int cellsize){

	int numCells = 62;
	int fgfCells = 26;
	int wntCells = numCells - fgfCells;
	int cellindex = 0;

	int y_middle = latticeHeight/2;
	int y_index = 0;
	int y_low = 0;
	int y_high = 0;
	
	// Dimension checks
	int h_lattice_sites = (1 + 3 + 2 + 3 + 4 + 1) * cellsize;
	int v_lattice_sites = 6 * cellsize;
	if (h_lattice_sites >= latticeWidth || v_lattice_sites >= latticeHeight){
		printf("Error: Lattice dimensions too small\n");
		return 1;
	}
	
	// Print FGF cells
	int x_index = xstart;
	while ( (x_index - xstart) < (1 + 3 + 2) * cellsize ){
		
		if (cellindex < 2*1){
			y_index = -1;
			while (y_index < 1){
				y_low = y_middle + y_index*cellsize;
				y_high = y_middle + y_index*cellsize + cellsize;
				printf("%d Back %d %d %d %d 0 0\n",cellindex++,x_index,x_index+cellsize,y_low,y_high);
				y_index++;
			}
		} else if (cellindex < 2*1 + 4*3){
			y_index = -2;
			while (y_index < 2){
				y_low = y_middle + y_index*cellsize;
				y_high = y_middle + y_index*cellsize + cellsize;
				printf("%d Back %d %d %d %d 0 0\n",cellindex++,x_index,x_index+cellsize,y_low,y_high);
				y_index++;
			}
		} else if (cellindex < fgfCells){
			y_index = -3;
			while (y_index < 3){
				y_low = y_middle + y_index*cellsize;
				y_high = y_middle + y_index*cellsize + cellsize;
				printf("%d Back %d %d %d %d 0 0\n",cellindex++,x_index,x_index+cellsize,y_low,y_high);
				y_index++;
			}
		}
		printf("\n");
		x_index = x_index + cellsize;
	}
	
	// Print Wnt cells
	x_index += split;
	while ( (x_index - xstart - split) < (1 + 3 + 2 + 3 + 4 + 1) * cellsize ){
		
		if (cellindex < fgfCells + 6*3){
			y_index = -3;
			while (y_index < 3){
				y_low = y_middle + y_index*cellsize;
				y_high = y_middle + y_index*cellsize + cellsize;
				printf("%d Front %d %d %d %d 0 0\n",cellindex++,x_index,x_index+cellsize,y_low,y_high);
				y_index++;
			}
		} else if (cellindex < fgfCells + 6*3 + 4*4){
			y_index = -2;
			while (y_index < 2){
				y_low = y_middle + y_index*cellsize;
				y_high = y_middle + y_index*cellsize + cellsize;
				printf("%d Front %d %d %d %d 0 0\n",cellindex++,x_index,x_index+cellsize,y_low,y_high);
				y_index++;
			}
		} else if (cellindex < numCells){
			y_index = -1;
			while (y_index < 1){
				y_low = y_middle + y_index*cellsize;
				y_high = y_middle + y_index*cellsize + cellsize;
				printf("%d Front %d %d %d %d 0 0\n",cellindex++,x_index,x_index+cellsize,y_low,y_high);
				y_index++;
			}
		}
		printf("\n");
		x_index = x_index + cellsize;
	}
	if (cellindex == fgfCells + wntCells - 1){
		return 0;
	}
	return 1;
}

/*
	Print medium sized PLLP (version 1 - thick)
	Back configuration: 4 X 2 | 6 X 3 | 8 X 3 = 50
	Front configuration: 8 X 6 | 6 X 3 | 4 X 2 = 74
*/
int printMedPrimordiumV1(int latticeWidth, int latticeHeight, int xstart, int split, int cellsize){

	int numCells = 124;
	int fgfCells = 50;
	int wntCells = numCells - fgfCells;
	int cellindex = 0;

	int y_middle = latticeHeight/2;
	int y_index = 0;
	int y_low = 0;
	int y_high = 0;
	
	// Dimension checks
	int h_lattice_sites = (2 + 3 + 3 + 6 + 3 + 2) * cellsize;
	int v_lattice_sites = 8 * cellsize;
	if (h_lattice_sites >= latticeWidth || v_lattice_sites >= latticeHeight){
		printf("Error: Lattice dimensions too small\n");
		return 1;
	}
	
	// Print FGF cells
	int x_index = xstart;
	while ( (x_index - xstart) < (2 + 3 + 3) * cellsize ){
		
		if (cellindex < 4*2){
			y_index = -2;
			while (y_index < 2){
				y_low = y_middle + y_index*cellsize;
				y_high = y_middle + y_index*cellsize + cellsize;
				printf("%d Back %d %d %d %d 0 0\n",cellindex++,x_index,x_index+cellsize,y_low,y_high);
				y_index++;
			}
		} else if (cellindex < 4*2 + 6*3){
			y_index = -3;
			while (y_index < 3){
				y_low = y_middle + y_index*cellsize;
				y_high = y_middle + y_index*cellsize + cellsize;
				printf("%d Back %d %d %d %d 0 0\n",cellindex++,x_index,x_index+cellsize,y_low,y_high);
				y_index++;
			}
		} else if (cellindex < fgfCells){
			y_index = -4;
			while (y_index < 4){
				y_low = y_middle + y_index*cellsize;
				y_high = y_middle + y_index*cellsize + cellsize;
				printf("%d Back %d %d %d %d 0 0\n",cellindex++,x_index,x_index+cellsize,y_low,y_high);
				y_index++;
			}
		}
		printf("\n");
		x_index = x_index + cellsize;
	}
	
	// Print Wnt cells
	x_index += split;
	while ( (x_index - xstart - split) < (2 + 3 + 3 + 6 + 3 + 2) * cellsize ){
		
		if (cellindex < fgfCells + 8*6){
			y_index = -4;
			while (y_index < 4){
				y_low = y_middle + y_index*cellsize;
				y_high = y_middle + y_index*cellsize + cellsize;
				printf("%d Front %d %d %d %d 0 0\n",cellindex++,x_index,x_index+cellsize,y_low,y_high);
				y_index++;
			}
		} else if (cellindex < fgfCells + 8*6 + 6*3){
			y_index = -3;
			while (y_index < 3){
				y_low = y_middle + y_index*cellsize;
				y_high = y_middle + y_index*cellsize + cellsize;
				printf("%d Front %d %d %d %d 0 0\n",cellindex++,x_index,x_index+cellsize,y_low,y_high);
				y_index++;
			}
		} else if (cellindex < numCells){
			y_index = -2;
			while (y_index < 2){
				y_low = y_middle + y_index*cellsize;
				y_high = y_middle + y_index*cellsize + cellsize;
				printf("%d Front %d %d %d %d 0 0\n",cellindex++,x_index,x_index+cellsize,y_low,y_high);
				y_index++;
			}
		}
		printf("\n");
		x_index = x_index + cellsize;
	}
	if (cellindex == fgfCells + wntCells - 1){
		return 0;
	}
	return 1;
}


/*
	Print medium sized PLLP (version 2 - thin)
	Back configuration: 4 X 5 | 6 X 5 = 50
	Front configuration: 6 X 7 | 4 X 8 = 74
*/
int printMedPrimordiumV2(int latticeWidth, int latticeHeight, int xstart, int split, int cellsize){

	int numCells = 124;
	int fgfCells = 50;
	int wntCells = numCells - fgfCells;
	int cellindex = 0;

	int y_middle = latticeHeight/2;
	int y_index = 0;
	int y_low = 0;
	int y_high = 0;
	
	// Dimension checks
	int h_lattice_sites = (5 + 5 + 7 + 8) * cellsize;
	int v_lattice_sites = 6 * cellsize;
	if (h_lattice_sites >= latticeWidth || v_lattice_sites >= latticeHeight){
		printf("Error: Lattice dimensions too small\n");
		return 1;
	}
	
	// Print FGF cells
	int x_index = xstart;
	while ( (x_index - xstart) < (5 + 5) * cellsize ){
		
		if (cellindex < 4*5){
			y_index = -2;
			while (y_index < 2){
				y_low = y_middle + y_index*cellsize;
				y_high = y_middle + y_index*cellsize + cellsize;
				printf("%d Back %d %d %d %d 0 0\n",cellindex++,x_index,x_index+cellsize,y_low,y_high);
				y_index++;
			}
		} else if (cellindex < fgfCells){
			y_index = -3;
			while (y_index < 3){
				y_low = y_middle + y_index*cellsize;
				y_high = y_middle + y_index*cellsize + cellsize;
				printf("%d Back %d %d %d %d 0 0\n",cellindex++,x_index,x_index+cellsize,y_low,y_high);
				y_index++;
			}
		}
		printf("\n");
		x_index = x_index + cellsize;
	}
	
	// Print Wnt cells
	x_index += split;
	while ( (x_index - xstart - split) < (5 + 5 + 7 + 8) * cellsize ){
		
		if (cellindex < fgfCells + 6*7){
			y_index = -3;
			while (y_index < 3){
				y_low = y_middle + y_index*cellsize;
				y_high = y_middle + y_index*cellsize + cellsize;
				printf("%d Front %d %d %d %d 0 0\n",cellindex++,x_index,x_index+cellsize,y_low,y_high);
				y_index++;
			}
		} else if (cellindex < numCells){
			y_index = -2;
			while (y_index < 2){
				y_low = y_middle + y_index*cellsize;
				y_high = y_middle + y_index*cellsize + cellsize;
				printf("%d Front %d %d %d %d 0 0\n",cellindex++,x_index,x_index+cellsize,y_low,y_high);
				y_index++;
			}
		}
		printf("\n");
		x_index = x_index + cellsize;
	}
	if (cellindex == fgfCells + wntCells - 1){
		return 0;
	}
	return 1;
}
