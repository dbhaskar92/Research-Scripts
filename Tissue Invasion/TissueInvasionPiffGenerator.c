#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int intcmp(const void *aa, const void *bb){
	const int *a = aa;
	const int *b = bb;
	return (*a < *b) ? -1 : (*a > *b);
}

int main(){

	int i,j,k;
	int randomnum;
	int flag;
	int counter = 1;

	int numberofcells = 400;	// Default: 400
	int numleaders = 160;		// Default: 80
	int leaders[numleaders];
	int width = 40;			// Default: 40
	int height = 10;		// Default: 10
	int cellsize = 20;		// Default: 5

	int y_max = cellsize*height - cellsize;
	int x_max = cellsize*width - cellsize;

	srand(time(NULL));

	// randomly select leader cells

	i = 0;
	while(i<numleaders){
		randomnum = (rand() % numberofcells) + 1;
		flag = 0;
		for(k=0; k<i; k++){
			if (leaders[k] == randomnum){
				flag = 1;
			}
		}
		if (flag == 0){
			leaders[i] = randomnum;
			i++;
		}
	}

	qsort(leaders,numleaders,sizeof(int),intcmp);

	printf("Selected leader cells: \n");
	for(i=0; i<numleaders; i++){
		printf("%d ",leaders[i]);
	}
	printf("\n\n");

	// generate piff output

	for(j=0; j<=y_max; j=j+cellsize){
		for(i=0; i<=x_max; i=i+cellsize){	
			flag = 0;
		
			for (k=0; k<numleaders; k++){
				if(leaders[k] == counter){
					flag = 1;
				}
			}
			if (flag == 1){
				printf("%d Red %d %d %d %d 0 0\n",counter,i, i+cellsize, j, j+cellsize);
			} else {
				printf("%d Green %d %d %d %d 0 0\n",counter,i, i+cellsize, j, j+cellsize);
			}
			counter++;
	
		}
		printf("\n");
	}

	return 0;

}
