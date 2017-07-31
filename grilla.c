#include <stdio.h>
#include <stdlib.h>
#define F 4
#define C 5
void llenar(int* v);
void imprimir(int* v);

int main(){
	int* v;
	v = malloc(F*C*sizeof(int)); 
	llenar(v);
	imprimir(v);
	return 0;
}

void llenar(int* v){
	int i;
	for(i=0;i<(F*C);i=i+1){
		v[i] = rand() % 2;
	}
}

void imprimir(int* v){
	int i, j;
	for (i=0; i<F; i=i+1){
		for(j=0; j<C; j=j+1){	

     		printf("%d   ",v[i*F+j]);
		if(j==C-1){ printf("\n"); }
		}	
	} 
}	
/*
void imprimir(int* v){
	int r, h;
		for(r=0;r<(F*C);r=r+1){
			for(h=1;h<F;h=h+1){
			printf("%d", v[r]);
			if(r==(h*C-1)){ printf("\n"); }
			}
		}
}
*/
