#include <iostream>
using namespace std;

void printArray(float arr[5][5], int N)
{
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			cout << arr[i][j]<<" ";
		}
		cout << endl;
	}
}
void printVector(float* arr, int N)
{
	for (int j = 0; j < N; ++j) {
		cout << arr[j]<<" ";
	}
	cout << endl;
}

float* backward_substitution(float input[5][5], float vector[5], int N){
    static float out[5]; //geht das so?
    for(int x = N-1; x >= 0; x--){
        float subtraction = 0;
        for(int i = x+1; i < N; i++){
            subtraction += input[x][i]*out[i]; 
        }
        out[x] = (vector[x]-subtraction)/input[x][x];
    }
    return out;
}

void GaussJordan_wo_pivot(float input[5][5], float vector[5], int N){
    for(int a=0; a < N-1; a++){
        float factor = input[a][a];
        for(int i=0;i<=N-1; i++){
            input[a][i] = input[a][i] / factor;
        }
        vector[a] = vector[a] / factor;
        for(int c=a+1; c<N; c++){            
            float subtraction_factor = input[c][a];
            for(int i=0;i<N; i++){
                input[c][i] = input[c][i]-subtraction_factor*input[a][i];
            }
            vector[c] = vector[c]-subtraction_factor*vector[a];
        }
    }
}

int index_largest_abs(float vector[5][5], int a){
    int index = a;
    float max = 0;
    for(int i=a;i<5;i++){
        if(vector[i][a] < -max){
            index = i;
            max = -vector[i][a];
        }
        if(vector[i][a] > max){
            index = i;
            max = vector[i][a];
        }
    }
    return index;
}

void GaussJordan_w_pivot(float input[5][5], float vector[5], int N){
    for(int a=0; a < N-1; a++){
        // pivoting
        int maxindex = index_largest_abs(input, a);
        // cout<<"maxindex: "<<maxindex<<"\n";
        float tmprow;
        for(int l=0;l<N;l++){
            tmprow = input[a][l];
            input[a][l] = input[maxindex][l];
            input[maxindex][l] = tmprow;
        }
        float tmpvec = vector[a];
        vector[a] = vector[maxindex];
        vector[maxindex] = tmpvec;
        // GaussJordan
        float factor = input[a][a];
        for(int i=0;i<=N-1; i++){ 
            input[a][i] = input[a][i] / factor;
        }
        vector[a] = vector[a] / factor;
        for(int c=a+1; c<N; c++){            
            float subtraction_factor = input[c][a];
            for(int i=0;i<N; i++){
                input[c][i] = input[c][i]-subtraction_factor*input[a][i];
            }
            vector[c] = vector[c]-subtraction_factor*vector[a];
        }
    }
}

float* matrix_vector_mult(float input[5][5], float vector[5]){
    static float out[5] = {0,0,0,0,0};
    for(int i=0;i<5;i++){
        for(int j=0;j<5;j++){
            out[i] += input[i][j]*vector[j];
        }
    }
    return out;
}


int main(){
    //input Daten
    float aufgabenmatrix[5][5] = {3, -1, 1, 4, 1,
                                1, 0, 1, 1, 1,
                                -1, -1, -1, -1,-1,
                                0, 0, 0, 7, 2,
                                3, 4, -3, 4, 1};
    float b0[] = {1,-1,1,-1,1};
    float b1[] = {0,1,2,3,4};
    const int N = 5;
    GaussJordan_w_pivot(aufgabenmatrix, b0, N);
    printArray(aufgabenmatrix, N);
    cout<<"\n";
    float* output = backward_substitution(aufgabenmatrix, b0, N);
    printVector(output, N);
    cout<<"\n";
    float aufgabenmatrix2[5][5] = {3, -1, 1, 4, 1, 1, 0, 1, 1, 1,-1, -1, -1, -1,-1,0, 0, 0, 7, 2,3, 4, -3, 4, 1};
    float* testergebnis = matrix_vector_mult(aufgabenmatrix2, output);
    printVector(testergebnis, N);
    while(1){}
    return 0;
}