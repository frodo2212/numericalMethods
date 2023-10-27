#include <iostream>
using namespace std;


// printet mir halt nen array am pointer
void printArray(float** arr, int N)
{
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			cout << arr[i][j]<<" ";
		}
		cout << endl;
	}
}

void printArray(float arr[5][5], int N)
{
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			cout << arr[i][j]<<" ";
		}
		cout << endl;
	}
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

int index_largest_abs(float vector[5]){
    int index = 0;
    float max = 0;
    for(int i=0;i<5;i++){
        if(vector[i] < -max || vector[i]>max){
            index = i;
            max = vector[i];
        }
    }
    return index;
}

void GaussJordan_w_pivot(float input[5][5], float vector[5], int N){
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


int main(){
    //input Daten
    float aufgabenmatrix[5][5] = {3, -1, 1, 4, 1,
                                1, 0, 1, 1, 1,
                                -1, -1, -1, -1,-1,
                                0, 0, 0, 7, 2,
                                3, 4, -3, 4, 1};
    float testmatrix[5][5] = {3, -1, 1, 4, 1,1, 1, 1, 3, 1,0, -1, -1, -1, -1,0, 0, 0, 7, 2,1, 4, -3, 4, 1};

    float b0[] = {1,-1,1,-1,1};
    float b1[] = {0,1,2,3,4};
    const int N = 5;
    // rechnungen
    printArray(testmatrix, N);
    cout<<"\n\n";
    GaussJordan_wo_pivot(testmatrix, b0, N);
    printArray(testmatrix, N);
    // float* output = backward_substitution(aufgabenmatrix, b0, N);
    // das hier plottet den Output
    //for(int i=0 ; i<5; i++)
	//	cout<<output[i]<<"\t";

}