#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sim_bp.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdbool.h>
using namespace std;

/*  argc holds the number of command line arguments
    argv[] holds the commands themselves

    Example:-
    sim bimodal 6 gcc_trace.txt
    argc = 4
    argv[0] = sim
    argv[1] = bimodal
    argv[2] = 6
    ... and so on
*/
// ****************************************************************************************************************************//
unsigned char *prediction_table_b = NULL;
unsigned int prediction, misprediction; 
void initialize_prediction_table_b(int m) {
	int size = 1 << m; 
    prediction_table_b = (unsigned char *)malloc(size * sizeof(unsigned char));
    
    for (int i = 0; i < size; i++) {
        prediction_table_b[i] = 2; // 2 in binary is 10
    }
}

//***************************** BIMODAL ************************************************//
int bimodal_func(int m, unsigned long int PC_addr, char actual_outcome) {
    
    //step 1
    unsigned long int index = (PC_addr >> 2) & ((1 << m) - 1);
    // step2
    unsigned char pred_val = prediction_table_b[index];

    char predicted_outcome;
    if (pred_val >= 2) {
        predicted_outcome = 't';  // Predicted taken
    } else 
        predicted_outcome = 'n';  // Predicted not taken

    //step 3
    if ((actual_outcome == 't') && (pred_val < 3)) prediction_table_b[index] += 1;
    else if((actual_outcome == 'n') && (pred_val > 0)) prediction_table_b[index] -= 1;

    return (predicted_outcome != actual_outcome) ? 1 : 0;

}

//*************************************** GSHARE ***************************************************//

unsigned char *prediction_table_g = NULL;
void initialize_prediction_table_g(int m) {
    int size = 1 << m;
    prediction_table_g = (unsigned char *)malloc(size * sizeof(unsigned char));
    if (prediction_table_g == NULL) {
        printf("Memory allocation failed for prediction_table_g\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < size; i++) {
        prediction_table_g[i] = 2;  // Initialize to "weakly taken"
    }
}


bool isNBitLong(unsigned int GHR, int n) {
    // Create a mask with n bits set to 1. For example, if n = 4, mask = 0b1111 = 15.
    unsigned int mask = (1 << n) - 1;

    // Check if GHR fits within n bits by using the mask.
    // If (GHR & ~mask) is 0, GHR is n-bit long.
    return (GHR & ~mask) == 0;
}

int gshare_func(int m, unsigned long int PC_addr, int n, unsigned int &globalBranchHistoryreg, char actual_outcome) {
        
        uint32_t pred_val;
        PC_addr = PC_addr >> 2;

        uint32_t full_m_mask = (1 << (m)) - 1; 
        uint32_t m_bits_PC = PC_addr & full_m_mask;

        uint32_t mask = (1 << (m-n)) - 1;
        uint32_t m_n_bits = m_bits_PC & mask;

        uint32_t gbhr_index = globalBranchHistoryreg & ((1 << n) - 1); 
        uint32_t n_outof_m = m_bits_PC >> (m-n);

        long int Xor = n_outof_m ^ gbhr_index;
    
        Xor = Xor << (m-n);
        uint32_t index = Xor | m_n_bits;
      //  cout<<"index = "<<index<<endl;

/*
        unsigned int upper_n_bits = (PC_addr >> (m-n)) & ((1 << n) - 1);
        unsigned int xor_bhr_upperN = globalBranchHistoryreg ^ upper_n_bits; 
        unsigned int lower_mn_bits = PC_addr & ((1 << (m-n)) - 1);
        unsigned int index =  (xor_bhr_upperN << (m-n)) | lower_mn_bits;
*/
    // Debug: Check the index value
    // printf("Calculated index: %u\n", index);

    // Access prediction table and determine predicted outcome
    pred_val = prediction_table_g[index];
    char predicted_outcome = (pred_val >= 2) ? 't' : 'n';

    // Check if the prediction was a misprediction
    int is_misprediction = (predicted_outcome != actual_outcome);

    // Update prediction table based on actual outcome
    if (actual_outcome == 't') {
        if (pred_val < 3) prediction_table_g[index]++;
    } else {
        if (pred_val > 0) prediction_table_g[index]--;
    }

    // Return whether it was a misprediction
    return is_misprediction;
}

void updateReg(uint32_t &globalBranchHistoryreg, char actual_outcome, int n){
    globalBranchHistoryreg = (globalBranchHistoryreg >> 1) | ((actual_outcome == 't' ? 1 : 0) << (n - 1));
}


//****************************************** HYBRID *****************************************************/
unsigned char *chooser_table;

void initialize_chooser_table(int k) {
	int size = 1 << k;
    chooser_table = (unsigned char *)malloc(size * sizeof(unsigned char));
    for (int i = 0; i < size; i++) {
        chooser_table[i] = 1; 
    }
}
int hybrid_func(int M2, int M1, int k, int n, unsigned long int PC_addr, unsigned int globalBranchHistoryreg, char actual_outcome) {
    uint32_t pred_val_g, pred_val_b;
    PC_addr = PC_addr >> 2;

    // Gshare index calculation
    uint32_t full_m_mask = (1 << M1) - 1;
    uint32_t m_bits_PC = PC_addr & full_m_mask;
    uint32_t mask = (1 << (M1 - n)) - 1;
    uint32_t m_n_bits = m_bits_PC & mask;
    uint32_t gbhr_index = globalBranchHistoryreg & ((1 << n) - 1);
    uint32_t n_outof_m = m_bits_PC >> (M1 - n);
    uint32_t Xor = (n_outof_m ^ gbhr_index) << (M1 - n);
    uint32_t index_g = Xor | m_n_bits;

    // Bimodal index calculation
    unsigned long int index_b = (PC_addr) & ((1 << M2) - 1);

    // Get chooser table index
    unsigned long int index_c = (PC_addr) & ((1 << k) - 1);

    // Predictions from both predictors
    pred_val_g = prediction_table_g[index_g];
    pred_val_b = prediction_table_b[index_b];
    char predicted_outcome_g = (pred_val_g >= 2) ? 't' : 'n';
    char predicted_outcome_b = (pred_val_b >= 2) ? 't' : 'n';

    // Chooser table outcome to select the predictor
    unsigned char chosen_val = chooser_table[index_c];
    char predicted_outcome = (chosen_val >= 2) ? predicted_outcome_g : predicted_outcome_b;

    // Update the predictor that was chosen
    if (chosen_val >= 2) {
        if (actual_outcome == 't' && pred_val_g < 3) {
            prediction_table_g[index_g]++;
        } else if (actual_outcome == 'n' && pred_val_g > 0) {
            prediction_table_g[index_g]--;
        }
    } else {
        if (actual_outcome == 't' && pred_val_b < 3) {
            prediction_table_b[index_b]++;
        } else if (actual_outcome == 'n' && pred_val_b > 0) {
            prediction_table_b[index_b]--;
        }
    }

    // Check for misprediction
    int misprediction = (predicted_outcome != actual_outcome);

    // Update chooser table based on the accuracy of each predictor
    if ((actual_outcome == predicted_outcome_g) && (actual_outcome != predicted_outcome_b)) {
        if (chosen_val < 3) chooser_table[index_c]++;
    } else if ((actual_outcome != predicted_outcome_g) && (actual_outcome == predicted_outcome_b)) {
        if (chosen_val > 0) chooser_table[index_c]--;
    }

    return misprediction;
}
int main (int argc, char* argv[])
{
    FILE *FP;               // File handler
    char *trace_file;       // Variable that holds trace file name
    bp_params params;       // Struct to hold parameters
    params.M1 = 0;
    params.M2 = 0;
    params.N = 0;
    params.K = 0;
    char outcome;           // Variable to hold branch outcome
    unsigned long int addr; // Variable to hold address read from input file
    unsigned int globalBranchHistoryreg = 0;
    int total_predictions = 0;
    int total_mispredictions = 0;

    if (!(argc == 4 || argc == 5 || argc == 7)) {
        printf("Error: Wrong number of inputs:%d\n", argc - 1);
        exit(EXIT_FAILURE);
    }

    params.bp_name = argv[1];

    if (strcmp(params.bp_name, "bimodal") == 0) {
        if (argc != 4) {
            printf("Error: %s wrong number of inputs:%d\n", params.bp_name, argc - 1);
            exit(EXIT_FAILURE);
        }
        params.M2 = strtoul(argv[2], NULL, 10);
        trace_file = argv[3];
        printf("COMMAND\n%s %s %lu %s\n", argv[0], params.bp_name, params.M2, trace_file);
    } else if (strcmp(params.bp_name, "gshare") == 0) {
        if (argc != 5) {
            printf("Error: %s wrong number of inputs:%d\n", params.bp_name, argc - 1);
            exit(EXIT_FAILURE);
        }
        params.M1 = strtoul(argv[2], NULL, 10);
        params.N = strtoul(argv[3], NULL, 10);
        trace_file = argv[4];
        printf("COMMAND\n %s %s %lu %lu %s\n", argv[0], params.bp_name, params.M1, params.N, trace_file);
    } else if (strcmp(params.bp_name, "hybrid") == 0) {
        if (argc != 7) {
            printf("Error: %s wrong number of inputs:%d\n", params.bp_name, argc - 1);
            exit(EXIT_FAILURE);
        }
        params.K = strtoul(argv[2], NULL, 10);
        params.M1 = strtoul(argv[3], NULL, 10);
        params.N = strtoul(argv[4], NULL, 10);
        params.M2 = strtoul(argv[5], NULL, 10);
        trace_file = argv[6];
        printf("COMMAND\n %s %s %lu %lu %lu %lu %s\n", argv[0], params.bp_name, params.K, params.M1, params.N, params.M2, trace_file);
    } else {
        printf("Error: Wrong branch predictor name:%s\n", params.bp_name);
        exit(EXIT_FAILURE);
    }

    FP = fopen(trace_file, "r");
    if (FP == NULL) {
        printf("Error: Unable to open file %s\n", trace_file);
        exit(EXIT_FAILURE);
    }

    initialize_prediction_table_b(params.M2);
    initialize_prediction_table_g(params.M1);
    initialize_chooser_table(params.K);

    char str[2];
    while (fscanf(FP, "%lx %s", &addr, str) != EOF) {
        outcome = str[0];
        

        //int misprediction = 0;
        if (strcmp(params.bp_name, "bimodal") == 0) {
		total_predictions++;
            total_mispredictions += bimodal_func(params.M2, addr, outcome);
        } else if (strcmp(params.bp_name, "gshare") == 0) {
		total_predictions++;
            total_mispredictions += gshare_func(params.M1, addr, params.N, globalBranchHistoryreg, outcome);
            updateReg(globalBranchHistoryreg, outcome, params.N);
           // cout<<"ghr = "<<globalBranchHistoryreg;
        } else if (strcmp(params.bp_name, "hybrid") == 0) {
		total_predictions++;
            total_mispredictions += hybrid_func(params.M2, params.M1, params.K, params.N, addr, globalBranchHistoryreg, outcome);
            updateReg(globalBranchHistoryreg, outcome, params.N);
        }

    }
    fclose(FP);

    float misprediction_rate = (float)total_mispredictions / total_predictions * 100;
    printf("OUTPUT\n");
    printf(" number of predictions:    %d\n", total_predictions);
    printf(" number of mispredictions: %d\n", total_mispredictions);
    printf(" misprediction rate:       %.2f%%\n", misprediction_rate);

    if (strcmp(params.bp_name, "bimodal") == 0) {
        printf("FINAL BIMODAL CONTENTS\n");
        for (int i = 0; i < (1 << params.M2); i++) {
            printf(" %d\t%d\n", i, prediction_table_b[i]);
        }
    } else if (strcmp(params.bp_name, "gshare") == 0) {
    	printf("FINAL GSHARE CONTENTS\n");
    	for (int i = 0; i < (1 << params.M1); i++) {
        	printf(" %d\t%d\n", i, prediction_table_g[i]);
      }
    }else if (strcmp(params.bp_name, "hybrid") == 0) {
    printf("FINAL CHOOSER CONTENTS\n");
    for (int i = 0; i < (1 << params.K); i++) {
        printf(" %d\t%d\n", i, chooser_table[i]);
    }
    printf("FINAL GSHARE CONTENTS\n");
    for (int i = 0; i < (1 << params.M1); i++) {
    printf(" %d\t%d\n", i, prediction_table_g[i]); 
    }
    printf("FINAL BIMODAL CONTENTS\n");
    for (int i = 0; i < (1 << params.M2); i++) {
    printf(" %d\t%d\n", i, prediction_table_b[i]);
    }
    
}
    free(prediction_table_b);
    free(prediction_table_g);
    free(chooser_table);

    return 0;
}
