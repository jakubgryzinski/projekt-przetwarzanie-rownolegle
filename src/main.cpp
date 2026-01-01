#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

typedef struct {
    float time_seconds;
    long long operations;
} ComputeResult;

ComputeResult cpu_sequential(float* input, float* output, int N, int R);
ComputeResult gpu_efficient_global(float* input, float* output, int N, int R, int k);
ComputeResult gpu_inefficient_global(float* input, float* output, int N, int R, int k);
ComputeResult gpu_efficient_shared(float* input, float* output, int N, int R, int k);
ComputeResult gpu_inefficient_shared(float* input, float* output, int N, int R, int k);

void generate_random_data(float* data, int size);

void run_experiment_1();
void run_experiment_2(int N_large);
void run_experiment(int N, int R, int k, int BS, int experiment_num);

int verify_results(float* reference, float* test, int size, float tolerance);
void save_results_to_file(const char* filename, int experiment_num, int R, int BS, int k, int N, float* input, float* output, int file_index);
void save_performance_to_csv(int experiment_num, int N, int R, int BS, int k, ComputeResult results[5]);

const char* filenames[] = {
    "outputs/cpu_sequential_results.txt",
    "outputs/gpu_efficient_global_results.txt",
    "outputs/gpu_inefficient_global_results.txt",
    "outputs/gpu_efficient_shared_results.txt",
    "outputs/gpu_inefficient_shared_results.txt"
};
int num_filenames = 5;

const char* performance_csv = "outputs/performance.csv";

int main(int argc, char* argv[]) {
    srand(42);
    
    struct stat st = {0};
    if (stat("outputs", &st) == -1) {
        mkdir("outputs", 0755);
    }
    
    for (int i = 0; i < num_filenames; i++) {
        remove(filenames[i]);
    }
    remove(performance_csv);
    
    if (argc > 1) { // N_large odczytujemy z wyników pierwszego eksperymentu
        int N_large = atoi(argv[1]);
        run_experiment_2(N_large);
    } else {
        run_experiment_1();
    }
    
    // run_experiment_1();
    // run_experiment_2(2048);
    
    return 0;
}

void generate_random_data(float* data, int size) {
    for (int i = 0; i < size; i++) {
        data[i] = ((float)rand() / RAND_MAX) * 100.0f;
    }
}

void run_experiment_1() {
    printf("Eksperyment 1\n");
    int BS_values[] = {8, 16, 32};
    int num_BS = 3;
    
    int R_values[] = {6, 40};
    int num_R = 2;
    
    int N_values[] = {256, 512, 1024, 2048, 4096, 8192, 16384, 32768};
    int num_N = 8;
    
    // int N_values[] = {256, 512, 1024, 2048, 4096, 8192};
    // int num_N = 6;
    
    for (int r_idx = 0; r_idx < num_R; r_idx++) {
        int R = R_values[r_idx];
        
        for (int bs_idx = 0; bs_idx < num_BS; bs_idx++) {
            int BS = BS_values[bs_idx];
            
            for (int n_idx = 0; n_idx < num_N; n_idx++) {
                int N = N_values[n_idx];
                
                if (N <= 2 * R) {
                    continue;
                }
                
                printf("R=%d (%d/%d), BS=%d (%d/%d), N=%d (%d/%d)\n", R, r_idx + 1, num_R, BS, bs_idx + 1, num_BS, N, n_idx + 1, num_N);
                run_experiment(N, R, 1, BS, 1);
            }
        }
    }
}

void run_experiment_2(int N_large) {
    printf("Eksperyment 2\n");
    int N = 2 * N_large; // N_wys
    
    int BS_values[] = {8, 16, 32};
    int num_BS = 3;
    
    int R_values[] = {6, 40};
    int num_R = 2;
    
    int k_values[] = {1, 2, 8};
    int num_k = 3;
    
    for (int r_idx = 0; r_idx < num_R; r_idx++) {
        int R = R_values[r_idx];
        
        for (int bs_idx = 0; bs_idx < num_BS; bs_idx++) {
            int BS = BS_values[bs_idx];
            
            for (int k_idx = 0; k_idx < num_k; k_idx++) {
                int k = k_values[k_idx];
                printf("R=%d (%d/%d), BS=%d (%d/%d), k=%d (%d/%d), N=%d\n", R, r_idx + 1, num_R, BS, bs_idx + 1, num_BS, k, k_idx + 1, num_k, N);
                run_experiment(N, R, k, BS, 2);
            }
        }
    }
}

void run_experiment(int N, int R, int k, int BS, int experiment_num) {
    int input_size = N * N;
    int output_size = (N - 2*R) * (N - 2*R);
    float* input = (float*)malloc(input_size * sizeof(float));
    float* outputs[5];
    
    for (int i = 0; i < 5; i++) {
        outputs[i] = (float*)malloc(output_size * sizeof(float));
        if (!outputs[i]) {
            fprintf(stderr, "Memory allocation failed\n");
            for (int j = 0; j < i; j++) free(outputs[j]);
            free(input);
            return;
        }
    }
    
    generate_random_data(input, input_size);
    
    ComputeResult results[5];
    results[0] = cpu_sequential(input, outputs[0], N, R);
    results[1] = gpu_efficient_global(input, outputs[1], N, R, k);
    results[2] = gpu_inefficient_global(input, outputs[2], N, R, k);
    results[3] = gpu_efficient_shared(input, outputs[3], N, R, k);
    results[4] = gpu_inefficient_shared(input, outputs[4], N, R, k);
    
    save_performance_to_csv(experiment_num, N, R, BS, k, results);
    
    float tolerance = 1e-5;
    for (int i = 1; i < 5; i++) {
        if (!verify_results(outputs[0], outputs[i], output_size, tolerance)) {
            fprintf(stderr, "ERROR: Verification failed! N=%d, R=%d, k=%d, BS=%d\n", N, R, k, BS);
            exit(1);
        }
    }
    
    for (int i = 0; i < num_filenames; i++) {
        save_results_to_file(filenames[i], experiment_num, R, BS, k, N, input, outputs[i], i);
    }
    
    free(input);
    for (int i = 0; i < 5; i++) {
        free(outputs[i]);
    }
}

int verify_results(float* reference, float* test, int size, float tolerance) {
    for (int i = 0; i < size; i++) {
        float diff = fabsf(reference[i] - test[i]);
        float max_val = fmaxf(fabsf(reference[i]), fabsf(test[i]));
        float relative_error = diff / (max_val + 1e-10);
        
        if (relative_error > tolerance) {
            return 0;
        }
    }
    return 1;
}

void compute_stats(float* data, int size, float* min, float* max, float* mean, float* sum) {
    *min = data[0];
    *max = data[0];
    *sum = 0.0f;
    
    for (int i = 0; i < size; i++) {
        if (data[i] < *min) *min = data[i];
        if (data[i] > *max) *max = data[i];
        *sum += data[i];
    }
    
    *mean = *sum / size;
}

void print_separator(FILE* fp, const char* label, int value) {
    fprintf(fp, "------------------------------\n");
    fprintf(fp, "%s: %d\n", label, value);
    fprintf(fp, "------------------------------\n");
}

void save_array_info(FILE* fp, const char* label, float* data, int size) {
    float min, max, mean, sum;
    compute_stats(data, size, &min, &max, &mean, &sum);
    
    fprintf(fp, "%s size: %d elements\n", label, size);
    fprintf(fp, "%s stats: min=%.6f, max=%.6f, mean=%.6f, sum=%.6f\n", label, min, max, mean, sum);
    fprintf(fp, "%s sample (first 20):", label);
    
    int count = (size < 20) ? size : 20;
    for (int i = 0; i < count; i++) {
        fprintf(fp, " %.6f", data[i]);
    }
    fprintf(fp, "\n\n");
}

void save_results_header(FILE* fp, int experiment_num, int R, int BS, int k, int N, int file_index) {
    static int last_exp[5] = {-1, -1, -1, -1, -1};
    static int last_R[5] = {-1, -1, -1, -1, -1};
    static int last_BS[5] = {-1, -1, -1, -1, -1};
    static int last_k[5] = {-1, -1, -1, -1, -1};
    
    if (last_exp[file_index] != experiment_num) {
        print_separator(fp, "Eksperyment", experiment_num);
        last_exp[file_index] = experiment_num;
        last_R[file_index] = last_BS[file_index] = last_k[file_index] = -1;
    }
    if (last_R[file_index] != R) {
        print_separator(fp, "R", R);
        last_R[file_index] = R;
        last_BS[file_index] = last_k[file_index] = -1;
    }
    if (last_BS[file_index] != BS) {
        print_separator(fp, "BS", BS);
        last_BS[file_index] = BS;
        last_k[file_index] = -1;
    }
    if (last_k[file_index] != k) {
        print_separator(fp, "k", k);
        last_k[file_index] = k;
    }
    print_separator(fp, "N", N);
}

void save_results_to_file(const char* filename, int experiment_num, int R, int BS, int k, int N, float* input, float* output, int file_index) {
    FILE* fp = fopen(filename, "a");
    save_results_header(fp, experiment_num, R, BS, k, N, file_index);
    
    int input_size = N * N;
    int output_size = (N - 2*R) * (N - 2*R);
    
    save_array_info(fp, "Input", input, input_size);
    save_array_info(fp, "Output", output, output_size);
    
    fclose(fp);
}

void save_performance_to_csv(int experiment_num, int N, int R, int BS, int k, ComputeResult results[5]) {
    FILE* fp = fopen(performance_csv, "a");
    
    fseek(fp, 0, SEEK_END);
    long size = ftell(fp);
    
    if (size == 0) {
        fprintf(fp, "Experiment,N,R,BS,k,CPU_time,CPU_ops,GPU_Eff_Global_time,GPU_Eff_Global_ops,GPU_Ineff_Global_time,GPU_Ineff_Global_ops,GPU_Eff_Shared_time,GPU_Eff_Shared_ops,GPU_Ineff_Shared_time,GPU_Ineff_Shared_ops\n");
    }
    
    fprintf(fp, "%d,%d,%d,%d,%d,%.9f,%lld,%.9f,%lld,%.9f,%lld,%.9f,%lld,%.9f,%lld\n",
            experiment_num, N, R, BS, k,
            results[0].time_seconds, results[0].operations,
            results[1].time_seconds, results[1].operations,
            results[2].time_seconds, results[2].operations,
            results[3].time_seconds, results[3].operations,
            results[4].time_seconds, results[4].operations);
    
    fclose(fp);
}

ComputeResult cpu_sequential(float* input, float* output, int N, int R) {
    ComputeResult result;
    clock_t start = clock();
    
     //TODO
    int output_size = (N - 2*R) * (N - 2*R);
    memset(output, 0, output_size * sizeof(float));
     //TODO
    
    clock_t end = clock();
    result.time_seconds = (float)(end - start) / CLOCKS_PER_SEC;
    result.operations = 0; //TODO
    
    return result;
}

ComputeResult gpu_efficient_global(float* input, float* output, int N, int R, int k) {
    ComputeResult result;
    clock_t start = clock();
    
    //TODO
    int output_size = (N - 2*R) * (N - 2*R);
    memset(output, 0, output_size * sizeof(float));
     //TODO
    
    clock_t end = clock();
    result.time_seconds = (float)(end - start) / CLOCKS_PER_SEC;
    result.operations = 0;  //TODO
    
    return result;
}

ComputeResult gpu_inefficient_global(float* input, float* output, int N, int R, int k) {
    ComputeResult result;
    clock_t start = clock();
    
     //TODO
    int output_size = (N - 2*R) * (N - 2*R);
    memset(output, 0, output_size * sizeof(float));
     //TODO
    
    clock_t end = clock();
    result.time_seconds = (float)(end - start) / CLOCKS_PER_SEC;
    result.operations = 0; //TODO
    
    return result;
}

ComputeResult gpu_efficient_shared(float* input, float* output, int N, int R, int k) {
    ComputeResult result;
    clock_t start = clock();
    
     //TODO
    int output_size = (N - 2*R) * (N - 2*R);
    memset(output, 0, output_size * sizeof(float));
     //TODO
    
    clock_t end = clock();
    result.time_seconds = (float)(end - start) / CLOCKS_PER_SEC;
    result.operations = 0; //TODO
    
    return result;
}

ComputeResult gpu_inefficient_shared(float* input, float* output, int N, int R, int k) {
    ComputeResult result;
    clock_t start = clock();
    
     //TODO
    int output_size = (N - 2*R) * (N - 2*R);
    memset(output, 0, output_size * sizeof(float));
     //TODO
    
    clock_t end = clock();
    result.time_seconds = (float)(end - start) / CLOCKS_PER_SEC;
    result.operations = 0;  //TODO
    
    return result;
}