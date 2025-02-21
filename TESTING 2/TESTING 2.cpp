#include <iostream>
#include <omp.h>
#include <vector>
#include <chrono>

int main() {
    // Set array size
    const int array_size = 1000000000;
    std::vector<int> arr(array_size); // Fill array with 1s

    omp_set_num_threads(16); 

    // Parallel sum using OpenMP
    auto start = std::chrono::high_resolution_clock::now();
    long long parallel_sum = 0;

#pragma omp parallel for reduction(+:parallel_sum)
    for (int i = 0; i < array_size; ++i) {
        parallel_sum += arr[i];
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto parallel_duration = std::chrono::duration<double>(end - start).count();

    // Output result
    std::cout << "Time taken (parallel): " << parallel_duration << " seconds." << std::endl;

    return 0;
}
