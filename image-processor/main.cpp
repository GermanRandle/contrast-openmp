#include <omp.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>

using namespace std;

int main(int argc, char** argv) {
    // Reading and checking input
    if (argc != 5) {
        cout << "Program must be run with 4 arguments only" << endl;
        return 0;
    }
    int threads;
    try {
        threads = stoi(argv[1]);
    } catch (invalid_argument &e) {
        cout << "Number of threads must be an integer" << endl;
        return 0;
    } catch (out_of_range &e) {
        cout << "The number of threads is out of range" << endl;
        return 0;
    }
    if (threads < 0) {
        cout << "Number of threads mustn't be negative" << endl;
        return 0;
    }
    #ifdef _OPENMP
        if (threads != 0) {
            omp_set_num_threads(threads);
        }
    #endif
    string input_file_name = argv[2];
    string output_file_name = argv[3];
    double coeff;
    try {
        coeff = stod(argv[4]);
    } catch (invalid_argument &e) {
        cout << "Coefficient must be a real number" << endl;
        return 0;
    } catch (out_of_range &e) {
        cout << "Coefficient is out of range" << endl;
        return 0;
    }
    if (coeff < 0 || coeff >= 0.5) {
        cout << "Coefficient must belong to semi-interval [0; 0.5)" << endl;
        return 0;
    }
    // Reading header of the input file and checking
    ifstream input;
    input.open(input_file_name, ios::binary);
    if (!input.is_open()) {
        cout << "File cannot be opened :(" << endl;
        return 0;
    }
    string file_type;
    int width, height, max_color;
    try {
        input >> file_type >> width >> height >> max_color;
    } catch (ios_base::failure &e) {
        cout << "File header cannot be read :(" << endl;
        input.close();
        return 0;
    }
    if (file_type != "P5" && file_type != "P6") {
        cout << "Unsupported file format" << endl;
        input.close();
        return 0;
    }
    if (width <= 0 || height <= 0) {
        cout << "Incorrect width/height of the image" << endl;
        input.close();
        return 0;
    }
    if (max_color != 255) {
        cout << "Incorrect max color value" << endl;
        input.close();
        return 0;
    }
    if (file_type == "P5") {
        vector<unsigned char> pic(height * width);
        try {
            input.get(); // The first symbol is always garbage
            for (int i = 0; i < height * width; i++) {
                pic[i] = input.get();
            }
        } catch (ios::failure &e) {
            cout << "Error while reading the image :(" << endl;
            input.close();
            return 0;
        }
        input.close();
        // Algorithm
        vector<int> freq(max_color + 1, 0);
        vector<unsigned char> stretched(max_color + 1);
        #ifdef _OPENMP
            double before_algo_time = omp_get_wtime();
            #pragma omp parallel default(none) shared(height, width, max_color, pic, freq, coeff, stretched)
            {
                vector<int> thread_freq(max_color + 1, 0);
                #pragma omp for schedule(static)
                for (int i = 0; i < height * width; i++) {
                    thread_freq[pic[i]]++;
                }
                #pragma omp critical
                {
                    for (int i = 0; i <= max_color; i++) {
                        freq[i] += thread_freq[i];
                    }
                }
                #pragma omp barrier
                #pragma omp single
                {
                    int psum = freq[0];
                    unsigned char left_border = max_color;
                    unsigned char right_border = 0;
                    for (int i = 1; i <= max_color; i++) {
                        psum += freq[i];
                        if (float(psum) / float(height * width) >= coeff) {
                            left_border = i - 1;
                            break;
                        }
                    }
                    psum = freq[max_color];
                    for (int i = max_color - 1; i >= 0; i--) {
                        psum += freq[i];
                        if (float(psum) / float(height * width) >= coeff) {
                            right_border = i + 1;
                            break;
                        }
                    }
                    for (int i = 0; i <= max_color; i++) {
                        if (i < left_border) {
                            stretched[i] = 0;
                        } else if (i > right_border) {
                            stretched[i] = max_color;
                        } else {
                            float res = float(i - left_border) / float(right_border - left_border) * max_color;
                            stretched[i] = (unsigned char)(res);
                        }
                    }
                }
                #pragma omp barrier
                #pragma omp for schedule(static)
                for (int i = 0; i < height * width; i++) {
                    pic[i] = stretched[pic[i]];
                }
            }
            double after_algo_time = omp_get_wtime();
            printf("Time (%i thread(s)): %g ms\n", threads, (after_algo_time - before_algo_time) * 1000);
        #else
            auto start = std::chrono::system_clock::now();
            for (int i = 0; i < height * width; i++) {
                freq[pic[i]]++;
            }
            int psum = freq[0];
            unsigned char left_border = max_color;
            unsigned char right_border = 0;
            for (int i = 1; i <= max_color; i++) {
                psum += freq[i];
                if (float(psum) / float(height * width) >= coeff) {
                    left_border = i - 1;
                    break;
                }
            }
            psum = freq[max_color];
            for (int i = max_color - 1; i >= 0; i--) {
                psum += freq[i];
                if (float(psum) / float(height * width) >= coeff) {
                    right_border = i + 1;
                    break;
                }
            }
            for (int i = 0; i <= max_color; i++) {
                if (i < left_border) {
                    stretched[i] = 0;
                } else if (i > right_border) {
                    stretched[i] = max_color;
                } else {
                    float res = float(i - left_border) / float(right_border - left_border) * max_color;
                    stretched[i] = (unsigned char)(res);
                }
            }
            for (int i = 0; i < height * width; i++) {
                pic[i] = stretched[pic[i]];
            }
            auto end = std::chrono::system_clock::now();
            std::chrono::duration<double> diff = end - start;
            printf("Time (no threads): %g ms\n", diff.count() * 1000);
        #endif
        // End of algorithm
        // Output
        ofstream output;
        output.open(output_file_name);
        output << "P5\n" << width << " " << height << "\n255\n";
        for (int i = 0; i < width * height; i++) {
            output.put(char(pic[i]));
        }
        output.close();
        return 0;
    }
    vector<unsigned char> pic(height * width * 3);
    try {
        input.get(); // The first symbol is always garbage
        for (int i = 0; i < height * width * 3; i++) {
            pic[i] = input.get();
        }
    } catch (ios::failure &e) {
        cout << "Error while reading the image :(" << endl;
        input.close();
        return 0;
    }
    input.close();
    // Algorithm
    vector<int> freq_r(max_color + 1, 0);
    vector<int> freq_g(max_color + 1, 0);
    vector<int> freq_b(max_color + 1, 0);
    vector<unsigned char> stretched(max_color + 1);
    #ifdef _OPENMP
        double before_algo_time = omp_get_wtime();
        #pragma omp parallel default(none) shared(height, width, max_color, pic, freq_r, freq_b, freq_g, coeff, stretched)
        {
            vector<int> thread_freq_r(max_color + 1, 0);
            vector<int> thread_freq_g(max_color + 1, 0);
            vector<int> thread_freq_b(max_color + 1, 0);
            #pragma omp for schedule(static)
            for (int i = 0; i < height * width * 3; i += 3) {
                thread_freq_r[pic[i]]++;
                thread_freq_g[pic[i + 1]]++;
                thread_freq_b[pic[i + 2]]++;
            }
            #pragma omp critical
            {
                for (int i = 0; i <= max_color; i++) {
                    freq_r[i] += thread_freq_r[i];
                    freq_g[i] += thread_freq_g[i];
                    freq_b[i] += thread_freq_b[i];
                }
            }
            #pragma omp barrier
            #pragma omp single
            {
                int psum = freq_r[0];
                unsigned char left_border_r = max_color;
                unsigned char right_border_r = 0;
                for (int i = 1; i <= max_color; i++) {
                    psum += freq_r[i];
                    if (float(psum) / float(height * width) >= coeff) {
                        left_border_r = i - 1;
                        break;
                    }
                }
                psum = freq_r[max_color];
                for (int i = max_color - 1; i >= 0; i--) {
                    psum += freq_r[i];
                    if (float(psum) / float(height * width) >= coeff) {
                        right_border_r = i + 1;
                        break;
                    }
                }

                psum = freq_g[0];
                unsigned char left_border_g = max_color;
                unsigned char right_border_g = 0;
                for (int i = 1; i <= max_color; i++) {
                    psum += freq_g[i];
                    if (float(psum) / float(height * width) >= coeff) {
                        left_border_g = i - 1;
                        break;
                    }
                }
                psum = freq_g[max_color];
                for (int i = max_color - 1; i >= 0; i--) {
                    psum += freq_g[i];
                    if (float(psum) / float(height * width) >= coeff) {
                        right_border_g = i + 1;
                        break;
                    }
                }

                psum = freq_b[0];
                unsigned char left_border_b = max_color;
                unsigned char right_border_b = 0;
                for (int i = 1; i <= max_color; i++) {
                    psum += freq_b[i];
                    if (float(psum) / float(height * width) >= coeff) {
                        left_border_b = i - 1;
                        break;
                    }
                }
                psum = freq_b[max_color];
                for (int i = max_color - 1; i >= 0; i--) {
                    psum += freq_b[i];
                    if (float(psum) / float(height * width) >= coeff) {
                        right_border_b = i + 1;
                        break;
                    }
                }
                unsigned char left_border = min(left_border_b, min(left_border_g, left_border_r));
                unsigned char right_border = max(right_border_b, max(right_border_g, right_border_r));
                for (int i = 0; i <= max_color; i++) {
                    if (i < left_border) {
                        stretched[i] = 0;
                    } else if (i > right_border) {
                        stretched[i] = max_color;
                    } else {
                        float res = float(i - left_border) / float(right_border - left_border) * max_color;
                        stretched[i] = (unsigned char)(res);
                    }
                }
            }
            #pragma omp barrier
            #pragma omp for schedule(static)
            for (int i = 0; i < height * width * 3; i++) {
                pic[i] = stretched[pic[i]];
            }
        }
        double after_algo_time = omp_get_wtime();
        printf("Time (%i thread(s)): %g ms\n", threads, (after_algo_time - before_algo_time) * 1000);
    #else
        auto start = std::chrono::system_clock::now();
        for (int i = 0; i < height * width * 3; i += 3) {
            freq_r[pic[i]]++;
            freq_g[pic[i + 1]]++;
            freq_b[pic[i + 2]]++;
        }
        int psum = freq_r[0];
        unsigned char left_border_r = max_color;
        unsigned char right_border_r = 0;
        for (int i = 1; i <= max_color; i++) {
            psum += freq_r[i];
            if (float(psum) / float(height * width) >= coeff) {
                left_border_r = i - 1;
                break;
            }
        }
        psum = freq_r[max_color];
        for (int i = max_color - 1; i >= 0; i--) {
            psum += freq_r[i];
            if (float(psum) / float(height * width) >= coeff) {
                right_border_r = i + 1;
                break;
            }
        }

        psum = freq_g[0];
        unsigned char left_border_g = max_color;
        unsigned char right_border_g = 0;
        for (int i = 1; i <= max_color; i++) {
            psum += freq_g[i];
            if (float(psum) / float(height * width) >= coeff) {
                left_border_g = i - 1;
                break;
            }
        }
        psum = freq_g[max_color];
        for (int i = max_color - 1; i >= 0; i--) {
            psum += freq_g[i];
            if (float(psum) / float(height * width) >= coeff) {
                right_border_g = i + 1;
                break;
            }
        }

        psum = freq_b[0];
        unsigned char left_border_b = max_color;
        unsigned char right_border_b = 0;
        for (int i = 1; i <= max_color; i++) {
            psum += freq_b[i];
            if (float(psum) / float(height * width) >= coeff) {
                left_border_b = i - 1;
                break;
            }
        }
        psum = freq_b[max_color];
        for (int i = max_color - 1; i >= 0; i--) {
            psum += freq_b[i];
            if (float(psum) / float(height * width) >= coeff) {
                right_border_b = i + 1;
                break;
            }
        }
        unsigned char left_border = min(left_border_b, min(left_border_g, left_border_r));
        unsigned char right_border = max(right_border_b, max(right_border_g, right_border_r));
        for (int i = 0; i <= max_color; i++) {
            if (i < left_border) {
                stretched[i] = 0;
            } else if (i > right_border) {
                stretched[i] = max_color;
            } else {
                float res = float(i - left_border) / float(right_border - left_border) * max_color;
                stretched[i] = (unsigned char)(res);
            }
        }
        for (int i = 0; i < height * width * 3; i++) {
            pic[i] = stretched[pic[i]];
        }
        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> diff = end - start;
        printf("Time (no threads): %g ms\n", diff.count() * 1000);
    #endif
    // End of algorithm
    // Output
    ofstream output;
    output.open(output_file_name, ios::binary);
    output << "P6\n" << width << " " << height << "\n255\n";
    for (int i = 0; i < width * height * 3; i++) {
        output.put(char(pic[i]));
    }
    output.close();
    return 0;
}
