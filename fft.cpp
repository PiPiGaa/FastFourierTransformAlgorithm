#define _USE_MATH_DEFINES
#include <iostream>
#include <complex>
#include <math.h>
#include <string>
#include <vector>
#include <complex>

using namespace std;

typedef complex<double> ftype;
ftype i(0, -1);

//Class_For_FastFourierTransform
class fft {
public:
    //Constructor
    fft(vector<ftype> seq_,double n_, int h_) {
        h = h_ < log2(n_) ? h_ : log2(n_);
        N = n_;
        seq = seq_;
    }

    //Direct transform
    void calcDirFFT() {
        coeff.clear();
        for (int k = 0; k < N; k++) {
            pair<ftype, ftype> s = calc_C(seq, N, h, k);
            coeff.push_back(s.first + s.second);
        }
    }

    //Reverse transform
    void calcRevFFT() {
        newFFTseq.clear();
        for (int k = 0; k < N; k++) {
            pair<ftype, ftype> s = calc_X(coeff, N, h, k);
            newFFTseq.push_back((s.first + s.second)/N);
        }
    }

    //Check for error
    double Check() {
        double error = 0;
        for (int y = 0; y < coeff.size(); y++)
            error += abs(seq[y]) - abs(newFFTseq[y]);
        return error;
    }

private:
    int h; //conversion depth
    double N; //sequence size
    vector<ftype> seq; //input sequence
    vector<ftype> coeff; //frequency coefficients
    vector<ftype> newFFTseq; //output sequence

    //Recursive division into an even and odd sum of a given depth for direct transform
    pair<ftype, ftype> calc_C(vector<ftype> buf, double n_, int h, int k) {
        pair<ftype, ftype> e, o;
        vector<ftype> S_mas_even;
        vector<ftype> S_mas_odd;
        if (h > 1) {
            for (int n = 0; n < n_ / 2; n++) {
                S_mas_even.push_back(buf[n * 2]);
                S_mas_odd.push_back(buf[n * 2 + 1]);
            }
            e = calc_C(S_mas_even, n_ / 2, h - 1, k);
            o = calc_C(S_mas_odd, n_ / 2, h - 1, k);
            pair<ftype, ftype> a(e.first + e.second, exp(i * (2 * M_PI * k / n_)) * (o.first + o.second));
            return a;
        }
        else {
            ftype S_k_even = 0, S_k_odd = 0;
            for (int n = 0; n < n_ / 2; n++) {
                S_k_even += buf[2 * n] * exp(i * (2 * M_PI * 2 * n * k / n_));
                S_k_odd += exp(i * (2 * M_PI * k / n_)) * buf[2 * n + 1] * exp(i * (2 * M_PI * 2 * n * k / n_));  
            }
            pair<ftype, ftype> a(S_k_even, S_k_odd);
            return a;
        }
    }

    //Recursive division into an even and odd sum of a given depth for reverse transform
    pair<ftype, ftype> calc_X(vector<ftype> buf, double n_, int h, int k) {
        pair<ftype, ftype> e, o;
        vector<ftype> S_mas_even;
        vector<ftype> S_mas_odd;
        if (h > 1) {
            for (int n = 0; n < n_ / 2; n++) {
                S_mas_even.push_back(buf[n * 2]);
                S_mas_odd.push_back(buf[n * 2 + 1]);
            }
            e = calc_X(S_mas_even, n_ / 2, h - 1, k);
            o = calc_X(S_mas_odd, n_ / 2, h - 1, k);
            pair<ftype, ftype> a(e.first + e.second, exp(i * (-2 * M_PI * k / n_)) * (o.first + o.second));
            return a;
        }
        else {
            ftype S_k_even = 0, S_k_odd = 0;
            for (int n = 0; n < n_ / 2; n++) {
                S_k_even += buf[2 * n] * exp(i * (-2 * M_PI * 2 * n * k / n_));
                S_k_odd += exp(i * (-2 * M_PI * k / n_)) * buf[2 * n + 1] * exp(i * (-2 * M_PI * 2 * n * k / n_));
            }
            pair<ftype, ftype> a(S_k_even, S_k_odd);
            return a;
        }
    }
};

int main()
{
    //Generating a sequence of random complex numbers
    srand(time(NULL));
    double n = pow(2, 5);
    int h = 1; //any in range 1..n
    vector<ftype> seq;
    for (int i = 0; i < n; i++) {
        ftype item(rand() % 1000, rand() % 1000);
        seq.push_back(item);
    }
    fft gg(seq, n, h);

    //Running fourier transforms
    gg.calcDirFFT();
    gg.calcRevFFT();
    cout << "error: " << gg.Check();

    return 0;
}

