#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <functional>

using namespace std;

constexpr double PI = atan(1) * 4;

template <typename FunTip>
double RombergIntegration(FunTip f, double a, double b, double eps = 1e-8,
int nmax = 1000000, int nmin = 50);

template <typename FunTip>
double AdaptiveIntegration(FunTip f, double a, double b, double eps = 1e-10,
int maxdepth = 30, int nmin = 1);

class ChebyshevApproximation {
    vector<double> c;
    int m;
    double xmin, xmax;
    ChebyshevApproximation(vector<double> c, double xmin, double xmax)
    : c(c), m(c.size() - 1), xmin(xmin), xmax(xmax) {}
public:
    template <typename FunTip>
    ChebyshevApproximation(FunTip f, double xmin, double xmax, int n);
    void set_m(int m);
    void trunc(double eps);
    double operator()(double x) const;
    double derivative(double x) const;
    ChebyshevApproximation derivative() const;
    ChebyshevApproximation antiderivative() const;
    double integrate(double a, double b) const;
    double integrate() const;
};

int main()
{

    try {
        ChebyshevApproximation ch([](double x) { return x; }, 0, -1, 10);
    } catch(exception &e) {
        cout << e.what() << endl;
    }
    try {
        ChebyshevApproximation ch([](double x) { return x; }, 0, 1, 0);
    } catch(exception &e) {
        cout << e.what() << endl;
    }
    ChebyshevApproximation ch([](double x) { return 3 * x * x * x + 2 * x * x + 5; }, 0, 10, 10);
    try {
        ch.set_m(11);
    } catch(exception &e) {
        cout << e.what() << endl;
    }
    ch.set_m(9);
    try {
        cout << ch(-1) << endl;
    } catch(exception &e) {
        cout << e.what() << endl;
    }
    cout << ch(1) << endl;
    cout << ch.derivative(1) << endl;
    cout << ch.derivative()(1) << endl;
    ChebyshevApproximation ch2([](double x) { return 1. / x; }, 1, 2.71, 16);
    ch2.trunc(0.0001);

    cout << ch2(1) << endl;
    cout << ch2.antiderivative()(2.71) - ch2.antiderivative()(1) << endl;
    cout << ch2.integrate(1, 2.71) << endl;
    cout << ch2.integrate() << endl;

    ChebyshevApproximation ch3([](double x) { return sin(x); }, 0, PI, 16);

    cout << ch3(PI / 2) << endl;
    cout << ch3.antiderivative()(PI / 2) - ch3.antiderivative()(0) << endl;
    cout << ch3.integrate(PI / 4, PI / 2) << endl;
    cout << ch3.integrate() << endl;

    try {
        cout << ch3.integrate(-PI / 2, PI) << endl;
    } catch(exception &e) {
        cout << e.what() << endl;
    }
    cout << RombergIntegration([](double x) { return sin(x); }, 0, PI / 2) << endl;
    cout << RombergIntegration([](double x) { return exp(x); }, 0, 1) << endl;
    cout << RombergIntegration([](double x) { return 3 * x * x; }, 1, 3) << endl;

    cout << AdaptiveIntegration([](double x) { return sin(x); }, 0, PI / 2) << endl;
    cout << AdaptiveIntegration([](double x) { return exp(x); }, 0, 1) << endl;
    cout << AdaptiveIntegration([](double x) { return 3 * x * x; }, 1, 3) << endl;

    return 0;
}

template <typename FunTip>
ChebyshevApproximation::ChebyshevApproximation(FunTip f, double xmin, double xmax, int n)
: c(n + 1), m(n), xmin(xmin), xmax(xmax) {
    if(xmin >= xmax || n < 1)
        throw domain_error("Bad parameters");
    vector<double> w(n + 1), v(n + 1);
    for(int i = 0; i <= n; i++) {
        w[i] = PI * (2 * i + 1) / (2 * n + 2);
        v[i] = f((xmin + xmax  + (xmax - xmin) * cos(w[i])) / 2);
    }
    for(int k = 0; k <= n; k++) {
        double s(0);
        for(int i = 0; i <= n; i++)
            s += v[i] * cos(k * w[i]);
        c[k] = 2 * s / (n + 1);
    }
}

void ChebyshevApproximation::set_m(int m) {
    if(m > c.size() - 1 || m < 1)
        throw domain_error("Bad order");
    this->m = m;
}

void ChebyshevApproximation::trunc(double eps) {
    for(int i = m; i >= 0; i--)
        if(fabs(c[i]) > eps) {
            m = i;
            break;
        }
}

double ChebyshevApproximation::operator()(double x) const {
    if(x < xmin || x > xmax)
        throw domain_error("Bad argument");
    double t((2 * x - xmin - xmax) / (xmax - xmin)),
    p(1), q(t), s(c[0] / 2 + c[1] * t), r;
    for(int k = 2; k <= m; k++) {
        r = 2 * t * q - p;
        s += c[k] * r;
        p = q;
        q = r;
    }
    return s;
}

double ChebyshevApproximation::derivative(double x) const {
    double t((2 * x - xmin - xmax) / (xmax - xmin)),
    p(1), q (4 * t), s(c[1] + 4 * c[2] * t), r;
    for(int k = 3; k <= m; k++) {
        r = k * (2 * t * q / (k - 1) - p / (k - 2));
        s += c[k] * r;
        p = q;
        q = r;
    }
    return 2 * s / (xmax - xmin);
}

ChebyshevApproximation ChebyshevApproximation::derivative() const {
    double mi(4. / (xmax - xmin));
    vector<double> c2(c.size());
    c2[m - 1] = mi * m * c[m];
    c2[m - 2] = mi * (m - 1) * c[m - 1];
    for(int k = m - 3; k >= 0; k--)
        c2[k] = c2[k + 2] + mi * (k + 1) * c[k + 1];
    c2.resize(m - 1);
    return ChebyshevApproximation(c2, xmin, xmax);
}

ChebyshevApproximation ChebyshevApproximation::antiderivative() const {
    double mi((xmax - xmin) / 4);
    vector<double> c2(m + 1);
    c2[0] = 0;
    for(int k = 1; k <= m - 1; k++)
        c2[k] = mi / k * (c[k - 1] - c[k + 1]);
    c2[m] = mi / m * c[m - 1];
    c2[m + 1] = mi / (m + 1) * c[m];
    return ChebyshevApproximation(c2, xmin, xmax);
}

double ChebyshevApproximation::integrate(double a, double b) const {
    if(a < xmin || a > xmax || b < xmin || b > xmax)
        throw domain_error("Bad interval");
    auto F(antiderivative());
    return F(b) - F(a);
}

double ChebyshevApproximation::integrate() const {
    double s(0);
    for(int k = 1; k <= (m - 1) / 2; k++)
        s += 2 * c[2 * k] / (1 - 4 * k * k);
    s *= (xmax - xmin) / 2;
    s += c[0] * (xmax - xmin) / 2;
    return s;
}

template <typename FunTip>
double RombergIntegration(FunTip f, double a, double b, double eps = 1e-8,
int nmax = 1000000, int nmin = 50) {
    int N(2);
    double h((b - a) / N), s((f(a) + f(b)) / 2), Iold(s);
    vector<double> I;
    for(int i = 1; N <= nmax; i++) {
        for(int j = 1; j <= N / 2; j++)
            s += f(a + (2 * j - 1) * h);
        I.push_back(h * s);
        double p(4);
        for(int k = I.size() - 2; k >= 0; k--) {
            I[k] = (p * I[k + 1] - I[k]) / (p - 1);
            p *= 4;
        }
        if(N >= nmin && fabs(I[0] - Iold) <= eps)
            return I[0];
        Iold = I[0];
        h /= 2;
        N *= 2;
    }
    return Iold;
}

template <typename FunTip>
double AdaptiveAux(FunTip f, double a, double b, double eps, double f1, double f2, double f3, int R) {
    double c((a + b) / 2), I1((b - a) * (f1 + 4 * f3 + f2) / 6), f4(f((a + c) / 2)),
    f5(f((c + b) / 2)), I2((b - a) * (f1 + 4 * f4 + 2 * f3 + 4 * f5 + f2) / 12);
    if(R <= 0 || fabs(I1 - I2) <= eps)
        return I2;
    return AdaptiveAux(f, a, c, eps, f1, f3, f4, R - 1) + AdaptiveAux(f, c, b, eps, f3, f2, f5, R - 1);
}

template <typename FunTip>
double AdaptiveIntegration(FunTip f, double a, double b, double eps = 1e-10,
int maxdepth = 30, int nmin = 1) {
    double s(0), h((b - a) / nmin);
    for(int i = 1; i <= nmin; i++) {
        s += AdaptiveAux(f, a, a + h, eps, f(a), f(a + h), f(a + h / 2), maxdepth);
        a += h;
    }
    return s;
}


