#include <iostream>
#include <algorithm>
#include <vector>
#include <limits>
#include <stdexcept>
#include <cmath>
#include <functional>
#include <fstream>
#include <string>
#include <sstream>
#include <limits>
#define PI 4 * atan(1.)

using namespace std;

template <typename FpType>
bool fpequal(FpType x, FpType y);

template <typename FpType>
bool between(FpType x1, FpType x2, FpType x3);

template <typename Tip>
string to_string(Tip x);

vector<double> linspace(double xmin, double xmax, int n);
int find_range(const vector<pair<double, double>> &data, int cached_index, double x);
vector<double> feval(const vector<double>& x, function<double(double)> f);
void print_vector_to_file(const vector<double> &v, ofstream &o);
void create_scilab(const vector<double>&, const vector<double>&, const vector<double>&,
                   const vector<double>&, const vector<double>&, string);
void test(vector<double> xdata, vector<double> xinterp, function<double(double)> f,
          string type = "linear", int bary_order = -1);

class Comparator {
public:
    bool operator()(const pair<double, double> &lhs, const pair<double, double> &rhs) { return lhs.first < rhs.first; }
};

class LinearInterpolator {
    mutable int cached_index;
    vector<pair<double, double>> data;
public:
    LinearInterpolator(vector<pair<double, double>> _data);
    double operator()(double x) const;
};

class PolynomialInterpolator {
    vector<double> x, y, yn;
public:
    PolynomialInterpolator(std::vector<std::pair<double, double>> data);
    double operator()(double x) const;
    void AddPoint(std::pair<double, double> p);
    std::vector<double> GetCoefficients() const;
};

class SplineInterpolator {
    mutable int cached_index;
    vector<pair<double, double>> data;
    vector<double> r, q, s;
public:
    SplineInterpolator(std::vector<std::pair<double, double>> data);
    double operator()(double x) const;
};

class BarycentricInterpolator {
    vector<pair<double, double>> data;
    vector<double> w;
public:
    BarycentricInterpolator(std::vector<std::pair<double, double>> data, int order);
    double operator()(double x) const;
    std::vector<double> GetWeights() const { return w; };
};

template <typename FunTip>
double Limit(FunTip f, double x0, double eps = 1e-8, double nmax = 20);

int main()
{
    /*
     *  Pošto je u tekstu zadaæe naèin testiranja dat kao prijedlog
     *  i nije obavezan, testiranje je ovdje vršeno na drugi naèin,
     *  generišuæi scilab skripte za svaki poziv funkcije test
     *  koje iscrtavaju grafove funkcija, interpoliranih verzija
     *  tih funkcija i, onda kada je to moguæe,
     *  grafove interpolacija dobivenih korištenjem
     *  scilabove funkcije interp. Skripte imaju naziv u formatu
     *  "interp_i.sci" i generišu se u istom folderu u kojem je
     *  ovaj fajl.
     */

    if(false) { //Staviti true prilikom ruènog testiranja
        typedef pair<double, double> Point;

        //Linearna interpolacija
        //Skripte 1 - 4

        //Ponavljanje istog x
        try {
           LinearInterpolator(vector<Point>{Point(0, 0), Point(1, 1), Point(0, 1)});
        } catch(exception& e) {
           cout << e.what() << endl;
        }

        test(linspace(-1, 1, 2), linspace(-1, 1, 100), [](double x) { return 3 * x + 2; });
        test(linspace(-1, 1, 10), linspace(-1, 1, 100), [](double x) { return pow(x, 3) + 3 * pow(x, 2); });
        test(linspace(-2 * PI, 2 * PI , 20), linspace(-2 * PI, 2 * PI , 100), [](double x) { return sin(x); });
        test(linspace(-1, 1, 2), linspace(-1, 1, 100), [](double x) { return 3 * x + 2; });

        //Polinomijalna interpolacija
        //Skritpe 5 - 8

        //Dijeljenje nulom
        try {
           PolynomialInterpolator(vector<Point>{Point(0, 0), Point(1, 1), Point(0, 1)});
        } catch(exception& e) {
           cout << e.what() << endl;
        }

        test(linspace(-3, 3, 5), linspace(-3, 3, 100), [](double x) { return
         4 * pow(x, 4) + 3 * pow(x, 3) + 3 * pow(x, 2) + 4 * pow(x, 1); }, "poly");
        test(linspace(-3, 3, 6), linspace(-3, 3, 100), [](double x) { return
         5 * pow(x, 5) - 20 * pow(x, 4) + 5 * pow(x, 3) + 50 * pow(x, 2) - 20 * pow(x, 1) - 40; }, "poly");
        test(linspace(-2, 2, 6), linspace(-2, 2, 100), [](double x) { return
         - 1 * pow(x, 5) + 4 * pow(x, 3); }, "poly");
        test(linspace(-PI, PI, 10), linspace(-PI, PI, 100), [](double x) { return
         tan(x); }, "poly");

        //Test dodavanje nove taèke i raèunanja koeficijenata
        auto f([](double x) { return 7 * pow(x, 6) - 1 * pow(x, 5) + 4 * pow(x, 3); });
        vector<pair<double, double>> data;
        vector<double> yinterp;
        for(auto x : linspace(-2, 1.5, 6))
            data.push_back(pair<double, double>(x, f(x)));
        PolynomialInterpolator pi(data);
        pi.AddPoint(Point(2, f(2)));
        auto c(pi.GetCoefficients());
        double eps = 0.0001;
        bool found_non_zero(false);
        for(int i = 0; i < c.size(); i++) {
            if(fabs(c[i]) > eps) {
                if(i != 0 && found_non_zero)
                    cout << (c[i] > 0 ? " + " : " - ");
                cout << to_string(fabs(c[i]))  + string(" * x^") + to_string(i);
                found_non_zero = true;
            }
            if(i == c.size() - 1)
                cout << endl;
        }

        //Kubni splajn
        //Skripte 9 - 12

        /*
         *  Prilikom testiranja uoèene su dvije potencijalne greške koje se nalaze u pseudokodu iz predavanja.
         *
         *  Prva greška je u kodu za raèunanje vrijednosti vektora q. U drugoj verziji koda u predavanju,
         *  u onoj gdje se raèunaju vrijednosti odreðenog elementa vektora q samo onda kada je to potrebno,
         *  jedan od sabiraka se dijeli sa 3, dok to nije sluèaj u prvoj verziji koda, kada se svi elementi raèunaju
         *  zajedno. Ovdje je korištena druga verzija jer dovodi do ispravnijih rezultata
         *  u što se može uvjeriti mijenjajuæi kod u liniji 321, pokreæuæi bilo koju od skripti za crtanje splajnova
         *  i analizirajuæi graf.
         *
         *  Druga potencijalna greška je indeks do kojeg ide prva for petlja u pseudokodu za konstrukciju splajna.
         *  Ako se upotrijebi granica navedena u predavanju tj. indeks n - 2 (ovdje u kodu n - 3 zbog offseta), splajn
         *  pravi vidno odstupanje i od prave funkcije i od scilabove interpolacije. Meðutim, ako se ta granica
         *  pomjeri na n - 1 (ovdje u kodu n - 2), greška se otklanja. U ovo se može uvjeriti na isti naèin
         *  mijenjajuæi kod u liniji 306 i pokreèuæi skriptu "interp_9.sci"
         *
         *
         */

        test(linspace(-2 * PI, 2 * PI , 10), linspace(-2 * PI, 2 * PI , 100), [](double x) { return sin(x); }, "spline");
        test(linspace(-2 * PI, 2 * PI , 20), linspace(-2 * PI, 2 * PI , 100), [](double x) { return sinh(x); }, "spline");
        test(linspace(-2 * PI, 2 * PI , 20), linspace(-2 * PI, 2 * PI , 100), [](double x) { return x * sin(1 / x); }, "spline");
        test(linspace(-2 * PI, 2 * PI , 20), linspace(-2 * PI, 2 * PI , 100), [](double x) { return sin(1 / x); }, "spline");

        //Baricentrièna interpolacija
        //Skripte 13 - 16

        test(linspace(-1, 1, 10), linspace(-1, 1, 50), [](double x) { return exp(x); }, "barycentric", 2);
        test(linspace(-1, 1, 20), linspace(-1, 1, 50), [](double x) { return 10 * x; }, "barycentric", 5);
        test(linspace(-2 * PI, 2 * PI , 20), linspace(-2 * PI, 2 * PI , 100), [](double x) { return x * sin(1 / x); }, "barycentric", 0);
        test(linspace(0, 2 * PI , 15), linspace(0, 2 * PI , 100), [](double x) { return sin(3.5 * x + 2); }, "barycentric", 3);

        auto f2([](double x) { return cos(x); });
        data.clear(); yinterp.clear();
        for(auto x : linspace(0, 2 * PI , 10))
            data.push_back(pair<double, double>(x, f(x)));

        //Dijeljeje nulom
        try {
            BarycentricInterpolator bi(vector<Point>{Point(0, 0), Point(0, 1), Point(0, 2)}, 2);
        }  catch(exception& e) {
            cout << e.what() << endl;
        }

        //Neispravan red
        try {
            BarycentricInterpolator bi(data, -10);
        }  catch(exception& e) {
            cout << e.what() << endl;
        }

        //Tezine
        BarycentricInterpolator bi(data, 3);
        for(auto w : bi.GetWeights())
            cout << w << " ";
        cout << endl;

        //Limesi

        //e
        cout << Limit([](double x) { return  pow((1 +  x), 1 / x); }, 0) << endl;
        //1
        cout << Limit([](double x) { return  sin(x) / x; }, 0) << endl;
        //ln2
        cout << Limit([](double x) { return  (pow(2, x) - 1) / x; }, 0) << endl;
        //Ne postoji limes
        try {
            cout << Limit([](double x) { return  sin(1 / x); }, 0) << endl;
        } catch (exception& e) {
            cout << e.what() << endl;
        }
    }
    return 0;
}

LinearInterpolator::LinearInterpolator(vector<pair<double, double>> _data) : cached_index(-1), data(_data) {
    sort(data.begin(), data.end(), Comparator());
    for(int i = 0; i < data.size() - 1; i++)
        if(fpequal(data[i].first, data[i + 1].first))
            throw domain_error("Invalid data set");
}

double LinearInterpolator::operator()(double x) const {
    int i = cached_index = find_range(data, cached_index, x);
    double x1(data[i].first), x2(data[i + 1].first), y1(data[i].second), y2(data[i + 1].second);
    return (x2 - x) / (x2 - x1) * y1 + (x - x1) / (x2 - x1) * y2;
}

PolynomialInterpolator::PolynomialInterpolator(std::vector<std::pair<double, double>> data)
: x(data.size()), y(data.size()), yn(data.size()) {
    int n(data.size());
    for(int i = 0; i < n; i++) {
        x[i] = data[i].first;
        y[i] = data[i].second;
    }
    yn[0] = y[n - 1];
    for(int j = 1; j <= n - 1; j++) {
        for(int i = n; i >= j + 1; i--) {
            if(fpequal(x[i - 1], x[i - j - 1]))
                throw domain_error("Invalid data set");
            y[i - 1] = (y[i - 1] - y[i - 2]) / (x[i - 1] - x[i - j - 1]);
        }
        yn[j] = y[n - 1];
    }
}

double PolynomialInterpolator::operator()(double x) const {
    double s(y[0]), f(1);
    int n(y.size());
    for(int j = 1; j <= n - 1; j++) {
        f *= (x - this->x[j - 1]);
        s += y[j] * f;
    }
    return s;
}

void PolynomialInterpolator::AddPoint(std::pair<double, double> p) {
    double x_new(p.first), y_new(p.second);
    int n(x.size() + 1);
    x.push_back(x_new);
    y.push_back(y_new);
    yn.resize(n);
    for(int i = 1; i <= n - 1; i++) {
        double yn_temp(yn[i - 1]);
        yn[i - 1] = y[n - 1];
        y[n - 1] = (y[n - 1] - yn_temp) / (x[n - 1] - x[n - i - 1]);
    }
    yn[n - 1] = y[n - 1];
}


std::vector<double> PolynomialInterpolator::GetCoefficients() const {
    int n(x.size());
	std::vector<double> w(n + 1), p(n + 1);
	w[0] = 1;
	for (int i = 1; i <= n; i++)
	{
		for (int j = 0; j <= i; j++)
			p[j] += y[i - 1] * w[j];
		w[i] = w[i - 1];
		for (int j = i - 1; j > 0; j--)
			w[j] = w[j - 1] - x[i - 1] * w[j];
		w[0] *= -x[i - 1];
	}

	return p;
}


SplineInterpolator::SplineInterpolator(std::vector<std::pair<double, double>> _data)
: cached_index(-1), data(_data), r(data.size()), s(data.size()), q(data.size()) {
    sort(data.begin(), data.end(), Comparator());
    int n(data.size());
    for(int i = 1; i <= n - 2; i++) { //U pseudokodu u predavanju ide do n - 3 (nakon sto se oduzme offset)
        s[i] = 2 * (data[i + 1].first - data[i - 1].first);
        r[i] = 3 * ((data[i + 1].second - data[i].second) / (data[i + 1].first - data[i].first) - (data[i].second - data[i - 1].second) / (data[i].first - data[i - 1].first));
    }
    for(int i = 1; i <= n - 3; i++) {
        double mi((data[i].first - data[i - 1].first) / s[i]);
        s[i + 1] -= mi * (data[i].first - data[i - 1].first);
        r[i + 1] -= mi * r[i];
    }
    r[n - 2] /= s[n - 2];
    for(int i = n - 3; i >= 1; i--)
        r[i] = (r[i] - (data[i].first - data[i - 1].first) * r[i + 1]) / s[i];
    for(int i = 0; i <= n - 2; i++) {
        double dx(data[i + 1].first - data[i].first);
        s[i] = (r[i + 1] - r[i]) / (3 * dx);
        q[i] = (data[i + 1].second - data[i].second) / dx - dx * (r[i + 1] + 2 * r[i]) / 3;
    }
}

double SplineInterpolator::operator()(double x) const {
    int i = cached_index = find_range(data, cached_index, x);
    double t(x - data[i].first), y(data[i].second);
    return y + t * (q[i] + t * (r[i] + s[i] * t));
}

BarycentricInterpolator::BarycentricInterpolator(std::vector<std::pair<double, double>> _data, int order)
:data(_data), w(data.size()) {
    sort(data.begin(), data.end(), Comparator());
    int n(data.size());
    if(!(order >= 0 && order <= n))
        throw domain_error("Invalid order");
    for(int i = 1; i <= n; i++) {
        double p;
        for(int k = max(1, i - order); k <= min(i, n - order); k++) {
            p = 1;
            for(int j = k; j <= k + order; j++)
                if(j != i) {
                    auto xi(data[i - 1].first), xj(data[j - 1].first);
                    if(fpequal(xi, xj))
                        throw domain_error("Invalid data set");
                    p /= (xi - xj);
                }
            if(k % 2 == 0)
                p = -p;
        }
        w[i - 1] += p;
    }
}

double BarycentricInterpolator::operator()(double x) const {
    double p(0), q(0);
    for(int i = 1; i <= data.size(); i++) {
        if(fpequal(x, data[i - 1].first))
            return data[i - 1].second;
        double u(w[i - 1] / (x - data[i - 1].first));
        p += u * data[i - 1].second;
        q += u;
    }
    return p / q;
}

template <typename FunTip>
double Limit(FunTip f, double x0, double eps = 1e-8, double nmax = 20) {
    double p, yold(numeric_limits<double>::infinity()), h(1);
    vector<double> y(nmax);
    for(int i = 1; i <= nmax; i++) {
        y[i - 1] = f(x0 + h);
        p = 2;
        for(int k = i - 1; k >= 1; k--) {
            y[k - 1] = (p * y[k] - y[k - 1]) / (p - 1);
            p *= 2;
        }
        if(fabs(y[0] - yold) < eps)
            return y[0];
        yold = y[0];
        h /= 2;
    }
    throw logic_error("Accuracy goal is not achieved");
}

template <typename FpType>
bool fpequal(FpType x, FpType y) {
     FpType eps = 10 * std::numeric_limits<FpType>::epsilon()
        * (std::fabs(x) + std::fabs(y));
     return std::fabs(x - y) <= eps;
}

template <typename FpType>
bool between(FpType x1, FpType x2, FpType x3) {
     return ((x1 < x2) || fpequal(x1 , x2)) && ((x2 < x3) || fpequal(x1 , x2));
}

template <typename Tip>
string to_string(Tip x) {
    return static_cast<std::ostringstream&>((std::ostringstream() << std::dec << x)).str();
}


int find_range(const vector<pair<double, double>> &data, int cached_index, double x) {
    if(cached_index != -1) {
        if(between(data[cached_index].first, x, data[cached_index + 1].first))
            return cached_index;
        if(cached_index > 0 && between(data[cached_index - 1].first, x, data[cached_index].first))
            return cached_index - 1;
        if(cached_index < data.size() - 2 && between(data[cached_index + 1].first, x, data[cached_index + 2].first))
            return cached_index + 1;
    }
    if(x <= data[1].first)
        return 0;
    if(x >= data[data.size() - 2].first)
        return data.size() - 2;
    int first(0), last(data.size() - 2), middle;
    while(first <= last) {
        middle = (first + last) / 2;
        if(between(data[middle].first, x, data[middle + 1].first))
            return middle;
        if(data[middle].first > x)
            last = middle;
        else
            first = middle;
    }
}

vector<double> linspace(double xmin, double xmax, int n) {
    vector<double> x;
    double step((xmax - xmin) / (n - 1));
    for(double i = xmin; i < xmax || fpequal(i, xmax); i += step)
        x.push_back(i);
    return x;
}

vector<double> feval(const vector<double>& x, function<double(double)> f) {
    vector<double> y;
    for(double xi : x)
        y.push_back(f(xi));
    return y;
}

void print_vector_to_file(const vector<double> &v, ofstream &o) {
    o << " [ ";
    for(auto t : v)
        o << t << " ";
    o << "]";
}

void create_scilab(const vector<double>& xdata, const vector<double>& ydata,
                   const vector<double>& xinterp, const vector<double>& yinterp,
                   const vector<double>& yactual, string draw_type) {
    static int counter(1);
    ofstream output("interp_" + to_string(counter++) + ".sci");
    output << "xdata = ";
    print_vector_to_file(xdata, output);
    output << ";" << endl;
    output << "ydata = ";
    print_vector_to_file(ydata, output);
    output << ";" << endl;
    output << "xinterp = ";
    print_vector_to_file(xinterp, output);
    output << ";" << endl;
    output << "yinterp = ";
    print_vector_to_file(yinterp, output);
    output << ";" << endl;
    output << "yinterpactual = ";
    print_vector_to_file(yactual, output);
    output << ";" << endl;
    if(draw_type != "actual") {
        if(draw_type == "linear")
            output << "[yinterpsci] = interp1(xdata, ydata, xinterp, \"" << draw_type << "\");" << endl;
        else if(draw_type == "spline")
            output << "d = splin(xdata, ydata, \"natural\");\n[yinterpsci] = interp(xinterp, xdata, ydata, d);" << endl;
    }
    output << "dif = yinterpactual - yinterp;\n[n, m] = size(dif);"
    "\nmaxd = 0;\nfor i = 1:m\nif(abs(dif(1, i)) > maxd)\nmaxd = abs(dif(1, i));\n"
    "end\nend\ndisp(\"Max abs error: \" + string(maxd));\nplot(xdata, ydata, \'x\');\n"
    "plot(xinterp, yinterpactual);\n" <<
    ((draw_type != "actual") ? "plot(xinterp, yinterpsci, 'g');\n" : "") <<
    "plot(xinterp, yinterp, 'r');\n" <<
    "legend(\"ydata\", \"actual\"," <<
    (draw_type != "actual" ? " \"scilab interp\"," : "") <<
    " \"my interp\");";
    output.close();
}

void test(vector<double> xdata, vector<double> xinterp, function<double(double)> f, string type, int bary_order) {
    vector<pair<double, double>> data;
    vector<double> ydata;
    for(auto x : xdata) {
        auto y(f(x));
        ydata.push_back(y);
        data.push_back(pair<double, double>(x, y));
    }
    vector<double> yinterp;
    string draw_type;
    if(type == "linear") {
        LinearInterpolator li(data);
        draw_type = type;
        for(auto x : xinterp)
            yinterp.push_back(li(x));
    } else if(type == "poly") {
        PolynomialInterpolator pi(data);
        draw_type = "actual";
        for(auto x : xinterp)
            yinterp.push_back(pi(x));
    } else if(type == "spline") {
        SplineInterpolator si(data);
        draw_type = type;
        for(auto x : xinterp)
            yinterp.push_back(si(x));
    } else if(type == "barycentric") {
        if(bary_order == -1)
            bary_order = xdata.size();
        BarycentricInterpolator bi(data, bary_order);
        draw_type = "actual";
        for(auto x : xinterp)
            yinterp.push_back(bi(x));
    }
    create_scilab(xdata, ydata, xinterp, yinterp, feval(xinterp, f), draw_type);
    double max_dy(0), max_x(0);
    for(int i = 0; i < xinterp.size(); i++)
        if(max_dy < fabs(f(xinterp[i]) - yinterp[i])) {
                max_x = xinterp[i];
                max_dy = max(max_dy, fabs(f(xinterp[i]) - yinterp[i]));
        }
    //cout << max_dy << " at x = " << max_x << endl;
}
