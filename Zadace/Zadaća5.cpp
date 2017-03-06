#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include <limits>
#include <fstream>
#include <algorithm>
#include <map>
#include <complex>
#include <ctime>
#include <cstdlib>

#define SIGN(a, b) b < 0 ? -fabs(a) : fabs(a)
#define INF -9.99e99

using namespace std;

int DEBUG = 0;
int GENERATE_SCILAB = 0;

template <typename FpType>
bool fpequal(FpType x, FpType y) {
	FpType eps = 10 * std::numeric_limits<FpType>::epsilon()
		* (std::fabs(x) + std::fabs(y));
	return std::fabs(x - y) <= eps;
}

template <typename FpType>
bool is_zero(FpType x, double eps = 1e-10) {
	return abs(x) <= eps;
}

template <typename FpType>
bool lequal(FpType x, double eps = 1e-10) {
	return x < 0 || is_zero(x, eps);
}

template <typename FpType>
bool sgn(FpType x) {
	if (is_zero(x))
		return 0;
	else if (x < 0)
		return -1;
	else if (x > 0)
		return 1;
}

void DebugInfo(int i) {
	if (DEBUG)
		cout << "Iteracija: " << i + 1 << endl;
}

template <typename FunTip>
bool BracketRoot(FunTip f, double x0, double &a, double &b,
	double hinit = 1e-5, double hmax = 1e10, double lambda = 1.4) {
	if (hinit <= 0 || hmax <= 0 || lambda <= 0)
		throw domain_error("Invalid parameters");
	double a2(x0), f1(f(a2)), h(hinit);
	do {
		while (fabs(h) < hmax) {
			double b2(a2 + h), f2(f(b2));
			if (lequal(f1 * f2)) {
				a = a2;
				b = b2;
				return true;
			}
			h = lambda * h;
			a2 = b2;
			f1 = f2;
		}
		h *= -1;
	} while (h < 0);
	return false;
}

template <typename FunTip>
double ModifiedRegulaFalsi(FunTip f, double a, double b, double eps = 1e-10,
	int maxiter = 100) {
	if (a >= b)
		throw range_error("Bad range");
	if (f(a) * f(b) > 0)
		throw range_error("Root must be bracketed");
	if (eps < 0 || maxiter < 0)
		throw domain_error("Invalid parameters");
	if (is_zero(f(a), 0))
		return a;
	if (is_zero(f(b), 0))
		return b;
	f = [f](double x) { return f(x) / (1 + abs(f(x))); };
	double f1(f(a)), f2(f(b));
	double c(a), cold(b);
	for (int i = 0; i <= maxiter && abs(c - cold) > eps; i++) {
		DebugInfo(i);
		if (i == maxiter)
			throw logic_error("Given accuracy has not achieved"); //Valjda been achieved
		cold = c;
		c = (a * f2 - b * f1) / (f2 - f1);
		double f3(f(c));
		if (is_zero(f3, 0))
			return c;
		if (lequal(f1 * f3, 0))
			b = c,
			f2 = f3;
		else
			a = c,
			f1 = f3;
	}
	return c;
}

//Kod uzet it Numerical Recipes in C++

template <typename FunTip>
double Ridders2(FunTip f, double a, double b, double eps = 1e-10,
	int maxiter = 100) {
	if (a >= b)
		throw range_error("Bad range");
	if (eps < 0 || maxiter < 0)
		throw domain_error("Invalid parameters");
	double fl(f(a)), fh(f(b));
	if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
		double xl(a), xh(b), ans(INF);
		for (int j = 0; j < maxiter; j++) {
			DebugInfo(j);
			double xm(0.5 * (xl + xh)), fm(f(xm)),
				s(sqrt(fm * fm - fl * fh));
			if (is_zero(s, 0)) return ans;
			double xnew(xm + (xm - xl) * ((fl >= fh) ? 1.0 : -1.0) * fm / s);
			if (abs(xnew - ans) <= eps) return ans;
			ans = xnew;
			double fnew(f(ans));
			if (is_zero(fnew, 0)) return ans;
			if (SIGN(fm, fnew) != fm) {
				xl = xm;
				fl = fm;
				xh = ans;
				fh = fnew;
			}
			else if (SIGN(fl, fnew) != fl) {
				xh = ans;
				fh = fnew;
			}
			else if (SIGN(fh, fnew) != fh) {
				xl = ans;
				fl = fnew;
			}
			if (abs(xh - xl) <= eps) return ans;
		}
		throw logic_error("Given accuracy has not achieved");
	}
	else {
		if (fl == 0.0) return a;
		if (fh == 0.0) return b;
		throw range_error("Root must be bracketed");
	}
}

//Kod iz predavanja

template <typename FunTip>
double Ridders(FunTip f, double a, double b, double eps = 1e-10,
	int maxiter = 100) {
	if (a >= b)
		throw range_error("Bad range");
	if (f(a) * f(b) > 0)
		throw range_error("Root must be bracketed");
	if (eps < 0 || maxiter < 0)
		throw domain_error("Invalid parameters");
	auto f1(f(a)), f2(f(b));
	for (int i = 0; i <= maxiter && fabs(b - a) > eps; i++) {
		DebugInfo(i);
		if (i == maxiter)
			throw logic_error("Given accuracy has not achieved");
		auto c((a + b) / 2), f3(f(c));
		if (f3 == 0)
			return c;
		double sign((f1 - f2) < 0 ? -1 : 1);
		auto d(c + f3 * (c - a) * (sign / sqrt(f3 * f3 - f1 * f2)));
		auto f4(f(d));
		if (f4 == 0)
			return d;
		if (f3 * f4 <= 0) {
			a = c;
			b = d;
			f1 = f3;
			f2 = f4;
		}
		else if (f1 * f4 <= 0) {
			b = d;
			f2 = f4;
		}
		else {
			a = d;
			f1 = f4;
		}
	}
	return (a + b) / 2;
}


template <typename FunTip>
double NewtonRaphson(FunTip f, FunTip fprim, double x0, double eps = 1e-10,
	int maxiter = 100) {
	if (eps < 0 || maxiter < 0)
		throw domain_error("Invalid parameters");
	double dx(INF);
	for (int i = 0; i <= maxiter && fabs(dx) > eps; i++) {
		DebugInfo(i);
		if (i == maxiter)
			throw logic_error("Convergence has not achieved");
		double v(f(x0));
		if (is_zero(fabs(v), eps))
			return x0;
		if (is_zero(fprim(x0), eps))
			throw logic_error("Convergence has not achieved");
		dx = v / fprim(x0);
		x0 -= dx;
	}
	return x0;
}


complex<double> Laguerre(vector<complex<double>> p, int n, complex<double> x,
	bool& c, double eps = 1e-10, double maxiter = 100) {
	complex<double> dx(INF);
	int k(1);
	if (DEBUG) {
		cout << "Aproksimacija pocinje od (" << x.real() << ", " << x.imag() << ")." << endl;
	}
	while (abs(dx) > eps && k < maxiter) {
		DebugInfo(k - 1);
		complex<double> f(p[n]), d(0), s(0);
		for (int i = n - 1; i >= 0; i--) { //Indeks treba ici do 0, ne 2
			s = s * x + 2.0 * d;
			d = d * x + f;
			f = f * x + p[i];
		}

		if (is_zero(f, 0)) {
			c = true;
			return x;
		}
		complex<double> r(sqrt(double(n - 1) * (double(n - 1) * d * d - double(n) * f * s)));
		if (abs(d + r) > abs(d - r))
			dx = double(n) * f / (d + r);
		else
			dx = double(n) * f / (d - r);
		if (DEBUG) {
			cout << "d : " << d << ", r: " << r << " dx: " << dx << endl;
		}
		x -= dx;
		if (DEBUG) {
			cout << "Dobivena vrijednost x : (" << x.real() << ", " << x.imag() << ")." << endl;
		}
		k++;
	}
	if (abs(dx) <= eps) {
		c = true;
		return x;
	}
	c = false;
	return x;
}

complex<double> Laguerre(vector<double> p, int n, complex<double> x,
	bool& c, double eps = 1e-10, double maxiter = 100) {
	vector<complex<double>> p2(p.size());
	for (int i = 0; i < p2.size(); i++)
		p2[i] = complex<double>(p[i]);
	return Laguerre(p2, n, x, c, eps, maxiter);
}

double fRand(double fMin, double fMax)
{
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}

complex<double> RandomComplex(double rmin, double rmax, double imin, double imax) {
	return complex<double>(fRand(rmin, rmax), fRand(imin, imax));
}

vector<complex<double>> PolyRoots(vector<complex<double>> p,
	double eps = 1e-10, int maxiters = 100, int tmax = 10) {
	if (eps <= 0 || maxiters <= 0 || tmax <= 0)
		throw domain_error("Invalid parameters");
	vector<complex<double>> roots(p.size()),
		porig(p);
	for (int i = p.size() - 1; i >= 1; i--) {
		int t(1);
		bool c(false);
		complex<double> x;
		while (!c && t < tmax) {
			x = Laguerre(p, i, RandomComplex(-10, 10, -10, 10),
				c, eps, maxiters);
			t++;
		}
		if (!c)
			throw logic_error("Convergence has not achieved");
		//Poliranje
		auto xpol(Laguerre(porig, porig.size() - 1, x, c, eps, maxiters));
		if (c)
			x = xpol;
		if (abs(x.imag()) <= eps)
			x = real(x);
		if (abs(x.real()) <= eps)
			x = { 0, x.imag() };
		roots[i] = x;
		auto v(p[i]);
		for (int j = i - 1; j >= 0; j--) {
			auto w(p[j]);
			p[j] = v;
			v = w + x * v;
		}
	}
	roots.erase(roots.begin());
	return roots;
}


vector<complex<double>> PolyRoots(vector<double> p,
	double eps = 1e-10, int maxiters = 100, int tmax = 10) {
	if (eps <= 0 || maxiters <= 0 || tmax <= 0)
		throw domain_error("Invalid parameters");
	vector<complex<double>> roots(p.size());
	vector<double> porig(p);
	int i(p.size() - 1);
	while (i >= 1) {
		int t(1);
		bool c(false);
		complex<double> x;
		while (!c && (t < tmax)) {
			x = Laguerre(p, i, RandomComplex(-10, 10, -10, 10),
				c, eps, maxiters);
			t++;
		}
		if (!c)
			throw logic_error("Convergence has not achieved");

		//Poliranje
		if (DEBUG)
			cout << "Poliranje " << endl;
		auto xpol(Laguerre(porig, porig.size() - 1, x, c, eps, maxiters));
		if (c)
			x = xpol;
		if (DEBUG)
			if (c)
				cout << "Poliranje uspjelo: " << xpol << endl;
		if (abs(x.imag()) <= eps) {
			x = x.real();
			if (DEBUG) {
				cout << "Odbacujem imag " << x << endl;
			}
			roots[i] = x;
			double v(p[i]);
			for (int j = i - 1; j >= 0; j--) {
				double w(p[j]);
				p[j] = v;
				v = w + x.real() * v;
			}
			if (DEBUG) {
				cout << "Podijeljeni polinom : ";
				for (int l = i - 1; l >= 0; l--)
					cout << p[l] << ", ";
				cout << endl;
			}
			i--;
		}
		else {
			if (DEBUG)
				cout << "Ne odbacujem imaginarni dio" << endl;
			if (abs(x.real()) <= eps) {
				x = { 0, x.imag() };
				if (DEBUG)
					cout << "Odbacujem imag " << x << endl;
			}
			roots[i] = x;
			roots[i - 1] = conj(x);
			if (DEBUG)
				cout << "Dobivene nule : " << roots[i] << ", " << roots[i - 1] << endl;
			auto alpha(2 * x.real());
			auto beta(abs(x) * abs(x));
			auto u(p[i]);
			auto v(p[i - 1] + alpha * u);
			for (int j = i - 2; j >= 0; j--) {
				auto w(p[j]);
				p[j] = u;
				u = v;
				v = w + alpha * v - beta * p[j];
			}
			if (DEBUG) {
				cout << "Podijeljeni polinom : ";
				for (int l = i - 2; l >= 0; l--)
					cout << p[l] << ", ";
				cout << endl;
			}
			i -= 2;
		}
	}
	roots.erase(roots.begin());
	return roots;
}

template <typename FunTip>
double BracketMinimum(FunTip f, double x0, double& xfinal,
	double eps = 1e-8, double hinit = 1e-5, double hmax = 1e10,
	double lambda = 1.4) {
	double h(hinit), x(x0);
	while (abs(h) < hmax) {
		if (f(x - h) < f(x))
			h *= -1;
		else if (!(f(x + h) < f(x))) {
			xfinal = x;
			return ((h < 0) ? -h : h);
		}
		x += h;
		h *= lambda;
	}
	throw logic_error("Minimum has not found");
}

template <typename FunTip>
double FindMinimum(FunTip f, double xinit, double eps = 1e-8,
	double hinit = 1e-5, double hmax = 1e10, double lambda = 1.4) {
	if (eps <= 0 || hinit <= 0 || hmax <= 0 || lambda <= 0)
		throw domain_error("Invalid parameters");
	double xfinal, hfinal(BracketMinimum(f, xinit, xfinal, eps, hinit, hmax, lambda));
	double a(xfinal - hfinal), c(xfinal), b(xfinal + hfinal), d;
	if (DEBUG)
		cout << "Ogradjen sa (" << a << ", " << c << ", " << b << ")" << endl;
	double R(0.61803399), C(1.0 - R);
	double x1, x2, x0(a), x3(b);
	if (abs(b - c) > abs(c - a)) {
		x1 = c;
		x2 = c + C * (b - c);
	}
	else {
		x2 = c;
		x1 = c - C * (c - a);
	}
	double f1(f(x1)), f2(f(x2));
	while (abs(x0 - x3) > eps) {
		if (f2 < f1) {
			x0 = x1;
			x1 = x2;
			x2 = R * x2 + C * x3;
			f1 = f2;
			f2 = f(x2);
		}
		else {
			x3 = x2;
			x2 = x1;
			x1 = R * x1 + C * x0;
			f2 = f1;
			f1 = f(x1);
		}
	}
	return (x0 + x3) / 2;
}

template <typename FunTip>
double RK4Step(FunTip f, double x, double y, double h) {
	double K1(f(x, y));
	double K2(f(x + h / 2, y + h * K1 / 2));
	double K3(f(x + h / 2, y + h * K2 / 2));
	double K4(f(x + h, y + h * K3));
	return y + h * (K1 + 2 * K2 + 2 * K3 + K4) / 6;
}

template <typename FunTip>
std::vector<std::pair<double, double>> RK4Integrator(FunTip f, double x0,
	double y0, double xmax, double h, double eps = 1e-8, bool adaptive = false) {
	if (eps <= 0 || h == 0)
		throw domain_error("Invalid parameters");
	if ((xmax < x0 && h > 0) || (xmax > x0 && h < 0))
		return{ { x0, y0 } };
	double x(x0), y(y0);
	vector<pair<double, double>> points;
	if (adaptive) {
		points.push_back({ x, y });
		while (h < 0 ? x >= xmax : x <= xmax) {
			double u = RK4Step(f, x, y, h / 2);
			double v = RK4Step(f, x + h / 2, u, h / 2);
			double w = RK4Step(f, x, y, h);
			double sigma = abs(w - v) / abs(h);
			if (sigma < eps)
			{
				x += h;
				y = v;
				points.push_back({ x, y });
			}
			//Greska u predavanju
			//Mada meni ni ova formula ne radi
			h *= min(5.0, 0.9 * pow((eps / sigma), 1. / 4));
		}
		if (points.back().first != xmax) {
			double x1(points[points.size() - 2].first),
				x2(points[points.size() - 1].first),
				y1(points[points.size() - 2].second),
				y2(points[points.size() - 1].second),
				ymax(y1 + (xmax - x1) * (y2 - y1) / (x2 - x1));
			points[points.size() - 1] = { xmax, ymax }; //Linearna interpolacija
		}
	}
	else {
		while (h < 0 ? (x >= xmax - eps) : x <= xmax + eps) {
			points.push_back({ x, y });
			double w(RK4Step(f, x, y, h));
			y = w;
			x += h;
			if (h < 0 ? x >= xmax - eps : x <= xmax + eps)
				points.push_back({ x, y });
		}
	}
	return points;
}

template <typename Tip>
string to_string(Tip x) {
	return static_cast<std::ostringstream&>((std::ostringstream() << std::dec << x)).str();
}

void print_vector_to_file(const vector<double> &v, ofstream &o) {
	o << " [ ";
	for (auto t : v)
		o << t << " ";
	o << "]";
}

void create_scilab(const vector<double>& xdata, const vector<double>& yinitial, const vector<double>& ysolved, string function) {
	static int counter(1);
	ofstream output("rk4_" + to_string(counter++) + ".sci");
	output << "t = ";
	print_vector_to_file(xdata, output);
	output << ";" << endl;
	output << "t0 = " << xdata[0] << "; " << endl;
	output << "y0 = ";
	print_vector_to_file(yinitial, output);
	output << ";" << endl;
	output << "ysolved = ";
	print_vector_to_file(ysolved, output);
	output << ";" << endl;
	output << "deff(\"yprim = f(x, y)\", \"yprim = " + function + "\");" << endl;
	output << "plot(t, ysolved, 'x');\nplot(t, ode(\"rk\", y0,t0,t,f), 'r');" << endl
		<< "legend([\"my integration\", \"scilab integration\"]);" << endl;
	output.close();
}

void test(double xinit, double xmax, double yinit, double h, function<double(double, double)> fxy,
	string fxystring, bool adaptive = false) {

	auto d(RK4Integrator(fxy, xinit, yinit, xmax, h, 1e-8, adaptive));

	vector<double> xdata, ysolved;
	for (auto point : d) {
		xdata.push_back(point.first);
		ysolved.push_back(point.second);
	}

	create_scilab(xdata, { yinit }, ysolved, fxystring);
}



int main() {

	//1. Brojevi iteracija potrebnih da se izvrši pronalazak nule se ne poklapaju u potpunosti sa 
	//	onim koji su navedeni u predavanju. Pretpostavljam da je dio razloga za ovo naèin na 
	//	koji se vrši poreðenje sa nulom. Ovdje su napisane funkcije is_zero() i lequal()
	//	koje porede brojeve sa nulom koristeæi proslijeðeni epsilon. Morao sam štimati epsilon
	//	u razlièitim funkcijama da bih postigao taèan broj iteracija. Recimo, za navedenu funkciju 
	//	iz predavanja, pod datim uslovima, taèno 54 iteracije sam uspio postiæi poredivši direktno 
	//	sa nulom, umjesto sa proslijeðenim epsilonom, dok sam bacanje izuzetka zbog 
	//	prvog izvoda u NewtonRaphson-u uspio dobiti tek poreðenjem sa epsilonom koji je reda 1e-15.
	//	Koristeæi i kod iz predavanja i kod iz "Numerical Recipes", nisam uspio dostiæi broj iteracija
	//	manjih od 9 za istu funkciju. 
	//2. U kodu za raèunanje vrijednosti polinoma, kao i njegovih izvoda u Laguerre funkciji, 
	//	mislim da donja granica petlje treba da bude 0, ne 2. Tako je navedeno u kodu u "Numerical Recipes", a
	//	meni ne radi bez te izmjene.
	//3. Druga potencijalna greška je u kodu za RK4Integrator funkciju, taènije za njenu adaptivnu verziju. 
	//	Ako sam dobro shvatio teoretsku analizu prije pseuodokoda u predavanju, h bi se trebalo raèunati po formuli
	//	h = h * min{5, 0.9 * pow((eps / sigma), 1. / 4) }, a ne h = h * min{ 5, 0.9 * h * pow((eps / sigma), 1. / 4) }.
	//	Meðutim, i u ispravljenoj verziji sam imao problema da se h konstantno smanjuje bez obzira da li je uvjet 
	//	taènosti zadovoljen ili ne.
	//4. Poliranje nula je implementovano.
	//5. Testiranje RK4 integratora se vrši generisanje scilab skripti u kojima se dobiveno rješenje poredi sa
	//	rezultatom dobivenim scilabovom ode() funkcijom sa "rk" parametrom. Da bi se omoguæilo generisanje
	//	treba globalnu varijablu GENERATE_SCILAB postaviti na true
        //6. Kod za FindMinimum je iz "Numerical Recipes". Zasniva se na istom algoritmu ali su malo preuredjene formule 
        //      za racunanje c i d i iz nekog razloga mi je davao bolje rezultate.

	const double e(exp(1)), pi(atan(1) * 4);
	double a, b;
	vector<function<double(double)>> f{
		[](double x) { return x * x;  },
		[](double x) { return sin(x);  },
		[](double x) { return cos(x); },
		[e](double x) { return 0.05 * (pow(e, 10 * (x - 3)) - 1); },
		[pi](double x) { return 3 * x * x + 1 / pow(pi, 4) * pow(log(pi - x), 2); },
		[](double x) { return x * x * x - 2 * x * x - 11 * x + 12;  },
		[](double x) { return 3 * x * x - 4 * x - 11;  },
		[](double x) { return -x * x;  },
		[](double x) { return 0.00000001 * x * x;  }, 
		[](double x) { return (x - 7) * (x - 11) * (x - 69);  },
		[](double x) { return 3 * x * x - 174 * x + 1319;  }
	};

	if (BracketRoot(f[0], -1, a, b))
		cout << a << " " << b << endl;
	else
		cout << "Nula nije pronadjena" << endl; // Nula nije pronadjena

	if (BracketRoot(f[4], pi + 0.1, a, b))
		cout << a << " " << b << endl;
	else
		cout << "Nula nije pronadjena" << endl; // Nula nije pronadjena

	if (BracketRoot(f[1], -1, a, b))
		cout << a << " " << b << endl; // a < 0 < b
	else
		cout << "Nula nije pronadjena" << endl;

	try {
		ModifiedRegulaFalsi(f[0], -1, 1);
	}
	catch (exception& e) {
		cout << e.what() << endl; // Root must be bracketed
	}

	try {
		ModifiedRegulaFalsi(f[1], 2, 1);
	}
	catch (exception& e) {
		cout << e.what() << endl; // Bad range
	}

	try {
		cout << ModifiedRegulaFalsi(f[1], a, b) << endl; // Oko 0
	}
	catch (exception& e) {
		cout << e.what() << endl;
	}

	//Test problematiène funkcije za regulu falsi
	cout << ModifiedRegulaFalsi(f[3], 1, 4) << endl;

	//Za datu funkciju transformacije je potrebno 54 iteracije
	try {
		cout << ModifiedRegulaFalsi(f[3], 1, 4, 1e-10, 54) << endl; //3
	}
	catch (exception& e) {
		cout << e.what() << endl;
	}

	try {
		cout << ModifiedRegulaFalsi(f[3], 1, 4, 1e-10, 53) << endl;
	}
	catch (exception& e) {
		cout << e.what() << endl; //Given accuracy not achieved
	}

	try {
		Ridders(f[0], -1, 1);
	}
	catch (exception& e) {
		cout << e.what() << endl; // Root must be bracketed
	}

	try {
		Ridders(f[1], 2, 1);
	}
	catch (exception& e) {
		cout << e.what() << endl; // Bad range
	}

	try {
		cout << Ridders(f[1], a, b) << endl; // Oko 0
	}
	catch (exception& e) {
		cout << e.what() << endl;
	}

	//Test problematiène funkcije za regulu falsi
	cout << Ridders(f[3], 1, 4) << endl; // Oko 3

	cout << "Poreðenje metoda Ridders2(Numerical Recipes) i Ridders(Predavanja) za navedene uvjete: " << endl;
	try {
		//Test problematiène funkcije za regulu falsi
		Ridders2(f[3], 1, 4, 1e-10, 9);
		cout << "Ridders(Recipes) stize za 9" << endl;
	}
	catch (...) {
		cout << "Ridders(Recipes) ne stize za 9" << endl;
	}

	try {
		//Test problematiène funkcije za regulu falsi
		Ridders2(f[3], 1, 4, 1e-10, 10);
		cout << "Ridders(Recipes) stize za 10" << endl;
	}
	catch (...) {
		cout << "Ridders(Recipes) ne stize za 10" << endl;
	}

	try {
		//Test problematiène funkcije za regulu falsi
		Ridders(f[3], 1, 4, 1e-10, 8);
		cout << "Ridders(Predavanja) stize za 8" << endl;
	}
	catch (...) {
		cout << "Ridders(Predavanja) ne stize za 8" << endl;
	}

	try {
		//Test problematiène funkcije za regulu falsi
		Ridders(f[3], 1, 4, 1e-10, 9);
		cout << "Ridders(Predavanja) stize za 9" << endl;
	}
	catch (...) {
		cout << "Ridders(Predavanja) ne stize za 9" << endl;
	}


	cout << NewtonRaphson(f[5], f[6], 2.35283) << endl; // Oko 1

	cout << NewtonRaphson(f[5], f[6], 2.35284) << endl; // Oko -3

	cout << NewtonRaphson(f[5], f[6], 2.35285) << endl; // Oko 4

	cout << NewtonRaphson(f[9], f[10], 8) << endl; // Oko 7

	cout << NewtonRaphson(f[9], f[10], 50) << endl; // Oko 69

	try {
		cout << NewtonRaphson(f[5], f[6], (2 - sqrt(37)) / 3) << endl;
	}
	catch (exception& e) {
		cout << e.what() << endl; //Convergence has not achieved
	}

	vector<complex<double>> x;

	x = PolyRoots(vector<double>{ 1, 1, 1 });

	for (auto xi : x)
		cout << xi << ", "; //1/2 +- sqrt(3)/2
	cout << endl;

	x = PolyRoots(vector<double>{ -6, 11, -6, 1 });


	for (auto xi : x)
		cout << xi << ", "; //1, 2 i 3
	cout << endl;
	x = PolyRoots(vector<double>{ -21, -25, -45, -28, -23, -3, 1 });

	for (auto xi : x)
		cout << xi << ", "; //+-i, 1/2 +- sqrt(3)/2, -3, 7
	cout << endl;

	vector<complex<double>> x2;
	x2 = PolyRoots(vector<complex<double>>{ -18, { 3, -12 }, { -3, 2 }, { 1, 2 }, 1});

	for (auto xi : x2)
		cout << xi << ", "; // -3, 2, i, -3i
	cout << endl;
	try {
		x2 = PolyRoots(vector<complex<double>>{ -18, { 3, -12 }, { -3, 2 }, { 1, 2 }, 1}, 1e-10, 3);
	}
	catch (exception& e) {
		cout << e.what() << endl; // Convergence has not achieved
	}


	//x^2
	cout << FindMinimum(f[0], 10) << endl; // 0

	cout << FindMinimum(f[0], -10) << endl; // 0

											//-x^2
	try {
		cout << FindMinimum(f[7], 5) << endl;
	}
	catch (exception& e) {
		cout << e.what() << endl; //Minimum has not found
	}

	//-x^2
	try {
		cout << FindMinimum(f[8], 100) << endl;
	}
	catch (exception& e) {
		cout << e.what() << endl; //Minimum has not found
	}

	if (GENERATE_SCILAB) {

		//rk4_1.sci
		test(0, 10, 0, 0.1, [](double x, double y) { return 1;  }, "1"); //y = x;

		//rk4_2.sci
		test(0, 10, 0, 0.1, [](double x, double y) { return 3 * x + 2 * y;  }, "3 * x + 2 * y"); //e ^ 2x - 3 * x / 2 - 3 / 4

		//rk4_3.sci
		test(0, -10, 0, -0.1, [](double x, double y) { return 3 * x + 2 * y;  }, "3 * x + 2 * y"); //e ^ 2x - 3 * x / 2 - 3 / 4 (unazad)

		/*
		//Ne radi
		//rk4_4.sci
		test(0, 10, 0, 0.1, [](double x, double y) { return 3 * x + 2 * y;  }, "3 * x + 2 * y", true); //e ^ 2x - 3 * x / 2 - 3 / 4 (adaptive)
		*/
	}
	return 0;
}