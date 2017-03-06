#include <iostream>
#include <vector>
#include <stdexcept>
#include <iomanip>
#include <algorithm>

using namespace std;

class Vector;
class Matrix;

double Epsilon(Matrix m);
double Epsilon(Vector v);

class Vector : public vector<double> {
    vector<double> niz;
    static void JednakeDuzine(const Vector& v1, const Vector& v2) {
        if(v1.NElems() != v2.NElems())
            throw domain_error("Incompatible formats");
    }
public:
    explicit Vector(int n) {
        if(n < 0)
            throw range_error("Bad dimension");
        niz = vector<double>(n);
    };
    Vector(std::initializer_list<double> l) : niz(l) {
        if(l.size() == 0)
            throw range_error("Bad dimension");
        niz = vector<double>(l);
    }
    int NElems() const { return niz.size(); }
    double &operator[](int i) { return niz[i]; }
    double operator[](int i) const { return niz[i]; }
    double &operator()(int i) {
        if(i < 1 || i > niz.size())
            throw range_error("Invalid index");
        return niz[i - 1];
    }
    double operator()(int i) const {
        if(i < 1 || i > niz.size())
            throw range_error("Invalid index");
        return niz[i - 1];
    }
    void Print(char separator = '\n') const {
        for(int i = 0; i < niz.size(); i++) {
            cout << niz[i];
            if((i < niz.size() - 1) || separator == '\n')
                cout << separator;
        }
    }
    friend Vector operator +(const Vector &v1, const Vector &v2);
    Vector &operator +=(const Vector &v) {
        JednakeDuzine(*this, v);
        for(int i = 0; i < niz.size(); i++)
            niz[i] += v.niz[i];
        return *this;
    }
    friend Vector operator -(const Vector &v1, const Vector &v2);
    Vector &operator -=(const Vector &v) {
        JednakeDuzine(*this, v);
        for(int i = 0; i < niz.size(); i++)
            niz[i] -= v.niz[i];
        return *this;
    }
    friend Vector operator *(double s, const Vector &v);
    friend Vector operator *(const Vector &v, double s);
    Vector &operator *=(double s) {
        for(int i = 0; i < niz.size(); i++)
            niz[i] *= s;
        return *this;
    }
    friend double operator *(const Vector &v1, const Vector &v2);
    friend Vector operator /(const Vector &v, double s);
    Vector &operator /=(double s) {
        if(s == 0) throw domain_error("Division by zero");
        for(int i = 0; i < niz.size(); i++)
            niz[i] /= s;
        return *this;
    }
};

Vector operator +(const Vector &v1, const Vector &v2) {
    auto v3(v1);
    return v3 += v2;
}

Vector operator -(const Vector &v1, const Vector &v2) {
    auto v3(v1);
    return v3 -= v2;
}

Vector operator *(double s, const Vector &v) {
    auto  v2(v);
    return v2 *= s;
}
Vector operator *(const Vector &v, double s) {
    return s * v;
}

double operator *(const Vector &v1, const Vector &v2) {
    Vector::JednakeDuzine(v1, v2);
    int sum(0);
    for(int i = 0; i < v1.niz.size(); i++)
        sum += v1.niz[i] * v2.niz[i]; // Kahanovo mnozenje?
    return sum;
}

Vector operator /(const Vector &v, double s) {
    auto v2(v);
    return v2 /= s;
}

class Matrix {
    vector<vector<double>> m;
    void TestKvadratneMatrice() const {
        if(NRows() != NCols())
            throw domain_error("Divisor matrix is not square");
    }
    int RREFAndRank();
public:
    Matrix(int m, int n);
    Matrix(const Vector &v) : m(v.NElems(), vector<double>(1)) {
        for(int i = 0; i < m.size(); i++)
            m[i][0] = v[i];
    }
    Matrix(initializer_list<std::vector<double>> l);
    int NRows() const { return m.size(); }
    int NCols() const { return m[0].size(); }
    double *operator[](int i) { return &m[i][0]; }
    const double *operator[](int i) const { return &m[i][0]; };
    double &operator()(int i, int j) {
        if(i < 1 || j < 1 || i > m.size() || j > m[0].size())
            throw range_error("Invalid index");
        return m[i - 1][j - 1];
    }
    double operator()(int i, int j) const {
        if(i < 1 || j < 1 || i > m.size() || j > m[0].size())
            throw range_error("Invalid index");
        return m[i - 1][j - 1];
    }
    void Print(int width = 10) const {
        for(int i = 0; i < NRows(); i++) {
            for(int j = 0; j < NCols(); j++)
                cout << setw(width) << m[i][j];
            cout << ((i == NRows() - 1) ? "" : "\n");
        }
    }
    friend Matrix operator +(const Matrix &m1, const Matrix &m2);
    Matrix &operator +=(const Matrix &m);
    friend Matrix operator -(const Matrix &m1, const Matrix &m2);
    Matrix &operator -=(const Matrix &m);
    friend Matrix operator *(double s, const Matrix &m);
    friend Matrix operator *(const Matrix &m, double s);
    Matrix &operator *=(double s);
    friend Matrix operator *(const Matrix &m1, const Matrix &m2);
    Matrix &operator *=(const Matrix &m);
    friend Vector operator *(const Matrix &m, const Vector &v);
    friend Matrix Transpose(const Matrix &m);
    void Transpose();

    friend Matrix LeftDiv(Matrix m1, Matrix m2);
    friend Vector LeftDiv(Matrix m, Vector v);
    friend Matrix operator /(const Matrix &m, double s);
    Matrix& operator /=(double s);
    friend Matrix operator /(Matrix m1, Matrix m2);
    Matrix& operator /=(Matrix m);
    double Det() const;
    friend double Det(Matrix m);
    void Invert();
    friend Matrix Inverse(Matrix m);
    void ReduceToRREF();
    friend Matrix RREF(Matrix m);
    int Rank() const;
    friend int Rank(Matrix m);

};

Matrix LeftDiv(Matrix m1, Matrix m2) {
    m1.TestKvadratneMatrice();
    if(m1.NRows() != m2.NRows())
        throw domain_error("Incompatible formats");
    int n(m2.NRows()), m(m2.NCols());
    double epsilon(Epsilon(m1));
    for(int k = 0; k < n; k++) {
        int p = k;
        for(int i = k + 1; i < n; i++)
            if(fabs(m1[i][k]) > fabs(m1[p][k]))
                p = i;
        if(fabs(m1[p][k]) < epsilon)
            throw domain_error("Divisor matrix is singular");
        if(p != k) {
            swap(m1.m[k], m1.m[p]);
            swap(m2.m[k], m2.m[p]);
        }

        for(int i = k + 1; i < n; i++) {
            auto u(m1[i][k] / m1[k][k]);
            for(int j = k + 1; j < n; j++)
                m1[i][j] -= u * m1[k][j];
            for(int j = 0; j < m; j++)
                m2[i][j] -= u * m2[k][j];
        }
    }
    Matrix x(n, m);
    for(int k = 0; k < m; k++) {
        for(int i = n - 1; i >= 0; i--) {
            auto s(m2[i][k]);
            for(int j = i + 1; j < n; j++)
                s -= m1[i][j] * x[j][k];
            x[i][k] = s / m1[i][i];
        }
    }
    return x;
}

Vector LeftDiv(Matrix m, Vector v) {
    auto temp(LeftDiv(m, Matrix(v)));
    Vector ret(temp.NRows());
    for(int i = 0; i < ret.NElems(); i++)
        ret[i] = temp[i][0];
    return ret;
}

Matrix operator / (const Matrix &m, double s) {
    Matrix ret(m);
    ret /= s;
    return ret;
}

Matrix& Matrix::operator /=(double s) {
    if(s == 0)
        throw domain_error("Division by zero");
    for(int i = 0; i < NCols(); i++)
        for(int j = 0; j < NRows(); j++)
            (*this)[i][j] /= s;
    return *this;
}

Matrix operator /(Matrix m1, Matrix m2) {
    return m1 /= m2;
}
Matrix& Matrix::operator /=(Matrix m1) {
    Matrix& m2(*this);
    m1.TestKvadratneMatrice();
    if(m2.NCols() != m1.NRows())
        throw domain_error("Incompatible formats");
    int n(m2.NRows()), m(m2.NCols());
    double epsilon(Epsilon(m2));
    for(int k = 0; k < n; k++) {
        int p = k;
        for(int i = k + 1; i < n; i++)
            if(fabs(m1[k][i]) > fabs(m1[k][p]))
                p = i;
        if(fabs(m1[k][p]) < epsilon)
            throw domain_error("Divisor matrix is singular");
        if(p != k) {
            for(int l = 0; l < m1.NRows(); l++)
                swap(m1[l][k], m1[l][p]);
            for(int l = 0; l < m2.NRows(); l++)
                swap(m2[l][k], m2[l][p]);
        }
        for(int i = k + 1; i < n; i++) {
            auto u(m1[k][i] / m1[k][k]);
            for(int j = k + 1; j < n; j++)
                m1[j][i] -= u * m1[j][k];
            for(int j = 0; j < m; j++)
                m2[j][i] -= u * m2[j][k];
        }
    }
    for(int k = 0; k < m; k++) {
        for(int i = n - 1; i >= 0; i--) {
            auto s(m2[k][i]);
            for(int j = i + 1; j < n; j++)
                s -= m1[j][i] * m2[k][j];
            m2[k][i] = s / m1[i][i];
        }
    }
    return *this;
}

double Matrix::Det() const {
    if(NCols() != NRows())
        throw domain_error("Matrix is not square");
    Matrix m(*this);
    int n(NCols());
    double d(1), epsilon(Epsilon(m));
    for(int k = 0; k < n; k++) {
        int p(k);
        for(int i = k + 1; i < n; i++)
            if(fabs(m[i][k]) > m[p][k])
                p = i;
        if(fabs(m[p][k]) < epsilon)
            return 0;
        if(p != k)
            swap(m.m[p], m.m[k]),
            d *= -1;
        for(int i = k + 1; i < n; i++) {
            auto u(m[i][k] / m[k][k]);
            for(int j = k + 1; j < n; j++)
                m[i][j] -= u * m[k][j];
        }
    }
    for(int i = 0; i < n; i++)
        d *= m[i][i];
    return d;
}

double Det(Matrix m) {
    return m.Det();
}

void Matrix::Invert() {
    if(Det() == 0)
        throw domain_error("Matrix is singular"); //Problematicno
    if(NCols() != NRows())
        throw domain_error("Matrix is not square");
    int n(NCols());
    Matrix& m(*this);
    for(int k = 0; k < n; k++) {
        auto u(m[k][k]);
        m[k][k] = 1;
        for(int j = 0; j < n; j++)
            m[k][j] /= u;
        for(int i = 0; i < n; i++)
            if(i != k) {
                u = m[i][k];
                m[i][k] = 0;
                for(int j = 0; j < n; j++)
                    m[i][j] -= u * m[k][j];
            }
    }
}

Matrix Inverse(Matrix m) {
    m.Invert();
    return m;
}

int Matrix::RREFAndRank() {
    Matrix& a(*this);
    int k(-1), l(-1), n(NRows()), m(NCols()), p;
    double epsilon(Epsilon(a));
    while(k < m && l < n) {
        l++;
        k++;
        int v(0);
        while(v < epsilon && l < n) {
            p = k;
            for(int i = k; i < m; i++)
                if(fabs(a[i][l]) > v) {
                    v = fabs(a[i][l]);
                    p = i;
                }
            if(v < epsilon)
                l++;
        }
        if(l < n) {
            if(p != k)
                swap(a.m[k], a.m[p]);
            auto u(a[k][l]);
            for(int j = l; j < n; j++)
                a[k][j] /= u;
            for(int i = 0; i < m; i++)
                if(i != k) {
                    u = a[i][l];
                    for(int j = l; j < n; j++) {
                        a[i][j] -= u * a[k][j];
                    //    if(a[i][j] < std::numeric_limits<double>::epsilon())
                     //       a[i][j] = 0;
                    }
                }
        }
    }
    return k;
}

void Matrix::ReduceToRREF() {
    RREFAndRank();
}

Matrix RREF(Matrix m) {
    m.ReduceToRREF();
    return m;
}

int Matrix::Rank() const {
    Matrix temp(*this);
    return temp.RREFAndRank();
}

int Rank(Matrix m) {
    return m.Rank();
}

Matrix::Matrix(int m, int n) {
    if(m <= 0 || n <= 0)
        throw range_error("Bad dimension");
    this->m = vector<vector<double>>(m, vector<double>(n));
}


Matrix::Matrix(initializer_list<std::vector<double>> l) {
    if(l.size() == 0)
        throw range_error("Bad dimension");
    for(auto it = l.begin(); it != l.end(); it++)
        if(it->size() != l.begin()->size())
            throw logic_error("Bad matrix");
    for(auto it = l.begin(); it != l.end(); it++)
        if(it->size() == 0)
            throw range_error("Bad dimension");
    m = vector<vector<double>>(l.size());
    auto it = l.begin();
    for(int i = 0; i < m.size(); i++, it++)
        m[i] = vector<double>(*it);
}

Matrix operator +(const Matrix &m1, const Matrix &m2) {
    auto ret(m1);
    return ret += m2;
}

Matrix& Matrix::operator +=(const Matrix &m) {
    if(NRows() != m.NRows() || NCols() != m.NCols())
        throw domain_error("Incompatible formats");
    for(int i = 0; i < NRows(); i++)
        for(int j = 0; j < NCols(); j++)
            this->m[i][j] += m[i][j];
    return *this;
}

Matrix operator -(const Matrix &m1, const Matrix &m2) {
    auto ret(m1);
    return ret -= m2;
}

Matrix& Matrix::operator -=(const Matrix &m) {
    if(NRows() != m.NRows() || NCols() != m.NCols())
        throw domain_error("Incompatible formats");
    for(int i = 0; i < NRows(); i++)
        for(int j = 0; j < NCols(); j++)
            this->m[i][j] -= m[i][j];
    return *this;
}

Matrix operator *(double s, const Matrix &m) {
    return m * s;
}

Matrix operator *(const Matrix &m, double s) {
    auto ret(m);
    return ret *= s;
}

Matrix& Matrix::operator *=(double s) {
    for(auto& vd : m)
        for(auto& d : vd)
            d *= s;
    return *this;
}

Matrix operator *(const Matrix &m1, const Matrix &m2) {
    auto ret(m1);
    return ret *= m2;
}

Matrix& Matrix::operator *=(const Matrix &m) {
    if(NCols() != m.NRows())
        throw domain_error("Incompatible formats");
    Matrix ret(NRows(), m.NCols());
    for(int i = 0; i < NRows(); i++)
        for(int j = 0; j < m.NCols(); j++)
            for(int k = 0; k < NCols(); k++)
                ret[i][j] += this->m[i][k] * m[k][j];
    this->m = ret.m;
    return *this;
}

Vector operator *(const Matrix &m, const Vector &v) {
    auto temp(m);
    temp *= v;
    Vector ret(temp.NRows());
    for(int i = 0; i < temp.NRows(); i++)
        ret[i] = temp[i][0];
    return ret;
}

Matrix Transpose(const Matrix &m) {
    auto ret(m);
    ret.Transpose();
    return ret;
}

void Matrix::Transpose() {
    if(NCols() == NRows())
        for(int i = 0; i < NRows() - 1; i++)
            for(int j = i + 1; j < NCols(); j++)
                swap(m[i][j], m[j][i]);
    else {
        vector<vector<double>> temp(NCols(), vector<double>(NRows()));
        for(int i = 0; i < NRows(); i++)
            for(int j = 0; j < NCols(); j++)
                temp[j][i] = m[i][j];
        m = temp;
    }
}

class LUDecomposer {
    Matrix a;
    Vector w;
    double l(int i, int j) const {
        if(j > i) return 0;
        if(j == i) return 1;
        return a[i][j];
    }
    double u(int i, int j) const {
        if(j < i) return 0;
        return a[i][j];
    }
public:
    LUDecomposer(Matrix m);
    void Solve(const Vector &b, Vector &x) const;
    Vector Solve(Vector b) const;
    void Solve(Matrix &b, Matrix &x) const;
    Matrix Solve(Matrix b) const;
    Matrix GetCompactLU() const;
    Matrix GetL() const;
    Matrix GetU() const;
    Vector GetPermuation() const;
};

LUDecomposer::LUDecomposer(Matrix m) : a(m.NRows(), m.NCols()), w(m.NRows()) {
    if(m.NCols() != m.NRows())
        throw domain_error("Matrix is not square");
    int n(m.NCols()), p;
    double s, epsilon(Epsilon(m));
    a = m;
    for(int j = 0; j < n; j++) {
        for(int i = 0; i <= j; i++) {
           s = a[i][j];
           for(int k = 0; k <= i - 1; k++)
                s -= a[i][k] * a [k][j];
           a[i][j] = s;
        }
        p = j;
        for(int i = j + 1; i < n; i++) {
            s = a[i][j];
            for(int k = 0; k <= j - 1; k++)
                s -= a[i][k] * a[k][j];
            a[i][j] = s;
            if(fabs(s) > fabs(a[p][j]))
                p = i;
        }
        if(fabs(a[p][j]) < epsilon)
            throw domain_error("Matrix is singular");
        if(p != j)
            for(int i = 0; i < n; i++)
                swap(a[j][i], a[p][i]);
        w[j] = p;
        auto u(1 / a[j][j]);
        for(int i = j + 1; i < n; i++)
            a[i][j] = u * a[i][j];
    }
 }

void LUDecomposer::Solve(const Vector &b0, Vector &x) const{
    Vector b(b0);
    if(b.NElems() != a.NRows() || b.NElems() != x.NElems())
        throw domain_error("Incompatible formats");
    int n(b.NElems());
    for(int i = 0; i < n; i++) {
        int p(w[i]);
        double s(b[p]);
        b[p] = b[i];
        for(int j = 0; j <= i - 1; j++)
            s -= l(i, j) * x[j];
        x[i] = s;
    }
    for(int i = n - 1; i >= 0; i--) {
        double s(x[i]);
        for(int j = i + 1; j < n; j++)
            s -= u(i, j) * x[j];
        x[i] = s / u(i, i);
    }
}

Vector LUDecomposer::Solve(Vector b) const {
    Vector temp(b.NElems());
    Solve(b, temp);
    return temp;
}

void LUDecomposer::Solve(Matrix &b, Matrix &x) const {
    if(b.NCols() != x.NCols())
        throw domain_error("Incompatible formats");
    for(int i = 0; i < b.NCols(); i++) {
        Vector vb(b.NRows()), vx(b.NRows());
        for(int j = 0; j < b.NRows(); j++)
            vb[j] = b[j][i];
        Solve(vb, vx);
        for(int j = 0; j < b.NRows(); j++)
            x[j][i] = vx[j];
    }
}

Matrix LUDecomposer::Solve(Matrix b) const {
    Matrix temp(b.NRows(), b.NCols());
    Solve(b, temp);
    return temp;
}

Matrix LUDecomposer::GetCompactLU() const {
    return a;
}

Matrix LUDecomposer::GetL() const {
    Matrix ret(a.NRows(), a.NCols());
    for(int i = 0; i < a.NRows(); i++)
        for(int j = 0; j < a.NCols(); j++)
            ret[i][j] = l(i, j);
    return ret;
}

Matrix LUDecomposer::GetU() const {
    Matrix ret(a.NRows(), a.NCols());
    for(int i = 0; i < a.NRows(); i++)
        for(int j = 0; j < a.NCols(); j++)
            ret[i][j] = u(i, j);
    return ret;
}

Vector LUDecomposer::GetPermuation() const {
    return w;
}


class QRDecomposer {
    Matrix a;
    Vector d;
    double r(int i, int j) const {
        if(j > i)
            return a[i][j];
        if(j == i)
            return d[i];
        return 0;
    }
    Matrix MulQWithMatrix(Matrix b, bool transposed = false) const;
    Vector MulQWithVector(Vector v, bool transposed = false) const;
public:
    QRDecomposer(Matrix m);
    void Solve(const Vector &b, Vector &x) const;
    Vector Solve(Vector b) const;
    void Solve(Matrix &b, Matrix &x) const;
    Matrix Solve(Matrix b) const;
    Vector MulQWith(Vector v) const;
    Matrix MulQWith(Matrix m) const;
    Vector MulQTWith(Vector v) const;
    Matrix MulQTWith(Matrix m) const;
    Matrix GetQ() const;
    Matrix GetR() const;
};

//Singularnost
QRDecomposer::QRDecomposer(Matrix a) : a(a.NRows(), a.NCols()), d(a.NCols()) {
    if(a.NRows() < a.NCols())
        throw domain_error("Invalid matrix format");
    int m(a.NRows()), n(a.NCols());
    double epsilon(Epsilon(a));
    for(int k = 0; k < n; k++) {
        double s(0);
        for(int i = k; i < m; i++)
            s += a[i][k] * a[i][k];
        s = sqrt(s);
        double u(sqrt(s * (s + fabs(a[k][k]))));
        if(fabs(u) < epsilon)
            throw domain_error("Matrix is singular");
        if(a[k][k] < 0)
            s *= -1;
        a[k][k] = (a[k][k] + s) / u;
        for(int i = k + 1; i < m; i++)
            a[i][k] /= u;
        d[k] = -s;
        for(int j = k + 1; j < n; j++) {
            s = 0;
            for(int i = k; i < m; i++)
                s += a[i][k] * a[i][j];
            for(int i = k; i < m; i++)
                a[i][j] -= s * a[i][k];
        }
    }
    this->a = a;
}

void QRDecomposer::Solve(const Vector &b, Vector &x) const {
    if(a.NCols() != a.NRows())
        throw domain_error("Matrix is not square");
    if(b.NElems() != a.NRows() || b.NElems() != x.NElems())
       throw domain_error("Incompatible formats");
    int n(a.NCols());
    Vector b2(MulQTWith(b));
    for(int i = n - 1; i >= 0; i--) {
        double s(b2[i]);
        for(int j = i + 1; j < n; j++)
            s -= r(i, j) * x[j];
        x[i] = s / r(i, i);
   }
}

Vector QRDecomposer::Solve(Vector b) const {
    Vector temp(b.NElems());
    Solve(b, temp);
    return temp;
}

void QRDecomposer::Solve(Matrix &b, Matrix &x) const {
    if(b.NCols() != x.NCols())
        throw domain_error("Incompatible formats");
    for(int i = 0; i < b.NCols(); i++) {
        Vector vb(b.NRows()), vx(b.NRows());
        for(int j = 0; j < b.NRows(); j++)
            vb[j] = b[j][i];
        Solve(vb, vx);
        for(int j = 0; j < b.NRows(); j++)
            x[j][i] = vx[j];
    }
}

Matrix QRDecomposer::Solve(Matrix b) const {
    Matrix temp(b.NRows(), b.NCols());
    Solve(b, temp);
    return temp;
}

Matrix QRDecomposer::MulQWithMatrix(Matrix b, bool transposed) const {
    for(int i = 0; i < b.NCols(); i++) {
        Vector temp(a.NRows());
        for(int j = 0; j < temp.NElems(); j++)
            temp[j] = b[j][i];
        temp = MulQWithVector(temp, transposed);
        for(int j = 0; j < temp.NElems(); j++)
            b[j][i] = temp[j];
    }
    return b;
}

Vector QRDecomposer::MulQWithVector(Vector v, bool transposed) const {
    int m(a.NRows()), n(a.NCols());
    if(m != v.NElems())
        throw domain_error("Incompatible formats");
    int k;
    for(transposed ? k = 0 : k = n - 1; transposed ? k < n : k >= 0; transposed ? k++ : k--) {
        double s(0);
        for(int i = k; i < m; i++)
            s += a[i][k] * v[i];
        for(int i = k; i < m; i++)
            v[i] -= s * a[i][k];
    }
    double epsilon(Epsilon(v));
    for(int i = 0; i < m; i++)
        if(fabs(v[i]) < epsilon)
            v[i] = 0;
    return v;
}

Matrix QRDecomposer::MulQWith(Matrix b) const {
   return MulQWithMatrix(b);
}

Vector QRDecomposer::MulQWith(Vector v) const {
    return MulQWithVector(v);
}

Matrix QRDecomposer::MulQTWith(Matrix m) const {
    return MulQWithMatrix(m, true);
}

Vector QRDecomposer::MulQTWith(Vector v) const {
    return MulQWithVector(v, true);
}

Matrix QRDecomposer::GetQ() const {
    int m(a.NRows()), n(a.NCols());
    Matrix q(m, m);
    for(int j = 0; j < m; j++) {
        for(int i = 0; i < m; i++)
            q[i][j] = 0;
        q[j][j] = 1;
        for(int k = n - 1; k >= 0; k--) {
            double s(0);
            for(int i = k; i < m; i++)
                s += a[i][k] * q[i][j];
            for(int i = k; i < m; i++)
                q[i][j] -= s * a[i][k];
        }
    }
    return q;
}
Matrix QRDecomposer::GetR() const {
    Matrix ret(a.NRows(), a.NCols());
    for(int i = 0; i < ret.NRows(); i++)
        for(int j = 0; j < ret.NCols(); j++)
            ret[i][j] = r(i, j);
    return ret;
}

double Epsilon(Matrix m) {
    double sum(0);
    for(int i = 0; i < m.NRows(); i++)
        for(int j = 0; j < m.NCols(); j++)
            sum += fabs(m[i][j]);
    return sum * numeric_limits<double>::epsilon();
}

double Epsilon(Vector v) {
    double sum(0);
    for(int j = 0; j < v.NElems(); j++)
        sum += fabs(v[j]);
    return sum * numeric_limits<double>::epsilon();
}

int main()
{
    Matrix a{{2, 1, 3}, {2, 6, 8}, {6, 8, 18}}, b{{1}, {3}, {5}};
    cout << endl << "Dijeljenje nekvadratnom matricom: " << endl;
    try {
        LeftDiv(Matrix{{1, 2, 3}, {3, 4, 5}}, Matrix{{1, 2, 3}, {3, 4, 5}}).Print();
    } catch(exception& e) {
        cout << e.what() << endl;
    }
    cout << endl << "Neispravni formati: " << endl;
    try {
        LeftDiv(Matrix{{1, 2, 3}, {3, 4, 5}, {5, 6, 7}}, Matrix{{1, 2, 3}, {3, 4, 5}}).Print();
    } catch(exception& e) {
        cout << e.what() << endl;
    }
    cout << endl << "A / B: " << endl;
    LeftDiv(a, b).Print();
    Vector bv{1, 3, 5};
    cout << endl << "A / B(vektorski oblik): " << endl;
    LeftDiv(a, bv).Print();
    Matrix c{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, d{{7, 5, 6}, {2, 0, 8}, {5, 7, 1}};
    cout << endl << "Dijeljenje nulom: " << endl;
    try {
        c / 0;
    } catch(exception &e) {
        cout << e.what() << endl;
    }
    cout << endl << "C / 5: " << endl;
    (c / 5).Print();
    cout << endl << "C / D: " << endl;
    (c / d).Print();
    cout << endl << "Determinanta nekvadratne: " << endl;
    try {
        Matrix{{1, 2, 3}, {4, 5, 6}}.Det();
    } catch(exception &e) {
        cout << e.what() << endl;
    }

    cout << endl << "Det(C): " << c.Det() << ", Det(D): " << Det(d) << endl;
    cout << "Inverzija nekvadratne matrice: " << endl;
    try {
        c.Invert();
    } catch(exception &e) {
        cout << e.what() << endl;
    }
    cout << endl << "Inverse(D): " << endl;
    Inverse(d).Print();
    Matrix e{{3, 4, 18, 34, 0, 2, 31},
             {1, -3, -7, -6, 2, 4, 26},
             {2, 1, 7, 16, 3, -1, 27},
             {5, 11, 43, 74, 2, 0, 56},
             {3, -3, -3, 6, -1, 14, 55},
             {-2, 0, -4, -12, 1, 5, 6},
             {1, -6, -16, -18, 4, 4, 33}};
    cout << endl << "RREF(E): " << endl;
    RREF(e).Print();
    cout << endl << "Rang(A): " << a.Rank() << endl << "Rang(B): " << b.Rank() << endl;
    cout << "Rang(C): " << Rank(c) << endl << "Rang(D): " << Rank(d) << endl << "Rang(E): " << Rank(e) << endl;
    cout << "RREF(A): " << endl;
    RREF(a).Print();
    LUDecomposer lud(Matrix{{1, -2, -2, -3},
                            {3, -9, 0, -9},
                            {-1, 2, 4, 7},
                            {-3, -6, 26, 2}});
    cout << endl << "L: " << endl;
    lud.GetL().Print();
    cout << endl << "U: " << endl;
    lud.GetU().Print();
    cout << endl << "W: " << endl;
    lud.GetPermuation().Print();
    cout << endl << "Kompaktno LU: " << endl;
    lud.GetCompactLU().Print();
    cout << endl << "L * U: " << endl;
    (lud.GetL() * lud.GetU()).Print();
    cout << endl << "LU dekompozicija nekvadratne matrice:" << endl;
    try {
        LUDecomposer(Matrix{{1, 2, 3}, {4, 5, 6}});
    } catch(exception &e) {
        cout << e.what() << endl;
    }
    cout << endl << "LU dekompozicija singularne matrice:" << endl;
    try {
        LUDecomposer(Matrix{{1, 2, 3}, {4, 5, 6}, {0, 0, 0}});
    } catch(exception &e) {
        cout << e.what() << endl;
    }
    cout << endl << "A / B(lu dekompozicija): " << endl;
    LUDecomposer(a).Solve(b).Print(); //Kako ova metoda vodi ka nizu poziva u kojem su sve 4 solve metode pozvane
                                      //broji se kao test svih ostalih
    QRDecomposer qrd(Matrix{{1, 1, 1}, {1, 2, 3}, {1, 3, 6}});
    cout << endl << "Test LU pivotizacije: " << endl;
    Matrix f{{0, 1, 2}, {82, 10, 4}, {54, 9, 1}};
    lud = LUDecomposer(f);
    cout << endl << "F: " << endl;
    f.Print();
    cout << endl << "L * U: " << endl;
    (lud.GetL() * lud.GetU()).Print();
    cout << endl << "Paskalova matrica: " << endl;
    cout << endl << "R: " << endl;
    qrd.GetR().Print();
    cout << endl << "Q: " << endl;
    qrd.GetQ().Print();
    cout << endl << "Q * R: " << endl;
    (qrd.GetQ() * qrd.GetR()).Print();
    cout << endl << "Q * Q^T: " << endl;
    cout << endl << "QR dekompozicija singularne matrice:" << endl;
    try {
        QRDecomposer(Matrix{{1, 2, 3}, {4, 5, 6}, {0, 0, 0}});
    } catch(exception &e) {
        cout << e.what() << endl;
    }
    Matrix E(qrd.GetQ().NRows(), qrd.GetQ().NCols());
    for(int i = 0; i < E.NRows(); i++)
        E[i][i] = 1;
    qrd.MulQWith(qrd.MulQTWith(E)).Print();
    cout << endl << "Mnozenje nepodudarajucih formata: " << endl;
    try {
        qrd.MulQTWith(Vector(qrd.GetQ().NCols() + 1));
    } catch(exception &e) {
        cout << e.what() << endl;
    }
    cout << endl << "A / B(qr dekompozicija): " << endl;
    QRDecomposer(a).Solve(b).Print();
    cout << endl << "Test QR pivotizacije: " << endl;
    qrd = QRDecomposer(f);
    cout << endl << "F: " << endl;
    f.Print();
    cout << endl << "Q * R: " << endl;
    (qrd.GetQ() * qrd.GetR()).Print();
    return 0;
}
