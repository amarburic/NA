#include <iostream>
#include <vector>
#include <stdexcept>
#include <iomanip>
#include <algorithm>

using namespace std;

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
};

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

int main()
{
    Matrix m1(3, 5), m2{{1, 2, 3, 4, 5}, {6, 7, 8, 9, 10}, {11, 12, 13, 14, 15}}, m3(Vector({5, 4, 3, 2, 1}));
    try {
        Matrix(-1, 1);
    } catch(exception& e) {
        cout << e.what() << endl;
    }
    try {
        Matrix{{}, {}};
    } catch(exception& e) {
        cout << e.what() << endl;
    }
    try {
        Matrix{{}, {1, 2}};
    } catch(exception& e) {
        cout << e.what() << endl;
    }
    try {
        Matrix{{1, 2, 3}, {5, 6}};
    } catch(exception& e) {
        cout << e.what() << endl;
    }

    const Matrix& crefm1(m1);
    cout << crefm1.NRows() << " " << crefm1.NCols() << endl;
    m1[0][0] = 1;
    cout << crefm1[0][0] << endl;
    cout << m2(2, 4) << " " << crefm1(1, 1) << endl;
    try {
        cout << m2(1, 6) << endl;
    } catch(exception& e) {
        cout << e.what() << endl;
    }
    m1.Print(5);
    cout << endl;
    m2.Print(5);
    cout << endl;
    m3.Print(5);
    cout << endl;
    (m1 += m2).Print(5);
    cout << endl;
    (m1 *= 5).Print(5);
    cout << endl;
    m1.Transpose();
    m2.Transpose();
    (m1 -= m2).Print(5);
    cout << endl;
    m2.Transpose();
    (m2 *= m3).Print(5);
    cout << endl;
    Matrix m4{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    m4.Transpose();
    m4.Print(5);
    Vector v1(10), v2{1, 2, 3, 4, 5};
    try {
        Vector v3(-5);
    } catch(exception& e) {
        cout << e.what() << endl;
    }
    try {
        Vector v3{};
    } catch(exception& e) {
        cout << e.what() << endl;
    }
    const Vector v3(v2);
    for(int i = 0; i < v1.NElems(); i++)
        cout << v1[i] << " ";
    cout << endl;
    for(int i = 0; i < v2.NElems(); i++)
        cout << v2[i] << " ";
    cout << endl;
    v2[0] = 100;
    v2.Print();
    //v3[0] *= 100;
    //v3.Print(",");
    try {
        v2(5) = 6;
        v2(7) = 9;
    } catch(exception& e) {
        cout << e.what() << endl;
    }
    try {
        v2(-3) = 6;
    } catch(exception& e) {
        cout << e.what() << endl;
    }
    v2.Print(',');
    cout << endl;
    try {
        (v2 - v3).Print(',');
        cout << endl;
        (v1 + v2).Print(',');
        cout << endl;
    } catch(exception& e) {
        cout << e.what() << endl;
    }
    (3 * v3).Print(',');
    cout << endl;
    (v2 *= 2).Print(',');
    cout << endl;
    (v2 /= 5).Print(',');
    try {
        cout << v2 * v3 << endl;
        cout << v1 * v2 << endl;
    } catch(exception& e) {
        cout << e.what() << endl;
    }
    try {
        v1 / 0;
    } catch(exception& e) {
        cout << e.what() << endl;
    }
    return 0;
}
