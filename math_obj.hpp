#include <iostream>
#include <cmath>
#include <fstream>
#include <any>
#include <cstring>
#include <array>
#include <vector>

#define INF_VEC Vec3D(INFINITY, INFINITY, INFINITY)
#define ZERO_VEC Vec3D(0, 0, 0)
#define NAN_VEC Vec3D(NAN, NAN, NAN);
#define X_HAT Vec3D(1, 0, 0)
#define Y_HAT Vec3D(0, 1, 0)
#define Z_HAT Vec3D(0, 0, 1)
#define EPSILON_0 8.85418782e-12        // m^-3 kg^-1 s^4 A^2
#define MU_0 1.25663706e-6              // m kg s^-2 A^-2
#define PI M_PI                         // dimensionless
#define LIGHT_SPEED 299792458.0         // m s^-1
#define ETA LIGHT_SPEED * 4 * PI * 1e-7 // m^2 kg s^-3 A^-2
#define M_PROTON 1.67262192369e-27      // kg
#define M_ELECTRON 9.10938356e-31       // kg
#define M_NEUTRON 1.674927471e-27       // kg
#define E_CHARGE 1.60217662e-19         // s A
#define AVOGADRO 6.02214086e+23         // dimensionless
#define k_B 1.38064852e-23              // m^2 kg s^-2 K^-1

// A vector function class
class Vec3D
{
public:
    double x, y, z;
    Vec3D(double x_, double y_, double z_)
    {
        x = x_;
        y = y_;
        z = z_;
    }

    Vec3D()
    {
        x = 0;
        y = 0;
        z = 0;
    }

    Vec3D operator-=(Vec3D n)
    {
        this->x - n.x;
        this->y - n.y;
        this->z - n.z;
        return *this;
    }

    Vec3D operator+=(Vec3D n)
    {
        this->x + n.x;
        this->y + n.y;
        this->z + n.z;
        return *this;
    }

    Vec3D operator-(Vec3D n)
    {
        return add(n * -1);
    }

    Vec3D operator-()
    {
        return *this * -1;
    }

    Vec3D operator+(Vec3D n)
    {
        return add(n);
    }

    Vec3D operator*(Vec3D n)
    {
        return multiply(n);
    }

    Vec3D operator*=(Vec3D n)
    {
        this->x *n.x;
        this->y *n.y;
        this->z *n.z;
        return *this;
    }

    Vec3D operator*(double n)
    {
        return multiply(n);
    }

    Vec3D operator*=(double n)
    {
        this->x *n;
        this->y *n;
        this->z *n;
        return *this;
    }

    Vec3D operator/(double n)
    {
        return multiply(1.0 / n);
    }

    double operator%(Vec3D n)
    {
        return dot(n);
    }

    Vec3D operator()(double x, double y, double z)
    {
        Vec3D d = Vec3D(x, y, z);
        return d;
    }

    bool operator==(Vec3D n)
    {
        return n.x == x && n.y == y && n.z == z;
    }

    double dot(Vec3D n)
    {
        return x * n.x + y * n.y + z * n.z;
    }

    Vec3D cross(Vec3D n)
    {
        return Vec3D(y * n.z - z * n.y, z * n.x - x * n.z, x * n.y - y * n.x);
    }

    Vec3D add(Vec3D n)
    {
        return Vec3D(x + n.x, y + n.y, z + n.z);
    }

    Vec3D add(double n)
    {
        return Vec3D(x + n, y + n, z + n);
    }

    Vec3D multiply(Vec3D n)
    {
        return Vec3D(x * n.x, y * n.y, z * n.z);
    }

    Vec3D multiply(double n)
    {
        return Vec3D(x * n, y * n, z * n);
    }

    std::string cstr()
    {
        return "x: " + std::to_string(x) + "\r\ny: " + std::to_string(y) + "\r\nz: " + std::to_string(z);
    }

    std::string vstr()
    {
        return "<" + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) + ">";
    }

    std::string sstr()
    {
        return std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z);
    }

    std::string cstr(Vec3D n)
    {
        return "x: " + std::to_string(n.x) + "\r\ny: " + std::to_string(n.y) + "\r\nz: " + std::to_string(n.z);
    }

    std::string vstr(Vec3D n)
    {
        return "<" + std::to_string(n.x) + ", " + std::to_string(n.y) + ", " + std::to_string(n.z) + ">";
    }

    std::string sstr(Vec3D n)
    {
        return std::to_string(n.x) + " " + std::to_string(n.y) + " " + std::to_string(n.z);
    }

    Vec3D sinVec(Vec3D v)
    {
        return Vec3D(sin(v.x), sin(v.y), sin(v.z));
    }

    Vec3D cosVec(Vec3D v)
    {
        return Vec3D(cos(v.x), cos(v.y), cos(v.z));
    }

    Vec3D tanVec(Vec3D v)
    {
        return Vec3D(tan(v.x), tan(v.y), tan(v.z));
    }

    double norm(Vec3D n)
    {
        return sqrt(n.dot(n));
    }

    double norm()
    {
        return sqrt(x * x + y * y + z * z);
    }

    double normSq()
    {
        return x * x + y * y + z * z;
    }

    Vec3D normalize(Vec3D n)
    {
        return n / n.norm();
    }

    Vec3D normalize()
    {
        double l = norm();
        x /= l;
        y /= l;
        z /= l;
        return *this;
    }

    Vec3D transform(double matrix[3][3])
    {
        double x_ = matrix[0][0] * x + matrix[0][1] * y + matrix[0][2] * z;
        double y_ = matrix[1][0] * x + matrix[1][1] * y + matrix[1][2] * z;
        double z_ = matrix[2][0] * x + matrix[2][1] * y + matrix[2][2] * z;
        return Vec3D(x_, y_, z_);
    }

    Vec3D operator*(double matrix[3][3])
    {
        return transform(matrix);
    }

    Vec3D operator*=(double matrix[3][3])
    {
        return this->transform(matrix);
    }
};

// Holds the size (requires integers)
struct SizeTuple
{
    const int x, y, z;
};

// Represents a 3-D field of Scalars
template <uint a, uint b, uint c>
class ScalarField
{
public:
    unsigned const int dims = 3;
    const unsigned int rows = a;
    const unsigned int cols = b;
    const unsigned int depth = c;
    std::array<unsigned int, 3> matrixSize = {rows, cols, depth};
    double *data = (double*) malloc(sizeof(double) * a * b * c);

    ScalarField(){};

    int index(uint i, uint j, uint k)
    {
        return k + j * b + i * a * b;
    }

    double operator()(unsigned int i, unsigned int j, unsigned int k)
    {
        return data[index(i, j, k)];
    }

    void operator()(uint i, uint j, uint k, double v)
    {
        data[index(i, j, k)] = v;
    }

    ScalarField<a, b, c> operator+(double n)
    {
        ScalarField<a, b, c> N;
        for (int i = 0; i < matrixSize[0]; ++i)
        {
            for (int j = 0; j < matrixSize[1]; ++j)
            {
                for (int k = 0; k < matrixSize[2]; ++k)
                {
                    N(i, j, k, operator()(i, j, k) + n);
                }
            }
        }
        return N;
    }

    ScalarField<a, b, c> operator*(double n)
    {
        ScalarField<a, b, c> N;
        for (int i = 0; i < matrixSize[0]; ++i)
        {
            for (int j = 0; j < matrixSize[1]; ++j)
            {
                for (int k = 0; k < matrixSize[2]; ++k)
                {
                    N(i, j, k, operator()(i, j, k) * n);
                }
            }
        }
        return N;
    }

    ScalarField<a, b, c> operator+(ScalarField<a, b, c> A)
    {
        ScalarField<a, b, c> N;
        for (int i = 0; i < matrixSize[0]; ++i)
        {
            for (int j = 0; j < matrixSize[1]; ++j)
            {
                for (int k = 0; k < matrixSize[2]; ++k)
                {
                    N(i, j, k, operator()(i, j, k) + A(i, j, k));
                }
            }
        }
        return N;
    }

    ScalarField<a, b, c> operator*(ScalarField<a, b, c> A)
    {
        ScalarField<a, b, c> N;
        for (int i = 0; i < matrixSize[0]; ++i)
        {
            for (int j = 0; j < matrixSize[1]; ++j)
            {
                for (int k = 0; k < matrixSize[2]; ++k)
                {
                    N(i, j, k, operator()(i, j, k) * A(i, j, k));
                }
            }
        }
        return N;
    }

    ScalarField<a, b, c> operator+=(double n)
    {
        for (int i = 0; i < matrixSize[0]; ++i)
        {
            for (int j = 0; j < matrixSize[1]; ++j)
            {
                for (int k = 0; k < matrixSize[2]; ++k)
                {
                    data[index(i, j, k)] += n;
                }
            }
        }
        return *this;
    }

    ScalarField<a, b, c> operator*=(double n)
    {
        for (int i = 0; i < matrixSize[0]; ++i)
        {
            for (int j = 0; j < matrixSize[1]; ++j)
            {
                for (int k = 0; k < matrixSize[2]; ++k)
                {
                    data[index(i, j, k)] *= n;
                }
            }
        }
        return *this;
    }

    ScalarField<a, b, c> operator+=(ScalarField<a, b, c> A)
    {
        for (int i = 0; i < matrixSize[0]; ++i)
        {
            for (int j = 0; j < matrixSize[1]; ++j)
            {
                for (int k = 0; k < matrixSize[2]; ++k)
                {
                    data[index(i, j, k)] += A(i, j, k);
                }
            }
        }
        return *this;
    }

    ScalarField<a, b, c> operator*=(ScalarField<a, b, c> A)
    {
        for (int i = 0; i < matrixSize[0]; ++i)
        {
            for (int j = 0; j < matrixSize[1]; ++j)
            {
                for (int k = 0; k < matrixSize[2]; ++k)
                {
                    data[index(i, j, k)] *= A(i, j, k);
                }
            }
        }
        return *this;
    }

    ScalarField<a, b, c> operator^(double n)
    {
        ScalarField<a, b, c> N;
        for (int i = 0; i < matrixSize[0]; ++i)
        {
            for (int j = 0; j < matrixSize[1]; ++j)
            {
                for (int k = 0; k < matrixSize[2]; ++k)
                {
                    if(operator()(i, j, k) == 0)
                    {
                        N(i, j, k, 0);
                    }
                    else
                    {
                        N(i, j, k, pow(operator()(i, j, k), n));
                    }
                }
            }
        }
        return N;
    }

    ScalarField<a, b, c> operator<<(double dataIn[a][b][c])
    {
        for (int i = 0; i < matrixSize[0]; ++i)
        {
            for (int j = 0; j < matrixSize[1]; ++j)
            {
                for (int k = 0; k < matrixSize[2]; ++k)
                {
                    operator()(i, j, k, dataIn[i][j][k]);
                }
            }
        }
        return *this;
    }

    ScalarField<a, b, c> operator<<(std::array<std::array<std::array<double, c>, b>, a> dataIn)
    {
        for (int i = 0; i < matrixSize[0]; ++i)
        {
            for (int j = 0; j < matrixSize[1]; ++j)
            {
                for (int k = 0; k < matrixSize[2]; ++k)
                {
                    operator()(i, j, k, dataIn[i][j][k]);
                }
            }
        }
        return *this;
    }

    ScalarField<a, b, c> operator<<(ScalarField<a, b, c> dataIn)
    {
        for (int i = 0; i < matrixSize[0]; ++i)
        {
            for (int j = 0; j < matrixSize[1]; ++j)
            {
                for (int k = 0; k < matrixSize[2]; ++k)
                {
                    operator()(i, j, k, dataIn(i, j, k));
                }
            }
        }
        return *this;
    }

    std::array<std::array<std::array<Vec3D, c>, b>, a> gradient(Vec3D spacing)
    {
        std::array<std::array<std::array<Vec3D, c>, b>, a> Gradiented;
        double dFdx, dFdy, dFdz;
        for (int i = 0; i < matrixSize[0]; ++i)
        {
            for (int j = 0; j < matrixSize[1]; ++j)
            {
                for (int k = 0; k < matrixSize[2]; ++k)
                {
                    if (i > 0 && j > 0 && k > 0 && i + 1 < matrixSize[0] && j + 1 < matrixSize[1] && k + 1 < matrixSize[2])
                    {
                        dFdx = (operator()(i + 1, j, k) - operator()(i - 1, j, k)) / (2 * spacing.x);
                        dFdy = (operator()(i, j + 1, k) - operator()(i, j - 1, k)) / (2 * spacing.y);
                        dFdz = (operator()(i, j, k + 1) - operator()(i, j, k - 1)) / (2 * spacing.z);
                    }
                    if (i == 0)
                    {
                        dFdx = (operator()(i + 1, j, k) - operator()(i, j, k)) / (spacing.x);
                    }
                    if (i + 1 == matrixSize[0])
                    {
                        dFdx = (operator()(i, j, k) - operator()(i - 1, j, k)) / (spacing.x);
                    }
                    if (k == 0)
                    {
                        dFdz = (operator()(i, j, k + 1) - operator()(i, j, k)) / (spacing.z);
                    }
                    if (k + 1 == matrixSize[2])
                    {
                        dFdz = (operator()(i, j, k) - operator()(i, j, k - 1)) / (spacing.z);
                    }
                    if (j == 0)
                    {
                        dFdy = (operator()(i, j + 1, k) - operator()(i, j, k)) / (spacing.y);
                    }
                    if (j + 1 == matrixSize[1])
                    {
                        dFdy = (operator()(i, j, k) - operator()(i, j - 1, k)) / (spacing.y);
                    }
                    Gradiented[i][j][k] = Vec3D(dFdx, dFdy, dFdz);
                }
            }
        }
        return Gradiented;
    }

    void toFile(std::string filename)
    {
        std::ofstream myfile(filename);
        if (myfile.is_open())
        {
            myfile << matrixSize[0] << "\n";
            myfile << matrixSize[1] << "\n";
            myfile << matrixSize[2] << "\n";

            for (int i = 0; i < matrixSize[0]; ++i)
            {
                for (int j = 0; j < matrixSize[1]; ++j)
                {
                    for (int k = 0; k < matrixSize[2]; ++k)
                    {
                        myfile << i << " " << j << " " << k << " " << operator()(i, j, k) << "\n";
                    }
                }
            }
            myfile.close();
        }
        else
        {
            std::cout << filename + " was not openable" << std::endl;
        }
    }
};

// Represents a 3-D field of Vec3D
template <uint a, uint b, uint c>
class VectorField
{
public:
    unsigned const int dims = 3;
    const unsigned int rows = a;
    const unsigned int cols = b;
    const unsigned int depth = c;
    std::array<unsigned int, 3> matrixSize = {rows, cols, depth};
    Vec3D *data = (Vec3D*) malloc(sizeof(Vec3D) * a * b * c);

    VectorField(){};

    int index(uint i, uint j, uint k)
    {
        return k + j * b + i * a * b;
    }

    Vec3D operator()(unsigned int i, unsigned int j, unsigned int k)
    {
        return data[index(i, j, k)];
    }

    void operator()(uint i, uint j, uint k, Vec3D v)
    {
        data[index(i, j, k)] = v;
    }

    ScalarField<a, b, c> getX()
    {
        ScalarField<a, b, c> X;
        for (int i = 0; i < matrixSize[0]; ++i)
        {
            for (int j = 0; j < matrixSize[1]; ++j)
            {
                for (int k = 0; k < matrixSize[2]; ++k)
                {
                    X(i, j, k, operator()(i, j, k).x);
                }
            }
        }
        return X;
    }

    ScalarField<a, b, c> getY()
    {
        ScalarField<a, b, c> Y;
        for (int i = 0; i < matrixSize[0]; ++i)
        {
            for (int j = 0; j < matrixSize[1]; ++j)
            {
                for (int k = 0; k < matrixSize[2]; ++k)
                {
                    Y(i, j, k, operator()(i, j, k).y);
                }
            }
        }
        return Y;
    }

    ScalarField<a, b, c> getZ()
    {
        ScalarField<a, b, c> Z;
        for (int i = 0; i < matrixSize[0]; ++i)
        {
            for (int j = 0; j < matrixSize[1]; ++j)
            {
                for (int k = 0; k < matrixSize[2]; ++k)
                {
                    Z(i, j, k, operator()(i, j, k).z);
                }
            }
        }
        return Z;
    }

    VectorField<a, b, c> operator+(Vec3D n)
    {
        VectorField<a, b, c> N;
        for (int i = 0; i < matrixSize[0]; ++i)
        {
            for (int j = 0; j < matrixSize[1]; ++j)
            {
                for (int k = 0; k < matrixSize[2]; ++k)
                {
                    N(i, j, k, operator()(i, j, k) + n);
                }
            }
        }
        return N;
    }

    VectorField<a, b, c> operator*(double n)
    {
        VectorField<a, b, c> N;
        for (int i = 0; i < matrixSize[0]; ++i)
        {
            for (int j = 0; j < matrixSize[1]; ++j)
            {
                for (int k = 0; k < matrixSize[2]; ++k)
                {
                    N(i, j, k, operator()(i, j, k) * n);
                }
            }
        }
        return N;
    }

    VectorField<a, b, c> operator+(VectorField<a, b, c> A)
    {
        VectorField<a, b, c> N;
        for (int i = 0; i < matrixSize[0]; ++i)
        {
            for (int j = 0; j < matrixSize[1]; ++j)
            {
                for (int k = 0; k < matrixSize[2]; ++k)
                {
                    N(i, j, k, operator()(i, j, k) + A(i, j, k));
                }
            }
        }
        return N;
    }

    VectorField<a, b, c> operator*(ScalarField<a, b, c> A)
    {
        VectorField<a, b, c> N;
        for (int i = 0; i < matrixSize[0]; ++i)
        {
            for (int j = 0; j < matrixSize[1]; ++j)
            {
                for (int k = 0; k < matrixSize[2]; ++k)
                {
                    N(i, j, k, operator()(i, j, k) * A(i, j, k));
                }
            }
        }
        return N;
    }

    ScalarField<a, b, c> operator*(VectorField<a, b, c> A)
    {
        ScalarField<a, b, c> N;
        for (int i = 0; i < matrixSize[0]; ++i)
        {
            for (int j = 0; j < matrixSize[1]; ++j)
            {
                for (int k = 0; k < matrixSize[2]; ++k)
                {
                    N(i, j, k, operator()(i, j, k).dot(A(i, j, k)));
                }
            }
        }
        return N;
    }

    VectorField<a, b, c> operator+=(Vec3D n)
    {
        for (int i = 0; i < matrixSize[0]; ++i)
        {
            for (int j = 0; j < matrixSize[1]; ++j)
            {
                for (int k = 0; k < matrixSize[2]; ++k)
                {
                    data[index(i, j, k)] += n;
                }
            }
        }
        return *this;
    }

    VectorField<a, b, c> operator*=(double n)
    {
        for (int i = 0; i < matrixSize[0]; ++i)
        {
            for (int j = 0; j < matrixSize[1]; ++j)
            {
                for (int k = 0; k < matrixSize[2]; ++k)
                {
                    data[index(i, j, k)] *= n;
                }
            }
        }
        return *this;
    }

    VectorField<a, b, c> operator+=(VectorField<a, b, c> A)
    {
        for (int i = 0; i < matrixSize[0]; ++i)
        {
            for (int j = 0; j < matrixSize[1]; ++j)
            {
                for (int k = 0; k < matrixSize[2]; ++k)
                {
                    data[index(i, j, k)] += A(i, j, k);
                }
            }
        }
        return *this;
    }

    VectorField<a, b, c> operator*=(ScalarField<a, b, c> A)
    {
        for (int i = 0; i < matrixSize[0]; ++i)
        {
            for (int j = 0; j < matrixSize[1]; ++j)
            {
                for (int k = 0; k < matrixSize[2]; ++k)
                {
                    data[index(i, j, k)] *= A(i, j, k);
                }
            }
        }
        return *this;
    }

    VectorField<a, b, c> operator<<(Vec3D dataIn[a][b][c])
    {
        for (int i = 0; i < matrixSize[0]; ++i)
        {
            for (int j = 0; j < matrixSize[1]; ++j)
            {
                for (int k = 0; k < matrixSize[2]; ++k)
                {
                    operator()(i, j, k, dataIn[i][j][k]);
                }
            }
        }
        return *this;
    }

    VectorField<a, b, c> operator<<(std::array<std::array<std::array<Vec3D, c>, b>, a> dataIn)
    {
        for (int i = 0; i < matrixSize[0]; ++i)
        {
            for (int j = 0; j < matrixSize[1]; ++j)
            {
                for (int k = 0; k < matrixSize[2]; ++k)
                {
                    operator()(i, j, k, dataIn[i][j][k]);
                }
            }
        }
        return *this;
    }

    VectorField<a, b, c> operator<<(ScalarField<a, b, c> fields[3])
    {
        for (int i = 0; i < matrixSize[0]; ++i)
        {
            for (int j = 0; j < matrixSize[1]; ++j)
            {
                for (int k = 0; k < matrixSize[2]; ++k)
                {
                    operator()(i, j, k, Vec3D(fields[0](i, j, k), fields[1](i, j, k), fields[2](i, j, k)));
                }
            }
        }
        return *this;
    }

    VectorField<a, b, c> operator<<(VectorField<a, b, c> field)
    {
        for (int i = 0; i < matrixSize[0]; ++i)
        {
            for (int j = 0; j < matrixSize[1]; ++j)
            {
                for (int k = 0; k < matrixSize[2]; ++k)
                {
                    operator()(i, j, k, field(i, j, k));
                }
            }
        }
        return *this;
    }

    VectorField<a, b, c> cross(VectorField<a, b, c> A)
    {
        VectorField<a, b, c> N;
        for (int i = 0; i < matrixSize[0]; ++i)
        {
            for (int j = 0; j < matrixSize[1]; ++j)
            {
                for (int k = 0; k < matrixSize[2]; ++k)
                {
                    N(i, j, k, operator()(i, j, k).cross(A(i, j, k)));
                }
            }
        }
        return N;
    }

    VectorField<a, b, c> curl(Vec3D spacing)
    {
        double dXdy, dXdz, dYdx, dYdz, dZdx, dZdy;
        VectorField<a, b, c> Curled;
        for (int i = 0; i < matrixSize[0]; ++i)
        {
            for (int j = 0; j < matrixSize[1]; ++j)
            {
                for (int k = 0; k < matrixSize[2]; ++k)
                {
                    if (i > 0 && j > 0 && k > 0 && i + 1 < matrixSize[0] && j + 1 < matrixSize[1] && k + 1 < matrixSize[2])
                    {
                        dXdy = (operator()(i, j + 1, k).x - operator()(i, j - 1, k).x) / (2 * spacing.y);
                        dXdz = (operator()(i, j, k + 1).x - operator()(i, j, k - 1).x) / (2 * spacing.z);
                        dYdx = (operator()(i + 1, j, k).y - operator()(i - 1, j, k).y) / (2 * spacing.x);
                        dYdz = (operator()(i, j, k + 1).y - operator()(i, j, k - 1).y) / (2 * spacing.z);
                        dZdx = (operator()(i + 1, j, k).z - operator()(i - 1, j, k).z) / (2 * spacing.x);
                        dZdy = (operator()(i, j + 1, k).z - operator()(i, j - 1, k).z) / (2 * spacing.y);
                    }
                    if (i == 0)
                    {
                        dZdx = (operator()(i + 1, j, k).z - operator()(i, j, k).z) / (spacing.x);
                        dYdx = (operator()(i + 1, j, k).y - operator()(i, j, k).y) / (spacing.x);
                    }
                    if (i + 1 == matrixSize[0])
                    {
                        dZdx = (operator()(i, j, k).z - operator()(i - 1, j, k).z) / (spacing.x);
                        dYdx = (operator()(i, j, k).y - operator()(i - 1, j, k).y) / (spacing.x);
                    }
                    if (k == 0)
                    {
                        dXdz = (operator()(i, j, k + 1).x - operator()(i, j, k).x) / (spacing.z);
                        dYdz = (operator()(i, j, k + 1).y - operator()(i, j, k).y) / (spacing.z);
                    }
                    if (k + 1 == matrixSize[2])
                    {
                        dXdz = (operator()(i, j, k).x - operator()(i, j, k - 1).x) / (spacing.z);
                        dYdz = (operator()(i, j, k).y - operator()(i, j, k - 1).y) / (spacing.z);
                    }
                    if (j == 0)
                    {
                        dXdy = (operator()(i, j + 1, k).x - operator()(i, j, k).x) / (spacing.y);
                        dZdy = (operator()(i, j + 1, k).z - operator()(i, j, k).z) / (spacing.y);
                    }
                    if (j + 1 == matrixSize[1])
                    {
                        dXdy = (operator()(i, j, k).x - operator()(i, j - 1, k).x) / (spacing.y);
                        dZdy = (operator()(i, j, k).z - operator()(i, j - 1, k).z) / (spacing.y);
                    }
                    Curled(i, j, k, Vec3D(dZdy - dYdz, dXdz - dZdx, dYdx - dXdy));
                }
            }
        }
        return Curled;
    }

    ScalarField<a, b, c> divergence(Vec3D spacing)
    {
        double dXdx, dYdy, dZdz;
        ScalarField<a, b, c> Diverged;
        for (int i = 0; i < matrixSize[0]; ++i)
        {
            for (int j = 0; j < matrixSize[1]; ++j)
            {
                for (int k = 0; k < matrixSize[2]; ++k)
                {
                    if (i > 0 && j > 0 && k > 0 && i + 1 < matrixSize[0] && j + 1 < matrixSize[1] && k + 1 < matrixSize[2])
                    {
                        dXdx = (operator()(i + 1, j, k).x - operator()(i - 1, j, k).x) / (2 * spacing.x);
                        dYdy = (operator()(i, j + 1, k).y - operator()(i, j - 1, k).y) / (2 * spacing.y);
                        dZdz = (operator()(i, j, k + 1).z - operator()(i, j, k - 1).z) / (2 * spacing.z);
                    }
                    if (i == 0)
                    {
                        dXdx = (operator()(i + 1, j, k).x - operator()(i, j, k).x) / (spacing.x);
                    }
                    if (i + 1 == matrixSize[0])
                    {
                        dXdx = (operator()(i, j, k).x - operator()(i - 1, j, k).x) / (spacing.x);
                    }
                    if (k == 0)
                    {
                        dZdz = (operator()(i, j, k + 1).z - operator()(i, j, k).z) / (spacing.z);
                    }
                    if (k + 1 == matrixSize[2])
                    {
                        dZdz = (operator()(i, j, k).z - operator()(i, j, k - 1).z) / (spacing.z);
                    }
                    if (j == 0)
                    {
                        dYdy = (operator()(i, j + 1, k).y - operator()(i, j, k).y) / (spacing.y);
                    }
                    if (j + 1 == matrixSize[1])
                    {
                        dYdy = (operator()(i, j, k).y - operator()(i, j - 1, k).y) / (spacing.y);
                    }
                    Diverged(i, j, k, dXdx + dYdy + dZdz);
                }
            }
        }
        return Diverged;
    }
    void toFile(std::string filename)
    {
        std::ofstream myfile(filename);
        if (myfile.is_open())
        {
            myfile << matrixSize[0] << "\n";
            myfile << matrixSize[1] << "\n";
            myfile << matrixSize[2] << "\n";

            for (int i = 0; i < matrixSize[0]; ++i)
            {
                for (int j = 0; j < matrixSize[1]; ++j)
                {
                    for (int k = 0; k < matrixSize[2]; ++k)
                    {
                        std::cout << i << " " << j << " " << k << " " << operator()(i, j, k).sstr() << "\n";
                        myfile << i << " " << j << " " << k << " " << operator()(i, j, k).sstr() << "\n";
                    }
                }
            }
            myfile.close();
        }
        else
        {
            std::cout << filename + " was not openable" << std::endl;
        }
    }
};