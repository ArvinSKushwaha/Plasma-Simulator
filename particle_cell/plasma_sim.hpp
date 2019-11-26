#include <cmath>
#include "../math_obj.hpp"

double gammaFromV(Vec3D velocity)
{
    return 1.0 / sqrt(1.0 - velocity.normSq() / pow(LIGHT_SPEED, 2));
}

template <uint a, uint b, uint c>
class YeeLattice
{
public:
    VectorField<a + 1, b + 1, c + 1> electricField;
    VectorField<a + 1, b + 1, c + 1> currentDensity;
    VectorField<a, b, c> magneticField;

    Vec3D getElectricField(Vec3D position)
    {
        return ZERO_VEC;
    }

    Vec3D getCurrentDensity(Vec3D position)
    {
        return ZERO_VEC;
    }

    Vec3D getMagneticField(Vec3D position)
    {
        return ZERO_VEC;
    }
};

class Particle
{
public:
    Vec3D position, velocity;
    double charge, mass, gamma;

    Particle(Vec3D pos, Vec3D vel, double q, double m)
    {
        position = pos;
        velocity = vel;
        charge = q;
        mass = m;
        gamma = gammaFromV(vel);
    }

    template <uint a, uint b, uint c>
    void update(double dt, YeeLattice<a, b, c> fields)
    {
        velocity += (fields.getElectricField(position) * charge / mass) * (dt / gamma);
    }
};