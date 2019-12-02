#include <iostream>
#include <cmath>
#include "../math_obj.hpp"

// Fundamental Constants
#define dt 1e-12                            // s
#define REAL_SIZE Vec3D(1, 1, 1)            // m
constexpr SizeTuple3D SIM_RES = {20, 10, 100}; // dimensionless max of 13500 elements

// Spacing between elements on the matrix
Vec3D spacing = Vec3D(REAL_SIZE.x / SIM_RES.x, REAL_SIZE.y / SIM_RES.y, REAL_SIZE.z / SIM_RES.z);

// This function converts matrix indices to a real position
Vec3D getPosition(int i, int j, int k)
{
    return Vec3D(i, j, k) * spacing + spacing / 2.0;
}

double moles_H2 = 1;
double t = 0;
double adiabatic_index = 5.0 / 3.0;
double total_charge = 2 * moles_H2 * E_CHARGE * AVOGADRO;
double total_mass = moles_H2 * 2 * M_PROTON * AVOGADRO;

// Relations:
// momentum = mass_density * velocity
// pressure = 2 * k_B/m_P * mass_density * temperature

int main()
{

    // Initialize Fields
    VectorField<SIM_RES.x, SIM_RES.y, SIM_RES.z> magneticField;
    VectorField<SIM_RES.x, SIM_RES.y, SIM_RES.z> velocity;
    ScalarField<SIM_RES.x, SIM_RES.y, SIM_RES.z> massDensity;
    VectorField<SIM_RES.x, SIM_RES.y, SIM_RES.z> currentDensity;
    ScalarField<SIM_RES.x, SIM_RES.y, SIM_RES.z> temperature;
    ScalarField<SIM_RES.x, SIM_RES.y, SIM_RES.z> pressure;
    for(int i = 0; i < SIM_RES.x; ++i)
    {
        for(int j = 0; j < SIM_RES.y; ++j)
        {
            for(int k = 0; k < SIM_RES.z; ++k)
            {
                massDensity(i, j, k, total_mass/(SIM_RES.y * SIM_RES.z * SIM_RES.x));
                velocity(i, j, k, Vec3D(10000, 0, 0));
                temperature(i, j, k, 1000);
            }
        }
    }
    // Calculation Intermediaries
    VectorField<SIM_RES.x, SIM_RES.y, SIM_RES.z> pressureGrad;
    ScalarField<SIM_RES.x, SIM_RES.y, SIM_RES.z> udotgradus[3];
    VectorField<SIM_RES.x, SIM_RES.y, SIM_RES.z> udotgradu;

    pressure << temperature * massDensity * 2 * (k_B / M_PROTON);
    currentDensity << velocity * massDensity * (total_charge / total_mass);

    Vec3D r_prime;
    for(int i1 = 0; i1 < SIM_RES.x; ++i1)
    {
        for(int j1 = 0; j1 < SIM_RES.y; ++j1)
        {
            for(int k1 = 0; k1 < SIM_RES.z; ++k1)
            {
                for(int i2 = 0; i2 < SIM_RES.x; ++i2)
                {
                    for(int j2 = 0; j2 < SIM_RES.y; ++j2)
                    {
                        for(int k2 = 0; k2 < SIM_RES.z; ++k2)
                        {
                            r_prime = getPosition(i1, j1, k1) - getPosition(i2, j2, k2);
                            if(r_prime.normSq() != 0)
                            {
                                // std::cout << Vec3D(i1, j1, k1).sstr() << " " << (velocity(i2, j2, k2).cross(r_prime) * (E_CHARGE/M_PROTON) * (MU_0/(4 * PI * pow(r_prime.norm(), 3)))*massDensity(i2, j2, k2)).sstr() << std::endl;
                                magneticField.data[magneticField.index(i1, j1, k1)] = magneticField.data[magneticField.index(i1, j1, k1)] + velocity(i2, j2, k2).cross(r_prime) * (E_CHARGE/M_PROTON) * massDensity(i2, j2, k2) * (MU_0/(4 * PI * pow(r_prime.norm(), 3)));
                            }
                        }
                    }
                }
                std::cout << magneticField(i1, j1, k1).sstr() << std::endl;
            }
        }
    }
    // Total Bytes of Storage
    // std::cout << sizeof(magneticField) + sizeof(velocity) + sizeof(massDensity) + sizeof(currentDensity) + sizeof(temperature) + sizeof(pressure) + sizeof(pressureGrad) + sizeof(udotgradus) + sizeof(udotgradu) << std::endl;

    // Time Propagation
    // Equations come from pages 9 and 10 of "Computational Methods in Plasma Physics" written by Steven Jardin and published on June 2nd, 2010
    // Equations 1.60, 1.61, 1.62, 1.63
    do
    {
        massDensity += (velocity * massDensity).divergence(spacing) * -dt;
        magneticField += (velocity.cross(magneticField)).curl(spacing) * dt;
        udotgradu << velocity.getX().gradient(spacing);
        udotgradus[0] << velocity * udotgradu;
        udotgradu << velocity.getY().gradient(spacing);
        udotgradus[1] << velocity * udotgradu;
        udotgradu << velocity.getZ().gradient(spacing);
        udotgradus[2] << velocity * udotgradu;
        udotgradu << udotgradus;
        pressureGrad << pressure.gradient(spacing);
        currentDensity << velocity * massDensity * (total_charge / total_mass);
        velocity += ((((currentDensity.cross(magneticField)) + (pressureGrad * -1)) * (massDensity ^ (-1))) + (udotgradu * -1)) * dt;
        pressure += ((velocity * pressureGrad) + (velocity.divergence(spacing) * pressure * adiabatic_index)) * -dt;
        temperature << pressure * (massDensity ^ (-1)) * (M_PROTON / (2 * k_B));
        t += dt;
        if((int)(t/dt) % 100 == 0)
        {
            std::cout << t << std::endl;
        }
    } while (t < 1e-9);
}