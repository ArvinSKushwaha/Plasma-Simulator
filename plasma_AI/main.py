# Replace keras with tensorflow.keras for Tensorflow 2.0
from keras import layers, models, optimizers, losses, backend

total_charge = 100
total_mass = 100

res = (100, 100, 100)
scalar_res = (100, 100, 100, 1)
vector_res = (100, 100, 100, 3)
x_res, y_res, z_res = res

inputs = [
    layers.Input(shape=vector_res),  # Magnetic Field
    layers.Input(shape=vector_res),  # Velocity Field
    layers.Input(shape=scalar_res),  # Pressure Field
    layers.Input(shape=scalar_res),  # Mass Density
]

gradConvPressure = layers.Conv3D(3, 3, padding="same")(inputs[2])
gradConvMassDensity = layers.Conv3D(3, 3, padding="same")(inputs[3])
divConvMagField = layers.Conv3D(1, 3, padding="same")(inputs[0])
divConvVelField = layers.Conv3D(1, 3, padding="same")(inputs[1])
curlConvMagField = layers.Conv3D(3, 3, padding="same")(inputs[0])
curlConvVelField = layers.Conv3D(3, 3, padding="same")(inputs[1])

udotgradu = layers.Conv3D(1, 3, padding="same")(layers.Multiply(inputs[1], backend.repeat_elements(divConvVelField, 3, -1)))

currentDensity = layers.Multiply(inputs[1], inputs[3], backend.ones_like(inputs[3])*total_charge/total_mass)

velocityNew = layers.Subtract(layers.Multiply(layers.Subtract(layers.Conv3D(3, 3, padding="same")(layers.Concatenate(currentDensity, inputs[0])), gradConvPressure), inputs[4]**-1), udotgradu)

# model = models.Model(inputs=inputs, outputs=outputs)