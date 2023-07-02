import numpy as np
from scipy.fft import fft2, fftn, fft

# Sinal de entrada de 2 dimensões
signal = np.array([[1, 2, 3],
                   [4, 5, 6],
                   [7, 8, 9]])

# Calculando a DFT bidimensional usando fft2
dft2 = fft2(signal)

# Calculando a DFT em todas as dimensões usando fftn
dftn = fftn(signal)

print("DFT 2D:")
print(dft2)
print()
print("DFT N-Dimensional:")
print(dftn)

print(fft([1, 2, 3]))
print(fft([1, 4, 7]))