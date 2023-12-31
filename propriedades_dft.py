import numpy as np
from scipy.fft import fft, ifft, fftshift

def prop_simetria(): # função para demonstrar a propriedade da simetria
    x = np.array([1, 2, 3, 4, 5])

    # Cálculo da DFT
    X = fft(x, 10)
    X = fftshift(X)

    print("X[k] = X[-k]:", X)

def prop_linearidade(): # função para demonstrar a propriedade de linearidade
    # Sinais de entrada
    x1 = np.array([1, 2, 3, 4, 5])
    x2 = np.array([6, 7, 8, 9, 10])

    # Constantes de ponderação
    a = 2
    b = 3

    # Cálculo da DFT
    X1 = fft(a*x1 + b*x2)
    X2 = a * fft(x1) + b * fft(x2)

    print("DFT(a * x1[n] + b * x2[n]):")
    print(fftshift(X1))
    print()
    print("a * DFT(x1[n]) + b * DFT(x2[n]):")
    print(fftshift(X2))

prop_linearidade()
prop_simetria()
