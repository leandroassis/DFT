#### 1. A Transformada de Fourier no Tempo Discreto: Características e Propriedades

A Transformada de Fourier no Tempo Discreto (DTFT) permite analisar as características espectrais de sinais discretos no domínio da frequência. Dentre essas características, podemos citar a amplitude, fase, resposta em frequência e espectro de frequência de um sinal discreto.

![Sinal Discreto $x[n] = 0.8^{n}u[n]$](./sinaldiscreto.png)
![Diagrama de Módulo (a) e Fase (b) da DTFT do sinal discreto $x[n] = 0.8^{n}u[n]$](./espectograma.png)

Matematicamente, a DTFT é o somatório ponderado de exponenciais complexas, onde cada exponencial completa representa a fase e módulo para uma frequência específica.

A fórmula geral para a DTFT de um sinal discreto x[n] é:

$$ X(\Omega) = \sum_{n = -\infty}^{\infty} (x[n] \cdot e^{-j\Omega n})$$

Em outras palavras, a Transformada de Fourier no Tempo Discreto é a transformação linear de um sinal discreto para sua projeção em uma base de exponenciais complexas.

A análise de Fourier em tempo discreto pode ser dividida em três tipos:
* Série Discreta de Fourier (DFS): É a DTFT de forma exclusiva para sinais periódicos. Nesse caso particular, o somatório ponderado de exponenciais complexas ocorre apenas para o número de amostras dentro de um período $N_{0}$ do sinal discreto x[n] (0 a $N_{0}-1$). A saída representa o módulo e desvio de fase para cada harmônico (múltiplo inteiro) da frequência $\omega$ do sinal de entrada.
  [Adicionar exemplo DFS]
  
* Transformada de Fourier de Tempo Discreto (DTFT): É uma transformação contínua que mapeia um sinal discreto x[n] para um espectro contínuo na frequência X($\Omega$). Nesse caso, a DTFT é uma função contínua de frequência (definida para qualquer valor de $\Omega$), representando a amplitude e fase para cada frequência contínua do sinal de entrada. A DTFT é definida para sinais que se estendem de -$\infin$ a $\infin$ no tempo discreto.
  [Adicionar exemplo DTFT]

* Transformada de Fourier Discreta (DFT): Se trata de uma versão discretizada da DTFT. A DFT é aplicada a sequências discretas de comprimento finito (o somatório ocorre de 0 a $N - 1$, onde N é o número de amostras na sequência) e produz um espectro de frequência discreto e periódico. Para o cálculo da DFT é utilizado o algoritmo FFT (Fast Fourier Transform), alvo de análise deste trabalho. A DFT é a transformada utilizada para a determinação de espectros de frequência de forma computacional, uma vez que estes não conseguem lidar com um número infinito de valores. Note que se o número de amostras for suficientemente grande, pode-se aproximar com a precisão necessária a DTFT com a DFT.

  [Adicionar exemplo DFT]

A DFT possui as seguintes propriedades:
1. Linearidade;
    [Adicionar demonstração linearidade]
2. Distributividade;
    [Adicionar demonstração distributividade]
3. Simetria conjugada complexa para sinais reais;
   [Adicionar demonstração simetria conjugada complexa]
4. Deslocamento no tempo;
   [Adicionar demonstração deslocamento no tempo]
5. Convolução no tempo;
   [Adicionar demonstração convolução no tempo]

#### 2. A DFT no Python

Dentre as implementações da DFT disponíveis na biblioteca scipy, destacam-se as funções fft e stft, que implementam a DFT e a STFT (Short-Time Fourier Transform), respectivamente.

* scipy.fft.fft ou numpy.fft.fft: Calcula a DFT de um sinal discreto. A função fft retorna o espectro de frequência discreto e periódico de um sinal discreto. A função fft é uma implementação do algoritmo FFT (Fast Fourier Transform), que é um algoritmo eficiente para o cálculo da DFT. A função fft é utilizada para o cálculo da DFT de um sinal discreto, enquanto a função fftshift é utilizada para o deslocamento do espectro de frequência para o centro do gráfico.
    O algoritmo FFT consiste no seguinte
    [Adicionar demonstração FFT Wikipedia]

    Em particular, a função fft do scipy.fft, bem como o algoritmo FFT, são extremamente aprimorados para tamanhos de dados múltiplos de 2.

    A biblioteca scipy.fft possui outras funções para a determinação da DFT para dados de N dimensões usando o FFT.

    O algoritmo FFT também pode ser utilizado para realizar a convolução de sinais de forma eficiente por conta da propriedade de convolução no tempo (fftconvolve na biblioteca scipy.signal).

    De forma geral, as funções fft possuem as mesmas características no Scipy e no Numpy. A diferença entre as funções fft do Scipy e do Numpy é que a função fft do Scipy possui um parâmetro opcional para a normalização do sinal de saída, enquanto a função fft do Numpy não possui esse parâmetro. Além disso, existe uma diferença na performance apesar de ambas as implementações serem baseadas no algoritmo FFT demonstrado.

    A biblioteca Numpy possui outras funções para a determinação da DFT caso a entrada seja puramente real (utilizando a propriedade de simetria conjugada complexa para sinais reais). na prática essa implementação não difere da FFT, apenas não computa os valores negativos do espectro de frequência.

* scipy.signal.stft: Calcula a STFT de um sinal de entrada. A STFT trata-se da DFT calculada em uma janela deslizante de amostras de uma sequência. Dessa forma, aplicando a STFT em um sinal é possível analisar as mudanças espectrais no sinal ao longo do tempo.

Ambas as funções fft e stft possuem suas inversas que apenas utilizam a equação inversa de Fourier para a reconstrução do sinal de entrada.

$$x[n] = \frac{1}{N}\sum_{k = 0}^{N - 1} X(\Omega) \cdot e^{+j\Omega n}$$
Como a DFT (e consequentemente a STFT) de um sinal é uma representação discretizada da resposta em frequência do sinal x[n], então $\Omega_{k} = \frac{2\pi k}{N}$. Onde k são os índices de frequência (harmônicos) e N é o número de amostras do espectro X[$\Omega$].

#### 3. Cálculo da DFT e análise de um espectro de frequência utilizando Python

[Gerar espectograma de um sinal de áudio]
[Explicar leitura do espectograma]

[DFT usando scipy.fft.fft]
[DFT usando numpy.fft.fft]
[DFT usando numpy.fft.rfft]
[DFT usando scipy.signal.stft]
[Algoritmo DFT baseado na fórmula]

[Comparar performance das funções]
[Comparar saída das funções]

#### 4. Bibliografia
1. Sinais e Sistemas Lineares, 2ª Edição, B. P. Lathi, LTC, 2007.
2. [Transformada Rápida de Fourier](https://pt.wikipedia.org/wiki/Transformada_rápida_de_Fourier)
3. [Documentação scipy.fft](https://docs.scipy.org/doc/scipy/reference/fft.html)
4. [Documentação scipy.signal](https://docs.scipy.org/doc/scipy/reference/signal.html)
5. [Documentação numpy.fft](https://numpy.org/doc/stable/reference/routines.fft.html)