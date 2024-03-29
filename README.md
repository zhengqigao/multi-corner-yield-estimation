# Multi-Corner Yield Estimation

## 1. Introduction

Parametric yield estimation over multiple process corners plays an important role in robust circuit design. To guarantee the reliability of a given circuit, we must analyze the parametric yields over all process corners. The number of process corners might remarkably grow at the new technology nodes due to the great complexity of semiconductor manufactory process. Therefore, we propose a novel BI-BD algorithm for multi-corner yield estimation. For more details, please refer to our paper published on IEEE TCAD ([see here](https://ieeexplore.ieee.org/document/8832253)).

## 2. usage

Obviously, MATLAB is required for running the algorithm. All the source files are in the directory named 'src'. Note that the .mat files are the circuit simulation results. Since we simulate the circuit using Hspice provided a commerical Process Design Kit (PDK), we are not able to provide simulation details. Moreover, our focus should be put on the algorithm side, and so these files are also not important :-)

## 3. Propsoed Algorithm (Key Concept)

We introduce a Gaussian distribution as prior. Combining with the likelihood function (in the form of Bernoulli distribution), we obtain a posterior distribuion. The optimal yield estimation results can be extracted by Maximum-a-Posteriori (MAP). However, there is no closed form for the posterior, we adope IRLS method to calculate its mode. Furthermore, the hyper-parameters in the model are efficiently inferenced by Laplacian Approximation and Expectation-Maximization.

## 4. Contact

Please feel free to contact me at <18212020014@fudan.edu.cn> or <zhengqigao@163.com>.