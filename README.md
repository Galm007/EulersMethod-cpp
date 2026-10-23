# Eulers Method

A homework submission for my Differential Equations class to approximate the solution to an ODE-IVP using Euler's Method, plotted using [Mathplot++](https://github.com/alandefreitas/matplotplusplus).

<img width="560" height="420" alt="plots" src="https://github.com/user-attachments/assets/7edadf65-0ced-48fb-b5fc-fa983fe1fc11" />

## Error tables
ODE-IVP 1: dx/dt = t + x , x(0) = 1 , 0 <= t <= 1
| h         | Value at t = 1 | Error From Exact Solution
|-----------|----------------|----------------------------
| h = 0.1   | 3.18748 | 0.249079   |
| h = 0.01  | 3.40963 | 0.0269358  |
| h = 0.001 | 3.43385 | 0.00271562 |
| Exact     | 3.43656 | 0.0        |

ODE-IVP 2: dx/dt = (x/3)(8-x) , x(0) = 1 , 0 <= t <= 5
| n         | Value at t = 5 | Error From Exact Solution
|-----------|----------------|----------------------------
n = 5       | 5.30429 | 2.69562    |
n = 10      | 8.00024 | -0.0003    |
n = 20      | 8.00000 | -0.0001    |
| Exact     | 7.99991 | 0.0        |

## How to run (Linux):
```bash
# Clone the repository
git clone https://github.com/Galm007/EulersMethod-cpp && cd EulersMethod-cpp

# Compile the program
g++ main.cpp -lm -lmatplot -o eulers_method

# Run the program (Generates plots.svg and ode-ivp-1-table.txt)
./eulers_method
```
