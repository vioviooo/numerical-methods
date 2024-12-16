import numpy as np
import matplotlib.pyplot as plt

# Define the given function
def f_given(T):
    x = T * 1e-4
    return (247.8955 + 25.3032 * np.log(x) + 0.0008283 * (1.0 / (x * x))
            - 0.22124 * (1.0 / x) + 98.3806 * x 
            - 64.039 * x**2 + 23.304 * x**3)

# Define the predicted function
def F_predicted(T):
    x = T * 1e-4
    return (261.047321 + 30.6373513 * np.log(x) - 0.000281369 * (1.0 / (x * x))
            + 0.015643651 * (1.0 / x) + 63.174977293 * x
            - 24.0775730439 * x**2 + 4.1166014 * x**3)

# Генерируем значения температур от 100K до 1000K
T = np.linspace(100, 1000, 1000)

f_given_values = f_given(T)
F_predicted_values = F_predicted(T)

plt.figure(figsize=(10, 6))
plt.plot(T, f_given_values, label='Ф(Т) из справочника', color='blue', linewidth=2)
plt.plot(T, F_predicted_values, label='Ф(Т) выведенная', color='red', linestyle='--', linewidth=2)

plt.title('Сравнение выведенной формулы Ф(Т) со справочной [газ]', fontsize=16)
plt.xlabel('Температура, K', fontsize=14)
plt.ylabel('Ф(Т), Дж * K^(-1) * моль^(-1)', fontsize=14)
plt.legend(fontsize=12)
plt.grid(True)

# Show the plot
plt.tight_layout()
plt.show()