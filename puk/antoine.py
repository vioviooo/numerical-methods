import numpy as np
import matplotlib.pyplot as plt

# Константы для уравнения Антуана (для воды)
# Действительны для диапазона температур от 1°C до 100°C (приблизительно)
A = 8.07131
B = 1730.63
C = 233.426

# Диапазон температур (°C)
T = np.linspace(1, 100, 500)

# Переводим температуру в Кельвины для оси X
T_K = T + 273.15

# Вычисляем давление с использованием уравнения Антуана (мм рт. ст.)
P = 10 ** (A - B / (T + C))

# Переводим давление в кПа для универсальности (1 мм рт. ст. = 0.133322 кПа)
P_kPa = P * 0.133322

# Строим график с температурой в Кельвинах
plt.figure(figsize=(8, 6))
plt.plot(T_K, P_kPa, label="Кривая насыщения воды")
plt.title("Кривая насыщения воды (уравнение Антуана)", fontsize=14)
plt.xlabel("Температура, K", fontsize=12)
plt.ylabel("Давление насыщения, кПа", fontsize=12)
plt.grid(True)
plt.legend(fontsize=12)
plt.show()
