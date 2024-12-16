import numpy as np
import matplotlib.pyplot as plt

def F0T_water(T):  # Расчет приведенной энергии Гиббса для воды при температуре Т
    x = T
    return -61.53274 + 17.64285 * np.log(x) - 823.02681312 * (1.0 / (x * x)) + 351.948940 * (1.0 / x) \
           - 0.205265594 * x + 0.0006830665 * x * x - 0.0000005178652 * x * x * x

# Диапазон температур
T_range = np.linspace(5, 600, 100)  # Используем 0.1 вместо 0, чтобы избежать log(0)

# Вычисляем значения F0T_water для диапазона температур
F0T_values = F0T_water(T_range)

# Пользовательские точки (пример)
user_points = {
    5: 3.408,
    100: 8.745,
    200: 16.336,
    300: 25.708,
    400: 39.729,
    500: 52.041,
    600: 63.111
}

# Преобразуем пользовательские точки в массивы для удобства
user_T = np.array(list(user_points.keys()))
user_F = np.array(list(user_points.values()))

# Построение графиков
plt.figure(figsize=(10, 6))

# График F0T_water
plt.plot(T_range, F0T_values, label="Выведенная Ф(Т)", color="blue")

# График пользовательских точек
plt.scatter(user_T, user_F, label="Данные из справочника", color="red", zorder=5)

# Настройка графика
plt.title("Сравнение выведенной формулы Ф(Т) со справочными данными [жидкость]", fontsize=14)
plt.xlabel("Температура, K", fontsize=12)
plt.ylabel("Ф(Т), Дж * K^(-1) * моль^(-1)", fontsize=12)
plt.legend(fontsize=12)
plt.grid()

# Отображение графика
plt.show()
