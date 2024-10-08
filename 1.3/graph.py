import matplotlib.pyplot as plt

dots_simple_iter = [(-2, 7), (-4, 11), (-6, 14), (-8, 17), (-10, 20), (-12, 24), (-14, 27), (-16, 30), (-18, 30), (-20, 30)]
dots_seidel = [(-2, 5), (-4, 7), (-6, 9), (-8, 12), (-10, 13), (-12, 15), (-14, 18), (-16, 20), (-18, 20), (-20, 20)]

x_vals_simple_iter = [x for x, y in dots_simple_iter]
y_vals_simple_iter = [y for x, y in dots_simple_iter]

x_vals_seidel = [x for x, y in dots_seidel]
y_vals_seidel = [y for x, y in dots_seidel]

plt.plot(x_vals_simple_iter, y_vals_simple_iter, marker = 'o', linestyle = '-', color = 'b', label = 'Метод простых итераций')

plt.plot(x_vals_seidel, y_vals_seidel, marker = 'o', linestyle = '-', color = 'r', label = 'Метод Зейделя')

plt.xlabel('Ось X, Степень (1eX)')
plt.ylabel('Oсь Y, Количество итераций')
plt.title('График зависимости количества итераций от точности вычислений')

plt.grid(True)

plt.legend()

plt.show()

# -20 30 20
# -18 30 20
# -16 30 20
# -14 27 18
# -12 24 15
# -11 22 15
# -10 20 13
# -9 19 12
# -8 17 12
# -7 15 10
# -6 14 9
# -5 12 8
# -4 11 7
# -3 9 6
# -2 7 5
# -1 5 4
