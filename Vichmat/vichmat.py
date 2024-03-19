import matplotlib.pyplot as plt

x = []  # список для значений по оси x
y11 = []
y12 = []
y13 = []
y21 = []
y22 = []
y23 = []

# Считываем данные из файла
with open('./output1.txt', 'r') as file:
    for line in file:
        data = line.split()
        x.append(float(data[2]))  # данные для оси x
        y11.append(float(data[0]))
        y21.append(float(data[1]))

with open('./output2.txt', 'r') as file:
    for line in file:
        data = line.split()
        y12.append(float(data[0]))
        y22.append(float(data[1]))

with open('./output3.txt', 'r') as file:
    for line in file:
        data = line.split()
        y13.append(float(data[0]))
        y23.append(float(data[1]))

plt.plot(x, y11, label='Эйлера 1')
plt.plot(x, y12, label='Рунге 1')
plt.plot(x, y13, label='Предиктор корректор 1')

plt.plot(x, y21, label='Эйлера 2')
plt.plot(x, y22, label='Рунге 2')
plt.plot(x, y23, label='Предиктор корректор 2')

plt.xlabel('Значения оси x')
plt.ylabel('Значения оси y')
plt.title('Графики из файла')
plt.legend()

plt.show()
