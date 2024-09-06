#ПРИЛОЖЕНИЕ 1. КОД РЕАЛИЗАЦИИ СОЛИТОНОВ
import numpy as np
import matplotlib.pyplot as plt
import math
# Задание параметров
alfa = 6
beta = 1
h = 0.2
tau = 0.000174
T = 0.0522
L=40
Q = 4
x0 = 10
delta = math.sqrt(4*beta/Q)
A = 3*Q/alfa
t = 0
N_x = int(L / h)
print('N_x = M = ',N_x)
print('A = ',A)
print('delta = ', delta)
#задаем шаг движения волны(время)
shag = 0.5
N_t = int(2/shag+1)
print('N_t = K =',N_t)
# Создание сетки значений функции
u = np.zeros((N_t, N_x), dtype=np.float64)
t = -1
for j in range(N_t):
    t = t + 1
    for i in range(N_x):
        u[j][i] = A/(np.cosh((h*(i)-x0-Q*t*shag)/delta)**2)
plt.axis([5,20,-0.5,2.5]) #задаем интервалы для вывода графика
plt.grid(True)
# Построение графика численного решения
x=np.zeros(N_x)
for i in range(0, N_x-1):
    x[i] = i*h
x[N_x-1]=L
plt.plot(x,u[0])
plt.xlabel('x')
plt.ylabel('u(x, t)')
plt.title('График солитона, t=0')
plt.show()
Q1 = 4
x1 = 2
Q2 = 2
x2 = 10
A1 = 3*Q1/alfa
A2 = 3*Q2/alfa
delta1 = math.sqrt(4*beta/Q1)
delta2 = math.sqrt(4*beta/Q2)
#задаем шаг движения волны(время)
shag = 0.5
N_t = int(10/shag+1)
print('N_t = K =',N_t)
# Создание сетки значений функции
u = np.zeros((N_t, N_x), dtype=np.float64)
t = -1
for j in range(N_t):
    t = t + 1
    for i in range(N_x):
        u[j][i] = A1/(np.cosh((h*(i)-x1-Q1*t*shag)/delta1)**2)+A2/(np.cosh((h*(i)-x2-Q2*t*shag)/delta2)**2)
plt.axis([0,40,-0.5,3]) #задаем интервалы для вывода графика
plt.grid(True)
# Построение графика численного решения
x=np.zeros(N_x)
for i in range(0, N_x-1):
    x[i] = i*h
x[N_x-1]=L
plt.plot(x,u[11])
plt.xlabel('x')
plt.ylabel('u(x, t)')
plt.title('Движение двух солитонов, t=5.5')
plt.show()

#ПРИЛОЖЕНИЕ 2. КОД СХЕМЫ №1
import numpy as np
import matplotlib.pyplot as plt
import math
# Задание параметров
alfa = 6
beta = 1
h = 0.2
tau = 0.000174
T = 0.0522
L = 19.6
Q = 4
x0 = 10
delta = math.sqrt(4*beta/Q)
A = 3*Q/alfa
t = 0
# Вычисление количества узлов по пространству и времени
N_x = int(L / h)
N_t = int(T / tau)
print('N_t = K =',N_t)
print('N_x = M = ',N_x)
# Создание сетки значений функции
u = np.zeros((N_t, N_x+4), dtype=np.float64)
for i in range(2, N_x+2):
    u[0][i] = A/(np.cosh((h*(i-2)-x0-Q*t)/delta)**2)
# Задаем граничное условие (x=0)
for i in range(N_t-1):
    u[i][2] = 0
# Задаем граничное условие (x=L)
for i in range(N_t-1):
    u[i][N_x+1] = 0
# Условие периодичности
u[0][0]=u[0][N_x]
u[0][1]=u[0][N_x+1]
u[0][N_x+2]=u[0][2]
u[0][N_x+3]=u[0][3]
# матрица Численного решения уравнения Кортевега де Фриза
reshenie=np.zeros((N_t, N_x), dtype=float)
#Разностная схема
for n in range(N_t-1):
    for i in range(3, N_x+1):
        u[n+1][i]= u[n][i]-tau/h*u[n][i]*(u[n][i+1]-u[n][i-1])-beta*tau/h**3*(u[n][i+2]-2*u[n][i+1]+2*u[n][i-1]-u[n][i-2])
    #Периодичность
    u[n+1][N_x+2]=u[n+1][2]
    u[n+1][N_x+3]=u[n+1][3]
    u[n+1][0]=u[n+1][N_x]
    u[n+1][1]=u[n+1][N_x+1]
print('u ', u)
for i in range(N_t):
    for j in range(N_x):
        reshenie[i][j]=u[i][j+2]
print('Матрица решений для всех слоёв\n',reshenie)
# Построение графика численного решения
x=np.zeros(N_x)
for i in range(0, N_x-1):
    x[i] = i*h
x[N_x-1]=L
plt.plot(x,reshenie[299])
print('Выведем график решения')
plt.grid(True)
plt.xlabel('x')
plt.ylabel('u(x, t)')
plt.title('Численное решение KdV на 300 временном слое')
plt.show()
# Создание сетки точных значений функции
U = np.zeros((N_t, N_x), dtype=float)
for n in range(N_t):
    for i in range(N_x):
        U[n][i] = A/(np.cosh((h*(i)-x0-Q*t)/delta)**2)
tmp = np.zeros(N_x)
tmp=U[20]-reshenie[20]
err = np.zeros(N_x)
sum = 0
for i in range (N_x):
    err[i]=(tmp[i])**2
    sum = sum + err[i]
error=np.sqrt(sum)
print(error)
print("max", max(tmp))

#ПРИЛОЖЕНИЕ 3. КОД СХЕМЫ №2
import numpy as np
import matplotlib.pyplot as plt
import math
# Задание параметров
alfa = 6
beta = 1
h = 0.2
tau = 0.000174
T = 0.0522
L = 19.6
Q = 4
x0 = 10
delta = math.sqrt(4*beta/Q)
A = 3*Q/alfa
t = 0
# Вычисление количества узлов по пространству и времени
N_x = int(L / h)
N_t = int(T / tau)
print('N_t = K =',N_t)
print('N_x = M = ',N_x)
# Создание сетки значений функции
u = np.zeros((N_t, N_x+6), dtype=np.float64)
for i in range(3, N_x+3):
    u[0][i] = A/(np.cosh((h*(i-3)-x0-Q*t)/delta)**2)
# Задаем граничное условие (x=0)
for i in range(N_t-1):
    u[i][2] = 0
# Задаем граничное условие (x=L)
for i in range(N_t-1):
    u[i][N_x+1] = 0
# Условие периодичности
u[0][0]=u[0][N_x]
u[0][1]=u[0][N_x+1]
u[0][2]=u[0][N_x+2]
u[0][N_x+3]=u[0][3]
u[0][N_x+4]=u[0][4]
u[0][N_x+5]=u[0][5]
# матрица Численного решения уравнения Кортевега де Фриза
reshenie=np.zeros((N_t, N_x), dtype=float)
#Разностная схема
for n in range(N_t-1):
    for i in range(4, N_x+2):
        u[n+1][i]= u[n][i]-tau/(12*h)*u[n][i]*(u[n][i+2]-8*u[n][i+1]+8*u[n][i-1]-u[n][i-2])-beta*tau/(8*h**3)*(u[n][i+3]-8*u[n][i+2]+13*u[n][i+1]-13*u[n][i-1]+8*u[n][i-2]-u[n][i-3])
    #Периодичность
    u[n+1][N_x+5]=u[n+1][5]
    u[n+1][N_x+4]=u[n+1][4]
    u[n+1][N_x+3]=u[n+1][3]
    u[n+1][0]=u[n+1][N_x]
    u[n+1][1]=u[n+1][N_x+1]
    u[n+1][2]=u[n+1][N_x+2]
print('u ', u)
for i in range(N_t):
    for j in range(N_x-1):
        reshenie[i][j]=u[i][j+3]
print('Матрица решений для всех слоёв\n',reshenie)
U = np.zeros((N_t, N_x), dtype=float)
# Построение графика численного решения
x=np.zeros(N_x)
for i in range(0, N_x-1):
    x[i] = i*h
x[N_x-1]=L
plt.plot(x,reshenie[218])
print('Выведем график решения')
plt.xlabel('x')
plt.ylabel('u(x, t)')
plt.title('Численное решение KdV на 218 временном слое')
plt.show()
for n in range(N_t):
    for i in range(N_x):
        U[n][i] = A/(np.cosh((h*(i)-x0-Q*t)/delta)**2)
tmp = np.zeros(N_x)
tmp=U[20]-reshenie[0]
err = np.zeros(N_x)
sum = 0
for i in range (N_x):
    err[i]=(tmp[i])**2
    sum = sum + err[i]
error=np.sqrt(sum)
print(error)
print("max", max(tmp))

#ПРИЛОЖЕНИЕ 4. КОД СХЕМЫ №3
import numpy as np
import matplotlib.pyplot as plt
import math
# Задание параметров
alfa = 6
beta = 1
h = 0.2
tau = 0.000174
T = 0.0522
L = 19.6
Q = 4
x0 = 10
delta = math.sqrt(4*beta/Q)
A = 3*Q/alfa
t = 0
# Вычисление количества узлов по пространству и времени
N_x = int(L / h)
N_t = int(T / tau)
print('N_t = K =',N_t)
print('N_x = M = ',N_x)
# Создание сетки значений функции
u = np.zeros((N_t, N_x+4), dtype=np.float64)
for i in range(2, N_x+2):
    u[0][i] = A/(np.cosh((h*(i-2)-x0-Q*t)/delta)**2)
# Задаем граничное условие (x=0)
for i in range(N_t-1):
    u[i][2] = 0
# Задаем граничное условие (x=L)
for i in range(N_t-1):
    u[i][N_x+1] = 0
# Условие периодичности
u[0][0]=u[0][N_x]
u[0][1]=u[0][N_x+1]
u[0][N_x+2]=u[0][2]
u[0][N_x+3]=u[0][3]
# матрица Численного решения уравнения Кортевега де Фриза
reshenie=np.zeros((N_t, N_x), dtype=float)
#Разностная схема
for n in range(N_t-1):
    for i in range(2, N_x+2):
        u[n+1][i]= u[n][i]-tau/h*1/3*(u[n][i+1]+u[n][i]+u[n][i-1])*(u[n][i+1]-u[n][i-1])-beta*tau/h**3*(u[n][i+2]-2*u[n][i+1]+2*u[n][i-1]-u[n][i-2])
    #Периодичность
    u[n+1][N_x+2]=u[n+1][2]
    u[n+1][N_x+3]=u[n+1][3]
    u[n+1][0]=u[n+1][N_x]
    u[n+1][1]=u[n+1][N_x+1]
print('u ', u)
for i in range(N_t):
    for j in range(N_x):
        reshenie[i][j]=u[i][j+2]
print('Матрица решений для всех слоёв\n',reshenie)
# Построение графика численного решения
x=np.zeros(N_x)
for i in range(0, N_x-1):
    x[i] = i*h
x[N_x-1]=L
plt.plot(x,reshenie[115])
print('Выведем график решения')
plt.xlabel('x')
plt.ylabel('u(x, t)')
plt.title('Численное решение KdV на 115 временном слое')
plt.show()
# Создание сетки точных значений функции
U = np.zeros((N_t, N_x), dtype=float)
for n in range(N_t):
    for i in range(N_x):
        U[n][i] = A/(np.cosh((h*(i)-x0-Q*t)/delta)**2)
tmp = np.zeros(N_x)
tmp=U[0]-reshenie[20]
err = np.zeros(N_x)
sum = 0
for i in range (N_x):
    err[i]=(tmp[i])**2
    sum = sum + err[i]
error=np.sqrt(sum)
print(error)
print("max", max(tmp))

#ПРИЛОЖЕНИЕ 5. КОД СХЕМЫ №4
import numpy as np
import matplotlib.pyplot as plt
import math
# Задание параметров
alfa = 6
beta = 1
h = 0.2
tau = 0.000174
T = 0.0522
L = 19.6
Q = 4
x0 = 10
delta = math.sqrt(4*beta/Q)
A = 3*Q/alfa
t = 0
# Вычисление количества узлов по пространству и времени
N_x = int(L / h)
N_t = int(T / tau)
print('N_t = ',N_t)
print('N_x = ',N_x)
print(A)
print(delta)
# Создание сетки значений функции
u = np.zeros((N_t, N_x+4), dtype=np.float64)
# Задаем начальное условие (t=0)
for i in range(2, N_x+1):
    u[0][i] = A/(np.cosh((h*(i-2)-x0)/delta)**2)
# Задаем граничное условие (x=0)
for i in range(N_t-1):
    u[i][2] = 0
# Задаем граничное условие (x=L)
for i in range(N_t-1):
    u[i][N_x+1] = 0
# Условие периодичности
u[0][0]=u[0][N_x]
u[0][1]=u[0][N_x+1]
u[0][N_x+2]=u[0][2]
u[0][N_x+3]=u[0][3]
print(u)
# матрица Численного решения уравнения Кортевега де Фриза
reshenie=np.zeros((N_t, N_x), dtype=float)
#Разностная схема
for n in range(N_t-1):
    for i in range(3, N_x):
        u[n+1][i]= u[n][i]-alfa*tau/h*(((u[n][i+1]+u[n][i])**2)/8 - ((u[n][i-1]+u[n][i])**2)/8)- beta*tau/(2*(h**3))*(u[n][i+2]-2*u[n][i+1]+2*u[n][i-1]-u[n][i-2])
#Периодичность
    u[n+1][N_x+2]=u[n+1][2]
    u[n+1][N_x+3]=u[n+1][3]
    u[n+1][0]=u[n+1][N_x]
    u[n+1][1]=u[n+1][N_x+1]
print(u)
for i in range(N_t):
    for j in range(N_x):
        reshenie[i][j]=u[i][j+2]
print('Матрица решений для всех слоёв\n',reshenie)
# Построение графика численного решения
x=np.zeros(N_x)
for i in range(0, N_x-1):
    x[i] = i*L/(N_x)
x[N_x-1]=L
plt.plot(x,reshenie[299])
print('Выведем график решения')
plt.xlabel('x')
plt.ylabel('u(x, t)')
plt.title('Численное решение KdV на 300 временном слое')
plt.show()
# Создание сетки точных значений функции
U = np.zeros((N_t, N_x), dtype=float)
for n in range(N_t):
    for i in range(N_x):
        U[n][i] = A/(np.cosh((h*(i)-x0-Q*t)/delta)**2)
tmp = np.zeros(N_x)
tmp=U[20]-reshenie[20]
err = np.zeros(N_x)
sum = 0
for i in range (N_x):
    err[i]=(tmp[i])**2
    sum = sum + err[i]
error=np.sqrt(sum)
print(error)
print("max", max(tmp))