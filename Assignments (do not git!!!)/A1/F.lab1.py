#Q1(b)

import numpy as np
data = np.loadtxt('cdata.txt') #loading in data
print('Q1b')
#defining equation 1
def std1(data):
    x_avg = 1/data.size * np.sum(data) #mean
    return np.sqrt((1 / (np.size(data) - 1) * np.sum(data - x_avg)**2))

#defining equation 2
def std2(data):
    x_avg = 1/data.size * np.sum(data) #mean
    if (1 / (np.size(data) - 1)) * np.sum(np.square(data) - data.size * x_avg**2) < 0:
        print("2 pass negative value error")
        return np.nan
    else:
        return np.sqrt((1 / (np.size(data) - 1)) * (np.sum(data)**2 - data.size * x_avg**2))

#mod from 1d
def newstd2(data):
    n = 0
    mean = 0
    x2 = 0  # sum of squares of x_i - xbar placeholder

    for x in data:
        n += 1
        xd = x - mean
        mean += xd / n
        delta2 = x - mean
        x2 += xd * delta2

    return np.sqrt(x2 / (n - 1))

# calculating relative error with "correct" std
def rel_err(x,data):
    return (x-np.std(data,ddof=1))/np.std(data,ddof=1)

print(rel_err(std1(data),data), rel_err(std2(data),data),rel_err(newstd2(data),data))

#Q1c
print('Q1c')

#generating the 2 sequences
seq1 = np.random.normal(0,1,2000)
seq2 = np.random.normal(1.e-7,1,2000)

print(rel_err(std1(seq1),seq1), rel_err(std2(seq1),seq1),rel_err(newstd2(seq1),seq1))
print(rel_err(std1(seq2),seq2), rel_err(std2(seq2),seq2),rel_err(newstd2(seq2),seq2))


#Q2a
import matplotlib.pyplot as plt
#defining p(u) and q(u)
def func_p(u):
    return (1-u)**8

def func_q(u):
    return  1 - 8*u + 28*u**2 - 56*u**3 + 70*u**4 - 56*u**5 + 28*u**6 - 8*u**7 + u**8

u = np.linspace(0.98,1.02,500)

#plotting p(u) and q(u)
plt.plot(u,func_q(u),'g',label='q(u)')
plt.plot(u,func_p(u),'r',label='p(u)')
plt.title('Q2a plot of q(u) and p(u)')
plt.legend()
plt.savefig('Q2a.png')
plt.show()


#Q2b
print('Q2b')

#plotting difference of q(u) and p(u)
plt.plot(u,func_p(u)-func_q(u),'g')
plt.title('Q2b plot of p(u) - q(u)')
plt.savefig('Q2b.1.png')
plt.show()

#plotting difference of q(u) and p(u) as histogram
plt.hist(func_p(u) - func_q(u),bins=50,edgecolor='black')
plt.title('Q2b histogram of p(u) - q(u) values')
plt.savefig('Q2b.2.png')
plt.show()


# standard deviation of p(u)- q(u)
print(np.std(func_p(u)-func_q(u)))
#defining  eqn3 to calculate standard deviation
def eqn3 (data):
    x2_bar = np.mean(np.square(data))
    return (10e-16) * np.sqrt(data.size) * np.sqrt(x2_bar)
print(eqn3(func_p(u)-func_q(u)))

#Q2c
print('Q2c')

u_new = np.linspace(0.98,1.00,500)
def eqn4(data):
    x2_bar = np.mean(np.square(data))
    x_bar = np.mean(data)
    return (10e-16)/np.sqrt(data.size) * np.sqrt(x2_bar)/x_bar
print(eqn4(func_p(u_new)-func_q(u_new)))

def func_rel(p,q):
    return np.abs(p-q)/np.abs(p)

u_test = np.linspace(0.98,0.988,500)
plt.plot(u_test,func_rel(func_p(u_test),func_q(u_test)),'g')
plt.title('Q2c plot of |p(u)-q(u)| / |p(u)|')
plt.savefig('Q2c.png')
plt.show()


#Q2(d)
print('Q2d')

#defining function f
def func_f(data):
    return u**8/((u**4)*(u**4))

print(np.std(func_f(u)))

plt.plot(u,func_f(u)-1,'g')
plt.title('Q2d plot of f(u)-1')
plt.savefig('Q2d.png')
plt.show()



#Q3b
print('Q3b')

def f(x):
    return 4/(1+x**2)

def trapezoidal_rule(a, b, n): #def trapezoid
    h = (b - a) / n
    x = np.linspace(a, b, n+1)
    y = f(x)
    integral = (h/2) * (y[0] + 2 * np.sum(y[1:n]) + y[n])
    return integral

def simpson_rule(a,b,n): #def simpson
    h = (b - a) / n
    x = np.linspace(a, b, n + 1)
    y = f(x)
    integral = (h / 3) * (y[0] + 4 * np.sum(y[1:n:2]) + 2 * np.sum(y[2:n - 1:2]) + y[n])

    return integral

a = 0 #def intervals and Number of subdivisions
b = 1
n = 4

print(trapezoidal_rule(a,b,n))
print(simpson_rule(a,b,n))

#Q3
import time
print('Q3c')

#trial and error for trap rule
n=4096
trap_t_start = time.time()
trap_diff = trapezoidal_rule(a,b,n)-np.pi
trap_t_end = time.time()
print(trap_diff,trap_t_end-trap_t_start)
#trial and error for simp rule
n=16
simp_t_start = time.time()
simp_diff = simpson_rule(a,b,n)-np.pi
simp_t_end = time.time()
print(simp_diff,simp_t_end-simp_t_start)

#Q3d
print('Q3d')

#from eqn 5.28
def err(I1,I2):
    return 1/3 * (I2-I1)
n1=16
n2=32
print(err(trapezoidal_rule(a,b,n1),trapezoidal_rule(a,b,n2)))
