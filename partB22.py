'''
Alice Aidlin 206448326
Maayan Nadivi 208207068
Bar Sela 206902355
Bar Weizman 206492449
#part B
#Q.14
'''

import datetime
import sympy as sp
import math
from sympy.utilities.lambdify import lambdify

eps = 0.00001

def Calculate_derivative(f):
    f_prime = f.diff(x)
    return f_prime

x = sp.symbols('x')


def error(start, end):
    upup = eps / (end - start)
    up = (-1) * math.log(upup, math.e)
    return (up // math.log(2, math.e)) + 1

def Newton_Raphson(func, func_tag, start, end):
    xr= round((start+end)/2,5)
    f = lambdify(x, func)
    fn = lambdify(x, func_tag)
    xr_1=round(xr-(f(xr)/fn(xr)),5)
    counter=1
    print("Xr=", round(xr,5), " f(x)=", round(f(xr),5), " f'(x)=", round(fn(xr),5))
    while (abs(xr-xr_1)>eps):
        print("Xr=", round(xr_1,5), " f(x)=", round(f(xr_1),5), " f'(x)=", round(fn(xr_1),5))
        if counter > 100:
            return "The function does not converge"
        xr = round(xr_1,5)
        xr_1 = round(xr - (f(xr) / fn(xr)),5)
        counter+=1

    print("num of iteration {}".format(counter))
    return round(xr_1,5)

def secant_method(func, start, end):
    xr_minus1 = round(start,5)
    xr = end
    f = lambdify(x, func)
    xr_1 = round((xr_minus1*f(xr)-xr*f(xr_minus1))/(f(xr)-f(xr_minus1)),5)
    counter = 1
    print("Xr=", round(xr,5), "Xr+1", round(xr_1,5), " f(x)=", round(f(xr),5))
    while (abs(xr - xr_1) > eps):
        print("Xr=", round(xr,5), "Xr+1", round(xr_1,5), " f(x)=", round(f(xr),5))
        if counter > 100:
            return "The function does not converge"
        xr_minus1=round(xr,5)
        xr=round(xr_1,5)
        xr_1=round((xr_minus1*f(xr)-xr*f(xr_minus1))/(f(xr)-f(xr_minus1)),5)
        counter += 1

    print("num of iterition {}".format(counter))
    return round(xr_1,5)

def help(func, start, end,n):
    f = lambdify(x, func)
    func_tag = Calculate_derivative(func)
    fn = lambdify(x, func_tag)
    next = start + 0.1
    next = round(next, 2)
    while round(start,2) < end:
        if f(start) * f(next) < 0:
            if n == 1:
                print("x = {}00000131550".format(Newton_Raphson(func,func_tag, start, next)))
                print("")
            if n == 2:
                print("x = {}00000131550".format(secant_method(func, start, next)))
                print("")
        if f(start) == 0:
            print("num of iteration 0")
            print("x = {}00000131550".format(float(start)))
            print("")

        start += 0.1
        start=round(start, 2)
        next += 0.1
        next=round(next, 2)
    if f(end)==0:
        print("num of iteration 0")
        print("x = {}00000131550".format((end)))
        print("")


'''Integral'''

def I(func,a,b):
    f = lambdify(x, func)
    return 0.5*(b-a)*(f(a)+f(b))

def TrapezoidalRule(func, a, b, n):
    h = float(b - a) / n
    sum=0.0
    for i in range(n):
        sum+=I(func,a,a+h)
        a+=h
    return sum


def SimpsonRule(func, a, b, n):
    f = lambdify(x, func)
    if n % 2 != 0:
        return "n must be even"
    h = (b - a) / n
    integral = round((f(a) + f(b)),3)
    print("integral=f(a) + f(b)=",round(integral,3))
    for i in range(1,n):
        r = a + i * h
        if i % 2 == 0:
            integral += 2 * round(f(r),3)
            print("integral+= 2 *", round(f(r),3), "=", round(integral,3))
        else:
            integral += 4 * round(f(r),3)
            print("integral+= 4 *", round(f(r),3), "=", round(integral,3))
    integral *= (h/3)
    print("integral=", round(integral,3))
    return round(integral,3)



def RombergsMethod(func, a, b, n):
    mat = list(range(n))  # make it list
    for i in range(n):
        mat[i] = list(range(n))
    for k in range(n):
        mat[k][0] = round(TrapezoidalRule(func, a, b, 2 ** k),3)
        for j in range(0, k):
            mat[k][j + 1] = (4 ** (j + 1) * round(mat[k][j],3) - round(mat[k - 1][j],3)) / (4 ** (j + 1) - 1)

            print("R[",k,"][",j + 1,"] = ",round(mat[k][j + 1],3))
    return round(mat[n-1][n-1],3)



'''main'''
def main():

    f = (x*math.e**(-x**2+5*x))*(2*x**2-3*x-5)

    start = 0
    end = 3
    print("\n````````````````````````Newton Raphson````````````````````````\n")
    help(f, start, end, 1)
    print("\n````````````````````````Secant Method````````````````````````\n")
    help(f, start, end, 2)
    print("\n````````````````````````INTEGRAL````````````````````````\n")
    print("\n````````````````````````Simpson Rule, Integral Value in Countries[0,1]````````````````````````\n")
    print("S={}00000131550".format(round(SimpsonRule(f,0.5,1,4),3)))
    print("\n````````````````````````Rombergs Method, Integral Value in Countries[0,1]````````````````````````\n")
    print("S={}00000131550".format(round(RombergsMethod(f,0.5,1,4),3)))

main()
