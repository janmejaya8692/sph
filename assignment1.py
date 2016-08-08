# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 23:54:48 2016

@author: janme
"""
from matplotlib import pyplot as plt
import numpy as np
import math


def initialize():
    n = [50, 100, 150, 200, 300]
#    num = raw_input("Enter how many no. of points you want:")
#    print 'Enter numbers in array: '
#    for i in range(int(num)):
#        z = raw_input("num :")
#        n.append(int(z))
    X_min = -1.0
    X_max = 1.0
    dx = []
    for i in range(len(n)):
        dx.append((X_max-X_min)/n[i])
    h = np.array([0.7, 0.8, 0.9, 1.0, 1.1, 1.3, 1.5, 1.8, 2.0, 2.2])*dx[0]
    f = raw_input("Enter G for gaussian kernel and C for cubicspline kernel")
    return h, dx, n, f, X_min, X_max


def kernel_decis(g, f, x, y, h, n):
        if(f == 'C'):
            if(g == 4):
                W = cubicspline(x, y, h, n)
            else:
                W = diff_cubics(x, y, h, n)
        elif(f == 'G'):
            if(g == 4):
                W = gaussian(x, y, h, n)
            else:
                W = diff_gaussian(x, y, h, n)
        else:
            print 'entered wrong option'
        return W


def cubicspline(x, y, h, n):
    W_val = []
    for i in range(n):
        q = abs(y-x[i])/h
        if(q < 1.0):
            W = (1.0/h)*(2.0/3-q**2+(q**3)/2.0)
        elif(1.0 < q < 2.0):
            W = (1.0/(6*h))*(2-q)**3
        else:
            W = 0
        W_val.append(W)
#    if(y == x[10]):
#        plt.plot(x, W_val)
#        plt.show()
    return W_val


def gaussian(x, y, h, n):
    W_val = []
    for i in range(n):
        q = abs(y-x[i])/h
        W = 1/(np.sqrt(np.pi)*h)*(math.exp(-q**2))
        W_val.append(W)
#    if(y == x[10]):
#        plt.plot(x, W_val)
#        plt.show()
    return W_val


def diff_cubics(x, y, h, n):
    W_val = []
    for i in range(n):
        q = abs(y-x[i])/h
        if(q < 1.0):
            W = (-2*q/h**2)+(3*q**2)/(2*h**2)
        elif(1.0 < q < 2.0):
            W = (-1/(2*h**2))*(2-q)**2
        else:
            W = 0
        W_val.append(W)
#    if(y == x[10]):
#       plt.plot(x, W_val)
#       plt.show()
    return W_val


def diff_gaussian(x, y, h, n):
    W_val = []
    for i in range(n):
        q = abs(y-x[i])/h
        W = -2*q/(np.sqrt(np.pi)*h**2)*(math.exp(-q**2))
        W_val.append(W)
#    if(y == x[10]):
#        plt.plot(x, W_val)
#        plt.show()
    return W_val


def cal_rho(n, h, f, x):
    rho = []
    m_j = 1.0
    g = 4
    for i in range(n):
        a = 0
        W = kernel_decis(g, f, x, x[i], h, n)
        for j in range(n):
            a += (m_j*W[j])
        rho.append(a)
    return rho


def cal_fun_val(n, h, f, x):
    rho = cal_rho(n, h, f, x)
    m_j = 1.0
    fun = []
    g = 4
    for i in range(n):
        b = 0
        W = kernel_decis(g, f, x, x[i], h, n)
        for j in range(n):
            b += ((m_j/rho[j])*np.sin(2*np.pi*x[j])*W[j])
        fun.append(b)
#    plt.plot(x, np.sin(2*np.pi*x))
#    plt.plot(x, fun)
#    plt.show()
    return fun


def cal_fun_diff(n, h, f, x):
    rho = cal_rho(n, h, f, x)
    diff_fun = []
    z = []
    m_j = 1.0
    g = 5
    for i in range(n):
        l = 0
        W = kernel_decis(g, f, x, x[i], h, n)
        for j in range(n):
            if(i != j):
                l += 1*(m_j/rho[j])*np.sin(2*np.pi*x[j])*W[j] *\
                                (abs(x[i]-x[j])/(x[i]-x[j]))
        diff_fun.append(l)
        z.append(2*np.pi*np.cos(2*np.pi*x[i]))
#    plt.plot(x, z)
#    plt.plot(x, diff_fun)
#    plt.show()
    return diff_fun


def error(n, h, f, x):
    c = 0
    m = 0
    fun = cal_fun_val(n, h, f, x)
    diff_fun = cal_fun_diff(n, h, f, x)
    for i in range(n):
        c += (((fun[i]-np.sin(2*np.pi*x[i]))/1.0)**2)
        m += (((diff_fun[i]-2*np.pi*np.cos(2*np.pi*x[i]))/1.0)**2)
    return np.sqrt(c), np.sqrt(m)


def plot():
    p = 10
    (h, dx, n, f, X_min, X_max) = initialize()
    for i in range(len(n)):
        x = np.arange(X_min, X_max, dx[i])
        err = []
        diff_err = []
        for j in range(len(h)):
            (er, diff_er) = error(n[i], h[j], f, x)
            err.append(er)
            diff_err.append(diff_er)
        plt.figure(i)
        plt.plot(h/dx[0], err, label="for no. of points %s" % n[i])
        plt.xlabel("h value")
        plt.ylabel("error")
        plt.title("error vs h value for %s no of points using %s kernel"
                  % (n[i], f), fontsize=10)
        plt.legend(loc="best")
        plt.savefig('%s.png' % i)
        plt.figure(p)
        plt.plot(h/dx[0], diff_err, label="for no. of points %s" % n[i])
        plt.xlabel("h value")
        plt.ylabel("differential error")
        plt.title("differential error vs h value for %s no of points \
using %s kernel" % (n[i], f), fontsize=10)
        plt.legend(loc="best")
        plt.savefig('%s.png' % p)
        p += 1
plot()
