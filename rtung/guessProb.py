# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 22:36:18 2015

@author: Ray
"""
from datetime import datetime
import random

length = 10000
a = .1758
t = .1852
c = .3200
g = .3190
n = 1000000
pal = "ATAT"

def generateChar(a, c, t, g):
    int1 = random.random()
    if int1 < a:
        return "A"
    if int1 < a+c:
        return "C"
    if int1 < a+c+t:
        return "T"
    return "G"
    
def strand(n, a, c, t, g, pal):
    str1 = generateChar(a, c, t, g) + generateChar(a, c, t, g) + generateChar(a, c, t, g) + generateChar(a, c, t, g)
    if pal == str1:
        return 0
    for i in range(n - 4):
        str1 = str1[1:] + generateChar(a, c, t, g)
        if pal == str1:
            return 0
    return 1

print(datetime.now().time().isoformat())

count = 0
for i in range(n):
    count = count + strand(length, a, c, t, g, pal)
print(count)

print(datetime.now().time().isoformat())