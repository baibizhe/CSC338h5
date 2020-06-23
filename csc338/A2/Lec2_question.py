sum=0
for i in range(23):
    sum+=2**i
base = 0
for i in range(8):
    base +=2**i
print(2**127)
print(2**(-126))
print(2**(-149))

# print(base**sum)
##??? ocnvert 0[11100..0]111111111...11111 to base 10
e , eps = 1,1
while 1+e >1:
    eps = e
    e = e/10
print(eps)
from decimal import *
getcontext().prec = 6
a,b = Decimal(1.2345e12),Decimal(124)
print(a+b==a,a-a+b)

#