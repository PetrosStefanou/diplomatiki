
import my_test_functions as mtf 
import matplotlib.pyplot as plt

fig = plt.figure()
a = []
b = []
x = []
y = []
for i in range(10):
    a.append(i)
    b.append(i + 1)

    x.append(mtf.add(a[i], b[i]))

    y.append(mtf.square(a[i]))
    print('the sum is {:1.2f} and the square is {:1.2f}'.format(x[i], y[i]))

plt.plot(a, y)
plt.show()
   

