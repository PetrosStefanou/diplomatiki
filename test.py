import my_test_functions as mtf 
import matplotlib.pyplot as plt

for i in range(10):
    a = i
    b = i + 1

    x = mtf.add(a, b)

    y = mtf.square(a)

    fig = plt.figure()
    plt.plot(a, y)
    plt.show()

    print('the sum is {:1.2f} and the square is {:1.2f}'.format(x, y))

