import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

'''Τεστ για την ορθότητα του αναδρομικού υπολογισμού της παραγώγου στη διαφορική εξίσωση.

Φαίνεται να έχει πολύ μικρή απόκλιση με τον ορθό υπολογισμό ωστόσο είναι απαραίτητο ο αναδρομικός όρος να είναι κατάλληλα διαιρεμένος ώστε να είναι μικρός.

Σε διαφορετική περίπτωση κρασάρει.

Εξάλλου για μεγάλα νούμερα της ανεξάρτητης μεταβλητής η απόκλιση μεγαλώνει αλλά όχι σημαντικά
'''

def right_way(y1, t):

    dy1dt = 100.*(t + t**2 + t**3)/101.

    return dy1dt

global der

der = [0.]


def wrong_way(y2, t):
    
    dy2dt = (t + t**2 + t**3) - der[-1]/100.
    
    der.append(dy2dt)
    
    # print(der[-1])
    
    return dy2dt    



t = np.linspace(-1., 1., 100)

y0 = 0.

y1 = odeint(right_way, y0, t)
y2 = odeint(wrong_way, y0, t)


fig = plt.figure()
plt.plot(t, y1)
plt.plot(t, y2, '--')
fig1 = plt.figure()
plt.plot(t, np.abs(y1-y2), 'r')
plt.show()