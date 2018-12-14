print('dose dyo akeraious h exit gia exodo')

x = int(input('dose 1o akeraio: '))

z = int(input('dose 2o akeraio: '))


while type(x) != int or type(z) != int:
    print('dose akeraio h exit gia exodo')
    x = input('dose 1o akeraio: ')
    z = input('dose 2o akeraio: ')
    if x == 'exit' or z == 'exit':
        break
else:
    s = x+z
    p = x*z
    print('to athroisma einai {} kai to ginomeno {}'.format(s,p))
