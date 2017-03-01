import threading
from threading import Thread


def func1():
    i= 0
    while i < 5:
        print ('Working thread 1 ' + str(i))
        i += 1   
    

def func2():
    i= 0
    while i < 5:
        print ('Working thread 2 ' + str(i))
        i += 1

t1 = Thread(target = func1)
t1.daemon = True
t1.start()

#t1.join()

t2 = Thread(target = func2)
t2.daemon = True
t2.start()


#t2.join()

print("done!")

