import numpy as np

ms = [80.3786609, 125.003415, 6598.84567, 6598.88561, 6599.48329, 70000]
diff = np.diff(ms)
print(diff)
lss = []
ls = []
tol = 800
ls.append(ms[0])
for i, d in enumerate(diff):
    print("indexs")
    print(i)
    print(d)
    if d < tol:
        ls.append(ms[i+1])
        if(i == len(ms) - 2):
        	print("got in here")
        	lss.append(ls)
   

    else:
        lss.append(ls)
        ls = []
        ls.append(ms[i+1])
        if(i == len(ms) - 2):
        	lss.append(ls)



for i in lss:
	print(i)
print(lss)