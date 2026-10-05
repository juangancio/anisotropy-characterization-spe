import numpy as np
import math

def perm_indices(ts, wl, lag):
    m = len(ts) - (wl - 1)*lag
    indcs = np.zeros(m, dtype=int)
    for i in range(1,wl):
        st = ts[(i - 1)*lag : m + ((i - 1)*lag)]
        for j in range(i,wl):
            zipped=zip(st,ts[j*lag : m+j*lag])
            indcs += [x > y for (x, y) in zipped]
            #print(st > ts[j*lag : m+j*lag])
        indcs*= wl - i
    return indcs + 1

def entropy(probs):
    h=0
    for i in range(len(probs)):
        if probs[i]==0:
            continue
        else:
            h=h-probs[i]*np.log(probs[i])
    return h

def probabilities(code,L):
    get_indexes = lambda x, xs: [k for (y, k) in zip(xs, range(len(xs))) if x == y]
    probs=[]
    for i in range(1,math.factorial(L)+1):
 
        probs=probs + [len(get_indexes(i,code))/len(code)
                   ]
        
    return probs


def h_entropy(data,L,lag):
    
    code=[]
    for i in range(len(data[:,0])):
        #new_data=[data[i,j] for j in range(len(data[i,:])) if data[i,j]>-1000]
        #code.extend(perm_indices(new_data,L,lag))
        code.append(perm_indices(data[i,:],L,lag))
    #print(len(np.array(code).reshape(1,np.size(code))[0]))
    return entropy(probabilities(np.array(code).reshape(1,np.size(code))[0],L))/np.log(math.factorial(L))

def v_entropy(data,L,lag):
    
    code=[]
    for i in range(len(data[0,:])):
        #new_data=[data[j,i] for j in range(len(data[:,i])) if data[j,i]>-1000]
        #code.extend(perm_indices(new_data,L,lag))
        code.append(perm_indices(data[:,i],L,lag))
    
    return entropy(probabilities(np.array(code).reshape(1,np.size(code))[0],L))/np.log(math.factorial(L))


def diag1_entropy(data,L,lag):
    # This is designed for a square array of data
    code = []
    N = data.shape[0]
    for i in range(N-1):
        diag = [data[i+k,k] for k in range(N-i)]
        #print(len(diag))
        if len(diag)<L+(L-1)*(lag-1):continue
        else:code.extend(perm_indices(diag,L,lag))
    for i in range(1,N-1):
        diag = [data[k,i+k] for k in range(N-i)]
        if len(diag)<L+(L-1)*(lag-1):continue
        else:code.extend(perm_indices(diag,L,lag))

    return entropy(probabilities(np.array(code),L))/np.log(math.factorial(L))

def diag2_entropy(data,L,lag):
    # This is designed for a square array of data
    data = np.rot90(data,k=1)
    code = []
    N = data.shape[0]
    for i in range(N-1):
        diag = [data[i+k,k] for k in range(N-i)]
        if len(diag)<L+(L-1)*(lag-1):continue
        else:code.extend(perm_indices(diag,L,lag))
    for i in range(1,N-1):
        diag = [data[k,i+k] for k in range(N-i)]
        if len(diag)<L+(L-1)*(lag-1):continue
        else:code.extend(perm_indices(diag,L,lag))

    return entropy(probabilities(np.array(code),L))/np.log(math.factorial(L))


def h_probs(data,L,lag):
    
    code=[]
    for i in range(len(data[:,0])):
        code.append(perm_indices(data[i,:],L,lag))
    #print(len(np.array(code).reshape(1,np.size(code))[0]))
    return probabilities(np.array(code).reshape(1,np.size(code))[0],L)

def v_probs(data,L,lag):
    
    code=[]
    for i in range(len(data[0,:])):
        code.append(perm_indices(data[:,i],L,lag))
    #print(len(np.array(code).reshape(1,np.size(code))[0]))
    return probabilities(np.array(code).reshape(1,np.size(code))[0],L)

def h_code(data,L,lag):
    
    code=[]
    for i in range(len(data[:,0])):
        code.append(perm_indices(data[i,:],L,lag))
    #print(len(np.array(code).reshape(1,np.size(code))[0]))
    return probabilities(np.array(code).reshape(1,np.size(code))[0],L)

def v_code(data,L,lag):
    
    code=[]
    for i in range(len(data[0,:])):
        code.append(perm_indices(data[:,i],L,lag))
    #print(len(np.array(code).reshape(1,np.size(code))[0]))
    return np.array(code)