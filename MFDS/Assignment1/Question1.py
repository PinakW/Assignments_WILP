import numpy as np
import sys
import copy

def ret_sigbit(x,dimension=0):
  global d
  if dimension==0:
    try:
      return np.around(x,decimals=(d - int(np.floor(np.log10(abs(x)))) - 1))
    except:
      if x == 0.0:
        return 0.0
      else:
        sys.exit("exception in ret_sigbit" + str(x))
  elif dimension==1:
    length = len(x)
    for i in range(length):
      try:
        x[i] = np.around(x[i],decimals=(d - int(np.floor(np.log10(abs(x[i])))) - 1))
      except:
        if x[i] ==0.0:
          x[i]= 0.0
        else:
          sys.exit("EXCEPTION in ret_sigbit" + str(x[i]))
    return x
  elif dimension == 2:
    col = len(mat[0,:])
    row = len(mat[:,0])
    for i in range(row):
      for j in range(col):
        try:
          x[i,j] = np.around(x[i,j],decimals=(d - int(np.floor(np.log10(abs(x[i,j])))) - 1))
        except:
          if x[i,j] == 0.0:
            x[i,j] = 0.0
          else:
            sys.exit("EXCEPTION in ret_sigbit " + str(x[i,j]))
    return x
def swap(row_1, row_2):
	row_1 = row_1 - row_2
	row_2 = row_1 + row_2
	row_1 = row_2 - row_1
	return row_1, row_2
def partial_pivot(mat, curr_index):
  ncol = len(mat[0,:])
  nrow = len(mat[:,0])
  if(curr_index>nrow):
    sys.exit("Some miscalculations in partial_pivot")
  column_under_consideration = mat[curr_index:nrow,curr_index]
  max_idx = np.argmax(abs(column_under_consideration))
  if max_idx !=0:
    mat[curr_index],mat[(curr_index + max_idx)] = swap(mat[curr_index],mat[(curr_index + max_idx)])
  return mat
def GE(mat,pivot=False):
  #col = len(mat[0,:])
  row = len(mat[:,0])
  for i in range(row - 1):
    #Do Pivoting if it is specified
    if pivot==True:
      mat = partial_pivot(mat, i)
      print("Pivot result")
      print(mat)
    for j in range(i+1,row):
      try:
        scale_factor = mat[j][i]/mat[i][i]
        #important we need to round the scale_factor to significant values as well
        scale_factor = ret_sigbit(scale_factor)
      except:
        sys.exit("Illegal math operation at " + str(i) +", " + str(j) + "\n" +str(mat[j][i]) + " / "+ str(mat[i][i])+"\nThe output of pivot = " + str(pivot) + " is wrong, possibly due to pivot element = zero")
      temp_subtractor = scale_factor * mat[i]
      temp_subtractor = ret_sigbit(temp_subtractor,dimension=1)
      mat[j] = mat[j] - temp_subtractor
  mat = ret_sigbit(mat,dimension=2)
  print("After Gaussian Elimination")
  print(mat)
  return mat
def back_sub(arr):
  #Assumption mat is augmented matrix
  ncol = len(mat[0,:])
  size =  ncol - 1
  y = np.zeros(size)

  y[size-1] = arr[size-1][size]/arr[size-1][size-1]
  
  for i in range(size-2,-1,-1):
    y[i] = arr[i][size]
    for j in range(i+1,size):
      temp_prod = ret_sigbit(arr[i][j]*y[j])
      y[i] = y[i] - temp_prod
      
    y[i] = y[i]/arr[i][i]
    
  y = ret_sigbit(y,dimension=1)
  return y

def solve(mat):
  mat_orig = copy.deepcopy(mat)
  print("Matrix after significant digit changes:")
  print(mat_orig)
  solution = back_sub(GE(mat, pivot=True))
  print("Result without Pivoting:")
  print(solution)
  mat = copy.deepcopy(mat_orig)
  solution = back_sub(GE(mat))
  print("Result with Pivoting:")
  print(solution)

#main
size = 4
multiplier =  100
d = int(input("Enter significant digits: "))
#mat = np.array([[0, 2, 0, 1, 0],[2, 2, 3, 2, -2], [4, -3, 0, 1, -7], [6, 1, -6, -5, 6]])
mat = multiplier * np.random.rand(size,size +1)
mat = mat.astype(float)
#print("Generated random matrix of size " + str(size) + " x " +str(size) + " :")
#print(mat)
solve(ret_sigbit(mat,dimension=2))
