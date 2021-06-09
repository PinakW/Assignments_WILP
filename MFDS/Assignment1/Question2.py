import numpy as np
import sys
import copy
def gen_mat(size, multiplier,diag_dom=False):
  ncols = size + 1
  nrows = size
  Positive_number_guard = 5
  A_aug = multiplier * np.random.rand(nrows,ncols)
  A = copy.deepcopy(A_aug)
  if diag_dom == True:
    print("Generating diagonally dominant matrix")
  elif diag_dom == False:
    print("NOT generating a diagonally dominant matrix")
  for i in range(nrows):
    if A[i][i] == 0.0:
      sys.exit("Diagonal element was generated to be 0.0 by the random generator\nCaught in gen_mat\nPlease re-run and hope for the best")
    abs_arr = np.abs(A[i])
    sum = np.sum(abs_arr)
    if sum > (2 * abs(A[i][i])):
      diff = sum - (2 * abs(A[i][i]))
        #This means that sum of non diagonal elements is greater than A[i][i]. This means the matrix is not diagonally dominant
      if diag_dom == True:
        #corrective action for making this row diagonally domminant
        #if A[i][i] is -ve number, we need to decrease it further by diff, so that its absolute value is increased
        if A[i][i] < 0:
          A[i][i] = A[i][i] - diff - Positive_number_guard
        elif A[i][i] > 0:
          A[i][i] = A[i][i] + diff + Positive_number_guard
      #else:
        #diagonal dominance is not requested, hence, no need to take any corrective action.
    elif sum <= (2 * abs(A[i][i])):
      diff =  (2 * abs(A[i][i])) - sum
      #This means abs(A[i][i]) is greater than non diagonal element addition.
      #if we need diagonal dominance then no corrective action needed
      if diag_dom == False:
        #if diagonal dominance is not requested, corrective action is needed
        #if A[i][i] is -ve number, we need to increase it further by diff, so that its absolute value is decreased
        if A[i][i] < 0:
          A[i][i] = A[i][i] + diff + Positive_number_guard
        elif A[i][i] > 0:
          A[i][i] = A[i][i] - diff - Positive_number_guard
  print("Generated A_aug")
  print(A)
  return A

def diag_1(A):
  global round_factor
  nrow = len(A[:,0])
  for i in range(nrow):
    A[i] = np.round(A[i]/A[i][i],decimals=round_factor)
  return A
def decompose(A):
  #A = diag_1(A)
  nrow = len(A[:,0])
  L = np.tril(A)
  U = np.triu(A)
  for i in range(nrow):
    L[i][i] = 0.0
    U[i][i] = 0.0
  return L,np.identity(nrow),U
def check_conv(x_new, x, tolerance):
  length = len(x)
  x_diff = np.zeros(length)
  for i in range(length):
    x_diff[i] = abs(x_new[i] - x[i])
  x_max_idx = abs(x_diff).argmax()
  xmax = x_diff.flat[x_max_idx]
  if xmax < (tolerance * x_new[x_max_idx]):
    return True
  else:
    return False
def Gauss_Seidel(A_aug, x, tolerance, N):
  global round_factor
  A_aug = diag_1(A_aug)
  #split A_aug int A and b
  ncol = len(A_aug[0,:])
  #b = A_aug[:,(ncol-1)].reshape(nrows,1)           #Reshape only for Gauss Jacobi
  b = A_aug[:,(ncol-1)]
  A = A_aug[:,0:(ncol - 1)]
  L, I, U = decompose(A)
  nrow = len(A[:,0])
  ncol = len(A[0,:])
  size = min(nrow,ncol)
  #Create x_new
  #temp = np.zeros(size).reshape(size,1)
  print("Matrix A")
  print(A)
  print("Vector b")
  print(b)
  temp = np.zeros(size)
  x_new = copy.deepcopy(x)
  length = len(x)
  convergence = False
  iter = 0
  while ((convergence == False) and (iter <= N)):
    for j in range(length):
      #print(str(b[j]) + " - " + str(L[j]) +" * " + str(x_new) + "-" + str(U[j]) + "*" + str(x))
      temp[j] = b[j] - np.dot(L[j], x_new) - np.dot(U[j],x)
      #Update immediately for Seibel
      x_new[j] = temp[j]
    x_new = np.around(x_new, decimals=round_factor)
    print("x after iteration " + str(iter))
    #print(x)
    print(x_new)
    print("")
    convergence = check_conv(x_new, x, tolerance)
    x = copy.deepcopy(x_new)
    iter +=1
  if convergence == False:
    print("The Gauss_Seidel Algorithm did not converge")
  else:
    print("iter = " + str(iter))
    return x_new
def Gauss_Jacobi(A_aug,x, tolerance, N):
  global round_factor
  A_aug = diag_1(A_aug)
  ncol = len(A_aug[0,:])
  #b = A_aug[:,(ncol-1)].reshape(nrows,1)           #Reshape only for Gauss Jacobi
  b = A_aug[:,(ncol-1)]
  A = A_aug[:,0:(ncol - 1)]
  print("Input for Gauss Jacobi")
  print(A)
  nrow = len(A[:,0])
  ncol = len(A[0,:])
  I = np.identity(nrow)
  size = min(nrow,ncol)
  temp = np.zeros(size)
  x_new = copy.deepcopy(x)
  length = len(x)
  convergence = False
  iter = 0
  mat_diff = copy.deepcopy(I - A)
  while ((convergence == False) and (iter <= N)):
    for j in range(length):
      #print(str(b[j]) + " - " + str(L[j]) +" * " + str(x_new) + "-" + str(U[j]) + "*" + str(x))
      temp[j] = b[j] + np.dot(mat_diff[j],x)
    #Update after one iteration for Jacobi
    x_new = copy.deepcopy(temp)
    x_new = np.round(x_new,decimals=round_factor)
    print("x after iteration " + str(iter))
    #print(x)
    print(x_new)
    print("")
    convergence = check_conv(x_new, x, tolerance)
    x = copy.deepcopy(x_new)
    iter +=1
  if convergence == False:
    print("The Gauss_Jacobi Algorithm did not converge")
  else:
    print("iter = " + str(iter))
    return x_new

#MAIN
size = 3
round_factor = 3
#We will use augmented matrix
ncols = size + 1
nrows = size
multiplier = 10
#np.array([[6, -2, 1, 11],[-2, 7, 2, 5],[1, 2, -5, -1]])#
A_aug = gen_mat(size, multiplier,diag_dom=True)#np.array([[1.0, -0.25, -0.25, 0, 50],[-0.25,1, 0, -0.25, 50],[-0.25, 0, 1, -0.25, 25],[0, -0.25, -0.25, 1, 25]])#multiplier * np.random.rand(nrows,ncols)
A_aug = A_aug.astype(float)
#nrow = len(A_aug[:,0])
#ncol = len(A_aug[0,:])
#for i in range(nrow):
#  A_aug[i] = np.round(A_aug[i],decimals=round_factor)

#round it off
#b = multiplier * np.random.rand(ncols).reshape(ncols,1)
#x = np.zeros(ncols).reshape(ncols,1)
#x_new = np.zeros(ncols).reshape(ncols,1)

entries = list(map(float, input("Enter a space separated list of " + str(size) + " elements: ").split()))
try:
  #x is a column vector having number of rows = number of columns of A, number of columns = 1 (cause column vector)
  #x = np.array(entries).reshape(size, 1)   #reshaping will be required for Gauss_jacobi to speed it up. As we have direct matrix multiplication there
  x = np.array(entries) #For Gauss Seibel, do not reshape to a column
except:
  sys.exit("Input array length does not match " + str(size) + "\nExitting...")
print("The initial input values are:")
print(x)
tolerance = float(input("Enter the % tolerance value: "))
tolerance = tolerance/100
N = int(input("Enter the maximum iteration count: "))
#m = [1, 2, 3]
#m = np.array(m)
#print("m")
#print(m)
#print("np.matmul(m,x)")
#print(np.matmul(m,x))
A_aug_orig = copy.deepcopy(A_aug)
x_orig = copy.deepcopy(x)
y = Gauss_Seidel(A_aug,x, tolerance, N)
print("\n\nOutput of Gauss Siedel is: ")
print(y)
A_aug = copy.deepcopy(A_aug_orig)
x = copy.deepcopy(x_orig)
y_jac = Gauss_Jacobi(A_aug,x, tolerance, N)
print("\n\nOutput of Gauss Jacobi is: ")
print(y_jac)
