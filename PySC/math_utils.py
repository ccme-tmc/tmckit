

def fixeddensity(thelen, density):
  if thelen == density:
    yield [1]*thelen
  elif density == 0:
    yield [0]*thelen
  else:
    for leftlist in fixeddensity(thelen-1, density):
      yield leftlist + [0]
    for leftlist in fixeddensity(thelen-1, density-1): \
      yield leftlist + [1]

def Math_ksubset(n,k):
  """
  Generate all possible k-th order subsets of {1,2,3,...,n}
  """
  s = range(1,n+1)

  ksub=[]
  for it in fixeddensity(n,k):
    ksub.append([s[i] for i in range(n) if it[i]==1])

  return ksub 
