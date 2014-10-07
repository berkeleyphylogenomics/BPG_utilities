
def remove_duplicates(seq): 
  # order preserving
  noDupes = []
  [noDupes.append(i) for i in seq if not noDupes.count(i)]
  return noDupes

