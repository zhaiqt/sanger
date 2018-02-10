def count_singlekey(inputDict, keyword):
   # sample input
   # inputDict = {
   #   abName1: { dna: 'atgc', protein: 'x' }
   #   abName2: { dna: 'ctga', protein: 'y' }
   # }

   countDict = {}
   for abName, abInfo in inputDict.iteritems():
     if countDict.has_key(abInfo[keyword]):
       countDict[abInfo[keyword]][1] += 1
     else:
       countDict[abInfo[keyword]] = [abName, 1]
   return countDict


def count_multikey(inputDict, keywords):
  # sample input
  # inputDict = {
  #   abName1: { dna: 'atgc', protein: 'x' }
  #   abName2: { dna: 'ctga', protein: 'y' }
  # }
  #keywords = list(keywords)
  keywords.sort()
  keywords = tuple(keywords)
  countDict = {}
  for abName, abInfo in inputDict.iteritems():
    combinedKey = []
    for k in keywords:
      combinedKey.append(abInfo[k])
    combinedKey = tuple(combinedKey)
    if countDict.has_key(combinedKey):
      countDict[combinedKey][1] += 1
    else:
      countDict[combinedKey] = [abName, 1]
  return countDict



